import collections
import logging
import re

from PySide6 import QtGui
from PySide6.QtCore import Qt, QThread, Signal
from PySide6.QtWidgets import QDialog, QMessageBox, QProgressBar, QPushButton, QTextEdit, QVBoxLayout

from flumut.core.logger import LOGGER
from flumut.core.options import FluMutOptions
from flumut.core.workflows import whole_workflow
from flumut.gui.globals import ICON_PATH


class StdIO:
    def __init__(self) -> None:
        self.buffer = collections.deque()

    def readline(self):
        if len(self.buffer) == 0:
            return None
        return self.buffer.popleft()

    def write(self, value):
        val = value.strip()
        if val:
            self.buffer.append(val)

    def flush(self):
        pass

    readable = writable = lambda self: True


class FluMutWorker(QThread):
    started = Signal(int)
    error = Signal(Exception)
    ended = Signal(int)

    def __init__(self, options: FluMutOptions, stderr_stream):
        QThread.__init__(self)
        self.options = options
        self.stderr_stream = stderr_stream

    def __del__(self):
        self.wait()

    def terminate(self) -> None:
        self.ended.emit(2)
        return super().terminate()

    def run(self):
        total_sequences = 0
        for fasta in self.options.input.fasta_files:
            total_sequences += len(re.findall(r'^>.+', fasta.read(), re.MULTILINE))
            fasta.seek(0)

        self.started.emit(total_sequences)

        handler = logging.StreamHandler(self.stderr_stream)
        handler.setFormatter(logging.Formatter('[%(levelname)s] %(message)s'))
        LOGGER.addHandler(handler)
        LOGGER.propagate = False
        LOGGER.setLevel(logging.DEBUG)

        try:
            whole_workflow(self.options)
            self.ended.emit(0)
        except Exception as e:
            self.ended.emit(1)
            self.error.emit(e)
        finally:
            LOGGER.removeHandler(handler)
            LOGGER.propagate = True


class FluMutOutputReader(QThread):
    stderr = Signal(str)

    def __init__(self, stderr_stream):
        QThread.__init__(self)
        self._stop = False
        self.stderr_stream: StdIO = stderr_stream

    def __del__(self):
        self.wait()

    def run(self):
        while not self._stop:
            line = self.stderr_stream.readline()
            if line:
                self.stderr.emit(line)
            else:
                self.msleep(100)

    def stop(self):
        """End the loop, and block until it has ended.

        The thread must not outlive the window it reports to: closing that
        window frees the C++ object this thread is still running on. Setting
        the flag alone leaves it up to one `msleep` still running.
        """
        self._stop = True
        self.wait()


class ProgressWindow(QDialog):
    def __init__(self, options: FluMutOptions) -> None:
        super().__init__()
        self.init_ui()
        self.setModal(True)
        self.setWindowFlag(Qt.WindowType.WindowMinimizeButtonHint)
        self.setWindowFlag(Qt.WindowType.WindowMaximizeButtonHint)
        self.setWindowFlag(Qt.WindowType.WindowCloseButtonHint)

        self.flumut_options = options
        self.start_flumut()

    def init_ui(self):
        layout = QVBoxLayout()

        self.setLayout(layout)
        self.setWindowTitle('Executing FluMut')
        self.setWindowIcon(QtGui.QIcon(ICON_PATH))
        self.setMinimumWidth(450)
        self.setMinimumHeight(300)

        self.progress_bar = QProgressBar()
        self.progress_bar.setTextVisible(False)
        layout.addWidget(self.progress_bar)

        self.log_txt = QTextEdit()
        self.log_txt.setReadOnly(True)
        layout.addWidget(self.log_txt)

        self.cancel_btn = QPushButton('Cancel')
        layout.addWidget(self.cancel_btn)
        self.cancel_btn.clicked.connect(self.cancel_flumut)

    def start_flumut(self):
        def handle_start(total_sequences):
            self.progress_bar.setRange(0, 5)
            self.progress_bar.setValue(0)
            # self.log_txt.append('Starting FluMut analysis...')
            # self.log_txt.append(f'Detected {total_sequences} sequences.')

        def handle_end(exit_code):
            self.logger_thread.stop()
            while (line := self.stderr_stream.readline()) is not None:
                log_stderr(line)

            if exit_code == 0:
                self.log_txt.setTextColor(self.log_txt.palette().color(QtGui.QPalette.ColorRole.Text))
                self.log_txt.append('FluMut terminated successfully.')
                self.set_progress_bar_color(98, 201, 27)
                self.progress_bar.setValue(self.progress_bar.maximum())
            elif exit_code == 1:
                self.log_txt.setTextColor(Qt.GlobalColor.red)
                self.log_txt.append('FluMut terminated with errors.')
                self.set_progress_bar_color(238, 1, 1)
            elif exit_code == 2:
                self.log_txt.setTextColor(Qt.GlobalColor.red)
                self.log_txt.append('FluMut terminated by the user.')
                self.set_progress_bar_color(238, 1, 1)
            else:
                self.log_txt.setTextColor(Qt.GlobalColor.red)
                self.log_txt.append('FluMut terminated with unknown exit code.')
            self.cancel_btn.setText('Close')

        def handle_error(error):
            self.set_progress_bar_color(238, 1, 1)
            self.log_txt.setTextColor(Qt.GlobalColor.red)
            self.log_txt.append(f'{error.__class__.__name__}: {error!s}')
            QMessageBox.warning(self, error.__class__.__name__, str(error))

        def log_stderr(line: str):
            if line.startswith(('[WARNING]', '[ERROR]', '[CRITICAL]')):
                self.log_txt.setTextColor(Qt.GlobalColor.red)
            else:
                self.log_txt.setTextColor(self.log_txt.palette().color(QtGui.QPalette.ColorRole.Text))
            if line.startswith('[INFO]'):
                self.progress_bar.setValue(self.progress_bar.value() + 1)
            self.log_txt.append(line)

        self.stderr_stream = StdIO()

        self.flumut_thread = FluMutWorker(self.flumut_options, self.stderr_stream)
        self.flumut_thread.started.connect(handle_start)
        self.flumut_thread.ended.connect(handle_end)
        self.flumut_thread.error.connect(handle_error)

        self.logger_thread = FluMutOutputReader(self.stderr_stream)
        self.logger_thread.stderr.connect(log_stderr)

        self.logger_thread.start()
        self.flumut_thread.start()

    def cancel_flumut(self):
        if self.cancel_btn.text() == 'Close':
            self.close()
        else:
            self.log_txt.setTextColor(self.log_txt.palette().color(QtGui.QPalette.ColorRole.Text))
            self.log_txt.append('Stopping FluMut analysis...')
            self.logger_thread.terminate()
            self.flumut_thread.terminate()

    def set_progress_bar_color(self, r, g, b):
        custom_palette = QtGui.QPalette()
        custom_palette.setColor(QtGui.QPalette.ColorRole.Highlight, QtGui.QColor(r, g, b))

        self.progress_bar.setPalette(custom_palette)
