def main():
    try:
        from flumut.gui.launcher_window import launch_gui
    except ImportError:
        raise SystemExit('GUI dependencies are not installed.\nRun:  pip install flumut[gui]')
    launch_gui()
