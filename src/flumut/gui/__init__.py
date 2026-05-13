def main():
    try:
        from flumut.gui.LauncherWindow import launch_gui
    except ImportError:
        raise SystemExit('GUI dependencies are not installed.\nRun:  pip install flumut[gui]')
    launch_gui()
