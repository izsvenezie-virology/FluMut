from importlib.resources import files

from flumut import gui

ICON_PATH: str = str(files(gui) / 'resources' / 'flumut_icon.ico')
