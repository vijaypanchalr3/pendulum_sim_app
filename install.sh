#!/bin/bash

~/.local/bin/pyinstaller --onefile application.py

# {
    # wine ~/.wine/drive_c/python37/Scripts/pyinstaller.exe --onefile application.py
# }||{
    sudo apt install wine
    wget https://www.python.org/ftp/python-3.7.9.amd64.msi
    wine msiexec -i python-3.7.9.amd64.msi /qb
    cd ~/.wine/drive_c/python37
    wine python37.exe Scripts/pip.exe install pyinstaller
    wine ~/.wine/drive_c/python37/Scripts/pyinstaller.exe --onefile application.py
# }

