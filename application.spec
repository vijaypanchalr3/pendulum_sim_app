# -*- mode: python ; coding: utf-8 -*-


block_cipher = None


a = Analysis(
    ['application.py'],
    pathex=["/home/vijay/Projects/pendulum_sim_app"],
    binaries=[],
    datas=[],
    hiddenimports=[],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)
a.datas += [('Lato-BoldItalic.ttf','/home/vijay/Projects/pendulum_sim_app/Lato-BoldItalic.ttf', "DATA")]
a.datas += [('bitmap1.png','/home/vijay/Projects/pendulum_sim_app/bitmap1.png', "DATA")]
a.datas += [('bitmap2.png','/home/vijay/Projects/pendulum_sim_app/bitmap2.png', "DATA")]

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.zipfiles,
    a.datas,
    [],
    name='pensim',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
