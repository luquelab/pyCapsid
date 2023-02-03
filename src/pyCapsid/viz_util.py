"""Functions for visualizing results in chimeraX

"""

def chimeraxLaunchTest(chimerax_path='"C:\\Program Files\\ChimeraX\\bin"', file_path='.' ):
    import os
    chimerax_exe = chimerax_path + '\\ChimeraX.exe'
    os.system('"' + chimerax_exe + '"' + '--script ' + './chimerax_script.py')
