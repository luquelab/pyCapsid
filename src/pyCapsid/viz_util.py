"""Functions for visualizing results in chimeraX

"""

def chimeraxLaunchTest(labels, chimerax_path='C:\\Program Files\\ChimeraX\\bin', pdb_path='.', save_path='.',
                       script_path='../src/pyCapsid/scripts/chimerax_script.py', labels_path='.'):
    import os
    from numpy import save
    # from tempfile import NamedTemporaryFile
    save(labels_path, labels)
    print(labels_path)

    chimerax_exe = chimerax_path + '\\ChimeraX.exe'
    print(chimerax_exe)
    cmd_string = f'""{chimerax_exe}" --script "{script_path} {labels_path} {save_path} {pdb_path}""'
    print(cmd_string)
    os.system(cmd_string)
    print('???')

    # results = subprocess.run(['"' + chimerax_exe + '"', '--script',
    #                           f'"src/pyCapsid/scripts/chimerax_script.py {labels_path} {save_path} {pdb_path}"'],
    #                          stdout=subprocess.PIPE)
    # print(results.stdout)
