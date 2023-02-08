"""Functions for visualizing results in chimeraX

"""

def chimeraxViz(labels, pdb, remote=True, chimerax_path='C:\\Program Files\\ChimeraX\\bin', pdb_path='.', save_path='.',
                       script_path='../src/pyCapsid/scripts/chimerax_script.py'):
    """

    :param labels:
    :param chimerax_path:
    :param pdb_path:
    :param save_path:
    :param script_path:
    :param labels_path:
    :return:
    """

    import os
    from numpy import save
    from tempfile import NamedTemporaryFile
    with NamedTemporaryFile(suffix='.npy', delete=False) as temp_file:
        save(temp_file, labels)
        labels_path = temp_file.name.replace('\\','/') #os.path.abspath(temp_file.name)


    chimerax_exe = chimerax_path + '\\ChimeraX.exe'
    cmd_string = f'""{chimerax_exe}" --script "{script_path} {labels_path} {save_path} {pdb_path} {str(remote)} {pdb}""'
    print(cmd_string)
    os.system(cmd_string)
    temp_file.close()
    os.unlink(temp_file.name)

    # results = subprocess.run(['"' + chimerax_exe + '"', '--script',
    #                           f'"src/pyCapsid/scripts/chimerax_script.py {labels_path} {save_path} {pdb_path}"'],
    #                          stdout=subprocess.PIPE)
    # print(results.stdout)


# Adapted from py3dmol tutorial
# https://william-dawson.github.io/using-py3dmol.html
class Atom(dict):
    def __init__(self, line):
        self["type"] = line[0:6].strip()
        self["idx"] = line[6:11].strip()
        self["name"] = line[12:16].strip()
        self["resname"] = line[17:20].strip()
        self["resid"] = int(int(line[22:26]))
        self["x"] = float(line[30:38])
        self["y"] = float(line[38:46])
        self["z"] = float(line[46:54])
        self["sym"] = line[76:78].strip()

    def __str__(self):
        line = list(" " * 80)

        line[0:6] = self["type"].ljust(6)
        line[6:11] = self["idx"].ljust(5)
        line[12:16] = self["name"].ljust(4)
        line[17:20] = self["resname"].ljust(3)
        line[22:26] = str(self["resid"]).ljust(4)
        line[30:38] = str(self["x"]).rjust(8)
        line[38:46] = str(self["y"]).rjust(8)
        line[46:54] = str(self["z"]).rjust(8)
        line[76:78] = self["sym"].rjust(2)
        return "".join(line) + "\n"


class Molecule(list):
    def __init__(self, file):
        for line in file:
            if "ATOM" in line or "HETATM" in line:
                self.append(Atom(line))

    def __str__(self):
        outstr = ""
        for at in self:
            outstr += str(at)

        return outstr
    
# Colormap generation taken from stackexchange
# https://stackoverflow.com/questions/42697933/colormap-with-maximum-distinguishable-colours
def generate_colormap(number_of_distinct_colors: int = 80):
    import math

    import numpy as np
    from matplotlib.colors import ListedColormap
    from matplotlib.cm import hsv

    if number_of_distinct_colors == 0:
        number_of_distinct_colors = 80

    number_of_shades = 7
    number_of_distinct_colors_with_multiply_of_shades = int(math.ceil(number_of_distinct_colors / number_of_shades) * number_of_shades)

    # Create an array with uniformly drawn floats taken from <0, 1) partition
    linearly_distributed_nums = np.arange(number_of_distinct_colors_with_multiply_of_shades) / number_of_distinct_colors_with_multiply_of_shades

    # We are going to reorganise monotonically growing numbers in such way that there will be single array with saw-like pattern
    #     but each saw tooth is slightly higher than the one before
    # First divide linearly_distributed_nums into number_of_shades sub-arrays containing linearly distributed numbers
    arr_by_shade_rows = linearly_distributed_nums.reshape(number_of_shades, number_of_distinct_colors_with_multiply_of_shades // number_of_shades)

    # Transpose the above matrix (columns become rows) - as a result each row contains saw tooth with values slightly higher than row above
    arr_by_shade_columns = arr_by_shade_rows.T
    
     # Keep number of saw teeth for later
    number_of_partitions = arr_by_shade_columns.shape[0]

    # Flatten the above matrix - join each row into single array
    nums_distributed_like_rising_saw = arr_by_shade_columns.reshape(-1)

    # HSV colour map is cyclic (https://matplotlib.org/tutorials/colors/colormaps.html#cyclic), we'll use this property
    initial_cm = hsv(nums_distributed_like_rising_saw)

    lower_partitions_half = number_of_partitions // 2
    upper_partitions_half = number_of_partitions - lower_partitions_half

    # Modify lower half in such way that colours towards beginning of partition are darker
    # First colours are affected more, colours closer to the middle are affected less
    lower_half = lower_partitions_half * number_of_shades
    for i in range(3):
        initial_cm[0:lower_half, i] *= np.arange(0.2, 1, 0.8/lower_half)

    # Modify second half in such way that colours towards end of partition are less intense and brighter
    # Colours closer to the middle are affected less, colours closer to the end are affected more
    for i in range(3):
        for j in range(upper_partitions_half):
            modifier = np.ones(number_of_shades) - initial_cm[lower_half + j * number_of_shades: lower_half + (j + 1) * number_of_shades, i]
            modifier = j * modifier / upper_partitions_half
            initial_cm[lower_half + j * number_of_shades: lower_half + (j + 1) * number_of_shades, i] += modifier

    return ListedColormap(initial_cm)
  

import numpy as np
def open_pdb(pdb):
    with open(pdb + "_capsid.pdb") as ifile:
        mol = Molecule(ifile)
    return mol

    

  
  
def clusters_colormap_hexcolor(clusters):
    import matplotlib as mpl
    norm = mpl.colors.Normalize(vmin=np.min(clusters), vmax=np.max(clusters))
    cmap = generate_colormap(int(np.max(clusters)))
    print(cmap)
    rgba = cmap(norm(clusters))
    print(rgba*255)
    hexcolor = []
    for c in rgba:
      hexcolor.append(mpl.colors.rgb2hex(c))

    return hexcolor
      
def cluster_scheme(mol, hexcolor, clusters):
  r0 = mol[0]['resid']
  c0 = hexcolor[0]
  clust_scheme = []
  select = '@'
  print(r0, c0)
  i = 1
  j = 0
  for at in mol:
      r = at['resid']
      if r == r0:
          select += str(i) + ','
          i += 1
      else:
          clust_scheme.append([c0, select[:-1]])
          select = '@' + str(i) + ','
          r0 = r
          j +=1
          i +=1
          c0 = hexcolor[j]
          l0 = clusters[j]
  clust_scheme.append([c0, select[:-1]])

  return clust_scheme
  
def view_pdb_ngl(pdb, capsid, labels):
    mol = open_pdb(pdb)
    hexcolor = clusters_colormap_hexcolor(labels)
    clust_scheme = cluster_scheme(mol, hexcolor, labels)

    import nglview as ngl
    color_scheme = ngl.color._ColorScheme(clust_scheme, label="scheme_regions")
    view = ngl.show_prody(capsid, gui=False)
    view.clear_representations()


    view.add_representation("spacefill",  color=color_scheme)
    view._remote_call("setSize", target='Widget', args=['1000px','1000px'])
    return view

    
