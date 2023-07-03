"""Functions for visualizing results in chimeraX

"""


def visualizeResults(pdb, capsid, labels, method='chimerax', **kwargs):
    if method == 'chimerax':
        chimeraxViz(labels, pdb, **kwargs)
    elif method == 'nglview':
        return view_pdb_ngl(pdb, capsid, labels, **kwargs)
    else:
        print('method must be chimerax or nglview')

def chimeraxViz(labels, pdb, remote=True, chimerax_path=None, pdb_path='.', save_path='.',
                rwb_scale=False, clust_mask=None, **kwargs):
    """Launches ChimeraX and runs a script that visualizes the results.

    :param labels:
    :param pdb:
    :param remote:
    :param chimerax_path:
    :param pdb_path:
    :param save_path:
    :param script_path:
    :param labels_path:
    :return:
    """

    import os
    import platform
    if chimerax_path is None:
        print('No chimerax path specified, checking in default locations for your OS')
        if platform.system()=='Linux':
            chimerax_path = 'usr/bin/chimerax'
        elif platform.system()=='Windows':
            chimerax_path = 'C:\\Program Files\\ChimeraX\\bin\\ChimeraX.exe'
        elif platform.system()=='Darwin':
            chimerax_path = '/Applications/ChimeraX'
        else:
            print('No chimerax path is given and cannot check default locations')
            chimerax_path = ''

    import matplotlib as mpl
    norm = mpl.colors.Normalize(vmin=np.min(labels), vmax=np.max(labels))

    if rwb_scale:
        import matplotlib.pyplot as plt
        cmap = plt.get_cmap('coolwarm_r')
        rgba = cmap(norm(labels)) * 255
    else:
        cmap = generate_colormap(int(np.max(labels)))
        rgba = cmap(norm(labels)) * 255

    if clust_mask is not None:
        rgba = addUnclusteredColor(rgba, clust_mask)

    from numpy import save
    from tempfile import NamedTemporaryFile
    with NamedTemporaryFile(suffix='.npy', delete=False) as temp_file:
        save(temp_file, rgba)
        rgba_path = os.path.normpath(temp_file.name)
        rgba_path = rgba_path.replace('\\', '/') #os.path.abspath(temp_file.name)
        print(rgba_path)

    # get path to chimerax script
    import pyCapsid.scripts as scpath
    import os

    script_path = os.path.abspath(scpath.__file__)
    script_path = script_path.replace('__init__', 'chimerax_script')

    chimerax_exe = chimerax_path #+ '\\ChimeraX.exe'
    if platform.system()=='Windows':
        cmd_string = f'""{chimerax_exe}" --script "{script_path} {rgba_path} {save_path} {pdb_path} {str(remote)} {pdb} {str(rwb_scale)}""'
    elif platform.system()=='Linux':
        cmd_string = f'{chimerax_exe} --script "{script_path} {rgba_path} {save_path} {pdb_path} {str(remote)} {pdb} {str(rwb_scale)}"'
    else:
        cmd_string = f'{chimerax_exe} --script "{script_path} {rgba_path} {save_path} {pdb_path} {str(remote)} {pdb} {str(rwb_scale)}"'
    print(cmd_string)
    os.system(cmd_string)
    temp_file.close()
    os.unlink(temp_file.name)

def visualizeSavedResults(pdb, results_file, n_cluster=None, method='chimerax', **kwargs):
    import numpy as np
    results = np.load(results_file)
    labels = results['labels']

    nc_range = results['nc_range']
    scores = results['score']

    if n_cluster is None:
        print('Defaulting to highest score clustering')
        ind = np.argmax(scores)
        n_c = nc_range[ind]
    else:
        ind = np.argwhere(nc_range == n_cluster)[0][0]
        n_c = n_cluster

    print(f'Visualizing cluster results of {pdb} for {n_c} clusters')
    labels = labels[ind]
    print(labels)
    if method=='chimerax':
        chimeraxViz(labels, pdb, **kwargs)
    elif method=='nglview':
        view = createCapsidView(pdb)
        return labels, view
    else:
        print('Method must be one of chimerax or nglview')


def addUnclusteredColor(rgba, clust_mask):
    from numpy import insert
    print(rgba.shape)
    print(rgba)
    black = np.array([0,0,0,255])
    for i,val in enumerate(clust_mask):
        if not val:
            rgba = insert(rgba, i, black, axis=0)
    print(rgba.shape)
    print(rgba)
    return rgba



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
    if number_of_distinct_colors < (number_of_shades-1):
        number_of_shades = number_of_distinct_colors - 2

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

  
def clusters_colormap_hexcolor(clusters, rwb_scale):
    import matplotlib as mpl
    norm = mpl.colors.Normalize(vmin=np.min(clusters), vmax=np.max(clusters))

    if rwb_scale:
        cmap = mpl.pyplot.get_cmap('coolwarm_r')
        rgba = cmap(norm(clusters))
    else:
        cmap = generate_colormap(int(np.max(clusters)))
        rgba = cmap(norm(clusters))
    hexcolor = []
    for c in rgba:
      hexcolor.append(mpl.colors.rgb2hex(c))

    return hexcolor, cmap
      
def cluster_scheme(mol, hexcolor, clusters):
  r0 = mol[0]['resid']
  c0 = hexcolor[0]
  clust_scheme = []
  select = '@'
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
  

def view_pdb_ngl(pdb, capsid, labels, rwb_scale=False):
    import biotite.structure.io as strucio
    strucio.save_structure(pdb + '_capsid.pdb', capsid, hybrid36=True)

    mol = open_pdb(pdb)
    hexcolor, cmap = clusters_colormap_hexcolor(labels, rwb_scale)
    clust_scheme = cluster_scheme(mol, hexcolor, labels)

    import nglview as ngl
    color_scheme = ngl.color._ColorScheme(clust_scheme, label="scheme_regions")
    view = ngl.show_structure_file(pdb + '_capsid.pdb', gui=False, default_representation=False)
    view._remote_call("setSize", target='Widget', args=['800px', '800px'])
    view.add_representation("spacefill", color=color_scheme)


    if rwb_scale:
        print('Each atom in this structure is colored according to the clustering quality score of its residue.')
        import matplotlib.colorbar as colorbar
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(10,0.5))
        cb = colorbar.ColorbarBase(ax, orientation='horizontal',
                                   cmap=cmap, norm=plt.Normalize(np.min(labels), np.max(labels)))
        plt.show()
        code = """
                    var $text = $("<div></div>")
                                .css("position", "absolute")
                                .css("top", "3%")
                                .css("left", "30%")
                                .css("padding", "2px 5px 2px 5px")
                                .css("opacity", "1.0")
                                .css("font-size", "30px")
                                .css("color", "black")
                                .appendTo(this.$container);

                    $text.text("{0} residue cluster scores")
                    """
        view._execute_js_code(code.format(pdb))
    else:
        print('Each atom in this structure has the same color as other atoms in the same cluster.')
        code = """
            var $text = $("<div></div>")
                        .css("position", "absolute")
                        .css("top", "3%")
                        .css("left", "40%")
                        .css("padding", "2px 5px 2px 5px")
                        .css("opacity", "1.0")
                        .css("font-size", "30px")
                        .css("color", "black")
                        .appendTo(this.$container);

            $text.text("{0} ({1} clusters)")
            """
        n_clusters = str(int(np.max(labels) + 1))
        view._execute_js_code(code.format(pdb, n_clusters))


    view.center()
    print('rendering')
    view.render_image()

    return view

def createCapsidView(pdb, capsid=None):

    if capsid is None:
        print('No capsid structure provided, getting capsid.')
        from pyCapsid.PDB import getCapsid
        capsid, _, _, _, _, _, _ = getCapsid(pdb)

    import biotite.structure.io as strucio
    strucio.save_structure(pdb + '_capsid.pdb', capsid, hybrid36=True)

    import nglview as ngl
    view = ngl.NGLWidget()
    from nglview.adaptor import FileStructure

    view.add_component(FileStructure(pdb + '_capsid.pdb'))
    #view.clear_representations()

    return view

def showColorBar(labels):
    print('Each atom in this structure is colored according to the clustering quality score of its residue.')
    import matplotlib.colorbar as colorbar
    import matplotlib.pyplot as plt
    hexcolor, cmap = clusters_colormap_hexcolor(labels, rwb_scale=True)
    fig, ax = plt.subplots(figsize=(10, 0.5))
    cb = colorbar.ColorbarBase(ax, orientation='horizontal',
                               cmap=cmap, norm=plt.Normalize(np.min(labels), np.max(labels)))
    plt.show()

def createClusterRepresentation(pdb, labels, view, rwb_scale=False, rep_type='cartoon'):

    mol = open_pdb(pdb)
    hexcolor, cmap = clusters_colormap_hexcolor(labels, rwb_scale)
    clust_scheme = cluster_scheme(mol, hexcolor, labels)

    view.clear_representations()
    import nglview as ngl
    color_scheme = ngl.color._ColorScheme(clust_scheme, label="scheme_regions")
    view._remote_call("setSize", target='Widget', args=['800px', '800px'])
    view.add_representation(rep_type, color=color_scheme)

    if rwb_scale:
        print('Each atom in this structure is colored according to the clustering quality score of its residue.')
        import matplotlib.colorbar as colorbar
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(10, 0.5))
        cb = colorbar.ColorbarBase(ax, orientation='horizontal',
                                   cmap=cmap, norm=plt.Normalize(np.min(labels), np.max(labels)))
        plt.show()
        code = """
                        var $text = $("<div></div>")
                                    .css("position", "absolute")
                                    .css("top", "3%")
                                    .css("left", "30%")
                                    .css("padding", "2px 5px 2px 5px")
                                    .css("opacity", "1.0")
                                    .css("font-size", "30px")
                                    .css("color", "black")
                                    .appendTo(this.$container);

                        $text.text("{0} residue cluster scores")
                        """
        view._execute_js_code(code.format(pdb))
    else:
        print('Each atom in this structure has the same color as other atoms in the same cluster.')
        code = """
                var $text = $("<div></div>")
                            .css("position", "absolute")
                            .css("top", "3%")
                            .css("left", "40%")
                            .css("padding", "2px 5px 2px 5px")
                            .css("opacity", "1.0")
                            .css("font-size", "30px")
                            .css("color", "black")
                            .appendTo(this.$container);

                $text.text("{0} ({1} clusters)")
                """
        n_clusters = str(int(np.max(labels) + 1))
        view._execute_js_code(code.format(pdb, n_clusters))

    view.center()
    print('rendering')
    #view.render_image()