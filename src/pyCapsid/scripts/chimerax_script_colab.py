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
    if number_of_distinct_colors < (number_of_shades - 1):
        number_of_shades = number_of_distinct_colors - 2

    number_of_distinct_colors_with_multiply_of_shades = int(
        math.ceil(number_of_distinct_colors / number_of_shades) * number_of_shades)

    # Create an array with uniformly drawn floats taken from <0, 1) partition
    linearly_distributed_nums = np.arange(
        number_of_distinct_colors_with_multiply_of_shades) / number_of_distinct_colors_with_multiply_of_shades

    # We are going to reorganise monotonically growing numbers in such way that there will be single array with saw-like pattern
    #     but each saw tooth is slightly higher than the one before
    # First divide linearly_distributed_nums into number_of_shades sub-arrays containing linearly distributed numbers
    arr_by_shade_rows = linearly_distributed_nums.reshape(number_of_shades,
                                                          number_of_distinct_colors_with_multiply_of_shades // number_of_shades)

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
        initial_cm[0:lower_half, i] *= np.arange(0.2, 1, 0.8 / lower_half)

    # Modify second half in such way that colours towards end of partition are less intense and brighter
    # Colours closer to the middle are affected less, colours closer to the end are affected more
    for i in range(3):
        for j in range(upper_partitions_half):
            modifier = np.ones(number_of_shades) - initial_cm[lower_half + j * number_of_shades: lower_half + (
                        j + 1) * number_of_shades, i]
            modifier = j * modifier / upper_partitions_half
            initial_cm[lower_half + j * number_of_shades: lower_half + (j + 1) * number_of_shades, i] += modifier

    return ListedColormap(initial_cm)

def getRGBA(values, rwb_scale):
    import matplotlib as mpl
    import numpy as np
    norm = mpl.colors.Normalize(vmin=np.min(values), vmax=np.max(values))
    if rwb_scale == 'True':
        import matplotlib.pyplot as plt
        cmap = plt.get_cmap('coolwarm_r')
        rgba = cmap(norm(values)) * 255
    else:
        cmap = generate_colormap(int(np.max(values)))
        rgba = cmap(norm(values)) * 255
    return rgba

from chimerax.core.commands import *
from chimerax.std_commands import *
from chimerax.core.commands import all_objects
from chimerax.core.commands import run
from chimerax.atomic.molsurf import MolecularSurface
from chimerax.atomic import *
import numpy as np
import matplotlib as mpl
# Colormap taken from stackexchange

# from chimerax.core.commands import measure


from numpy.linalg import norm
import sys
print(sys.argv)
if len(sys.argv) == 2:
    report_dir = sys.argv[1]
    import os
    os.chdir(report_dir)
    # defaults
    rwb_scale = 'False'
    fs = [x for x in os.listdir('./') if x.endswith('.npz')]
    f = fs[0]
    pdb = f.split('_final')[0]
    n_clusters = np.load(f'{pdb}_final_results.npz')['nc']

    # params = np.loadtxt('chimerax_params.txt')

    files = [x for x in os.listdir('./') if (x.endswith('.pdb') or x.endswith('.cif'))]
    if len(files) == 0:
        remote = 'True'
    elif len(files) == 1:
        remote = 'False'
        pdb = files[0]
    else:
        raise Exception("More than one .pdb file in directory")
else:
    pdb, report_dir, remote, rwb_scale = (sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    #...

run(session, 'set bg white')
run(session, 'graphics silhouettes true')

if remote=='True':
    run(session, f'open {pdb}')
    run(session, 'view orient')
    run(session, 'lighting full')
    run(session, f'save ../figures/structures/{pdb}_asymmetric_unit.png')
    run(session, 'sym #1 assembly 1 copies True')
    run(session, 'close #1')
    run(session, 'view orient')
    run(session, f'save ../figures/structures/{pdb}_full_capsid.png')
else:
    run(session, f'open {pdb}')

# center of mass measurement from ChimeraX tutorial
run(session, 'sel protein')
run(session, 'del ~sel')
run(session, 'sel clear')
atoms = all_objects(session).atoms  # getting atom list
coords = atoms.scene_coords  # getting atom coords
modelCenter = coords.mean(axis=0)  # calculate center of mass
print("modelCenter", modelCenter)

if modelCenter.any:
    print('Aligning Models')
    modelX, modelY, modelZ = modelCenter
    run(session, 'move x ' + str(-1 * modelX) + ' coordinateSystem #2 models #2')  # adjust x coordinates
    run(session, 'move y ' + str(-1 * modelY) + ' coordinateSystem #2 models #2')  # adjust y coordinates
    run(session, 'move z ' + str(-1 * modelZ) + ' coordinateSystem #2 models #2')  # adjust z coordinates
else:
    print('Models Already Aligned')

coords = atoms.scene_coords
radius = norm(coords, axis=1).max()  # calculate the norms of the x coordinates, and then choose the maximum value
print("radius_est: ", radius)

results = np.load(f'{pdb}_final_results_full.npz')
labels = results['labels']
nc_range = results['nc_range']
scores = results['full_scores']

if n_clusters is None:
    print('Defaulting to highest score clustering')
    ind = np.argmax(scores)
    n_c = nc_range[ind]
else:
    ind = np.argwhere(nc_range == n_clusters)[0][0]
    n_c = n_clusters

print(f'Visualizing cluster results of {pdb} for {n_c} clusters')
labels = labels[ind]
score = scores[ind]

rgba_scores = getRGBA(score, rwb_scale='True')
rgba_clusters = getRGBA(labels, rwb_scale='False')
print(rgba_scores)
residues = atoms.unique_residues
cx_nr = len(residues)
enm_nr = rgba_scores.shape[0]

print('# of residues:', cx_nr)
print('# of ENM residues:', enm_nr)

for i in range(cx_nr):
    res = residues[i]
    res.ribbon_color = rgba_scores[i, :]
    ats = res.atoms
    for at in ats:
        at.color = rgba_scores[i, :]

# run(session, 'hkcage 1 0 alpha hexagonal-dual radius ' + str(radius) + ' spherefactor 0.2')
run(session, f'save ../figures/structures/{pdb}_residue_cluster_scores.png')

print('# of residues:', cx_nr)
print('# of ENM residues:', enm_nr)

for i in range(cx_nr):
    res = residues[i]
    res.ribbon_color = rgba_clusters[i, :]
    ats = res.atoms
    for at in ats:
        at.color = rgba_clusters[i, :]

# run(session, 'hkcage 1 0 alpha hexagonal-dual radius ' + str(radius) + ' spherefactor 0.2')
run(session, f'save ../figures/structures/{pdb}_highest_quality_clusters.png')

