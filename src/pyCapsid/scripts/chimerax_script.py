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


from chimerax.std_commands import split
from numpy.linalg import norm
import sys
print(sys.argv)
color_dir, save_dir, pdb_dir, remote, pdb, rwb_scale = (sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])

if remote=='True':
    run(session, f'open {pdb}')
    run(session, 'sym #1 assembly 1 copies True')
    run(session, 'close #1')
else:
    run(session, 'open ' + '\'' + pdb_dir + '\'')

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

rgba = np.load(color_dir)
# import matplotlib as mpl
# norm = mpl.colors.Normalize(vmin=np.min(labels), vmax=np.max(labels))
#
# if rwb_scale=='True':
#     import matplotlib.pyplot as plt
#     cmap = plt.get_cmap('coolwarm_r')
#     rgba = cmap(norm(labels))*255
# else:
#     cmap = generate_colormap(int(np.max(labels)))
#     rgba = cmap(norm(labels))*255

residues = atoms.unique_residues
cx_nr = len(residues)
enm_nr = rgba.shape[0]

print('# of residues:', cx_nr)
print('# of ENM residues:', enm_nr)

for i in range(cx_nr):
    res = residues[i]
    res.ribbon_color = rgba[i, :]
    ats = res.atoms
    for at in ats:
        at.color = rgba[i, :]

# run(session, 'hkcage 1 0 alpha hexagonal-dual radius ' + str(radius) + ' spherefactor 0.2')
run(session, 'view orient')
run(session, 'lighting full')
run(session, f'save {pdb}_capsid.png')

