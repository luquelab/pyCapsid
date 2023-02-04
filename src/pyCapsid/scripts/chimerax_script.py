from chimerax.core.commands import *
from chimerax.std_commands import *
from chimerax.core.commands import all_objects
from chimerax.core.commands import run
from chimerax.atomic.molsurf import MolecularSurface
from chimerax.atomic import *
import numpy as np
import matplotlib as mpl


# from chimerax.core.commands import measure


from chimerax.std_commands import split
from numpy.linalg import norm
import sys
print(sys.argv)
label_dir, save_dir, pdb_dir = (sys.argv[1], sys.argv[2], sys.argv[3])

run(session, 'open ' + '\'' + pdb_dir + '\'')
# center of mass measurement from ChimeraX tutorial
atoms = all_objects(session).atoms  # getting atom list
print(atoms)
coords = atoms.scene_coords  # getting atom coords
modelCenter = coords.mean(axis=0)  # calculate center of mass
print("modelCenter", modelCenter)
if modelCenter.any:
    print('Aligning Models')
    run(session, 'hkcage 1 0')
    modelX, modelY, modelZ = modelCenter
    run(session, 'move x ' + str(-1 * modelX) + ' coordinateSystem #1 models #1')  # adjust x coordinates
    run(session, 'move y ' + str(-1 * modelY) + ' coordinateSystem #1 models #1')  # adjust y coordinates
    run(session, 'move z ' + str(-1 * modelZ) + ' coordinateSystem #1 models #1')  # adjust z coordinates
    run(session, 'close #3')
else:
    print('Models Already Aligned')

radius = norm(coords, axis=1).max()  # calculate the norms of the x coordinates, and then choose the maximum value
print("radius_est: ", radius)

labels = np.load(label_dir)

import matplotlib as mpl
norm = mpl.colors.Normalize(vmin=np.min(labels), vmax=np.max(labels))
cmap = generate_colormap(int(np.max(labels)))
rgba = cmap(norm(labels))*255

residues = all_objects(session).residues
cx_nr = len(residues)
enm_nr = len(labels)

print('# of residues:', cx_nr)

for i in range(cx_nr):
    res = residues[i]
    res.ribbon_color = rgba[i % rgba.shape[0], :]
    ats = res.atoms
    for at in ats:
        at.color = rgba[i % rgba.shape[0], :]

run(session, 'hkcage 1 0 alpha hexagonal-dual radius ' + str(radius) + ' spherefactor 0.2')
run(session, 'view orient')
run(session, 'save ' + 'test.png')
