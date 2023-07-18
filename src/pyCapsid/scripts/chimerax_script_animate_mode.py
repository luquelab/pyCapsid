import numpy as np
def generateTrajFromModes(atoms, modes, evals, n_mode=0, max_amplitude=None, n_frames = 20):
    residues = atoms.unique_residues
    coords = atoms.scene_coords

    mode_vec = modes[:, n_mode].reshape((-1, 3))
    if max_amplitude is not None:
        vector_lengths = np.sqrt(np.sum(mode_vec**2, axis=-1))
        scale = max_amplitude / np.max(vector_lengths)
        mode_vec *= scale
    else:
        mode_vec *= 1/evals[n_mode]


    time = np.linspace(0, 2*np.pi, n_frames, endpoint=False)
    deviation = np.sin(time)[:, np.newaxis, np.newaxis] * mode_vec
    print(deviation.shape)
    oscillation = np.zeros((n_frames, coords.shape[0], 3))

    ind = 0
    for i in range(len(residues)):
        res = residues[i]
        n_atoms = res.num_atoms
        res_start = ind
        res_stop = ind + n_atoms - 1
        ind += n_atoms
        dev = deviation[:,i:i+1,:]
        oscillation[:,res_start:res_stop,:] = coords[res_start:res_stop,:] + dev
    return oscillation

from chimerax.core.commands import *
from chimerax.std_commands import *
from chimerax.core.commands import all_objects
from chimerax.core.commands import run
from chimerax.atomic.molsurf import MolecularSurface
from chimerax.atomic import *
import matplotlib as mpl
# Colormap taken from stackexchange

# from chimerax.core.commands import measure


from numpy.linalg import norm
import sys
print(sys.argv)

import argparse
import os

parser = argparse.ArgumentParser(description='ChimeraX script for visualization of pyCapsid results')
parser.add_argument('-report_dir', help='Location of pyCapsid report chimerax directory', default=None, required=False)
parser.add_argument('-nmode', help='Index of normal mode to visualize', default=None, required=False)
parser.add_argument('-remote', help='Whether to use a remote structure from the PDB database', default=None, required=False)
parser.add_argument('-pdb', help='If remote is True or none, PDBID of the target structure. Otherwise, the local filename of the target structure', default=None, required=False)
parser.add_argument('-amplitude', help='Maximum amplitude of motion along the mode for visualization purposes', default=None, required=False)
parser.add_argument('-frames', help='Number of frames to generate for the animation', default=20, required=False)
parser.add_argument('-save_frames', help='Whether to save the individual frames used to create the image', default=False, required=False)
args = vars(parser.parse_args())

if args['report_dir'] is None:
    report_dir = os.path.dirname(sys.argv[0])
    print(report_dir)
else:
    report_dir = args['report_dir']
print(report_dir)
os.chdir(report_dir)

if args['pdb'] is None:
    fs = [x for x in os.listdir('./') if x.endswith('_final_results.npz')]
    f = fs[0]
    pdb = f.split('_final')[0]
else:
    pdb = args['pdb']

if args['remote'] is None:
    files = [x for x in os.listdir('./') if (x.endswith('.pdb') or x.endswith('.cif'))]
    if len(files) == 0:
        remote = 'True'
    elif len(files) == 1:
        remote = 'False'
        pdb = files[0]
    else:
        raise Exception("More than one .pdb file in directory")
else:
    remote = args['remote']

if args['amplitude'] is None:
    max_amplitude = None # change to finding first non-degenerate mode
else:
    max_amplitude = float(args['amplitude'])

if args['nmode'] is None:
    modes = np.load(f'{pdb}_modes.npz')
    evals = modes['eigen_vals']
    uniques, inds, counts = np.unique(evals.round(decimals=8), return_index=True, return_counts=True)
    icoEvalInds = inds[counts == 1]
    print('Indices of non-degenerate eigenmodes: ', icoEvalInds)
    if len(icoEvalInds) == 0:
        print('No non-degenerate modes found')
        n_mode = 0
    else:
        print('Visualizing lowest frequency non-degenerate mode')
        n_mode=icoEvalInds[0]
else:
    n_mode = int(args['nmode'])

if args['frames'] is None:
    n_frames = 20
else:
    n_frames = int(args['frames'])

save_frames = args['save_frames']

run(session, 'set bg white')
run(session, 'graphics silhouettes true')

if remote=='True':
    run(session, f'open {pdb}')
    run(session, 'sym #1 assembly 1 copies True')
    run(session, 'close #1')
    run(session, 'view orient')
    run(session, 'hide atoms')
    run(session, 'show cartoons')
    run(session, 'lighting soft')
    run(session, f'save ../figures/structures/{pdb}_full_capsid.png')
    run(session, 'rename #2 id #1')
else:
    run(session, f'open {pdb}')
    run(session, 'view orient')
    run(session, 'lighting soft')
    run(session, 'hide atoms')
    run(session, 'show cartoons')
    run(session, f'save ../figures/structures/{pdb}_full_capsid.png')

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
    run(session, 'move x ' + str(-1 * modelX) + ' coordinateSystem #1 models #1')  # adjust x coordinates
    run(session, 'move y ' + str(-1 * modelY) + ' coordinateSystem #1 models #1')  # adjust y coordinates
    run(session, 'move z ' + str(-1 * modelZ) + ' coordinateSystem #1 models #1')  # adjust z coordinates
    run(session, 'view orient')
else:
    print('Models Already Aligned')

coords = atoms.scene_coords

modes = np.load(f'{pdb}_modes.npz')
evecs = modes['eigen_vecs']
evals = modes['eigen_vals']

atoms = all_objects(session).atoms
osc = generateTrajFromModes(atoms, evecs, evals, n_mode=n_mode, max_amplitude=max_amplitude, n_frames=n_frames)

print(f'Visualizing motion of {pdb} structure along mode # {n_mode}')
run(session, 'view orient pad -0.2')
run(session, 'color bychain')
run(session, 'movie record supersample 3 directory ./ format png')

from time import sleep
for i in range(n_frames):
    atoms.scene_coords = osc[i, :, :]
    run(session, 'wait 1')
    #run(session, f'save ../figures/structures/{pdb}_mode_{n_mode}_animation_frame_{i}.png')
if save_frames:
    run(session, f'movie encode quality high resetMode keep output ../figures/structures/{pdb}_mode_{n_mode}_animation.mp4')
else:
    run(session, f'movie encode quality high output ../figures/structures/{pdb}_mode_{n_mode}_animation.mp4')


