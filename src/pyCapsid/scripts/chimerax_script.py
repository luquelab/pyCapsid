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



run(session, 'hkcage 1 0 alpha hexagonal-dual radius 100 spherefactor 0.2')
run(session, 'view orient')
run(session, 'save ' + 'test.png')