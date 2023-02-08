"""Module with functions for downloading and dealing with PDB/PDBx files."""


def getCapsid(pdb, dir='.', pdbx=False, local=False, save=False, chains='', chains_clust=''):
    """Downloads and opens molecular data from a PDB entry or loads data from a local file.

    :param pdb: PDB id of entry to download. Can also be the name of a local file
    :param dir: Target directory where download will be placed
    :param pdbx: Whether the target structure should be acquired in pdbx/mmcif format
    :param local: Whether to instead load a local file
    :param save: Whether to save a copy of the complete assembly as pdb/pdbx. Necessary if visualizing in external software.
    :param chains: List of chains from the entry to include in the ENM model
    :param chains_clust: List of chains that will be assigned to quasi-rigid clusters. Must be a subset of 'chains'
    """
    if local:
        filename = dir + pdb
    else:
        filename = downloadPDB(pdb, dir, pdbx)

    if pdbx:
        capsid, calphas, coords, bfactors, chain_starts, title = loadPDBx(filename, pdb, save)
    else:
        capsid, calphas, coords, bfactors, chain_starts, title = loadPDB(filename, pdb, save)

    return capsid, calphas, coords, bfactors, chain_starts, title


def downloadPDB(pdb, dir='.', pdbx=False):
    """Downloads pdb and returns the filename

    :param pdb: PDB id of entry to download. Can also be the name of a local file
    :param dir: Target directory where download will be placed
    :param pdbx: Whether the target structure should be acquired in pdbx/mmcif format
    """
    from biotite.database.rcsb import fetch
    if pdbx:
        filename = fetch(pdb, target_path=dir, format='pdbx', overwrite=True, verbose=True)
    else:
        filename = fetch(pdb, target_path=dir, format='pdb', overwrite=True, verbose=True)

    return filename


def loadPDBx(filename, pdb, save):
    """Loads PDBx data from a file

    :param filename: Name of local file
    :param pdb: PDB id of entry
    :param save: Whether to save a copy of the complete assembly as pdb/pdbx.
        """
    import biotite.structure as struc
    import biotite.structure.io.pdbx as pdbx
    import biotite.structure.io as strucio
    import os

    pdbx_file = pdbx.PDBxFile()
    pdbx_file.read(filename)
    capsid = pdbx.get_assembly(pdbx_file, assembly_id="1", model=1, extra_fields=['b_factor'])

    title = pdb  # pdbx_file.get_category('pdbx_database_related')['details']

    print("Number of protein chains:", struc.get_chain_count(capsid))
    # capsid = capsid[capsid.filter_amino_acids()].copy()

    calphas = capsid[capsid.atom_name == 'CA']
    coords = calphas.coord

    if save:
        print('Saving complete capsid PDBx')
        strucio.save_structure(pdb + '_capsid.pdbx', capsid)

    chain_starts = struc.get_chain_starts(calphas)

    return capsid, calphas, coords, calphas.b_factor, chain_starts, title


def loadPDB(filename, pdb, save):
    """Loads PDBx data from a file

    :param filename: Name of local file
    :param pdb: PDB id of entry
    :param save: Whether to save a copy of the complete assembly as pdb/pdbx
    """
    from prody import parsePDB, writePDB
    import os

    capsid, header = parsePDB(filename, header=True, biomol=True, secondary=True, extend_biomol=True)
    asym_unit = parsePDB(filename)  # Maybe need a check to see if assembly

    if type(capsid) is list:
        print('list')
        capsid = capsid[0]

    ENM_capsid = capsid.select('protein').copy()
    calphas = ENM_capsid.select('calpha')

    ENM_capsid_asym = asym_unit.select('protein').copy()
    calphas_asym = ENM_capsid_asym.select('calpha').copy()

    chain_starts = getProdyChainStarts(calphas_asym)

    print('Number Of Residues: ', calphas.getCoords().shape[0])

    if save:
        print('Writing complete capsid PDB')
        #print('Number Of Saved Residues: ', ENM_capsid.numResidues())
        writePDB(pdb + '_capsid.pdb', ENM_capsid, hybrid36=True)

    if 'title' in header:
        title = header['title']
    else:
        title = pdb

    return ENM_capsid, calphas, calphas.getCoords(), calphas.getBetas(), chain_starts, title


def getProdyChainStarts(calphas_asym, n_units=60):
    """

    :param calphas_asym:
    :param n_units:
    :return:
    """
    import numpy as np
    n_asym = calphas_asym.numAtoms()
    n_chains = calphas_asym.numChains()
    chains = calphas_asym.getChids()
    chaindiff = np.where(chains[:-1] != chains[1:])[0]
    chain_starts = []
    for i in range(n_units):
        start = n_asym * i
        chain_starts.append(start)
        for c in chaindiff:
            chain_starts.append(start + c)

    chain_starts.append(n_asym * n_units)

    return np.array(chain_starts)

# def buildMassesCoords(atoms):
#     print(atoms[0])
#     bfs = []
#     masses = []
#     prot_ids = []
#
#
#     for pn, chain in enumerate(atoms.iterChains()):
#         print('Chain Id ' + pn + ':' + chain.getChid())
#         prot_ids.append(pn)
#         for res in chain.iterResidues():
#             mass = np.sum(res.getMasses())
#             masses.append(mass)
#
#             bfactor = np.mean(res.getBetas())
#             bfs.append(bfactor)
#
#     return np.asarray(bfs), np.asarray(masses), np.asarray(pn)
