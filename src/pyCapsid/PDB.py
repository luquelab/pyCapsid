"""Module with functions for downloading and dealing with PDB/PDBx files."""


def getCapsid(pdb, save_pdb_path='./', pdbx=False, local=False, save_full_pdb=False, chains=None, chains_clust=None, is_ico=False):
    """Downloads and opens molecular data from a PDB entry or loads data from a local file.

    :param pdb: PDB id of entry to download. Can also be the name of a local file
    :param save_pdb_path: Target directory where download will be placed
    :param pdbx: Whether the target structure should be acquired in pdbx/mmcif format
    :param local: Whether to instead load a local file
    :param save_full_pdb: Whether to save a copy of the complete assembly as pdb/pdbx. Useful for some visualizations.
    :param chains: List of chains from the entry to include in the ENM model. Not implemented
    :param chains_clust: List of chains that will be assigned to quasi-rigid clusters. Must be a subset of 'chains'. Not implemented
    """


    if isinstance(pdbx, str):
        pdbx = pdbx == 'true' or pdbx == 'True'
    if isinstance(local, str):
        local = local == 'true' or local == 'True'
    if isinstance(save_full_pdb, str):
        save_full_pdb = save_full_pdb == 'true' or save_full_pdb == 'True'

    if local:
        filename = save_pdb_path + pdb
    else:
        filename = downloadPDB(pdb, save_pdb_path, pdbx)

    if pdbx:
        capsid, calphas, asym, coords, bfactors, chain_starts, title = loadPDBx(filename, pdb, save_full_pdb)
    else:
        try:
            capsid, calphas, asym, coords, bfactors, chain_starts, title = loadPDB(filename, pdb, save_full_pdb)
        except:
            print('no .pdb file found. Checking for pdbx/mmcif')
            print('Add pdbx=True if you know you will be fetching a pdbx/mmcif file')
            capsid, calphas, asym, coords, bfactors, chain_starts, title = loadPDBx(filename, pdb, save_full_pdb)

    n_res = len(calphas)
    print(f'# of residues: {n_res}')

    return capsid, calphas, asym, coords, bfactors, chain_starts, title


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
        try:
            filename = fetch(pdb, target_path=dir, format='pdb', overwrite=True, verbose=True)
        except:
            print('no .pdb file found for this PDB id. Checking for pdbx/mmcif')
            print('Add pdbx=True if you know you will be fetching a pdbx/mmcif file')
            filename = fetch(pdb, target_path=dir, format='pdbx', overwrite=True, verbose=True)
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
    import numpy as np

    pdbx_file = pdbx.PDBxFile()
    pdbx_file.read(filename)
    capsid = pdbx.get_assembly(pdbx_file, assembly_id="1", model=1, extra_fields=['b_factor'])
    asym = pdbx.get_structure(pdbx_file, model=1, extra_fields=['b_factor'])

    capsid = capsid[struc.filter_amino_acids(capsid)]
    title = pdb  # pdbx_file.get_category('pdbx_database_related')['details']

    # capsid = capsid[capsid.filter_amino_acids()].copy()

    capsid = capsid[struc.filter_amino_acids(capsid)]
    asym = asym[struc.filter_amino_acids(asym)]
    print("Number of protein chains in full structure:", struc.get_chain_count(capsid))
    print("Number of protein chains in asymmetric unit:", struc.get_chain_count(asym))

    calphas = capsid[capsid.atom_name == 'CA']
    coords = calphas.coord

    if save:
        print('Saving complete capsid PDBx')
        strucio.save_structure(pdb + '_capsid.pdbx', capsid)

    chain_starts = struc.get_chain_starts(calphas)
    chain_starts = np.append(chain_starts, len(calphas))

    return capsid, calphas, asym, coords, calphas.b_factor, chain_starts, title

def loadPDB(filename, pdb_id, save):
    """Loads PDBx data from a file

    :param filename: Name of local file
    :param pdb_id: PDB id of entry
    :param save: Whether to save a copy of the complete assembly as pdb/pdbx.
        """
    import biotite.structure as struc
    import biotite.structure.io.pdb as pdb
    import biotite.structure.io as strucio
    import numpy as np

    pdb_file = pdb.PDBFile()
    pdb_file.read(filename)
    capsid = pdb.get_assembly(pdb_file, assembly_id="1", model=1, extra_fields=['b_factor'])
    asym = pdb.get_structure(pdb_file, model=1, extra_fields=['b_factor'])

    title = pdb_id  # pdbx_file.get_category('pdbx_database_related')['details']

    # capsid = capsid[capsid.filter_amino_acids()].copy()

    capsid = capsid[struc.filter_amino_acids(capsid)]
    asym = asym[struc.filter_amino_acids(asym)]
    print("Number of protein chains in full structure:", struc.get_chain_count(capsid))
    print("Number of protein chains in asymmetric unit:", struc.get_chain_count(asym))


    calphas = capsid[capsid.atom_name == 'CA']
    coords = calphas.coord

    if save:
        print('Saving complete capsid PDB')
        strucio.save_structure(pdb_id + '_capsid.pdb', capsid, hybrid36=True)

    chain_starts = struc.get_chain_starts(calphas)
    chain_starts = np.append(chain_starts, len(calphas))

    return capsid, calphas, asym, coords, calphas.b_factor, chain_starts, title


# def loadPDB(filename, pdb, save):
#     """Loads PDBx data from a file
#
#     :param filename: Name of local file
#     :param pdb: PDB id of entry
#     :param save: Whether to save a copy of the complete assembly as pdb/pdbx
#     """
#     from prody import parsePDB, writePDB
#
#     capsid, header = parsePDB(filename, header=True, biomol=True, secondary=True, extend_biomol=True)
#     asym_unit = parsePDB(filename)  # Maybe need a check to see if assembly
#
#
#     ENM_capsid = capsid.select('protein').copy()
#     calphas = ENM_capsid.select('calpha')
#
#     ENM_capsid_asym = asym_unit.select('protein').copy()
#     calphas_asym = ENM_capsid_asym.select('calpha').copy()
#
#     chain_starts = getProdyChainStarts(calphas_asym)
#
#     print('Number Of Residues: ', calphas.getCoords().shape[0])
#
#     if save:
#         print('Writing complete capsid PDB')
#         #print('Number Of Saved Residues: ', ENM_capsid.numResidues())
#         writePDB(pdb + '_capsid.pdb', ENM_capsid, hybrid36=True)
#
#     if 'title' in header:
#         title = header['title']
#     else:
#         title = pdb
#
#     return ENM_capsid, calphas, calphas.getCoords(), calphas.getBetas(), chain_starts, title


# def getProdyChainStarts(calphas_asym, n_units=60):
#     """
#
#     :param calphas_asym:
#     :param n_units:
#     :return:
#     """
#     import numpy as np
#     n_asym = calphas_asym.numAtoms()
#     n_chains = calphas_asym.numChains()
#     chains = calphas_asym.getChids()
#     chaindiff = np.where(chains[:-1] != chains[1:])[0]
#     chain_starts = []
#     for i in range(n_units):
#         start = n_asym * i
#         chain_starts.append(start)
#         for c in chaindiff:
#             chain_starts.append(start + c)
#
#     chain_starts.append(n_asym * n_units)
#
#     return np.array(chain_starts)

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
