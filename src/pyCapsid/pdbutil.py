def getCapsid(pdb, dir='.', pdbx=False, local=False, chains='', chains_clust=''):
    if local:
        filename = dir + pdb
    else:
        filename = downloadPDB(pdb, dir, pdbx)

    if pdbx:
        capsid, calphas, coords, bfactors, title = loadPDBx(filename, pdb)
    else:
        capsid, calphas, coords, bfactors, title = loadPDB(filename, pdb)

    return capsid, calphas, coords, bfactors, title


def downloadPDB(pdb, dir='.', pdbx=False):
    from biotite.database.rcsb import fetch
    if pdbx:
        filename = fetch(pdb, target_path=dir, format='pdbx', overwrite=False, verbose=True)
    else:
        filename = fetch(pdb, target_path=dir, format='pdb', overwrite=False, verbose=True)

    return filename


def loadPDBx(filename, pdb):
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

    if not os.path.exists(pdb + '_ca.pdbx'):
        print('Writing calphas PDBx')
        strucio.save_structure(pdb + '_ca.pdbx', calphas)
    if not os.path.exists(pdb + '_capsid.pdbx'):
        print('Writing complete capsid PDBx')
        strucio.save_structure(pdb + '_capsid.pdbx', capsid)

    chains = struc.get_chain_starts(calphas)

    return capsid, calphas, coords, calphas.b_factor, chains, title


def loadPDB(filename, pdb):
    from prody import parsePDB, writePDB
    import os

    capsid, header = parsePDB(filename, header=True, biomol=True, secondary=True, extend_biomol=True)

    if type(capsid) is list:
        capsid = capsid[0]
    ENM_capsid = capsid.select('protein').copy()

    calphas = ENM_capsid.select('calpha')
    print('Number Of Residues: ', calphas.getCoords().shape[0])

    if not os.path.exists(pdb + '_ca.pdb'):
        print('Writing calphas PDB')
        writePDB(pdb + '_ca.pdb', calphas, hybrid36=True)
    if not os.path.exists(pdb + '_capsid.pdb'):
        print('Writing complete capsid PDB')
        writePDB(pdb + '_capsid.pdb', ENM_capsid, hybrid36=True)

    if 'title' in header:
        title = header['title']
    else:
        title = pdb

    return ENM_capsid, calphas, calphas.getCoords(), calphas.getBetas(), title


def buildMassesCoords(atoms):
    print(atoms[0])
    bfs = []
    masses = []
    prot_ids = []

    # print('segments:', atoms.numSegments())

    for pn, chain in enumerate(atoms.iterChains()):
        print('Chain Id ' + pn + ':' + chain.getChid())
        prot_ids.append(pn)
        for res in chain.iterResidues():
            mass = np.sum(res.getMasses())
            masses.append(mass)

            bfactor = np.mean(res.getBetas())
            bfs.append(bfactor)

            # coord = res['CA'].getCoords()
            # coords.append(coord)

    return np.asarray(bfs), np.asarray(masses), np.asarray(pn)
