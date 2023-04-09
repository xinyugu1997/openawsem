from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import numpy as np
from Bio.PDB.PDBParser import PDBParser

def read_reference_structure_for_q_calculation_PBC(oa, pdb_file, min_seq_sep=3, max_seq_sep=np.inf, contact_threshold=0.95, Qflag=0):
    # default use all chains in pdb file.
    # this change use the canonical Qw/Qo calculation for reference Q
    # for Qw calculation is 0; Qo is 1;


    structure_interactions = []
    parser = PDBParser()
    structure = parser.get_structure('X', pdb_file)
    chains=structure[0]
    residues = []
    for chain in chains:
     for x in chain:
        residues.append(x)
    Lchain = len(chain)
    print("Monomer length:", Lchain, "res")
    print("Q value for ", len(residues), "res")
    for i, residue_i in enumerate(residues):
         for j, residue_j in enumerate(residues):
                    if (j-i) >= min_seq_sep and abs(i-j) <= max_seq_sep:  # taking the signed value to avoid double counting
                        ca_i = residue_i['CA']

                        ca_j = residue_j['CA']

                        r_ijN = abs(ca_i - ca_j)/10.0   #*nanometers # convert to nm
                        if Qflag ==1 and r_ijN >= contact_threshold: continue
                        if Qflag ==2 and (r_ijN >= contact_threshold or i//Lchain == j//Lchain): continue 
                        sigma_ij = 0.1*abs(i-j)**0.15 # 0.1 nm = 1 A
                        gamma_ij = 1.0
                        i_index = oa.ca[i]
                        j_index = oa.ca[j]
                        structure_interaction = [i_index, j_index, [gamma_ij, r_ijN, sigma_ij]]
                        structure_interactions.append(structure_interaction)
    print("read ref Qflag=", Qflag)
    print("normalization=", len(structure_interactions))
    return structure_interactions

def q_value_PBC(oa, reference_pdb_file, Qflag=0, min_seq_sep=3, max_seq_sep=np.inf, contact_threshold=0.95, PBC=True):
    # create bond force for q calculation
    qvalue = CustomBondForce("(1/normalization)*gamma_ij*exp(-(r-r_ijN)^2/(2*sigma_ij^2))")
    qvalue.addPerBondParameter("gamma_ij")
    qvalue.addPerBondParameter("r_ijN")
    qvalue.addPerBondParameter("sigma_ij")
    # create bonds
    structure_interactions = read_reference_structure_for_q_calculation_PBC(oa, reference_pdb_file, 
        min_seq_sep=min_seq_sep, max_seq_sep=max_seq_sep, contact_threshold=contact_threshold, Qflag=Qflag)
    #print(len(structure_interactions))
    #print(structure_interactions)
    qvalue.addGlobalParameter("normalization", len(structure_interactions))
    for structure_interaction in structure_interactions:
        qvalue.addBond(*structure_interaction)
    qvalue.setForceGroup(1)
    if PBC == True:
    	qvalue.setUsesPeriodicBoundaryConditions(True)
    return qvalue



def read_reference_structure_for_qh_calculation_PBC(oa, pdb_file, gamma_file=None, include_intraC=True, min_seq_sep=3, max_seq_sep=np.inf, contact_threshold=0.95, Qflag=0):
    structure_interactions = []
    parser = PDBParser()
    structure = parser.get_structure('X', pdb_file)
    chains = structure[0]
    Nchains = len(chains)
    chain_id = ["A", "B", "C", "D", "E", "F"]

    Mres = len(chains[chain_id[0]])
    rN = np.zeros((Nchains, Mres, Mres))
    print("Q value for ", Nchains, "chains, ", Mres, "resi per chain")

    for i in range(Mres):
        residue_i = chains[chain_id[0]][i+1]
        ca_i = residue_i['CA']
        for x in range(Nchains):
            for j in range(Mres):
               residue_j = chains[chain_id[x]][j+1]
               ca_j = residue_j['CA']
               rN[x][i][j] = abs(ca_i - ca_j)/10.0   #*nanometers # convert to nm

    np.save("./qh_rN.npy", rN)

    if gamma_file:
       gamma = np.loadtxt(gamma_file)
    else:
       gamma = np.ones(Nchains)


    Mchains = len(oa.chain_ends)
    #intra-chain pairs
    if include_intraC:
    	print('Include intra chain residue pairs')
    	for x in range(Mchains):
    	    for i in range(Mres):
    	       for j in range(i+min_seq_sep, Mres):
    	            sigma_ij = 0.1*abs(i-j)**0.15
    	            i_index = oa.ca[x*Mres+i]
    	            j_index = oa.ca[x*Mres+j]
    	            structure_interaction = [i_index, j_index, [gamma[0],0,0,0,0,0, rN[0][i][j],-1,-1,-1,-1,-1, sigma_ij]]
    	            structure_interactions.append(structure_interaction)
    	            print(structure_interaction)

    #inter-chain pairs
    for x in range(Mchains):
        for y in range(x+1, Mchains):
              for i in range(Mres):
                 for j in range(Mres):
                      i_index = oa.ca[x*Mres+i]
                      j_index = oa.ca[y*Mres+j]
                      structure_interaction = [i_index, j_index, [0, gamma[1], gamma[2], gamma[3], gamma[4], gamma[5], -1, rN[1][i][j], rN[2][i][j], rN[3][i][j], rN[4][i][j], rN[5][i][j], -1]]
                      structure_interactions.append(structure_interaction) 
                      print(structure_interaction)


    print("read ref Qflag=", Qflag)
    print("normalization=", len(structure_interactions))
    return structure_interactions

def qh_value_PBC(oa, reference_pdb_file, gamma_file=None, include_intraC=True, Qflag=0, min_seq_sep=3, max_seq_sep=np.inf, contact_threshold=0.95, sigma_sq=0.05):
    # create bond force for q calculation
    qvalue = CustomBondForce(f"(1/normalization)*(gamma0_ij*exp(-(r-r0_ijN)^2/(2*sigma_ij^2)) + gamma1_ij*exp(-(r-r1_ijN)^2/(2*{sigma_sq})) + gamma2_ij*exp(-(r-r2_ijN)^2/(2*{sigma_sq})) + gamma3_ij*exp(-(r-r3_ijN)^2/(2*{sigma_sq})) +gamma4_ij*exp(-(r-r4_ijN)^2/(2*{sigma_sq})) + gamma5_ij*exp(-(r-r5_ijN)^2/(2*{sigma_sq})))")
    qvalue.addPerBondParameter("gamma0_ij")
    qvalue.addPerBondParameter("gamma1_ij")
    qvalue.addPerBondParameter("gamma2_ij")
    qvalue.addPerBondParameter("gamma3_ij")
    qvalue.addPerBondParameter("gamma4_ij")
    qvalue.addPerBondParameter("gamma5_ij")
    qvalue.addPerBondParameter("r0_ijN")
    qvalue.addPerBondParameter("r1_ijN")
    qvalue.addPerBondParameter("r2_ijN")
    qvalue.addPerBondParameter("r3_ijN")
    qvalue.addPerBondParameter("r4_ijN")
    qvalue.addPerBondParameter("r5_ijN")
    qvalue.addPerBondParameter("sigma_ij")
    # create bonds
    structure_interactions = read_reference_structure_for_qh_calculation_PBC(oa, reference_pdb_file, gamma_file=gamma_file, include_intraC=include_intraC, 
        min_seq_sep=min_seq_sep, max_seq_sep=max_seq_sep, contact_threshold=contact_threshold, Qflag=Qflag)
    #print(len(structure_interactions))
    #print(structure_interactions)
    qvalue.addGlobalParameter("normalization", len(structure_interactions))
    for structure_interaction in structure_interactions:
        qvalue.addBond(*structure_interaction)
    qvalue.setForceGroup(1)
    qvalue.setUsesPeriodicBoundaryConditions(True)
    return qvalue


def qbias_PBC_term(oa, q0, reference_pdb_file, q_hexa=False, gamma_file=None, include_intraC=True, k_qbias=10000, Qflag=0, qbias_min_seq_sep=3, qbias_max_seq_sep=np.inf, qbias_contact_threshold=0.95, sigma_sq=0.05):
    qbias = CustomCVForce(f"0.5*{k_qbias}*(q-{q0})^2")
    if q_hexa:
       q = qh_value_PBC(oa, reference_pdb_file, gamma_file=gamma_file, include_intraC=include_intraC, Qflag=Qflag, min_seq_sep=qbias_min_seq_sep, max_seq_sep=qbias_max_seq_sep, contact_threshold=qbias_contact_threshold, sigma_sq=sigma_sq) 
    else:
       q = q_value_PBC(oa, reference_pdb_file, Qflag=Qflag, min_seq_sep=qbias_min_seq_sep, max_seq_sep=qbias_max_seq_sep, contact_threshold=qbias_contact_threshold)
    qbias.addCollectiveVariable("q", q)
    print("Qbias_PBC term ON, q0, Qflag= ", q0, Qflag)
    print("PBC: ", qbias.usesPeriodicBoundaryConditions())
    return qbias
