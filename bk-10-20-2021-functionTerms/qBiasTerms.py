from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import numpy as np
from Bio.PDB.PDBParser import PDBParser

def read_reference_structure_for_q_calculation_3(oa, pdb_file, inter_domain=False, link=0, min_seq_sep=3, max_seq_sep=np.inf, contact_threshold=0.95*nanometers, Qflag=0):
    # default use all chains in pdb file.
    # this change use the canonical Qw/Qo calculation for reference Q
    # for Qw calculation is 0; Qo is 1;
    structure_interactions = []
    parser = PDBParser()
    structure = parser.get_structure('X', pdb_file)
    model = structure[0]
    chain_start = 0
    count = 0
    for chain in model.get_chains():
        chain_start += count
        count = 0
        for i, residue_i in enumerate(chain.get_residues()):
            #  print(i, residue_i)
            count +=1
            for j, residue_j in enumerate(chain.get_residues()):
                    if abs(i-j) >= min_seq_sep and abs(i-j) <= max_seq_sep:  # taking the signed value to avoid double counting
                        if inter_domain and (i >= link or j < link): continue 
                        ca_i = residue_i['CA']

                        ca_j = residue_j['CA']

                        r_ijN = abs(ca_i - ca_j)/10.0   #*nanometers # convert to nm
                        if Qflag ==1 and r_ijN >= contact_threshold: continue
                        sigma_ij = 0.1*abs(i-j)**0.15 # 0.1 nm = 1 A
                        gamma_ij = 1.0
                        i_index = oa.ca[i+chain_start]
                        j_index = oa.ca[j+chain_start]
                        structure_interaction = [i_index, j_index, [gamma_ij, r_ijN, sigma_ij]]
                        structure_interactions.append(structure_interaction)
    print("read ref Qflag=", Qflag)
    return structure_interactions

def q_value(oa, reference_pdb_file, inter_domain=False, link=0, Qflag=0, min_seq_sep=3, max_seq_sep=np.inf, contact_threshold=0.95*nanometers):
    ### Modified by Mingchen to compute canonical QW/QO
    # create bond force for q calculation
    qvalue = CustomBondForce("(1/normalization)*gamma_ij*exp(-(r-r_ijN)^2/(2*sigma_ij^2))")
    qvalue.addPerBondParameter("gamma_ij")
    qvalue.addPerBondParameter("r_ijN")
    qvalue.addPerBondParameter("sigma_ij")
    # create bonds
    structure_interactions = read_reference_structure_for_q_calculation_3(oa, reference_pdb_file, inter_domain=inter_domain, link=link,
        min_seq_sep=min_seq_sep, max_seq_sep=max_seq_sep, contact_threshold=contact_threshold, Qflag=Qflag)
    #print(len(structure_interactions))
    #print(structure_interactions)
    qvalue.addGlobalParameter("normalization", len(structure_interactions))
    for structure_interaction in structure_interactions:
        qvalue.addBond(*structure_interaction)
    qvalue.setForceGroup(1)
    return qvalue


def qbias_term(oa, q0, reference_pdb_file, inter_domain=False, link=0, k_qbias=10000, Qflag=0, qbias_min_seq_sep=3, qbias_max_seq_sep=np.inf, qbias_contact_threshold=0.95*nanometers):
    qbias = CustomCVForce("0.5*k_qbias*(q-q0)^2")
    q = q_value(oa, reference_pdb_file, inter_domain=inter_domain, link=link, Qflag=Qflag, min_seq_sep=qbias_min_seq_sep, max_seq_sep=qbias_max_seq_sep, contact_threshold=qbias_contact_threshold)
    qbias.addCollectiveVariable("q", q)
    qbias.addGlobalParameter("k_qbias", k_qbias)
    qbias.addGlobalParameter("q0", q0)
    print("Qbias term ON, q0, Qflag= ", q0, Qflag)
    return qbias


def qc_sum(oa, structure_interactions, contact_threshold):
    ### Modified by Xinyu to compute Qc
    # create bond force for q calculation
    qsum = CustomBondForce(f"gamma_ij*exp(-(r-r_ijN)^2/(2*sigma_ij^2))*(1-(step(r-{contact_threshold}))*(step(r_ijN-{contact_threshold})))")
    qsum.addPerBondParameter("gamma_ij")
    qsum.addPerBondParameter("r_ijN")
    qsum.addPerBondParameter("sigma_ij")
    for structure_interaction in structure_interactions:
        qsum.addBond(*structure_interaction)
    return qsum

def qc_normalization(oa, structure_interactions, contact_threshold):
    normalization = CustomBondForce(f"1-step(r-{contact_threshold})*step(r_ijN-{contact_threshold})")
    normalization.addPerBondParameter("r_ijN")
    for structure_interaction in structure_interactions:
        tmp=[structure_interaction[0], structure_interaction[1], [structure_interaction[2][1]]]
        normalization.addBond(*tmp)
    print(structure_interaction, tmp)
    return normalization

def qc_value(oa, reference_pdb_file, min_seq_sep=3, max_seq_sep=np.inf, contact_threshold=0.95):
    qcvalue=CustomCVForce(f"qsum/normalization")
    structure_interactions = read_reference_structure_for_q_calculation_3(oa, reference_pdb_file,
              min_seq_sep=min_seq_sep, max_seq_sep=max_seq_sep, contact_threshold=contact_threshold, Qflag=0)
    qsum=qc_sum(oa, structure_interactions, contact_threshold=contact_threshold)
    normalization = qc_normalization(oa, structure_interactions, contact_threshold=contact_threshold)
    qcvalue.addCollectiveVariable("qsum", qsum)
    qcvalue.addCollectiveVariable("normalization", normalization)
    qcvalue.setForceGroup(1)
    return qcvalue

def qcbias_term(oa, q0, reference_pdb_file, k_qbias=10000, qbias_contact_threshold=0.95, qbias_min_seq_sep=3, qbias_max_seq_sep=np.inf):
    qbias = CustomCVForce(f"0.5*{k_qbias}*(q-{q0})^2")
    q = qc_value(oa, reference_pdb_file, min_seq_sep=qbias_min_seq_sep, max_seq_sep=qbias_max_seq_sep, contact_threshold=qbias_contact_threshold)
    qbias.addCollectiveVariable("q", q)
    print("Qcbias term ON")
    return qbias

