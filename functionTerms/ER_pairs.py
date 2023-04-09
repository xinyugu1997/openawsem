from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import numpy as np

def er_pairs_term(oa, er_pairs_file, k_er=4.184, er_well_width=0.1, forceGroup=25):
    import itertools
    k_er *= oa.k_awsem
    # create contact force
    er = CustomBondForce(f"-{k_er}*gamma_ij*exp(-(r-r_ijN)^2/(2*sigma_ij^2))")
    # # add global parameters
    er.addPerBondParameter("gamma_ij")
    er.addPerBondParameter("r_ijN")
    er.addPerBondParameter("sigma_ij")
    structure_interactions_er = []
    ### read in dat files from contact predictions;
    in_rnativeCACA = np.loadtxt(er_pairs_file)
    for p in range(len(in_rnativeCACA)):
                i = int(in_rnativeCACA[p][0])
                j = int(in_rnativeCACA[p][1])
                sigma_ij = er_well_width*abs(i-j)**0.15 # 0.1 nm = 1 A
                gamma_ij = 1.0
                r_ijN = in_rnativeCACA[p][2]/10.0*nanometers;
                structure_interactions_er.append([oa.ca[i], oa.ca[j], [gamma_ij, r_ijN, sigma_ij]])
    # create bonds
    for structure_interaction_er in structure_interactions_er:
        er.addBond(*structure_interaction_er)
        print(structure_interaction_er)
    er.setForceGroup(forceGroup)
    print("ER pair term is ON, k_er=", k_er)
    return er

def er_pairs_Ribbon_term(oa, er_pairs_file, k_er=0.01*4.184, rN=1.2, er_well_width=0.1, forceGroup=25):
    #rN=1.2*nanometers, er_well_width=0.1*nanometers
    import itertools
    k_er *= oa.k_awsem
    # create contact force
    er = CustomBondForce(f"-{k_er}*gamma_ij*exp(-max((r-{rN}),0)^2/(2*sigma_ij^2))")
    er.addPerBondParameter("gamma_ij")
    er.addPerBondParameter("sigma_ij")

    nres = len(oa.ca)
    nchain = len(oa.chain_ends)
    Mres = int(nres/nchain)
    print(Mres)

    structure_interactions_er = []
    ### read in dat files from contact predictions;
    in_rnativeCACA = np.loadtxt(er_pairs_file)
    for p in range(len(in_rnativeCACA)):
        for q in range(nchain):
                i = int(in_rnativeCACA[p][0]) + q*Mres
                j = int(in_rnativeCACA[p][1]) + q*Mres
                sigma_ij = er_well_width*abs(i-j)**0.15 # 0.1 nm = 1 A
                gamma_ij = in_rnativeCACA[p][2]
                structure_interactions_er.append([oa.ca[i], oa.ca[j], [gamma_ij, sigma_ij]])
    # create bonds
    for structure_interaction_er in structure_interactions_er:
        er.addBond(*structure_interaction_er)
        print(structure_interaction_er)
    er.setForceGroup(forceGroup)
    print("ER pair Ribbon term is ON, k_er=", k_er)
    return er
