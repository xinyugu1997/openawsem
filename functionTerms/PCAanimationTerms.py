from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import numpy as np


def PCA_fix_term(oa, index, in_rnativeCACA, k_er=4.184*200000, forceGroup=25):
    ### this is a structure prediction related term; Adapted from Sirovitz Schafer Wolynes 2017 Protein Science;
    ### See original papers for reference: Make AWSEM AWSEM-ER with Evolutionary restrictions
    ### ER restrictions can be obtained from multiple sources (RaptorX, deepcontact, and Gremlin)
    ### term modified from amh-go term, and the current strength seems to be high, and needs to be lowered somehow.
    ### amh-go normalization factor will be added soon. Based on Eastwood Wolynes 2000 JCP
    import itertools
    k_er *= oa.k_awsem
    # create contact force
    er = CustomBondForce(f"0.5*{k_er}*(r-r_ijN)^2")
    # # add global parameters
    er.addPerBondParameter("r_ijN")
    structure_interactions_er = []
    ### read in dat files from contact predictions;
    for i in range(len(index)):
                r_ijN = in_rnativeCACA[i]/10.0*nanometers;
                structure_interactions_er.append([oa.ca[index[i][0]], oa.ca[index[i][1]], [r_ijN]])
    # create bonds
    for structure_interaction_er in structure_interactions_er:
        er.addBond(*structure_interaction_er)
    er.setForceGroup(forceGroup)
    print("PCA fix term is ON")
    return er

def PCA_cg_fix_term(oa, groups, in_rnativeCACA, k_er=4.184*200000, forceGroup=25):
    ### this is a structure prediction related term; Adapted from Sirovitz Schafer Wolynes 2017 Protein Science;
    ### See original papers for reference: Make AWSEM AWSEM-ER with Evolutionary restrictions
    ### ER restrictions can be obtained from multiple sources (RaptorX, deepcontact, and Gremlin)
    ### term modified from amh-go term, and the current strength seems to be high, and needs to be lowered somehow.
    ### amh-go normalization factor will be added soon. Based on Eastwood Wolynes 2000 JCP
    print("PCA cg fix term is ON")
    import itertools
    from scipy.spatial import distance as sdist

    k_er *= oa.k_awsem
    # create contact force
    er = CustomCentroidBondForce(2, f"0.5*{k_er}*(distance(g1,g2)-r_ijN)^2")
    # # add global parameters
    er.addPerBondParameter("r_ijN")
    for g in range(len(groups)):
        er.addGroup([oa.ca[i] for i in groups[g]])            

    structure_interactions_er = []
    ### read in dat files from contact predictions;
    matrix=sdist.squareform(in_rnativeCACA)
    for i,j in itertools.combinations(range(len(groups)),2):
                r_ijN = matrix[i][j]/10.0*nanometers
                structure_interactions_er.append([[i,j], [r_ijN]])
    # create bonds
    for structure_interaction_er in structure_interactions_er:
        er.addBond(*structure_interaction_er)
    er.setForceGroup(forceGroup)
    return er

def pcvalue_cartesin(oa, pc_vector, mean_coords, forceGroup=3):
    pcvalue = CustomExternalForce("(x-x0)*vx+(y-y0)*vy+(z-z0)*vz")
    pcvalue.addPerParticleParameter("x0")
    pcvalue.addPerParticleParameter("y0")
    pcvalue.addPerParticleParameter("z0")
    pcvalue.addPerParticleParameter("vx")
    pcvalue.addPerParticleParameter("vy")
    pcvalue.addPerParticleParameter("vz")
    # add particles
    for i in range(int(len(pc_vector)/3)):
        pcvalue.addParticle(oa.ca[i], [mean_coords[3*i], mean_coords[3*i+1], mean_coords[3*i+2], pc_vector[3*i], pc_vector[3*i+1], pc_vector[3*i+2]])
    pcvalue.setForceGroup(forceGroup)
    return pcvalue


def PCA_cartesin_fix_term(oa, pc0, pc_vector, mean_coords, k_pc=4.184*2000, forceGroup=25):
    pcbias = CustomCVForce(f"0.5*{k_pc}*(pc-{pc0})^2")
    pc = pcvalue_cartesin(oa, pc_vector, mean_coords)
    pcbias.addCollectiveVariable("pc", pc)
    print("PC bias term ON, pc0,k_pc = ", pc0, k_pc)
    pcbias.setForceGroup(forceGroup)
    return pcbias

def cartesin_nail_term(oa, k_nail, nail_id, nail_coord):
    nail = CustomExternalForce(f"0.5*{k_nail}*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    nail.addPerParticleParameter("x0")
    nail.addPerParticleParameter("y0")
    nail.addPerParticleParameter("z0")
    # add particles
    for i in range(len(nail_id)):
        nail.addParticle(oa.ca[nail_id[i]], [nail_coord[3*i], nail_coord[3*i+1], nail_coord[3*i+2]])
    print("cartesin nail term ON")
    return nail
