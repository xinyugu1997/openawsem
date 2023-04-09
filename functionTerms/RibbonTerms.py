from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import numpy as np


def ribbon_term(oa, ribbon_group=False, group=[], ribbon_seq_sep=3, k_ribbon=4.184*100, Monomer=True, CB_flat=False, olig_num=1, forceGroup=26):
    import itertools
    k_ribbon *= oa.k_awsem
    # create contact force
    ribbon = CustomExternalForce(f"{k_ribbon}*(z-z_layer)^2")
    ribbon.addPerParticleParameter("z_layer")
    structure_interactions_ribbon = []
    #add interactions;
    if Monomer:
         if ribbon_group:
             Nres=int(len(oa.ca)/len(oa.chain_starts))
             print(Nres)
             for i in range(len(oa.chain_starts)):
                z_layer = i*5.0/10.0*nanometers
                print(i, z_layer)                
                for j in group:
                    structure_interactions_ribbon.append([oa.ca[j+i*Nres], [z_layer]])
                    if CB_flat and oa.cb[j+i*Nres] != -1:
                         structure_interactions_ribbon.append([oa.cb[j+i*Nres], [z_layer]])
                    #structure_interactions_ribbon.append([oa.n[j+i*Nres], [z_layer]])
                    #structure_interactions_ribbon.append([oa.c[j+i*Nres], [z_layer]])
         else:
             for i in range(len(oa.chain_starts)):
                z_layer = i*5.0/10.0*nanometers
                print(i, z_layer) 
                for j in range(oa.chain_starts[i]+1, oa.chain_ends[i]): 
                    if j%ribbon_seq_sep == 0:  
                       structure_interactions_ribbon.append([oa.ca[j], [z_layer]])
                       if CB_flat and oa.cb[j] != -1:
                          structure_interactions_ribbon.append([oa.cb[j], [z_layer]])
                       #structure_interactions_ribbon.append([oa.n[j], [z_layer]])
                       #structure_interactions_ribbon.append([oa.c[j], [z_layer]])
    else:
         Nmono = len(oa.chain_starts)/olig_num
         for i in range(len(oa.chain_starts)):
            z_layer = (i%Nmono)*5.0/10.0*nanometers
            print(i, z_layer)
            for j in range(oa.chain_starts[i]+1, oa.chain_ends[i]):
                if j%ribbon_seq_sep == 0:
                   structure_interactions_ribbon.append([oa.ca[j], [z_layer]])
                   #structure_interactions_ribbon.append([oa.n[j], [z_layer]])
                   #structure_interactions_ribbon.append([oa.c[j], [z_layer]])


    # create bonds
    for structure_interaction_ribbon in structure_interactions_ribbon:
        ribbon.addParticle(*structure_interaction_ribbon)
    ribbon.setForceGroup(forceGroup)
    print("Ribbon term is ON"+"\t fix "+str(len(structure_interactions_ribbon))+" atoms")
    return ribbon

