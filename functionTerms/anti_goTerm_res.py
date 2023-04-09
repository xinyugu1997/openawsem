from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import numpy as np
import pandas as pd

def lists_of_index_atoms(nres, atoms, index):
    atom_types=['n', 'h', 'ca', 'c', 'o', 'cb']
    atom_types_table = {'N':'n', 'H':'h', 'CA':'ca', 'C':'c', 'O':'o', 'CB':'cb'}
    index_atoms=[]
    atom_lists=dict(zip(atom_types,[[-1]*nres for i in range(len(atom_types))]))
    j = 0
    m = 0
    z = 0
    for atom in atoms:
        if m != atom.residue.index and z == 1:
           j = j+1
        if j >= len(index): break
        if atom.residue.index == index[j] -1 :
           atom_lists[atom_types_table[atom.name]][atom.residue.index] = atom.index
           index_atoms.append(atom.index)
           z = 1
        else: z = 0
        m = atom.residue.index
    return atom_lists, index_atoms

def lists_of_go_atoms(nres, atoms, go_chain_start):
    print("Go chain starts from " + str(go_chain_start))
    atom_types=['n', 'h', 'ca', 'c', 'o', 'cb']
    atom_types_table = {'N':'n', 'H':'h', 'CA':'ca', 'C':'c', 'O':'o', 'CB':'cb'}
    go_atoms=[]
    atom_lists=dict(zip(atom_types,[[-1]*nres for i in range(len(atom_types))]))
    for atom in atoms:
        if atom.residue.index >= go_chain_start:
           atom_lists[atom_types_table[atom.name]][atom.residue.index] = atom.index
           go_atoms.append(atom.index)
 
    return atom_lists, go_atoms


def additive_anti_go_term(oa, index, go_chain_id, k_antigo=8368, r_antigo=0.65): 
    # multiply interaction strength by overall scaling
    k_antigo *= oa.k_awsem
    # create contact force
    antigo = CustomNonbondedForce("k_antigo*step(r_0-r)*(r_0-r)^2")
    # # add global parameters
    antigo.addGlobalParameter("k_antigo", k_antigo)
    antigo.addGlobalParameter("r_0", r_antigo)
    for i in range(oa.natoms):
        antigo.addParticle()
    #create interaction groups
    oa.index_atom_lists, oa.index_atoms = lists_of_index_atoms(oa.nres, oa.pdb.topology.atoms(), index)
    oa.go_atom_lists, oa.go_atoms = lists_of_go_atoms(oa.nres, oa.pdb.topology.atoms(), go_chain_start = oa.chain_starts[go_chain_id]) 

#    oa.index_ca = oa.index_atom_lists['ca']
#    oa.index_cb = oa.index_atom_lists['cb'] 
#    oa.index_o = oa.index_atom_lists['o']
#    oa.go_ca = oa.go_atom_lists['ca']
#    oa.go_cb = oa.go_atom_lists['cb']
#    oa.go_o = oa.go_atom_lists['o']
#
    antigo.addInteractionGroup([x for x in oa.index_atoms if x >= 0], [x for x in oa.go_atoms if x >= 0])
    antigo.setCutoffDistance(r_antigo)
    antigo.setNonbondedMethod(antigo.CutoffNonPeriodic)
    antigo.createExclusionsFromBonds(oa.bonds, 1)
    antigo.setForceGroup(31)
    print("ANTI-GO term is ON")
    return antigo


def string_length_term(oa, inter, binder_chain_id):
    string_length = CustomCentroidBondForce(2, "distance(g1,g2)")
    oa.inter_atom_lists, oa.inter_atoms = lists_of_index_atoms(oa.nres, oa.pdb.topology.atoms(), inter)
    oa.binder_atom_lists, oa.binder_atoms = lists_of_go_atoms(oa.nres, oa.pdb.topology.atoms(), go_chain_start = oa.chain_starts[binder_chain_id])
    string_length.addGroup([x for x in oa.inter_atom_lists['ca'] if x >= 0])
    string_length.addGroup([x for x in oa.binder_atom_lists['ca'] if x >= 0])
    bondGroups=[]
    bondGroups.append(0)
    bondGroups.append(1)
    string_length.addBond(bondGroups)
    string_length.setForceGroup(2)
    return string_length


def string_term(oa, inter, binder_chain_id, r0, k_string=10*4.184):
    string = CustomCentroidBondForce(2, f"0.5*{k_string}*(distance(g1,g2)-{r0})^2")
#    string.addGlobalParameter("k_string",k_string)
#    string.addGlobalParameter("r0", r0)
    oa.inter_atom_lists, oa.inter_atoms = lists_of_index_atoms(oa.nres, oa.pdb.topology.atoms(), inter)
    oa.binder_atom_lists, oa.binder_atoms = lists_of_go_atoms(oa.nres, oa.pdb.topology.atoms(), go_chain_start = oa.chain_starts[binder_chain_id])
    string.addGroup([x for x in oa.inter_atom_lists['ca'] if x >= 0])
    string.addGroup([x for x in oa.binder_atom_lists['ca'] if x >= 0])
    bondGroups=[]
    bondGroups.append(0)
    bondGroups.append(1)

    print(string.getGroupParameters(0))
    print(string.getGroupParameters(1))

    string.addBond(bondGroups)
    print("String bias on: r0, k_string = ", r0, k_string)
    return string

def bond_two_res_term(oa, res1, res2, k_bond=0.415*4.184, k_bond1=0, r0=1.5):
    bond2res = CustomBondForce(f"-4*{k_bond}/r*exp(-r)+{k_bond1}*(r-{r0})^2")
    bond2res.addBond(oa.cb[res1], oa.cb[res2])
    print("Bond:", oa.cb[res1], oa.cb[res2], k_bond, k_bond1, r0)
    return bond2res


def attract_term(oa, group1, group2, r0, k_attract=0.0001*4.184):
    attract = CustomCentroidBondForce(2, f"step(distance(g1,g2)-{r0})*{k_attract}*(distance(g1,g2)-{r0})^2")
    oa.g1_atom_lists, oa.g1_atoms = lists_of_index_atoms(oa.nres, oa.pdb.topology.atoms(), group1)
    oa.g2_atom_lists, oa.g2_atoms = lists_of_index_atoms(oa.nres, oa.pdb.topology.atoms(), group2)
    attract.addGroup([x for x in oa.g1_atom_lists['ca'] if x >= 0])
    attract.addGroup([x for x in oa.g2_atom_lists['ca'] if x >= 0])
    bondGroups=[]
    bondGroups.append(0)
    bondGroups.append(1)
    attract.addBond(bondGroups)
    print("Attract bias on: r0, k_attract = ", r0, k_attract)
    return attract
