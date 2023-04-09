from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import numpy as np


def group_constraint_by_distance(oa, d0=0*angstrom, group1=[], group2=[], forceGroup=3, k=1*kilocalorie_per_mole/(nanometer**2)):
    # CustomCentroidBondForce only work with CUDA not OpenCL.
    # only CA, CB, O has mass. so the group have to include those.
    k = k.value_in_unit(kilojoule_per_mole/(nanometer**2))   # convert to kilojoule_per_mole, openMM default uses kilojoule_per_mole as energy.
    k_constraint = k * oa.k_awsem
    d0 = d0.value_in_unit(nanometer)   # convert to nm
    constraint = CustomCentroidBondForce(2, f"0.5*{k_constraint}*(distance(g1,g2)-{d0})^2")
    #constraint = CustomCentroidBondForce(2, f"0.5*step(distance(g1,g2)-{d0})*{k_constraint}*(distance(g1,g2)-{d0})^2")
    # example group set up group1=[oa.ca[7], oa.cb[7]] use the ca and cb of residue 8.
    g1 =[oa.ca[x-1] for x in group1]
    g2 =[oa.ca[x-1] for x in group2]
    constraint.addGroup(g1)    # group use particle index.
    constraint.addGroup(g2)
    constraint.addBond([0, 1])
    print(constraint.getGroupParameters(0))
    print(constraint.getGroupParameters(1))
    constraint.setForceGroup(forceGroup)
    print("Attract bias on: r0, k_attract = ", d0, k_constraint)
    return constraint

def group_distance(oa, group1=[], group2=[], forceGroup=2):
    dist = CustomCentroidBondForce(2, "distance(g1,g2)")
    # example group set up group1=[oa.ca[7], oa.cb[7]] use the ca and cb of residue 8.
    g1 =[oa.ca[x-1] for x in group1]
    g2 =[oa.ca[x-1] for x in group2]
    dist.addGroup(g1)    # group use particle index.
    dist.addGroup(g2)
    dist.addBond([0, 1])
    print(dist.getGroupParameters(0))
    print(dist.getGroupParameters(1))
    dist.setForceGroup(forceGroup)
    return dist

def group_constraint_by_distance_gaussian(oa, d0=0*angstrom, sigma=10*angstrom, group1=[], group2=[], forceGroup=3, k=1*kilocalorie_per_mole):
    # CustomCentroidBondForce only work with CUDA not OpenCL.
    # only CA, CB, O has mass. so the group have to include those.
    k = k.value_in_unit(kilojoule_per_mole)   # convert to kilojoule_per_mole, openMM default uses kilojoule_per_mole as energy.
    k_constraint = k * oa.k_awsem
    d0 = d0.value_in_unit(nanometer)   # convert to nm
    sigma = sigma.value_in_unit(nanometer)   # convert to nm
    constraint = CustomCentroidBondForce(2, f"-{k_constraint}*exp(-(distance(g1,g2)-{d0})^2/(2*{sigma}^2))")
    g1 =[oa.ca[x-1] for x in group1]
    g2 =[oa.ca[x-1] for x in group2]
    constraint.addGroup(g1)    # group use particle index.
    constraint.addGroup(g2)
    constraint.addBond([0, 1])
    print(constraint.getGroupParameters(0))
    print(constraint.getGroupParameters(1))
    constraint.setForceGroup(forceGroup)
    print("Gaussian Attract bias on: r0, sigma, k_attract = ", d0, sigma, k_constraint)
    return constraint
