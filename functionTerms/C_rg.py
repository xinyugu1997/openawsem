from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *

def cylindrical_rg_term(oa, atomGroup=-1, forceGroup=2):
    nres, ca = oa.nres, oa.ca
    if atomGroup == -1:
        group = list(range(nres))
    else:
        group = atomGroup          # atomGroup = [0, 1, 10, 12]  means include residue 1, 2, 11, 13.
    n = len(group)
    normalization = n * n
    rg_square = CustomCompoundBondForce(2, f"1/{normalization}*((x1-x2)^2+(y1-y2)^2)")

    for i in group:
        for j in group:
            if j <= i:
                continue
            rg_square.addBond([ca[i], ca[j]], [])

    rg = CustomCVForce(f"rg_square^0.5")
    rg.addCollectiveVariable("rg_square", rg_square)
    rg.setForceGroup(forceGroup)
    return rg
