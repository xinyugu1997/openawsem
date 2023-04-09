from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import numpy as np

def NB_debye_huckel_term(self, k_dh=4.15*4.184, forceGroup=30, screening_length=1.0, periodic=False, chargeFile=None):
        print("Nonbonded Debye Huckel term is ON")
        k_dh *= self.k_awsem*0.1
        k_screening = 1.0
        screening_length = 1.0  # (in the unit of nanometers)

        dh = CustomNonbondedForce(f"{k_dh}*charge1*charge2/r*exp(-{k_screening}*r/{screening_length})")
        dh.addPerParticleParameter("charge")
        charge_index = []
        if chargeFile is None:
             for i, atom in enumerate(self.pdb.topology.atoms()):
                   charge = 0.0
                   if atom.name == 'CB':
                       if self.seq[atom.residue.index] == "R" or self.seq[atom.residue.index]=="K":
                               charge = 1.0
                               charge_index.append(i)
                       if self.seq[atom.residue.index] == "D" or self.seq[atom.residue.index]=="E":
                              charge = -1.0
                              charge_index.append(i)
                   dh.addParticle([charge])
        else:
             print("Using charge file\n")
             chargeInfo = np.loadtxt(chargeFile)
             for i, atom in enumerate(self.pdb.topology.atoms()):
                   charge = 0.0
                   if atom.name == 'CB':
                       if chargeInfo[atom.residue.index] != 0:
                               charge = chargeInfo[atom.residue.index]
                               charge_index.append(i)
                   dh.addParticle([charge])

        dh.addInteractionGroup([x for x in charge_index], [x for x in charge_index])

        if periodic:
            dh.setCutoffDistance(1.2)
            dh.setNonbondedMethod(dh.CutoffPeriodic)
        else:
            dh.setNonbondedMethod(dh.NoCutoff)
        print("Contact cutoff ", dh.getCutoffDistance())
        print("NonbondedMethod: ", dh.getNonbondedMethod())

        dh.createExclusionsFromBonds(self.bonds, 1)
        dh.setForceGroup(forceGroup)
        return dh


