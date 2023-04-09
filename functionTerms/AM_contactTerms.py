from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import numpy as np
from simtk.unit.quantity import Quantity
import pandas as pd
from Bio.PDB.Polypeptide import three_to_one


gamma_se_map_1_letter = {   'A': 0,  'R': 1,  'N': 2,  'D': 3,  'C': 4,
                            'Q': 5,  'E': 6,  'G': 7,  'H': 8,  'I': 9,
                            'L': 10, 'K': 11, 'M': 12, 'F': 13, 'P': 14,
                            'S': 15, 'T': 16, 'W': 17, 'Y': 18, 'V': 19}


def read_gamma(gammaFile):
    data = np.loadtxt(gammaFile)
    gamma_direct = data[:210]
    gamma_mediated = data[210:]
    return gamma_direct, gamma_mediated


def inWhichChain(residueId, chain_ends):
    chain_table = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'
    for i, end_of_chain_resId in enumerate(chain_ends):
        if end_of_chain_resId < residueId:
            pass
        else:
            return chain_table[i]


def AM_contact_term(oa, k_contact=4.184, periodic=False, parametersLocation=".", burialPartOn=True, withExclusion=True, forceGroup=22,
                gammaName="gamma.dat", burialGammaName="burial_gamma.dat", membraneGammaName="membrane_gamma.dat", r_min=0.45):
    if isinstance(k_contact, float) or isinstance(k_contact, int):
        k_contact = k_contact * oa.k_awsem   # just for backward comptable
    elif isinstance(k_contact, Quantity):
        k_contact = k_contact.value_in_unit(kilojoule_per_mole)   # convert to kilojoule_per_mole, openMM default uses kilojoule_per_mole as energy.
        k_contact = k_contact * oa.k_awsem
    else:
        print(f"Unknown input, {k_contact}, {type(k_contact)}")
    # combine direct, burial, mediated.

    # r_min = .45
    r_max = .65
    r_minII = .65
    r_maxII = .95
    r_maxIII = .82
    eta = 50  # eta actually has unit of nm^-1.
    eta_sigma = 7.0
    rho_0 = 2.6
    min_sequence_separation = 3  # means j-i > 9
    min_sequence_separation_mem = 3
    nwell = 2
    eta_switching = 10
    gamma_ijm = np.zeros((nwell, 20, 20))
    water_gamma_ijm = np.zeros((nwell, 20, 20))
    protein_gamma_ijm = np.zeros((nwell, 20, 20))

    # read in seq data.
    seq = oa.seq
    # read in gamma info
    gamma_direct, gamma_mediated = read_gamma(os.path.join(parametersLocation, gammaName))

    burial_kappa = 4.0
    burial_ro_min = [0.0, 3.0, 6.0]
    burial_ro_max = [3.0, 6.0, 9.0]
    burial_gamma = np.loadtxt(os.path.join(parametersLocation, burialGammaName))

    contact = CustomGBForce()

    m = 0  # water environment
    count = 0
    for i in range(20):
        for j in range(i, 20):
            gamma_ijm[m][i][j] = gamma_direct[count][0]
            gamma_ijm[m][j][i] = gamma_direct[count][0]
            count += 1
    count = 0
    for i in range(20):
        for j in range(i, 20):
            water_gamma_ijm[m][i][j] = gamma_mediated[count][1]
            water_gamma_ijm[m][j][i] = gamma_mediated[count][1]
            count += 1
    count = 0
    for i in range(20):
        for j in range(i, 20):
            protein_gamma_ijm[m][i][j] = gamma_mediated[count][0]
            protein_gamma_ijm[m][j][i] = gamma_mediated[count][0]
            count += 1
    # residue interaction table (step(abs(resId1-resId2)-min_sequence_separation))
    res_table = np.zeros((nwell, oa.nres, oa.nres))
    for i in range(oa.nres):
        for j in range(oa.nres):
            resId1 = i
            chain1 = inWhichChain(resId1, oa.chain_ends)
            resId2 = j
            chain2 = inWhichChain(resId2, oa.chain_ends)
            if abs(resId1-resId2)-min_sequence_separation >= 0 or chain1 != chain2:
                res_table[0][i][j] = 1
            else:
                res_table[0][i][j] = 0

            if abs(resId1-resId2) >= 1 or chain1 != chain2:
                res_table[1][i][j] = 1
            else:
                res_table[1][i][j] = 0







    print("# of residue pairs:", np.sum(res_table))



    contact.addTabulatedFunction("gamma_ijm", Discrete3DFunction(nwell, 20, 20, gamma_ijm.T.flatten()))
    contact.addTabulatedFunction("water_gamma_ijm", Discrete3DFunction(nwell, 20, 20, water_gamma_ijm.T.flatten()))
    contact.addTabulatedFunction("protein_gamma_ijm", Discrete3DFunction(nwell, 20, 20, protein_gamma_ijm.T.flatten()))
    contact.addTabulatedFunction("burial_gamma_ij", Discrete2DFunction(20, 3, burial_gamma.T.flatten()))
    contact.addTabulatedFunction("res_table", Discrete3DFunction(nwell, oa.nres, oa.nres, res_table.T.flatten()))

    contact.addPerParticleParameter("resName")
    contact.addPerParticleParameter("resId")
    contact.addPerParticleParameter("isCb")

    contact.addComputedValue("rho", f"isCb1*isCb2*step(abs(resId1-resId2)-2)*0.25*(1+tanh({eta}*(r-{r_min})))*(1+tanh({eta}*({r_max}-r)))", CustomGBForce.ParticlePair)

    # replace cb with ca for GLY
    cb_fixed = [x if x > 0 else y for x,y in zip(oa.cb,oa.ca)]
    none_cb_fixed = [i for i in range(oa.natoms) if i not in cb_fixed]
    # print(oa.natoms, len(oa.resi), oa.resi, seq)
    for i in range(oa.natoms):
        contact.addParticle([gamma_se_map_1_letter[seq[oa.resi[i]]], oa.resi[i], int(i in cb_fixed)])


    intra = 0
    inter = 1
    inMembrane = 0
    contact.addEnergyTerm(f"-isCb1*isCb2*{k_contact}*(res_table({intra}, resId1, resId2)*gamma_ijm({inMembrane}, resName1, resName2)*theta+\
                            (thetaII*res_table({intra}, resId1, resId2)+2*thetaIII*res_table({inter}, resId1, resId2))*(sigma_water*water_gamma_ijm({inMembrane}, resName1, resName2)+\
                            sigma_protein*protein_gamma_ijm({inMembrane}, resName1, resName2)));\
                            sigma_protein=1-sigma_water;\
                            theta=0.25*(1+tanh({eta}*(r-{r_min})))*(1+tanh({eta}*({r_max}-r)));\
                            thetaII=0.25*(1+tanh({eta}*(r-{r_minII})))*(1+tanh({eta}*({r_maxII}-r)));\
                            thetaIII=0.25*(1+tanh({eta}*(r-{r_min})))*(1+tanh({eta}*({r_maxIII}-r)));\
                            sigma_water=0.25*(1-tanh({eta_sigma}*(rho1+2-{rho_0})))*(1-tanh({eta_sigma}*(rho2+2-{rho_0})))",
                            CustomGBForce.ParticlePair)


    if burialPartOn:
        # burial term
        for i in range(3):
            contact.addGlobalParameter(f"rho_min_{i}", burial_ro_min[i])
            contact.addGlobalParameter(f"rho_max_{i}", burial_ro_max[i])
        for i in range(3):
            contact.addEnergyTerm(f"-0.5*isCb*{k_contact}*burial_gamma_ij(resName, {i})*\
                                        (tanh({burial_kappa}*(rho+2-rho_min_{i}))+\
                                        tanh({burial_kappa}*(rho_max_{i}-rho-2)))", CustomGBForce.SingleParticle)

    print("Number of atom: ", oa.natoms, "Number of residue: ", len(cb_fixed))
    print("Inter-chain mediated term ON, r_maxIII:", r_maxIII)

    # withExclusion won't affect the result. But may speed up the calculation with CPU but slows down for GPU.
    if withExclusion:
        for e1 in none_cb_fixed:
            for e2 in none_cb_fixed:
                if e1 > e2:
                    continue
                contact.addExclusion(e1, e2)
        for e1 in none_cb_fixed:
            for e2 in cb_fixed:
                contact.addExclusion(e1, e2)

    # contact.setCutoffDistance(1.1)
    if periodic:
        contact.setNonbondedMethod(contact.CutoffPeriodic)
    else:
        contact.setNonbondedMethod(contact.CutoffNonPeriodic)
    print("Contact cutoff ", contact.getCutoffDistance())
    print("NonbondedMethod: ", contact.getNonbondedMethod())
    contact.setForceGroup(forceGroup)
    return contact

