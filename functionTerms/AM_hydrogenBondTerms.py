from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import numpy as np

se_map_1_letter = {'A': 0,  'P': 1,  'K': 2,  'N': 3,  'R': 4,
                   'F': 5,  'D': 6,  'Q': 7,  'E': 8,  'G': 9,
                   'I': 10, 'H': 11, 'L': 12, 'C': 13, 'M': 14,
                   'S': 15, 'T': 16, 'Y': 17, 'V': 18, 'W': 19}

def isChainStart(residueId, chain_starts, n=2):
    # return true if residue is near chain starts.
    # n=0 means always return False
    # n=1 means only return True if residue is the the first residue of a chain.
    # n=2 means return True if residue is the first or the one nearest to the first residue of a chain.
    atBegin = False
    for i in range(n):
        if (residueId-i) in chain_starts:
            atBegin = True
    return atBegin

def isChainEnd(residueId, chain_ends, n=2):
    # return true if residue is near chain ends.
    # n=0 means always return False
    # n=1 means only return True if residue is the the last residue of a chain.
    # n=2 means return True if residue is the last or the one nearest to the last residue of a chain.
    atEnd = False
    for i in range(n):
        if (residueId+i) in chain_ends:
            atEnd = True
    return atEnd
def isChainEdge(residueId, chain_starts, chain_ends, n=2):
    # n is how far away from the two ends count as in chain edge.
    return (isChainStart(residueId, chain_starts, n) or isChainEnd(residueId, chain_ends, n))
    # atBegin = False
    # atEnd = False
    # for i in range(n):
    #     if (residueId-i) in chain_starts:
    #         atBegin = True
    # for i in range(n):
    #     if (residueId+i) in chain_ends:
    #         atEnd = True
    # return (atBegin or atEnd)

def inWhichChain(residueId, chain_ends):
    chain_table = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'
    for i, end_of_chain_resId in enumerate(chain_ends):
        if end_of_chain_resId < residueId:
            pass
        else:
            return chain_table[i]

def read_beta_parameters():
    ### directly copied from Nick Schafer's
    # os.chdir(parameter_directory)
    in_anti_HB = open("anti_HB", 'r').readlines()
    in_anti_NHB = open("anti_NHB", 'r').readlines()
    in_para_HB = open("para_HB", 'r').readlines()
    in_para_one = open("para_one", 'r').readlines()
    in_anti_one = open("anti_one", 'r').readlines()

    p_par = np.zeros((20))
    p_anti = np.zeros((20))
    p_antihb = np.zeros((20,20,2))
    p_antinhb = np.zeros((20,20,2))
    p_parhb = np.zeros((20,20,2))

    for i in range(20):
        p_par[i] = float(in_para_one[i].strip())
        p_anti[i] = float(in_anti_one[i].strip())
        for j in range(20):
            p_antihb[i][j][0] = float(in_anti_HB[i].strip().split()[j])
            p_antinhb[i][j][0] = float(in_anti_NHB[i].strip().split()[j])
            p_parhb[i][j][0] = float(in_para_HB[i].strip().split()[j])

    for i in range(20):
        for j in range(20):
            p_antihb[i][j][1] = float(in_anti_HB[i+21].strip().split()[j])
            p_antinhb[i][j][1] = float(in_anti_NHB[i+21].strip().split()[j])
            p_parhb[i][j][1] = float(in_para_HB[i+21].strip().split()[j])
    return p_par, p_anti, p_antihb, p_antinhb, p_parhb


def get_lambda_by_index(i, j, lambda_i):


    lambda_table = [[1.37, 1.36, 1.17],
                    [3.89, 3.50, 3.52],
                    [0.00, 3.47, 3.62]]
    if abs(j-i) >= 4 and abs(j-i) < 18:
        return lambda_table[lambda_i][0]
    elif abs(j-i) >= 18 and abs(j-i) < 45:
        return lambda_table[lambda_i][1]
    elif abs(j-i) >= 45:
        return lambda_table[lambda_i][2]
    else:
        return 0

def get_alpha_by_index(i, j, alpha_i):
    alpha_table = [[1.30, 1.30, 1.30],
                    [1.32, 1.32, 1.32],
                    [1.22, 1.22, 1.22],
                    [0.00, 0.33, 0.33],
                    [0.00, 1.01, 1.01]]
    if abs(j-i) >= 4 and abs(j-i) < 18:
        return alpha_table[alpha_i][0]
    elif abs(j-i) >= 18 and abs(j-i) < 45:
        return alpha_table[alpha_i][1]
    elif abs(j-i) >= 45:
        return alpha_table[alpha_i][2]
    else:
        return 0

def get_pap_gamma_APH(donor_idx, acceptor_idx, chain_i, chain_j, gamma_APH):
    # if chain_i == chain_j and abs(j-i) < 13 or abs(j-i) > 16:
    # if abs(j-i) < 13 or abs(j-i) > 16:
    # if i-j < 13 or i-j > 16:
    # if (donor_idx - acceptor_idx >= 13 and donor_idx - acceptor_idx <= 16) or chain_i != chain_j:
    if (donor_idx - acceptor_idx >= 13 and donor_idx - acceptor_idx <= 16) and chain_i == chain_j:
        return gamma_APH
    else:
        return 0

def get_pap_gamma_AP(donor_idx, acceptor_idx, chain_i, chain_j, gamma_AP, ssweight):
    if ssweight[donor_idx][1] == 1 and ssweight[acceptor_idx][1] == 1:
        additional_scale = 1.5
    else:
        additional_scale = 1.0
    # if (donor_idx - acceptor_idx >= 17):
    if (donor_idx - acceptor_idx >= 17) or chain_i != chain_j:
        return additional_scale * gamma_AP
    else:
        return 0

def get_pap_gamma_P(donor_idx, acceptor_idx, chain_i, chain_j, gamma_P, ssweight):
    if ssweight[donor_idx][1] == 1 and ssweight[acceptor_idx][1] == 1:
        additional_scale = 1.5
    else:
        additional_scale = 1.0
    if (donor_idx - acceptor_idx >= 9) or chain_i != chain_j:
        return additional_scale * gamma_P
    else:
        return 0

def get_Lambda_2(i, j, p_par, p_anti, p_antihb, p_antinhb, p_parhb, a):
    Lambda = get_lambda_by_index(i, j, 1)
    Lambda += -0.5*get_alpha_by_index(i, j, 0)*p_antihb[a[i], a[j]][0]
    Lambda += -0.25*get_alpha_by_index(i, j, 1)*(p_antinhb[a[i+1], a[j-1]][0] + p_antinhb[a[i-1], a[j+1]][0])
    Lambda += -get_alpha_by_index(i, j, 2)*(p_anti[a[i]] + p_anti[a[j]])
    return Lambda

def get_Lambda_3(i, j, p_par, p_anti, p_antihb, p_antinhb, p_parhb, a):
    Lambda = get_lambda_by_index(i, j, 2)
    Lambda += -get_alpha_by_index(i, j, 3)*p_parhb[a[i+1], a[j]][0]
    Lambda += -get_alpha_by_index(i, j, 4)*p_par[a[i+1]]
    Lambda += -get_alpha_by_index(i, j, 3)*p_par[a[j]]
    return Lambda


# def beta_term_1(oa, k_beta=4.184):
#     print("beta_1 term ON")
#     nres, n, h, ca, o, res_type = oa.nres, oa.n, oa.h, oa.ca, oa.o, oa.res_type
#     # print(lambda_1)
#     r_ON = .298
#     sigma_NO = .068
#     r_OH = .206
#     sigma_HO = .076

#     lambda_1 = np.zeros((nres, nres))
#     for i in range(nres):
#         for j in range(nres):
#             lambda_1[i][j] = get_lambda_by_index(i, j, 0)
#     theta_ij = f"exp(-(r_Oi_Nj-{r_ON})^2/(2*{sigma_NO}^2)-(r_Oi_Hj-{r_OH})^2/(2*{sigma_HO}^2))"
#     beta_string_1 = f"-{k_beta}*lambda_1(res_i,res_j)*theta_ij;theta_ij={theta_ij};r_Oi_Nj=distance(a1,d1);r_Oi_Hj=distance(a1,d2);"
#     beta_1 = CustomHbondForce(beta_string_1)
#     beta_1.addPerDonorParameter("res_i")
#     beta_1.addPerAcceptorParameter("res_j")
#     beta_1.addTabulatedFunction("lambda_1", Discrete2DFunction(nres, nres, lambda_1.T.flatten()))
#     # print(lambda_1)
#     # print(len(oa.o), nres)
#     for i in range(nres):
#         if oa.o[i]!= -1:
#             beta_1.addAcceptor(oa.o[i], -1, -1, [i])
#         if oa.n[i]!=-1 and oa.h[i]!=-1:
#             beta_1.addDonor(oa.n[i], oa.h[i], -1, [i])
#     beta_1.setNonbondedMethod(CustomHbondForce.CutoffNonPeriodic)
#     beta_1.setCutoffDistance(1.0)
#     beta_1.setForceGroup(23)
#     # beta_2.setForceGroup(24)
#     # beta_3.setForceGroup(25)
#     return beta_1

def convert_units(k):
    if isinstance(k, float) or isinstance(k, int):
        k = k   # just for backward comptable
    elif isinstance(k, Quantity):
        k = k.value_in_unit(kilojoule_per_mole)   # convert to kilojoule_per_mole, openMM default uses kilojoule_per_mole as energy.
    else:
        print(f"Unknown input, {k}, {type(k)}")
    return k

def AM_beta_term_1(oa, k_beta=1*kilocalories_per_mole, periodic=False, forceGroup=27):
    print("AM beta_1 term ON")
    k_beta = convert_units(k_beta) * oa.k_awsem
    nres, n, h, ca, o, res_type = oa.nres, oa.n, oa.h, oa.ca, oa.o, oa.res_type
    # print(lambda_1)
    r_ON = .298
    sigma_NO = .068
    r_OH = .206
    sigma_HO = .076

    lambda_1 = np.zeros((nres, nres))
    for i in range(nres):
        for j in range(nres):
            lambda_1[i][j] = get_lambda_by_index(i, j, 0)
    theta_ij = f"exp(-(r_Oi_Nj-{r_ON})^2/(2*{sigma_NO}^2)-(r_Oi_Hj-{r_OH})^2/(2*{sigma_HO}^2))"
    mu_1 = 10  # nm^-1
    # mu_2 = 5   # nm^-1
    rcHB = 1.2  # in nm
    # v1i ensures the hydrogen bonding does not occur when five residue segment is shorter than 12 A
    # v1i = f"0.5*(1+tanh({mu_1}*(distance(a2,a3)-{rcHB})))"
    v1i = "1"
    beta_string_1 = f"-{k_beta}*lambda_1(res_i,res_j)*theta_ij*v1i;theta_ij={theta_ij};v1i={v1i};r_Oi_Nj=distance(a1,d1);r_Oi_Hj=distance(a1,d2);"
    beta_1 = CustomHbondForce(beta_string_1)
    beta_1.addPerDonorParameter("res_i")
    beta_1.addPerAcceptorParameter("res_j")
    beta_1.addTabulatedFunction("lambda_1", Discrete2DFunction(nres, nres, lambda_1.T.flatten()))
    # print(lambda_1)
    # print(len(oa.o), nres)
    for i in range(nres):
        if oa.o[i]!= -1:
            ca_i_minus_2 = oa.ca[0] if i <= 2 else oa.ca[i-2]
            ca_i_plus_2 = oa.ca[-1] if i+2 >= nres else oa.ca[i+2]
            # beta_1.addAcceptor(oa.o[i], ca_i_minus_2, ca_i_plus_2, [i])
            beta_1.addAcceptor(oa.o[i], -1, -1, [i])
        if oa.n[i]!=-1 and oa.h[i]!=-1:
            beta_1.addDonor(oa.n[i], oa.h[i], -1, [i])
#    beta_1.setNonbondedMethod(CustomHbondForce.CutoffNonPeriodic)
#    beta_1.setCutoffDistance(1.0)



    if periodic:
        beta_1.setNonbondedMethod(beta_1.CutoffPeriodic)
    else:
        beta_1.setNonbondedMethod(beta_1.CutoffNonPeriodic)
    beta_1.setCutoffDistance(1.0)
    print("Contact cutoff ", beta_1.getCutoffDistance())
    print("NonbondedMethod: ", beta_1.getNonbondedMethod())



    beta_1.setForceGroup(forceGroup)
    print('k_beta=' + str(k_beta) + '\n')
    return beta_1

def AM_beta_term_2(oa, k_beta=1*kilocalories_per_mole, periodic=False, forceGroup=27):
    print("AM beta_2 term ON")
    k_beta = convert_units(k_beta) * oa.k_awsem
    nres, n, h, ca, o, res_type = oa.nres, oa.n, oa.h, oa.ca, oa.o, oa.res_type
    # print(lambda_1)
    r_ON = .298
    sigma_NO = .068
    r_OH = .206
    sigma_HO = .076
    eta_beta_1 = 10.0
    eta_beta_2 = 5.0
    # r_HB_c = 0.4
    r_HB_c = 1.2
    p_par, p_anti, p_antihb, p_antinhb, p_parhb = read_beta_parameters()

    # for lookup table.
    a = []
    for ii in range(oa.nres):
        a.append(se_map_1_letter[oa.seq[ii]])

    lambda_2 = np.zeros((nres, nres))
    for i in range(nres):
        for j in range(nres):
            if isChainEdge(i, oa.chain_starts, oa.chain_ends, n=1) or \
                    isChainEdge(j, oa.chain_starts, oa.chain_ends, n=1):
                continue
            lambda_2[i][j] = get_Lambda_2(i, j, p_par, p_anti, p_antihb, p_antinhb, p_parhb, a)
    theta_ij = f"exp(-(r_Oi_Nj-{r_ON})^2/(2*{sigma_NO}^2)-(r_Oi_Hj-{r_OH})^2/(2*{sigma_HO}^2))"
    theta_ji = f"exp(-(r_Oj_Ni-{r_ON})^2/(2*{sigma_NO}^2)-(r_Oj_Hi-{r_OH})^2/(2*{sigma_HO}^2))"
    beta_string_2 = f"-{k_beta}*lambda_2(res_i,res_j)*theta_ij*theta_ji;\
                        theta_ij={theta_ij};r_Oi_Nj=distance(a1,d1);r_Oi_Hj=distance(a1,d2);\
                        theta_ji={theta_ji};r_Oj_Ni=distance(d3,a2);r_Oj_Hi=distance(d3,a3);"
    beta_2 = CustomHbondForce(beta_string_2)
    beta_2.addPerDonorParameter("res_i")
    beta_2.addPerAcceptorParameter("res_j")
    beta_2.addTabulatedFunction("lambda_2", Discrete2DFunction(nres, nres, lambda_2.T.flatten()))
    # print(lambda_1)
    # print(len(oa.o), nres)
    for i in range(nres):
        if o[i]!= -1 and n[i]!=-1 and h[i]!=-1:
            beta_2.addAcceptor(o[i], n[i], h[i], [i])
            beta_2.addDonor(n[i], h[i], o[i], [i])
#    beta_2.setNonbondedMethod(CustomHbondForce.CutoffNonPeriodic)
#    beta_2.setCutoffDistance(1.0)



    if periodic:
        beta_2.setNonbondedMethod(beta_2.CutoffPeriodic)
    else:
        beta_2.setNonbondedMethod(beta_2.CutoffNonPeriodic)
    beta_2.setCutoffDistance(1.0)
    print("Contact cutoff ", beta_2.getCutoffDistance())
    print("NonbondedMethod: ", beta_2.getNonbondedMethod())




    beta_2.setForceGroup(forceGroup)
    print('k_beta=' + str(k_beta) + '\n')
    return beta_2


def AM_beta_term_3(oa, k_beta=1*kilocalories_per_mole, periodic=False, forceGroup=27):
    print("AM beta_3 term ON")
    k_beta = convert_units(k_beta) * oa.k_awsem
    nres, n, h, ca, o, res_type = oa.nres, oa.n, oa.h, oa.ca, oa.o, oa.res_type
    # print(lambda_1)
    r_ON = .298
    sigma_NO = .068
    r_OH = .206
    sigma_HO = .076
    eta_beta_1 = 10.0
    eta_beta_2 = 5.0
    # r_HB_c = 0.4
    r_HB_c = 1.2
    p_par, p_anti, p_antihb, p_antinhb, p_parhb = read_beta_parameters()

    # for lookup table.
    a = []
    for ii in range(oa.nres):
        a.append(se_map_1_letter[oa.seq[ii]])

    lambda_3 = np.zeros((nres, nres))
    for i in range(nres):
        for j in range(nres):
            if isChainEdge(i, oa.chain_starts, oa.chain_ends, n=1) or \
                    isChainEdge(j, oa.chain_starts, oa.chain_ends, n=1):
                continue
            lambda_3[i][j] = get_Lambda_3(i, j, p_par, p_anti, p_antihb, p_antinhb, p_parhb, a)

    theta_ij = f"exp(-(r_Oi_Nj-{r_ON})^2/(2*{sigma_NO}^2)-(r_Oi_Hj-{r_OH})^2/(2*{sigma_HO}^2))"
    theta_jip2 = f"exp(-(r_Oj_Nip2-{r_ON})^2/(2*{sigma_NO}^2)-(r_Oj_Hip2-{r_OH})^2/(2*{sigma_HO}^2))"

    beta_string_3 = f"-{k_beta}*lambda_3(res_i,res_j)*theta_ij*theta_jip2;\
                        theta_ij={theta_ij};r_Oi_Nj=distance(a1,d1);r_Oi_Hj=distance(a1,d2);\
                        theta_jip2={theta_jip2};r_Oj_Nip2=distance(d3,a2);r_Oj_Hip2=distance(d3,a3);"
    beta_3 = CustomHbondForce(beta_string_3)

    beta_3.addPerDonorParameter("res_i")
    beta_3.addPerAcceptorParameter("res_j")
    beta_3.addTabulatedFunction("lambda_3", Discrete2DFunction(nres, nres, lambda_3.T.flatten()))
    # print(lambda_1)
    # print(len(oa.o), nres)
    for i in range(nres):
        if isChainEdge(i, oa.chain_starts, oa.chain_ends, n=2):
            continue
        if o[i] != -1 and n[i+2] !=-1 and h[i+2] !=-1:
            beta_3.addAcceptor(o[i], n[i+2], h[i+2], [i])
        if o[i] != -1 and n[i] !=-1 and h[i] !=-1:
            beta_3.addDonor(n[i], h[i], o[i], [i])
#    beta_3.setNonbondedMethod(CustomHbondForce.CutoffNonPeriodic)
#    beta_3.setCutoffDistance(1.0)




    if periodic:
        beta_3.setNonbondedMethod(beta_3.CutoffPeriodic)
    else:
        beta_3.setNonbondedMethod(beta_3.CutoffNonPeriodic)
    beta_3.setCutoffDistance(1.0)
    print("Contact cutoff ", beta_3.getCutoffDistance())
    print("NonbondedMethod: ", beta_3.getNonbondedMethod())



    beta_3.setForceGroup(forceGroup)
    print('k_beta=' + str(k_beta) + '\n')

    return beta_3


def AM_pap_term_1(oa, k_pap=1*kilocalories_per_mole, dis_i_to_i4=1.2, forceGroup=28, ssweightFileName="ssweight", periodic=False):
    print("AM pap_1 term ON")
    k_pap = convert_units(k_pap) * oa.k_awsem
    # dis_i_to_i4 should be in nm, it disfavor hydrogen bond when ca_i and ca_i+4 are 1.2 nm apart away.
    nres, ca = oa.nres, oa.ca
    # r0 = 2.0 # nm
    r0 = 0.8  # nm
    eta_pap = 70  # nm^-1
    eta_1 = 10 #nm^-1
    gamma_aph = 1.0
    gamma_ap = 0.4
    gamma_p = 0.4
    r_ON = .298
    sigma_NO = .068


    if not os.path.exists(ssweightFileName):
        print("No ssweight given, assume all zero")
        ssweight = np.zeros((nres, 2))
    else:
        ssweight = np.loadtxt(ssweightFileName)

    gamma_1 = np.zeros((nres, nres))
    gamma_2 = np.zeros((nres, nres))
    for i in range(nres):
        for j in range(nres):
            resId1 = i
            chain1 = inWhichChain(resId1, oa.chain_ends)
            resId2 = j
            chain2 = inWhichChain(resId2, oa.chain_ends)
            gamma_1[i][j] = get_pap_gamma_APH(i, j, chain1, chain2, gamma_aph)
            gamma_2[i][j] = get_pap_gamma_AP(i, j, chain1, chain2, gamma_ap, ssweight)

    constraint_i_and_i4 = f"0.5*(1+tanh({eta_1}*(distance(a1,a2)-{dis_i_to_i4})))"

#    pap_function = f"-{k_pap}*(gamma_1(donor_idx,acceptor_idx)+gamma_2(donor_idx,acceptor_idx))\
#                        *0.5*(1+tanh({eta_pap}*({r0}-distance(a1,d1))))\
#                        *0.5*(1+tanh({eta_pap}*({r0}-distance(a2,d2))))\
#                        *{constraint_i_and_i4}"

    theta_ij = f"exp(-(r_Oi_Nj-{r_ON})^2/(2*{sigma_NO}^2))"
    theta_ji = f"exp(-(r_Oj_Ni-{r_ON})^2/(2*{sigma_NO}^2))"
    pap_function = f"-{k_pap}*(gamma_1(donor_idx,acceptor_idx)+gamma_2(donor_idx,acceptor_idx))*theta_ij*theta_ji\
                        *{constraint_i_and_i4};\
                        theta_ij={theta_ij};r_Oi_Nj=distance(a1,d1);\
                        theta_ji={theta_ji};r_Oj_Ni=distance(a2,d2);"

    pap = CustomHbondForce(pap_function)
    pap.addPerDonorParameter("donor_idx")
    pap.addPerAcceptorParameter("acceptor_idx")
    pap.addTabulatedFunction("gamma_1", Discrete2DFunction(nres, nres, gamma_1.T.flatten()))
    pap.addTabulatedFunction("gamma_2", Discrete2DFunction(nres, nres, gamma_2.T.flatten()))
    # print(ca)
    # count = 0;
    i = 0

    for i in range(nres):
        if not isChainEnd(i, oa.chain_ends, n=4):
            if oa.o[i] != -1 and oa.o[i+4] != -1:
                pap.addAcceptor(oa.o[i], oa.o[i+4], -1, [i])

        if not isChainStart(i, oa.chain_starts, n=4):
            if oa.n[i] != -1 and oa.n[i-4] != -1:
                pap.addDonor(oa.n[i], oa.n[i-4], -1, [i])

#    pap.setNonbondedMethod(CustomHbondForce.CutoffNonPeriodic)
#    pap.setCutoffDistance(1.0)


    if periodic:
        pap.setNonbondedMethod(pap.CutoffPeriodic)
    else:
        pap.setNonbondedMethod(pap.CutoffNonPeriodic)
    pap.setCutoffDistance(1.0)
    print("Contact cutoff ", pap.getCutoffDistance())
    print("NonbondedMethod: ", pap.getNonbondedMethod())



    # print(count)
    pap.setForceGroup(forceGroup)
    print('k_pap=' + str(k_pap) + '\n')
    return pap

def AM_pap_term_2(oa, k_pap=1*kilocalories_per_mole, dis_i_to_i4=1.2, forceGroup=28, ssweightFileName="ssweight", periodic=False):
    print("AM pap_2 term ON")
    k_pap = convert_units(k_pap) * oa.k_awsem
    nres, ca = oa.nres, oa.ca
    # r0 = 2.0 # nm
    r0 = 0.8  # nm
    eta_pap = 70  # nm^-1
    eta_1 = 10 # nm^-1
    gamma_aph = 1.0
    gamma_ap = 0.4
    gamma_p = 0.4
    r_ON = .298
    sigma_NO = .068


    if not os.path.exists(ssweightFileName):
        print("No ssweight given, assume all zero")
        ssweight = np.zeros((nres, 2))
    else:
        ssweight = np.loadtxt(ssweightFileName)

    gamma_3 = np.zeros((nres, nres))
    for i in range(nres):
        for j in range(nres):
            resId1 = i
            chain1 = inWhichChain(resId1, oa.chain_ends)
            resId2 = j
            chain2 = inWhichChain(resId2, oa.chain_ends)
            gamma_3[i][j] = get_pap_gamma_P(i, j, chain1, chain2, gamma_p, ssweight)


    constraint_i_and_i4 = f"0.5*(1+tanh({eta_1}*(distance(a1,a2)-{dis_i_to_i4})))"
#    pap_function = f"-{k_pap}*gamma_3(donor_idx,acceptor_idx)\
#                        *0.5*(1+tanh({eta_pap}*({r0}-distance(a1,d1))))\
#                        *0.5*(1+tanh({eta_pap}*({r0}-distance(a2,d2))))\
#                        *{constraint_i_and_i4}"


    theta_ij = f"exp(-(r_Oi_Nj-{r_ON})^2/(2*{sigma_NO}^2))"
    theta_ji = f"exp(-(r_Oj_Ni-{r_ON})^2/(2*{sigma_NO}^2))"
    pap_function = f"-{k_pap}*gamma_3(donor_idx,acceptor_idx)*theta_ij*theta_ji\
                        *{constraint_i_and_i4};\
                        theta_ij={theta_ij};r_Oi_Nj=distance(a1,d1);\
                        theta_ji={theta_ji};r_Oj_Ni=distance(a2,d2);"



    pap = CustomHbondForce(pap_function)
    pap.addPerDonorParameter("donor_idx")
    pap.addPerAcceptorParameter("acceptor_idx")
    pap.addTabulatedFunction("gamma_3", Discrete2DFunction(nres, nres, gamma_3.T.flatten()))
    # print(oa.n)
    # count = 0;
    for i in range(nres):
        if not isChainEnd(i, oa.chain_ends, n=4):
            if oa.o[i] != -1 and oa.o[i+4] != -1:
                pap.addAcceptor(oa.o[i], oa.o[i+4], -1, [i])
            if oa.n[i] != -1 and oa.n[i+4] != -1:
                pap.addDonor(oa.n[i], oa.n[i+4], -1, [i])

#    pap.setNonbondedMethod(CustomHbondForce.CutoffNonPeriodic)
#    pap.setCutoffDistance(1.0)


    if periodic:
        pap.setNonbondedMethod(pap.CutoffPeriodic)
    else:
        pap.setNonbondedMethod(pap.CutoffNonPeriodic)
    pap.setCutoffDistance(1.0)
    print("Contact cutoff ", pap.getCutoffDistance())
    print("NonbondedMethod: ", pap.getNonbondedMethod())



    # print(count)
    pap.setForceGroup(forceGroup)
    print('k_pap=' + str(k_pap) + '\n')
    return pap


