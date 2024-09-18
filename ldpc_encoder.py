#LDPC encoder

import numpy as np
import random
from functools import reduce
from pyldpc import utils
from scipy.sparse import csr_matrix


def regular_coding_matrix_systematic(N, col_w, row_w, sparse=True):
    """Compute a coding matrix G in systematic format with an identity block.
    Parameters
    ----------
    H: array (n_equations, n_code). Parity-check matrix obtained from "MakeIrregularLDPCCode".
    sparse: (boolean, default True): if `True`, scipy.sparse is used
    to speed up computation if n_code > 1000.
    Returns
    -------
    H_new: (n_equations, n_code) array. Modified parity-check matrix given by a
        permutation of the columns of the provided H.
    G_systematic.T: Transposed Systematic Coding matrix associated to H_new.
    Notions
    -------
    H_new.shape may != H.shape
    """
    K = N - (N * col_w) / row_w
    '''
    while True:
        H = enc_H_gallager(N, col_w, row_w)
        if np.linalg.matrix_rank(H) == (N * col_w) / row_w:
            break
        else:
            print('re-generate H_matrix')
            continue
    '''
    H = enc_H_gallager(N, col_w, row_w)
    n_equations, n_code = H.shape

    if n_code > 1000 or sparse:
        sparse = True
    else:
        sparse = False

    P1 = np.identity(n_code, dtype=int)

    Hrowreduced = utils.gaussjordan(H)

    n_bits = n_code - sum([a.any() for a in Hrowreduced])

    # After this loop, Hrowreduced will have the form H_ss : | I_(n-k)  A |

    while(True):
        zeros = [i for i in range(min(n_equations, n_code))
                 if not Hrowreduced[i, i]]
        if len(zeros):
            indice_colonne_a = min(zeros)
        else:
            break
        list_ones = [j for j in range(indice_colonne_a + 1, n_code)
                     if Hrowreduced[indice_colonne_a, j]]
        if len(list_ones):
            indice_colonne_b = min(list_ones)
        else:
            break
        aux = Hrowreduced[:, indice_colonne_a].copy()
        Hrowreduced[:, indice_colonne_a] = Hrowreduced[:, indice_colonne_b]
        Hrowreduced[:, indice_colonne_b] = aux

        aux = P1[:, indice_colonne_a].copy()
        P1[:, indice_colonne_a] = P1[:, indice_colonne_b]
        P1[:, indice_colonne_b] = aux

    # Now, Hrowreduced has the form: | I_(n-k)  A | ,
    # the permutation above makes it look like :
    # |A  I_(n-k)|

    P1 = P1.T
    identity = list(range(n_code))
    sigma = identity[n_code - n_bits:] + identity[:n_code - n_bits]

    P2 = np.zeros(shape=(n_code, n_code), dtype=int)
    P2[identity, sigma] = np.ones(n_code)

    if sparse:
        P1 = csr_matrix(P1)
        P2 = csr_matrix(P2)
        H = csr_matrix(H)

    P = utils.binaryproduct(P2, P1)

    if sparse:
        P = csr_matrix(P)

    H_new = utils.binaryproduct(H, np.transpose(P))

    G_systematic = np.zeros((n_bits, n_code), dtype=int)
    G_systematic[:, :n_bits] = np.identity(n_bits)
    G_systematic[:, n_bits:] = \
        (Hrowreduced[:n_code - n_bits, n_code - n_bits:]).T

    return H_new.astype(int), G_systematic.T
def make_irregular_LDPC(N, K, var_degrees, check_degrees):
    """
    Parameters
    ----------
    N:               Length of the code.
    K:               Dimension of the code.
    check_degrees:    Hash indicating how many checks have 1 connection,
                     2 connections, 3 connections, etc.
    var_degrees:      Hash indicating how many vars have 1 connection,
                     2 connections, 3 connections, etc.

    Returns
    -------
    H_matrix : a parity-check matrix
    example: H_matrix = make_irregular_LDPC(10,3,{1:4,2:4,3:2},{2:3,3:4})

    """

    M = N - K
    num_checks = reduce(lambda x, y: x+y, list(check_degrees.values()))
    assert num_checks == M, (
        'Number of checks in check_degrees sums to ' + 'num_checks' +
        '!=' + 'M' + '.')
    assert reduce(lambda x, y: x+y, list(var_degrees.values())) == N, (
        'Number of vars in var_degrees does not sum to N.')


    result = [0]*M  # output an array with M dimension

    vars = []
    cur_var_index = 0
    for d in list(var_degrees.keys()):  # for each possible var degree
        for i in range(var_degrees[d]):  # for each of var with that degree
            vars.extend([cur_var_index]*d)
            cur_var_index = cur_var_index+1
    assert cur_var_index == N
    cur_check_index = 0
    for d in list(check_degrees.keys()):  # for each possible check degree
        for i in range(check_degrees[d]):  # for each check with that degree
            result[cur_check_index] = [None]*d
            for connectionForCheck in range(d):
                vIndex = random.randrange(len(vars))
                if (result[cur_check_index].count(vars[vIndex]) == 0):
                    result[cur_check_index][connectionForCheck] = vars.pop(
                        vIndex)
                else:
                    
                    vars.pop(vIndex)
            while (result[cur_check_index].count(None) > 0):
                result[cur_check_index].pop(result[cur_check_index].index(None))
            cur_check_index = cur_check_index + 1
    assert len(vars) == 0, 'vars should be empty but it is' + 'vars'

    # convert link array to matrix-type, H_matrix
    H_matrix = np.zeros((M, N))
    for i in range(M):
        for j in result[i]:
            H_matrix[i, j] = 1
    H_matrix = H_matrix.astype(int)
    return H_matrix


def compute_lam_rho(terms):
    """
    compute_lam_rho(terms):
    terms:   A hash table containing lambda or rho coefficients.
             For example, if lambda(x) = .4 x^2 + .6 x^3, then
             the input would be {2:.4, 3:.6}.
    Returns  The integral of the argument from 0 to 1.
    """
    sum = 0
    total = 0
    for i in list(terms.keys()):
        sum = sum + terms[i]
        total = total + terms[i]/float(i)
    assert(abs(sum-1.0) < .0001)
    return total


def compute_degree_seq(N, M, lam, rho):
    """
    N:     Block size.
    M:     Number of constraints.
    lam:   A hash table specifying the variable degress using the lambda
           notation.  Specifically, lam[i] = p denotes that the fraction
           of EDGES coming from a variable node of degree i is p.
    rho:   A hash table specifying the check degress using the rho
           notation.  Specifically, rho[i] = p denotes that the fraction
           of EDGES coming from a check node of degree i is p.

    Returns a pair of hash tables (var_degrees, check_degrees) where
    varDegress[i] = y indicates that there are y varialbes of degree i.
    """
    total_var_edges = float(N)/compute_lam_rho(lam)
    total_check_edges = float(M)/compute_lam_rho(rho)

    var_degrees = {}
    for key in list(lam.keys()):
        var_degrees[key] = int(round(total_var_edges*lam[key]/float(key)))

    check_degrees = {}
    for key in list(rho.keys()):
        check_degrees[key] = int(round(total_check_edges*rho[key]/float(key)))

    return (var_degrees, check_degrees)


def compute_edge_mismatch(var_degrees, check_degrees):
    edges_checks = np.dot(list(check_degrees.keys()),
                             list(check_degrees.values()))
    edges_vars = np.dot(list(var_degrees.keys()), list(var_degrees.values()))
    edge_mismatch = edges_checks-edges_vars
    return edge_mismatch


def make_irregular_LDPC_from_degree(N, K, lam, rho):
    """
    make_irregular_LDPC_from_degree(N,K,lam,rho):
    N:     Block size.
    K:     Dimension.
    lam:   A hash table specifying the variable degress using the lambda
           notation.  Specifically, lam[i] = p denotes that the fraction
           of EDGES coming from a variable node of degree i is p.
    rho:   A hash table specifying the check degress using the rho
           notation.  Specifically, rho[i] = p denotes that the fraction
           of EDGES coming from a check node of degree i is p.

    This function creates an irregular LDPC code for the desired parameters.
    """
    M = N - K
    total = 0

    (var_degrees, check_degrees) = compute_degree_seq(N, M, lam, rho)
    for key in list(check_degrees.keys()):
        total = total + check_degrees[key]
    cCleanupIndex = list(check_degrees.values()).index(
        max(check_degrees.values()))
    cCleanupIndex = list(check_degrees.keys())[cCleanupIndex]
    check_degrees[cCleanupIndex] = check_degrees[cCleanupIndex] - (total-M)
    assert check_degrees[cCleanupIndex] > 0

    total = 0
    for key in list(var_degrees.keys()):
        total = total + var_degrees[key]

    vCleanupIndex = list(var_degrees.values()).index(max(var_degrees.values()))
    vCleanupIndex = list(var_degrees.keys())[vCleanupIndex]
    var_degrees[vCleanupIndex] = var_degrees[vCleanupIndex] - (total-N)
    assert var_degrees[vCleanupIndex] > 0
    edge_mismatch = compute_edge_mismatch(var_degrees, check_degrees)



    k = list(var_degrees.keys())
    sorted(k)

    degreeDiff = k[1]-k[0]
    edgeDiff = edge_mismatch/degreeDiff + 1
    var_degrees[k[0]] = int(var_degrees[k[0]]-edgeDiff)
    var_degrees[k[1]] = int(var_degrees[k[1]]+edgeDiff)
    assert var_degrees[k[0]] > 0
    assert var_degrees[k[1]] > 0

    edge_mismatch = compute_edge_mismatch(var_degrees, check_degrees)
    k = list(check_degrees.keys())
    sorted(k)
    if (edge_mismatch < 0):
        check_degrees[k[0]] = check_degrees[k[0]] - 1
        edge_mismatch = edge_mismatch - k[0]
        if -edge_mismatch not in check_degrees:
            check_degrees[-edge_mismatch] = 0
        check_degrees[-edge_mismatch] = int(check_degrees[-edge_mismatch]+1)

    else:
        # haven't yet implemented this case
        assert 0

    return make_irregular_LDPC(N, K, var_degrees, check_degrees)


def coding_matrix_systematic(K, N, lam, rho, sparse=True):
    """Compute a coding matrix G in systematic format with an identity block.
    Parameters
    ----------
    H: array (n_equations, n_code). Parity-check matrix obtained from "MakeIrregularLDPCCode".
    sparse: (boolean, default True): if `True`, scipy.sparse is used
    to speed up computation if n_code > 1000.
    Returns
    -------
    H_new: (n_equations, n_code) array. Modified parity-check matrix given by a
        permutation of the columns of the provided H.
    G_systematic.T: Transposed Systematic Coding matrix associated to H_new.
    Notions
    -------
    H_new.shape may != H.shape
    """
    while True:
        H = make_irregular_LDPC_from_degree(N, K, lam, rho)
        if np.linalg.matrix_rank(H) == N - K:
            break
        else:

            continue

    n_equations, n_code = H.shape

    if n_code > 1000 or sparse:
        sparse = True
    else:
        sparse = False

    P1 = np.identity(n_code, dtype=int)

    Hrowreduced = utils.gaussjordan(H)

    n_bits = n_code - sum([a.any() for a in Hrowreduced])

    # After this loop, Hrowreduced will have the form H_ss : | I_(n-k)  A |

    while(True):
        zeros = [i for i in range(min(n_equations, n_code))
                 if not Hrowreduced[i, i]]
        if len(zeros):
            indice_colonne_a = min(zeros)
        else:
            break
        list_ones = [j for j in range(indice_colonne_a + 1, n_code)
                     if Hrowreduced[indice_colonne_a, j]]
        if len(list_ones):
            indice_colonne_b = min(list_ones)
        else:
            break
        aux = Hrowreduced[:, indice_colonne_a].copy()
        Hrowreduced[:, indice_colonne_a] = Hrowreduced[:, indice_colonne_b]
        Hrowreduced[:, indice_colonne_b] = aux

        aux = P1[:, indice_colonne_a].copy()
        P1[:, indice_colonne_a] = P1[:, indice_colonne_b]
        P1[:, indice_colonne_b] = aux

    # Now, Hrowreduced has the form: | I_(n-k)  A | ,
    # the permutation above makes it look like :
    # |A  I_(n-k)|

    P1 = P1.T
    identity = list(range(n_code))
    sigma = identity[n_code - n_bits:] + identity[:n_code - n_bits]

    P2 = np.zeros(shape=(n_code, n_code), dtype=int)
    P2[identity, sigma] = np.ones(n_code)

    if sparse:
        P1 = csr_matrix(P1)
        P2 = csr_matrix(P2)
        H = csr_matrix(H)

    P = utils.binaryproduct(P2, P1)

    if sparse:
        P = csr_matrix(P)

    H_new = utils.binaryproduct(H, np.transpose(P))

    G_systematic = np.zeros((n_bits, n_code), dtype=int)
    G_systematic[:, :n_bits] = np.identity(n_bits)
    G_systematic[:, n_bits:] = \
        (Hrowreduced[:n_code - n_bits, n_code - n_bits:]).T

    return H_new.astype(int), G_systematic.T


def ldpc_encoder(msg, G_matrix):
    """
    Parameters
    ----------
    msg : array(K)
        Message.
    generate_matrix : array(N,K)
        Transposed coding matrix obtained from "coding_matrix_systematic".

    Returns
    -------
    codeword : array(N)
        Coded message.

    """
    N, K = G_matrix.shape
    codeword = np.zeros(N)
    for i in range(N):
        codeword[i] = np.mod(np.sum(msg*G_matrix[i, :]), 2)
    codeword = codeword.astype(int)
    return codeword


def check_matrix(G_matrix, H_matrix):
    # check if G*H == 0
    N, K = G_matrix.shape
    for i in range(K):
        for j in range(N-K):
            check = np.mod(np.sum(G_matrix[:, i] * H_matrix[j, :]), 2)
            if not np.all(check == 0):
                print("Error in generating matrix and check matrix matching")
    
