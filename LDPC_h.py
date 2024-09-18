#LDPC codec in high err env


from pyldpc import utils
import copy
import math
from numpy import NaN, array, zeros, exp, log, empty, arctanh, tanh
import numpy as np
from functools import reduce

from scipy.sparse import csr_matrix
import random
random.seed(1)
np.random.seed(1)
#get regular Gallager
def enc_H_gallager(N, col_w, row_w):
    """
    enc_H_gallager(N,col_w,row_w):
    N:     Length of code.
    col_w: Column weight of code (i.e., how many 1's in a column).
    row_w: Row weight of code (i.e., how many 1's in a row).

    Create a regular LDPC code matrix using the construction in
    Gallager's book and return the result as a link array.  The
    return value, r, is a list of N lists where r[i][j] is the
    ith one in the jth row of H.
    The algorithm works by starting with a sub-matrix of n/row_w rows
    containing col_w 1's in the first positions in the first row,
    col_w 1's in the next set of positions in the second row, etc.

    H = [ 1 1 1 1 0 0 0 0 0 0 0 0 ...
          0 0 0 0 1 1 1 1 0 0 0 0 ...
          0 0 0 0 0 0 0 0 1 1 1 1 ...
          ...........................
    Next we create col_w-1 additional submatricies which are permuations of
    this and append them to the bottom of H.
    """
    num_rows = (N * col_w) / row_w
    num_sub_rows = int(num_rows / col_w)
    assert row_w*num_rows == N*col_w, 'N*col_w not divisible by row_w'
    assert (N/row_w)*row_w == N, 'N not divisible by row_w'

    H_sub = [0]*num_sub_rows
    for i in range(num_sub_rows):
        H_sub[i] = list(map(lambda x, y: x + y, [i*row_w]
                            * row_w, list(range(row_w))))
    H = copy.deepcopy(H_sub)

    for i in range(col_w-1):

        H_new_sub = [0]*num_sub_rows
        for m in range(num_sub_rows):
            H_new_sub[m] = [0]*row_w

        rand_perm = np.random.permutation(N)
        for j in range(num_sub_rows):
            for k in range(row_w):
                H_new_sub[j][k] = rand_perm[H_sub[j][k]]
        l = list(H_new_sub[j])
        l.sort()
        H_new_sub[j] = l
        H = H + copy.deepcopy(H_new_sub)
    #print(H)
    M = int(num_rows)
    H_matrix = np.zeros((M, N))
    for i in range(M):
        for j in H[i]:
            H_matrix[i, j] = 1
    H_matrix = H_matrix.astype(int)
    
    return H_matrix

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

#Generate non-regular LDPC codes

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


def hamming(x, y):

    len_ = x.size
    num_ = 0
    for i in range(len_):
        if x[i] != y[i]:
            num_ += 1
    # print(num_/len_)
    return num_/len_

def binaryproduct(X, Y):

    A = X.dot(Y)
    try:
        A = A.toarray()
    except AttributeError:
        pass
    return A % 2


def max_star_unit(x, y):
    output = max(x, y)+math.log(1+exp(-abs(x-y)))
    #output = max(x,y)                         
    return output


def max_star_fuc(nums_list):
    if len(nums_list) == 1:
        return nums_list[0]
    else:
        x = nums_list[0]
        y = nums_list[1]
        nums_list = nums_list[2:]
        output = max_star_unit(x, y)
        while True:
            if nums_list == []:
                break
            unit = nums_list[0]
            nums_list = nums_list[1:]
            output = max_star_unit(output, unit)
        return output


def _compute_M_N_Q_ex(H, n_equations, n_v, M_dict, N_dict, L_int, Q_metrics):
    for i in range(n_equations):
        for j in range(n_v):
            if H[i][j] == 1:
                M_dict[j].append(i)
                N_dict[i].append(j)
                Q_metrics[i, j] = L_int[j]



def _Iterative_R(H, n_equations, n_v, M_dict, N_dict, L_int, Q_metrics, R_metrics):
    for i in range(n_equations):
        for j in range(n_v):
            if H[i][j] == 1:
                unit = 1
                M_list = N_dict[i]
                for key in M_list:
                    if key == j:
                        continue
                    # print(unit,M_list,key,Q_metrics[i,key])
                    unit = unit * tanh(Q_metrics[i, key]/2)
                with np.errstate(divide='ignore', invalid='ignore'):
                    unit = max(min(unit,1),-1)
                    #R_value = max(min(2 * arctanh(unit),500),-500)
                    R_metrics[i, j] = max(min(2 * arctanh(unit),500),-500)

def _Isitzero(H, LLR_tatol):
    # print(LLR_tatol)
    code_list = []
    for i in LLR_tatol:
        if i > 0:
            code_list.append(0)
        else:
            code_list.append(1)
    # print(code_list)
    code = np.array(code_list)
    mat = np.dot(H, code)
    # print(mat)
    for i in mat:
        if i == 0 or i % 2 == 0:
            pass
        else:
            return False
    return [code, LLR_tatol]



def _compute_LLR(H, n_equations, n_v, M_dict, N_dict, L_int, Q_metrics, R_metrics, LLR_tatol, LLR_ex):
    for j in range(n_v):
        r_sum = 0
        for i in M_dict[j]:
            #print("i calculate")
            r_sum += R_metrics[i, j]
        LLR_tatol[j] = L_int[j] + r_sum
        # L_int or P
        LLR_ex[j] = r_sum
        
def _Iterative_Q(H, n_equations, n_v, M_dict, N_dict, L_int, Q_metrics, R_metrics):
    for i in range(n_equations):
        for j in range(n_v):
            if H[i][j] == 1:
                unit = 0
                for m in M_dict[j]:
                    if m == i:
                        continue
                    unit += R_metrics[m, j]
                Q_metrics[i, j] = L_int[j] + unit

def binaryproduct(X, Y):
    """Compute a matrix-matrix / vector product in Z/2Z."""
    A = X.dot(Y)
    try:
        A = A.toarray()
    except AttributeError:
        pass
    return A % 2

def transform_list(lst):
    return [1000 if x > 0 else -1000 for x in lst]
#########################################################################################################################
                
                
def log_LDPC_decode(H, G, L_int=[],round = 0):

    #L_int: External information transmitted from the line decoder

    n = H.shape[1]  # n
    k = n-H.shape[0]  # k

    len_code = n  # Original length of the code
    len_msg = int(len_code/n*k)  # Length of information bits

    L_metrics = zeros([len_code, 2])

    # Next is the flow of the sum-product algorithm

    n_equations = int(H.shape[0])

    n_v = int(H.shape[1])
    # print("n_equations,n_v:",n_equations,n_v)

    M_dict = {}  # Define M,as the number dimension of the calibration equation
    N_dict = {}  # Define N, the dimension of the number of code words

    Iter_nums = 0

    for i in range(n_v):
        M_dict[i] = []
    for i in range(n_equations):
        N_dict[i] = []

    Q_metrics = zeros([n_equations, n_v])

    R_metrics = zeros([n_equations, n_v])

    LLR_tatol = zeros(n)

    now_L_int = zeros(n)

    LLR_ex = zeros(n)

    fix_mode = 0

    # print("L_int1:",L_int)

    _compute_M_N_Q_ex(H, n_equations, n_v, M_dict,
                      N_dict, L_int, Q_metrics)  # Initialising Q messages

    #print("Q_metrics:", Q_metrics)
    _Iterative_R(H, n_equations, n_v, M_dict, N_dict,
                 L_int, Q_metrics, R_metrics)  # Calculating R information
    #print("R_metrics:", R_metrics)
    _compute_LLR(H, n_equations, n_v, M_dict, N_dict, L_int,
                 Q_metrics, R_metrics, LLR_tatol, LLR_ex)  # Calculate total_LLR
    # print("-------##############--")
    #print("LLR_tatol:", LLR_tatol)

    while True:
        #print(Iter_nums)
        if _Isitzero(H, LLR_tatol) != False:
            if round > 0 :
                LLR_ex = transform_list(LLR_tatol)
                LLR_tatol = transform_list(LLR_tatol)
            fix_mode = 1
            break

        if Iter_nums == 50:
            #print("stop_iteration_num = 10")
            break
        
        _Iterative_Q(H, n_equations, n_v, M_dict,
                     N_dict, L_int, Q_metrics, R_metrics)

        _Iterative_R(H, n_equations, n_v, M_dict,
                     N_dict, L_int, Q_metrics, R_metrics)

        _compute_LLR(H, n_equations, n_v, M_dict, N_dict, L_int,
                     Q_metrics, R_metrics, LLR_tatol, LLR_ex)
        Iter_nums += 1                
    # print("LLR_tatol:",LLR_tatol)
    fin_list=[]
    for i in range(n):
        if LLR_tatol[i] > 0 :
            fin_list.append(0)
        else:
            fin_list.append(1)
    #print("________________________________________")
    #print("L_int:",L_int)
    #print("tatol:",LLR_tatol)
    #print("ex:",LLR_ex)
    return [LLR_tatol, LLR_ex,fin_list,fix_mode]



#########################################################################################################################


def Column_Decode_LDPC(code_meritcs_list, H, G, multiple_nums, line_nums):
    """
    

    Parameters
    ----------
    code_meritcs_list : ndarray(block_index, row_index,[LLR值])
        从行解码获得的 LLR 值.
    H : matrix
        Calibration matrix.
    G : matrix
        Generate Matrix.
    multiple_nums : TYPE
        DESCRIPTION.
    line_nums : TYPE
        DESCRIPTION.

    Returns
    -------
    msg_meritcs : ndarray
        Decoded message.

    """
    
    len_code = G.shape[0]
    msg_L_meritcs = zeros([line_nums, len_code], dtype=float)
    msg_meritcs = zeros([line_nums, len_code], dtype=int)
    for code_list in code_meritcs_list:
        block_index = code_list[0]
        row_index = code_list[1]
        L_int = code_list[2]
        #print("L_int:", L_int)
        L_array = log_LDPC_decode(H, G, L_int)[0]
        msg_L_meritcs[row_index, :] = msg_L_meritcs[row_index, :] + L_array
        # print("msg_L_meritcs:",msg_L_meritcs)

    for i in range(msg_L_meritcs.shape[0]):
        for j in range(msg_L_meritcs.shape[1]):
            if msg_L_meritcs[i, j] >= 0:
                msg_meritcs[i, j] = 0
            else:
                msg_meritcs[i, j] = 1
    return msg_meritcs
