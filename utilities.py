
from unittest import result
from numpy import zeros
import numpy as np
import copy


#Decimalindex outputs a decimal 4
def int2quaternion(x):
    now_list=[]
    while x > 3:
        now_list.append(str(x%4))
        x = x // 4
    if x :
        now_list.append(str(x))
    x = "".join(reversed(now_list))
    return x 

#Quadratic to two binary,  0->00 1->01 2->10 3->11
def quaternion2binary_matrix(matrix) :
    binary_matrix_1 = copy.deepcopy(matrix)
    binary_matrix_2 = copy.deepcopy(matrix)
    binary_matrix_1[binary_matrix_1 <= 1] = 0
    binary_matrix_1[binary_matrix_1 > 1] = 1
    binary_matrix_2[binary_matrix_2 == 2] = 0
    binary_matrix_2[binary_matrix_2 > 0] = 1
    return [binary_matrix_1,binary_matrix_2]


#Aggregates two binary matrices into a single binary 4. The aggregation rules are 00->0,01->1,10->2,11->3
def binary2quaternion_matrix(binary_matrix_1,binary_matrix_2) :
    output_matrix = copy.deepcopy(binary_matrix_1)
    output_matrix[np.where((binary_matrix_1 ==1) & (binary_matrix_2 == 1))] = 3
    output_matrix[np.where((binary_matrix_1 ==1) & (binary_matrix_2 == 0))] = 2
    output_matrix[np.where((binary_matrix_1 ==0) & (binary_matrix_2 == 1))] = 1
    output_matrix[np.where((binary_matrix_1 ==0) & (binary_matrix_2 == 0))] = 0
    return output_matrix




#Decimal is directly converted to two binary bit arrays.
def get_b_list(x):
    arr_0 = []
    arr_1 = []
    if x == 0 :
        return [0],[0]
    else:
        while x > 0:
            arr_0.append(x%2)
            x = x >> 1
            arr_1.append(x%2)
            x = x >> 1
        return list(reversed(arr_1)),list(reversed(arr_0))


    

def base4todec(x):
    result = 0
    for i in x:
        result = result << 2 | i
    return result


def get_diff_nums(array1,array2):
    nums = 0
    for j in range(array1.shape[0]):
        if array1[j]*array2[j] <0 :
            nums +=1
    return nums


###########################################################################################################################
if __name__ == '__main__':
    a=99

    b,c = get_b_list(99)

    d =binary2quaternion_matrix(np.array(b),np.array(c))

