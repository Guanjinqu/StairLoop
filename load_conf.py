#Load encoder parameters

from numpy import array
import commpy.channelcoding.convcode as cc

from LDPC import coding_matrix_systematic,regular_coding_matrix_systematic

import random
random.seed(1)
import yaml

def generate_conv(index):
    #Determine the convolutional code structure used based on INDEX
    if index == 0:
        #Constructing convolutional codes with code rate 0.5
        memory1 = array([1])                     #Number of registers

        g_matrix1 = array([[1,2]])               #生成矩阵的样子
    
        fd1 = array([[3]])
        trellis = cc.Trellis(memory1,g_matrix1,feedback=fd1)
        trellis.output_table[0][0] = 1
        trellis.output_table[0][1] = 3
        trellis.output_table[1][0] = 0
        trellis.output_table[1][1] = 2
        trellis.next_state_table[0][0] = 0
        trellis.next_state_table[0][1] = 1
        trellis.next_state_table[1][0] = 0
        trellis.next_state_table[1][1] = 1

    if index == 1 :
        #Constructing codes with code rate 1/3 and GC-balance
        rate =1/3
        memory = array([1])                    
        g_matrix = array([[1,3,2]])               
        fd = array([[3]])
        trellis = cc.Trellis(memory,g_matrix,feedback=fd,code_type='rsc')
        #Modifications based on GC-balance
        trellis.output_table[0][0] = 1
        trellis.output_table[0][1] = 6
        trellis.output_table[1][0] = 2
        trellis.output_table[1][1] = 5
        trellis.next_state_table[0][0] = 1
        trellis.next_state_table[0][1] = 0
        trellis.next_state_table[1][0] = 0
        trellis.next_state_table[1][1] = 1

    return trellis

def generate_ldpc(index):
    #Determine the adopted LDPC code structure based on INDEX
    if index == 0:
        K = 82
        N = 100
        var_degrees, check_degrees = {3: 0.4052, 4: 0.3927,
                                    8: 0.1466, 9: 0.0555}, {23: 0.3109, 24: 0.6891}  # irregular_degree_sequence
        H_matrix, G_matrix = coding_matrix_systematic(K, N, var_degrees, check_degrees)

    if index == 1 :
        K = 64
        N = 100
        var_degrees, check_degrees = {2:0.4329, 3:0.1583, 
                                5:0.4088},{8:1} # irregular_degree_sequence
        H_matrix, G_matrix = coding_matrix_systematic(K, N, var_degrees, check_degrees) 

    if index == 2 :
        H_matrix, G_matrix = regular_coding_matrix_systematic(30,3,6)


    if index == 3 :
        
        K = 210
        N = 250
        var_degrees, check_degrees = {3: 0.4052, 4: 0.3927,8: 0.1466, 9: 0.0555}, {23: 0.3109, 24: 0.6891}
        H_matrix, G_matrix = coding_matrix_systematic(K, N, var_degrees, check_degrees) 

    return H_matrix, G_matrix


str_plot = """
  ______    __                __            __                                     
 /      \  /  |              /  |          /  |                                    
/$$$$$$  |_$$ |_     ______  $$/   ______  $$ |        ______    ______    ______  
$$ \__$$// $$   |   /      \ /  | /      \ $$ |       /      \  /      \  /      \ 
$$      \$$$$$$/    $$$$$$  |$$ |/$$$$$$  |$$ |      /$$$$$$  |/$$$$$$  |/$$$$$$  |
 $$$$$$  | $$ | __  /    $$ |$$ |$$ |  $$/ $$ |      $$ |  $$ |$$ |  $$ |$$ |  $$ |
/  \__$$ | $$ |/  |/$$$$$$$ |$$ |$$ |      $$ |_____ $$ \__$$ |$$ \__$$ |$$ |__$$ |
$$    $$/  $$  $$/ $$    $$ |$$ |$$ |      $$       |$$    $$/ $$    $$/ $$    $$/ 
 $$$$$$/    $$$$/   $$$$$$$/ $$/ $$/       $$$$$$$$/  $$$$$$/   $$$$$$/  $$$$$$$/  
                                                                         $$ |      
                                                                         $$ |      
                                                                         $$/       
                          Welcome to StairLoop!                                                                                                                                                                                            

"""

def load_config():
    with open('config.yaml',encoding='utf-8') as file1:
        data = yaml.load(file1,Loader=yaml.FullLoader)
        data["trellis"] = generate_conv(data["Conv_mode"])
        data['H_matrix'],data['G_matrix'] = generate_ldpc(data["LDPC_mode"])
    #print(str_plot)

    return data

def get_label():
    print(str_plot)

if __name__ == '__main__':
    get_label()