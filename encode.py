#encoding

from pickletools import uint1
import random
import commpy.channelcoding.convcode as cc
import numpy as np
from numpy import zeros,array
from utilities import quaternion2binary_matrix,binary2quaternion_matrix,int2quaternion,get_b_list
from LDPC import coding_matrix_systematic
from channel import channel_model_unit
import sys
import yaml
from ref2bin import file2ref
from load_conf import get_label
np.set_printoptions(threshold=sys.maxsize)
random.seed(1)

def binaryproduct(X, Y):
    """
    Computes the product of two matrices X and Y and returns the result of modelling the product on 2.
    
    Parameters.
    X: the first matrix.
    Y: the second matrix.
    
    Return value.
    The product matrix modulo 2.
    """
    A = X.dot(Y)
    try:
        A = A.toarray()
    except AttributeError:
        pass
    return A % 2

class Encode_unit :
    def __init__(

        self,
        G_matrix ,
        H_matrix,
        trellis,
        unit_name,
        msg_length = 150,
        index_length = 20 ,
        block_num = 4,
        input_file_name = None,
        input_mode = False
    ):
        self.G_matrix = G_matrix
        self.H_matrix = H_matrix
        self.trellis = trellis
        self.unit_name = unit_name
        self.msg_length = msg_length
        self.index_length = index_length
        self.block_num = block_num
        self.input_file_name = input_file_name
        self.input_mode = input_mode

        for i in range(150):
            if i**4 > self.G_matrix.shape[0] :
                self.line_index = i
                break
        self.block_index = self.index_length-self.line_index

    def generate_msg_matrix(self) :
        #Constructing the information matrix
        
        get_label()
        
        self.msg_metrics_list = [np.random.randint(0,1,size = [int(self.G_matrix.shape[0]/2),self.msg_length])]

        if self.input_mode == False:
            for i in range(int(self.block_num/2-1)):
                self.msg_metrics_list.append(np.random.randint(0,4,size = [self.G_matrix.shape[1],self.msg_length]))
        else:
            file2ref(self.input_file_name,"raw_data")
            unit_dict = {"A":0,"T":1,"G":2,"C":3}
            input_file = open("raw_data","r").read()
            lenth_file = len(input_file)
            with open('config.yaml', 'r') as file:
                config = yaml.safe_load(file)
            config['raw_length'] = lenth_file  
            with open('config.yaml', 'w') as file:
                yaml.dump(config, file)
            print("info_nums:",lenth_file)
            print("blocks_np_nums:",int(self.block_num/2-1)*self.G_matrix.shape[1]*self.msg_length)
            for i in range(int(self.block_num/2-1)):
                self.msg_metrics_list.append(np.random.randint(0,1,size = [self.G_matrix.shape[1],self.msg_length]))  
            for i in range(int(self.block_num/2-1)):
                for j in range(self.G_matrix.shape[1]):
                    for k in range(self.msg_length):
                        index = i*self.G_matrix.shape[1]*self.msg_length+j*self.msg_length+k
                        if index < lenth_file:
                            self.msg_metrics_list[i+1][j][k] = unit_dict[input_file[index]]

        self.msg_metrics_list.append(np.random.randint(0,1,size = [int(self.G_matrix.shape[0]/2),self.msg_length]))
        output_list = []
        for i in self.msg_metrics_list:
            output_list.append(i.tolist())
        output_name = self.unit_name + "_msg_matrix"
        output_file = open(output_name,"w")
        output_file.write(str(self.msg_metrics_list))
        output_file.close()
        print("The information matrix is ​​constructed")
        
    

    def generate_encode_matrix(self) :
        #Constructing the encoded matrix

        self.binary_msg_matrix_list_1 = []
        self.binary_msg_matrix_list_2 = []
        for i in range(len(self.msg_metrics_list)):
            msg_metrics = self.msg_metrics_list[i]
            binary_matrix_list = quaternion2binary_matrix(msg_metrics)
            self.binary_msg_matrix_list_1.append(binary_matrix_list[0])
            self.binary_msg_matrix_list_2.append(binary_matrix_list[1])
    
        #Column Encoding
        column_metrics_list_1 = [np.random.randint(0,1,size = [int(self.G_matrix.shape[0]/2),self.msg_length])]
        column_metrics_list_2 = [np.random.randint(0,1,size = [int(self.G_matrix.shape[0]/2),self.msg_length])]
        for i in range(int(self.block_num/2-1)) :
            msg_m_1 = self.binary_msg_matrix_list_1[i+1]
            msg_m_2 = self.binary_msg_matrix_list_2[i+1]
            column_m_1 = zeros([self.G_matrix.shape[0],self.msg_length],dtype= int)
            column_m_2 = zeros([self.G_matrix.shape[0],self.msg_length],dtype= int)
            for j in range(self.msg_length):
                column_m_1[:,j] = binaryproduct(self.G_matrix,msg_m_1[:,j])
                column_m_2[:,j] = binaryproduct(self.G_matrix,msg_m_2[:,j])
            column_metrics_list_1.append(column_m_1)    
            column_metrics_list_2.append(column_m_2) 
        column_metrics_list_1.append(np.random.randint(0,1,size = [int(self.G_matrix.shape[0]/2),self.msg_length]))
        column_metrics_list_2.append(np.random.randint(0,1,size = [int(self.G_matrix.shape[0]/2),self.msg_length]))

        #Row Encoding
        code_metrics_list_1 = []
        code_metrics_list_2 = []
        bf_nums = 0
        at_nums = 1
        now_block_index = 0
        now_line_index = 0
        for i in range(int(self.block_num/2)):
            bf_metrics_1 = column_metrics_list_1[bf_nums]
            bf_metrics_2 = column_metrics_list_2[bf_nums]
            at_metrics_1 = column_metrics_list_1[at_nums]
            at_metrics_2 = column_metrics_list_2[at_nums]
            bf_nums += 1
            at_nums += 1
            code_metrics_1 = zeros([int(self.G_matrix.shape[0]/2),int((2*self.msg_length+self.trellis.total_memory+self.index_length)/self.trellis.k*self.trellis.n)],dtype= int)
            code_metrics_2 = zeros([int(self.G_matrix.shape[0]/2),int((2*self.msg_length+self.trellis.total_memory+self.index_length)/self.trellis.k*self.trellis.n)],dtype= int)
            
            now_line_index = 0
            for j in range(int(self.G_matrix.shape[0]/2)):

                if i == 0:
                    bf_msg_1 = zeros(bf_metrics_1.shape[1])
                    bf_msg_2 = zeros(bf_metrics_2.shape[1])
                else:
                    bf_msg_1 = bf_metrics_1[int(self.G_matrix.shape[0]/2)+j]
                    bf_msg_2 = bf_metrics_2[int(self.G_matrix.shape[0]/2)+j]

                at_msg_1 = at_metrics_1[j]
                at_msg_2 = at_metrics_2[j]
                block_list_1,block_list_2 = get_b_list(now_block_index)
                ex_block_1 = zeros([int(self.block_index-len(block_list_1))])
                ex_block_2 = zeros([int(self.block_index-len(block_list_2))])
                block_index_1 = np.append(ex_block_1,np.array(block_list_1)).astype(int)
                block_index_2 = np.append(ex_block_2,np.array(block_list_2)).astype(int)
                
                line_list_1,line_list_2 = get_b_list(now_line_index)
                ex_line_1 = zeros([int(self.line_index-len(line_list_1))])
                ex_line_2 = zeros([int(self.line_index-len(line_list_2))])
                line_index_1 = np.append(ex_line_1,np.array(line_list_1)).astype(int)
                line_index_2 = np.append(ex_line_2,np.array(line_list_2)).astype(int)

                index_1 = np.append(block_index_1,line_index_1).astype(int)
                index_2 = np.append(block_index_2,line_index_2).astype(int)

                now_line_index += 1
                                                                                                                              

                msg_now_1 = np.append(bf_msg_1,at_msg_1).astype(int)
                msg_now_2 = np.append(bf_msg_2,at_msg_2).astype(int)

                fin_msg_1 = np.append(index_1,msg_now_1).astype(int)
                fin_msg_2 = np.append(index_2,msg_now_2).astype(int)

                code_now_1 = cc.conv_encode(fin_msg_1,self.trellis, termination='term')
                code_now_2 = cc.conv_encode(fin_msg_2,self.trellis, termination='term')
                for k in range(len(code_now_1)):
                    code_metrics_1[j,k] = code_now_1[k]
                    code_metrics_2[j,k] = code_now_2[k]
            code_metrics_list_1.append(code_metrics_1)
            code_metrics_list_2.append(code_metrics_2)
            now_block_index += 1
        self.quaternion_code_metrics_list = []
        
        for i in range(len(code_metrics_list_1)):
            binary_matrix_1 = code_metrics_list_1[i]
            binary_matrix_2 = code_metrics_list_2[i]
            self.quaternion_code_metrics_list.append(binary2quaternion_matrix(binary_matrix_1,binary_matrix_2))
        


        output_info = ""
        output_dict = {0:"A",1:"T",2:"G",3:"C"}
        output_name = self.unit_name + "_encode_info"
        output_file = open(output_name,"w")
        for block in self.quaternion_code_metrics_list:
            for line in block:
                output_info = ""
                for unit in line:
                    output_info += output_dict[unit]
                output_file.write(output_info)
                output_file.write('\n')


        output_file.close()     
        print("Channel coding completed")
        return self.quaternion_code_metrics_list
        
    def generate_channel_coding(self,extra_nums = 1,pr_dict= {
        "ps":0.0001,
        "pd":0.0001,
        "pi":0.0001,
        "column":0}):
        """
        Modelling channel errors

        This function simulates the DNA storage channel based on the error rate versus the sequencing depth and writes the results to a file

        Parameters.
        extra_nums (int): sequencing depth, default is 1
        pr_dict (dict): dictionary of channel parameters, contains the keys ps, pd, pi and column, default value is 0.0001

        Returns.
        No return value, but writes the encoding result to a file.
        """
        file_name = self.unit_name+"_coding"
        output_file = open(file_name,'w')
        ex_coding_nums = int(self.G_matrix.shape[0]/2*len(self.quaternion_code_metrics_list))
        self.quaternion_code = zeros([ex_coding_nums,int((2*self.msg_length+self.trellis.total_memory+self.index_length)/self.trellis.k*self.trellis.n)],dtype= int)
        now_index = 0
        for block in self.quaternion_code_metrics_list :
            for line in range(int(self.G_matrix.shape[0]/2)):
                self.quaternion_code[now_index,:] = block[line,:]
                now_index += 1
        for i in range(int(ex_coding_nums*extra_nums)):
            coding_index = random.randint(0,ex_coding_nums-1)
            now_coding = self.quaternion_code[coding_index,:]
            after_coding = channel_model_unit(now_coding,pr_dict).tolist()
            output_file.write(str(after_coding))
            output_file.write('\n')

    def generate_mpi_channel_coding(self,extra_nums = 1,mpi_nums = 2,pr_dict= {
        "ps":0.0001,
        "pd":0.0001,
        "pi":0.0001,
        "column":0}):
        """
        Modelling channel errors for MPI

        This function simulates the DNA storage channel based on the error rate versus the sequencing depth and writes the results to multiple file files

        Parameters.
        extra_nums (int): sequencing depth, default is 1
        mpi_nums (int): number of MPI processes.
        pr_dict (dict): dictionary of channel parameters, contains the keys ps, pd, pi and column, default value is 0.0001
        

        Returns.
        No return value, but writes the encoding result to a file.
        """
        ex_coding_nums = int(self.G_matrix.shape[0]/2*len(self.quaternion_code_metrics_list))
        unit_reads_nums = int(ex_coding_nums*extra_nums/mpi_nums)
        ex_name = self.unit_name+"_mpi_coding_"

        self.quaternion_code = zeros([ex_coding_nums,int((2*self.msg_length+self.trellis.total_memory+self.index_length)/self.trellis.k*self.trellis.n)],dtype= int)
        now_index = 0
        for block in self.quaternion_code_metrics_list :
            for line in range(int(self.G_matrix.shape[0]/2)):
                self.quaternion_code[now_index,:] = block[line,:]
                now_index += 1

        now_rxtra_nums = 0        
        for i in range(mpi_nums):
            file_name = ex_name +str(i)

            output_file = open(file_name,'w')
            for i in range(int(unit_reads_nums)):
                coding_index = random.randint(0,ex_coding_nums-1)
                now_coding = self.quaternion_code[coding_index,:]
                after_coding = channel_model_unit(now_coding,pr_dict).tolist()
                output_file.write(str(after_coding))
                output_file.write('\n')
            output_file.close()
        print("Simulated channel completed")

