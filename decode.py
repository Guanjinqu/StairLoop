

import numpy as np
from numpy import zeros,array
from BCJR import log_map_decode
from LDPC import log_LDPC_decode
from utilities import quaternion2binary_matrix,binary2quaternion_matrix,int2quaternion,get_b_list,base4todec,get_diff_nums
import copy
import tqdm
import random
random.seed(1)
class Decode_unit :
    def __init__(

        self,
        G_matrix ,
        H_matrix,
        trellis,
        msg_name,
        pr_dict = {
        "ps":0.0001,
        "pd":0.01,
        "pi":0.01,
        "column":0},
        msg_length = 150,
        index_length = 20 ,
        block_num = 4,
        iterations =2,
        mpi_mode = False,
        output_mode = False,
        acc_mode = False,
        fast_mode = False
    ):
        self.G_matrix = G_matrix
        self.H_matrix = H_matrix
        self.trellis = trellis
        self.msg_name = msg_name
        self.pr_dict = pr_dict
        self.msg_length = msg_length
        self.index_length = index_length
        self.block_num = block_num
        self.iterations = iterations
        self.mpi_mode = mpi_mode
        self.output_mode = output_mode
        self.acc_mode = acc_mode
        self.line_nums = int(self.G_matrix.shape[0]/2)
        for i in range(150):
            if i**4 > self.G_matrix.shape[0] :
                self.line_index = i
                break
        self.block_index = self.index_length-self.line_index
        self.pr_dict["pt"]=((1-self.pr_dict["pi"])*(1-self.pr_dict["pd"]))
        self.total_msg = int(self.index_length+2*msg_length)
        self.F_dict = {}
        self.fast_mode = fast_mode

    def get_fin_msg(self,list_1,list_2,mode = 1):
        """
        Processes message blocks according to the specified pattern and generates the final message block list
        
        Parameters.
        - list_1: first list of messages, containing multiple blocks
        - list_2: the second message list, contains multiple message blocks
        - mode: processing mode, decide whether or not to do special processing for message blocks, default is 1
        
        Returns: self.fin_msg_block
        - self.fin_msg_block_list: the final list of blocks after processing.
        """
        if mode == 2 :
            for i in range(len(list_1)):
                list_1[i][list_1[i] >= 0] = 0
                list_1[i][list_1[i] < 0] = 1
                list_2[i][list_2[i] >= 0] = 0
                list_2[i][list_2[i] < 0] = 1
        self.fin_msg_block_list = []
        for i in range(len(list_1)):
            fin_msg_block = binary2quaternion_matrix(list_1[i],list_2[i])
            self.fin_msg_block_list.append(fin_msg_block)
        if self.output_mode == True :
            result_file_name = self.msg_name + "_decode_result"
            result_file = open(result_file_name,'w')
            result_file.write(str(self.fin_msg_block_list))
            result_file.close()
        return self.fin_msg_block_list


    def fast_decode(self) :
        #StairLoop Decoding
        output_list = [self.pr_dict]

        n_row = self.trellis.n
        k_row = self.trellis.k
        memory_row = self.trellis.total_memory
        
        v_column = self.G_matrix.shape[0]
        msg_column = self.G_matrix.shape[1]  
        
        msg_block_nums = int((self.block_num-1)/2)

        before_row_index = -1
        before_LLR = 0
        self.fin_msg_block_list_1 = []
        self.fin_msg_block_list_2 = []
        msg_block_L_list_1 = []
        msg_block_L_list_2 = []  
        for i in range(msg_block_nums):

            msg_block_L_list_1.append(zeros([v_column,self.msg_length] ))
            msg_block_L_list_2.append(zeros([v_column,self.msg_length] ))
            self.fin_msg_block_list_1.append(zeros([v_column,self.msg_length] ))
            self.fin_msg_block_list_2.append(zeros([v_column,self.msg_length] ))
        coding_dict = {}
        check_dict = {}
        for block in range(int(self.block_num/2)):
            for line in range(self.line_nums):
                coding_dict[(block,line)] = []
                check_dict[(block,line)] = 0
        input_name = self.msg_name + "_coding"
        input_file = open(input_name,"r")

        print("First search:")
        total_reads = 0
        error_reads = 0
        for line in tqdm.tqdm(input_file.readlines() ):
            total_reads += 1
            now_code = np.array(eval(line))
            now_code_list = quaternion2binary_matrix(now_code)
            now_code_1 = now_code_list[0]
            now_code_2 = now_code_list[1]
            #print(self.F_dict)
            result_bcjr = log_map_decode(now_code_1,now_code_2,self.total_msg,self.trellis,self.pr_dict,D=2,mode="ex",F_dict = self.F_dict)

            b_code_1,L_array_1,fin_L_array_1 = result_bcjr[0][0],result_bcjr[0][1],result_bcjr[0][2]
            b_code_2,L_array_2,fin_L_array_2 = result_bcjr[1][0],result_bcjr[1][1],result_bcjr[1][2]
            q_code = binary2quaternion_matrix(b_code_1,b_code_2)
            q_block_index = q_code[:self.block_index]
            q_line_index = q_code[self.block_index:self.line_index+self.block_index]
            now_block_index = base4todec(q_block_index)#Changing from binary 4 to binary 10
            now_line_index = base4todec(q_line_index) 

            if (now_block_index,now_line_index) not in coding_dict:
                error_reads +=1
                #print("error here:",now_block_index,now_line_index,coding_dict)
                continue
            coding_dict[(now_block_index,now_line_index)].append(now_code)
            if now_block_index  == 0 :
                msg_block_L_list_1[0][now_line_index,:] += L_array_1[self.msg_length+self.index_length:-memory_row]
                msg_block_L_list_2[0][now_line_index,:] += L_array_2[self.msg_length+self.index_length:-memory_row]
                self.fin_msg_block_list_1[0][now_line_index,:] += fin_L_array_1[self.msg_length+self.index_length:-memory_row]
                self.fin_msg_block_list_2[0][now_line_index,:] += fin_L_array_2[self.msg_length+self.index_length:-memory_row]
            elif now_block_index  == msg_block_nums:
                msg_block_L_list_1[now_block_index-1][int(v_column/2+now_line_index),:] += L_array_1[self.index_length:self.msg_length+self.index_length]
                msg_block_L_list_2[now_block_index-1][int(v_column/2+now_line_index),:] += L_array_2[self.index_length:self.msg_length+self.index_length]
                self.fin_msg_block_list_1[now_block_index-1][int(v_column/2+now_line_index),:] += fin_L_array_1[self.index_length:self.msg_length+self.index_length]
                self.fin_msg_block_list_2[now_block_index-1][int(v_column/2+now_line_index),:] += fin_L_array_2[self.index_length:self.msg_length+self.index_length]
            else:

                msg_block_L_list_1[now_block_index-1][int(v_column/2+now_line_index),:] += L_array_1[self.index_length:self.msg_length+self.index_length]
                msg_block_L_list_1[now_block_index][now_line_index,:] += L_array_1[self.msg_length+self.index_length:-memory_row]
                msg_block_L_list_2[now_block_index-1][int(v_column/2+now_line_index),:] += L_array_2[self.index_length:self.msg_length+self.index_length]
                msg_block_L_list_2[now_block_index][now_line_index,:] += L_array_2[self.msg_length+self.index_length:-memory_row]

                self.fin_msg_block_list_1[now_block_index-1][int(v_column/2+now_line_index),:] += fin_L_array_1[self.index_length:self.msg_length+self.index_length]
                self.fin_msg_block_list_1[now_block_index][now_line_index,:] += fin_L_array_1[self.msg_length+self.index_length:-memory_row]
                self.fin_msg_block_list_2[now_block_index-1][int(v_column/2+now_line_index),:] += fin_L_array_2[self.index_length:self.msg_length+self.index_length]
                self.fin_msg_block_list_2[now_block_index][now_line_index,:] += fin_L_array_2[self.msg_length+self.index_length:-memory_row]
        for key in coding_dict:

            reads_nums = len(coding_dict[key])
            if reads_nums == 0 :
                reads_nums = 1
            now_block_index = key[0]
            now_line_index = key[1]
            if  not  now_block_index  == 0 :
                msg_block_L_list_1[now_block_index-1][int(v_column/2+now_line_index),:] = msg_block_L_list_1[now_block_index-1][int(v_column/2+now_line_index),:]/reads_nums
                msg_block_L_list_2[now_block_index-1][int(v_column/2+now_line_index),:] = msg_block_L_list_2[now_block_index-1][int(v_column/2+now_line_index),:]/reads_nums
            if not now_block_index  == msg_block_nums:
                msg_block_L_list_1[now_block_index][now_line_index,:]  = msg_block_L_list_1[now_block_index][now_line_index,:] /reads_nums
                msg_block_L_list_2[now_block_index][now_line_index,:]  = msg_block_L_list_2[now_block_index][now_line_index,:]/reads_nums  



        #Next perform a column decoding：
        for msg_block_index in range(len(msg_block_L_list_1)):
            msg_block_1 =msg_block_L_list_1[msg_block_index]
            msg_block_2 =msg_block_L_list_2[msg_block_index]
            for j in range(msg_block_1.shape[1]):
                codewords_1 = [msg_block_1[:,j]]
                codewords_2 = [msg_block_2[:,j]]

                fin_column_L_array_1  = log_LDPC_decode(self.H_matrix,self.G_matrix,L_int=msg_block_1[:,j])[0] 
                fin_column_L_array_2 = log_LDPC_decode(self.H_matrix,self.G_matrix,L_int=msg_block_2[:,j])[0]

                column_L_array_1  = log_LDPC_decode(self.H_matrix,self.G_matrix,L_int=msg_block_1[:,j])[1] 
                column_L_array_2 = log_LDPC_decode(self.H_matrix,self.G_matrix,L_int=msg_block_2[:,j])[1] 

                msg_block_L_list_1[msg_block_index][:,j] = column_L_array_1
                msg_block_L_list_2[msg_block_index][:,j] = column_L_array_2
                self.fin_msg_block_list_1[msg_block_index][:,j] = fin_column_L_array_1
                self.fin_msg_block_list_2[msg_block_index][:,j] = fin_column_L_array_2
        print("Round 1 accuracy:")

        self.get_fin_msg(self.fin_msg_block_list_1,self.fin_msg_block_list_2,2)
        self.Acc()

        for nothing in range(self.iterations-1):

            ex_msg_list_1 = copy.deepcopy(msg_block_L_list_1)
            ex_msg_list_2 = copy.deepcopy(msg_block_L_list_2)
            self.fin_msg_block_list_1 = []
            self.fin_msg_block_list_2 = []
            msg_block_L_list_1 = []
            msg_block_L_list_2 = []  
            for i in range(msg_block_nums):
                msg_block_L_list_1.append(zeros([v_column,self.msg_length] ))
                msg_block_L_list_2.append(zeros([v_column,self.msg_length] ))
                self.fin_msg_block_list_1.append(zeros([v_column,self.msg_length] ))
                self.fin_msg_block_list_2.append(zeros([v_column,self.msg_length] ))                
            for key in tqdm.tqdm(coding_dict) :
                if coding_dict[key] == [] :
                    continue
                block_index = key[0]
                row_index = key[1]
                #Calculate the L value of index
                block_list_1,block_list_2 = get_b_list(block_index)
                ex_block_1 = zeros([int(self.block_index-len(block_list_1))])
                ex_block_2 = zeros([int(self.block_index-len(block_list_2))])
                block_index_1 = np.append(ex_block_1,np.array(block_list_1)).astype(int)
                block_index_2 = np.append(ex_block_2,np.array(block_list_2)).astype(int)
                
                line_list_1,line_list_2 = get_b_list(row_index)
                ex_line_1 = zeros([int(self.line_index-len(line_list_1))])
                ex_line_2 = zeros([int(self.line_index-len(line_list_2))])
                line_index_1 = np.append(ex_line_1,np.array(line_list_1)).astype(int)
                line_index_2 = np.append(ex_line_2,np.array(line_list_2)).astype(int)

                index_1 = np.append(block_index_1,line_index_1).astype(int)
                index_2 = np.append(block_index_2,line_index_2).astype(int)

                index_L_1 = zeros(self.index_length)
                index_L_2 = zeros(self.index_length)

                for i in range(self.index_length) :
                    if index_1[i] == 1 :
                        index_L_1[i] = - 1000
                    else:
                        index_L_1[i] = 1000
                    if index_2[i] == 1 :
                        index_L_2[i] = -1000
                    else:
                        index_L_2[i] = 1000



                if block_index == 0:
                    L_int_f_1 =zeros([1,self.msg_length])
                    L_int_f_2 =zeros([1,self.msg_length])
                    for i in range(self.msg_length):
                        L_int_f_1[0,i]= 1000
                        L_int_f_2[0,i]= 1000
                else:
                    L_int_f_1 = ex_msg_list_1[block_index-1][int(v_column/2+row_index),:]
                    L_int_f_2 = ex_msg_list_2[block_index-1][int(v_column/2+row_index),:]
                if block_index  == msg_block_nums:
                    L_int_b_1 = zeros([1,self.msg_length])
                    L_int_b_2 = zeros([1,self.msg_length])
                    for i in range(self.msg_length):
                        L_int_b_1[0,i] = 1000
                        L_int_b_2[0,i] = 1000
                else:
                    L_int_b_1 = ex_msg_list_1[block_index][row_index,:]
                    L_int_b_2 = ex_msg_list_2[block_index][row_index,:]

                L_int_1 = np.append(L_int_f_1,L_int_b_1)
                L_int_2 = np.append(L_int_f_2,L_int_b_2)
                L_int_1 = np.append(L_int_1,zeros([1,memory_row]))
                L_int_2 = np.append(L_int_2,zeros([1,memory_row]))
                L_int_1 = np.append(index_L_1,L_int_1)
                L_int_2 = np.append(index_L_2,L_int_2)
                #print(block_index,row_index,L_int_2)

                if self.fast_decode == True and check_dict[key] ==3:
                    
                    if block_index == 0 :
                        msg_block_L_list_1[0][row_index,:] += L_int_1[self.msg_length+self.index_length:-memory_row]
                        msg_block_L_list_2[0][row_index,:] += L_int_2[self.msg_length+self.index_length:-memory_row]
                    elif block_index  == msg_block_nums:
                        msg_block_L_list_1[block_index-1][int(v_column/2+row_index),:] += L_int_1[self.index_length:self.msg_length+self.index_length]
                        msg_block_L_list_2[block_index-1][int(v_column/2+row_index),:] += L_int_2[self.index_length:self.msg_length+self.index_length]
                    else:
                        msg_block_L_list_1[block_index-1][int(v_column/2+row_index),:] += L_int_1[self.index_length:self.msg_length+self.index_length]
                        msg_block_L_list_1[block_index][row_index,:] += L_int_1[self.msg_length+self.index_length:-memory_row]
                        msg_block_L_list_2[block_index-1][int(v_column/2+row_index),:] += L_int_2[self.index_length:self.msg_length+self.index_length]
                        msg_block_L_list_2[block_index][row_index,:] += L_int_2[self.msg_length+self.index_length:-memory_row]   
                           
                else:
                    reads_nums = len(coding_dict[key])
                    last_read = zeros(self.total_msg)
                    fin_row_list = [zeros(self.total_msg),zeros(self.total_msg)]
                    
                    now_check_index = 0 
                    for line in coding_dict[key]:
                        codewords_1,codewords_2 = quaternion2binary_matrix(line)
                        if line.shape[0] == last_read.shape[0]:
                            if not (line == last_read).all() :
                                result_bcjr = log_map_decode(codewords_1,codewords_2,self.total_msg,self.trellis,self.pr_dict,D=2,mode="ex",L_int_1=L_int_1,L_int_2=L_int_2,F_dict = self.F_dict)
                            
                        else:
                            result_bcjr = log_map_decode(codewords_1,codewords_2,self.total_msg,self.trellis,self.pr_dict,D=2,mode="ex",L_int_1=L_int_1,L_int_2=L_int_2,F_dict = self.F_dict)
                        last_read = line

                        fin_L_array_1,L_array_1 =result_bcjr[0][2],result_bcjr[0][1]
                        fin_L_array_2,L_array_2 =result_bcjr[1][2],result_bcjr[1][1]
                        if (get_diff_nums(L_int_1,fin_L_array_1)+get_diff_nums(L_int_2,fin_L_array_2)) > self.total_msg :
                            del coding_dict[key][now_check_index]
                            now_line_index -= 1
                            continue
                        fin_row_list[0] += fin_L_array_1[:-1]
                        fin_row_list[1] += fin_L_array_2[:-1]     
                        now_check_index +=1
                        #Next, assign the LLR to the L matrix and then decode the columns.
                        if block_index  == 0 :
                            msg_block_L_list_1[0][row_index,:] += L_array_1[self.msg_length+self.index_length:-memory_row]
                            msg_block_L_list_2[0][row_index,:] += L_array_2[self.msg_length+self.index_length:-memory_row]
                            self.fin_msg_block_list_1[0][row_index,:] += fin_L_array_1[self.msg_length+self.index_length:-memory_row]
                            self.fin_msg_block_list_2[0][row_index,:] += fin_L_array_2[self.msg_length+self.index_length:-memory_row]
                        elif block_index  == msg_block_nums:
                            msg_block_L_list_1[block_index-1][int(v_column/2+row_index),:] += L_array_1[self.index_length:self.msg_length+self.index_length]
                            msg_block_L_list_2[block_index-1][int(v_column/2+row_index),:] += L_array_2[self.index_length:self.msg_length+self.index_length]
                            self.fin_msg_block_list_1[block_index-1][int(v_column/2+row_index),:] += fin_L_array_1[self.index_length:self.msg_length+self.index_length]
                            self.fin_msg_block_list_2[block_index-1][int(v_column/2+row_index),:] += fin_L_array_2[self.index_length:self.msg_length+self.index_length]
                        else:

                            msg_block_L_list_1[block_index-1][int(v_column/2+row_index),:] += L_array_1[self.index_length:self.msg_length+self.index_length]
                            msg_block_L_list_1[block_index][row_index,:] += L_array_1[self.msg_length+self.index_length:-memory_row]
                            msg_block_L_list_2[block_index-1][int(v_column/2+row_index),:] += L_array_2[self.index_length:self.msg_length+self.index_length]
                            msg_block_L_list_2[block_index][row_index,:] += L_array_2[self.msg_length+self.index_length:-memory_row]

                            self.fin_msg_block_list_1[block_index-1][int(v_column/2+row_index),:] += fin_L_array_1[self.index_length:self.msg_length+self.index_length]
                            self.fin_msg_block_list_1[block_index][row_index,:] += fin_L_array_1[self.msg_length+self.index_length:-memory_row]
                            self.fin_msg_block_list_2[block_index-1][int(v_column/2+row_index),:] += fin_L_array_2[self.index_length:self.msg_length+self.index_length]
                            self.fin_msg_block_list_2[block_index][row_index,:] += fin_L_array_2[self.msg_length+self.index_length:-memory_row]


                    #Dealing with the problem of coefficient multiplication caused by multiple reads
                    if not block_index  == 0 :
                        msg_block_L_list_1[block_index-1][int(v_column/2+row_index),:] = msg_block_L_list_1[block_index-1][int(v_column/2+row_index),:]/reads_nums
                        msg_block_L_list_2[block_index-1][int(v_column/2+row_index),:] = msg_block_L_list_2[block_index-1][int(v_column/2+row_index),:]/reads_nums
                    if not block_index  == msg_block_nums:
                        msg_block_L_list_1[block_index][row_index,:] = msg_block_L_list_1[block_index][row_index,:]/reads_nums
                        msg_block_L_list_2[block_index][row_index,:] = msg_block_L_list_2[block_index][row_index,:]/reads_nums
                    if (get_diff_nums(fin_row_list[0],L_array_1)+get_diff_nums(fin_row_list[1],L_array_2)) == 0:
                        check_dict[key] += 1
                    else:
                        check_dict[key] == 0



            for msg_block_index in range(len(msg_block_L_list_1)):
                msg_block_1 =msg_block_L_list_1[msg_block_index]
                msg_block_2 =msg_block_L_list_2[msg_block_index]

                for j in range(msg_block_1.shape[1]):
                    codewords_1 = [msg_block_1[:,j]]
                    codewords_2 = [msg_block_2[:,j]]
                    

                    fin_column_L_array_1  = log_LDPC_decode(self.H_matrix,self.G_matrix,L_int=msg_block_1[:,j])[0] 
                    fin_column_L_array_2 = log_LDPC_decode(self.H_matrix,self.G_matrix,L_int=msg_block_2[:,j])[0]
                    column_L_array_1  = log_LDPC_decode(self.H_matrix,self.G_matrix,L_int=msg_block_1[:,j])[1] 
                    column_L_array_2 = log_LDPC_decode(self.H_matrix,self.G_matrix,L_int=msg_block_2[:,j])[1] 

                    msg_block_L_list_1[msg_block_index][:,j] = column_L_array_1
                    msg_block_L_list_2[msg_block_index][:,j] = column_L_array_2
                    self.fin_msg_block_list_1[msg_block_index][:,j] = fin_column_L_array_1
                    self.fin_msg_block_list_2[msg_block_index][:,j] = fin_column_L_array_2

            print("Accuracy of round %d :"%(nothing+2))
            self.get_fin_msg(self.fin_msg_block_list_1,self.fin_msg_block_list_2,2)
            result = self.Acc()
            output_list.append(result)


    def Acc(self) :
        """
        Calculate the decoding accuracy.

        This method is mainly used to calculate the accuracy between the received message and the final message block.
        It reads the message matrix stored in the file and compares it to the final message block stored internally
        It reads the message matrix stored in the file and compares it to the internally stored final message block to determine how many bits of the message delivery were accurate.

        Returns.
            float: Ratio of accuracy, the number of error-free bits divided by the total number of bits.
        """
        input_name = self.msg_name + "_msg_matrix"
        input_file = open(input_name,"r").read()
        msg_list = eval(input_file)[1:-1]
        total_bit = 0
        error_bit = 0

        for k in range(len(msg_list)):
            msg_block = msg_list[k]
            fin_block = self.fin_msg_block_list[k]

            for i in range(msg_block.shape[0]):
                for j in range(msg_block.shape[1]):
                    total_bit +=1
                    if msg_block[i,j] != fin_block[i,j]:
                        #print("Error location：",k,i,j,msg_block[i,j],fin_block[i,j])
                        error_bit += 1
        print("Accuracy:",(total_bit-error_bit)/total_bit)
        return (total_bit-error_bit)/total_bit



