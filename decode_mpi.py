#mpi decoding


import numpy as np
from load_conf import load_config,get_label
from numpy import zeros,array
from BCJR import log_map_decode
from LDPC import log_LDPC_decode
from utilities import quaternion2binary_matrix,binary2quaternion_matrix,int2quaternion,get_b_list,base4todec,get_diff_nums
import copy
import tqdm
import random
from ref2bin import ref2file
from mpi4py import MPI
random.seed(1)
class Decode_mpi_unit :
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
        raw_length = 0,
        mpi_mode = False,
        output_mode = False,
        acc_mode = False,
        fast_mode = False,
        input_name = None,
        output_name = None

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
        self.input_name = input_name
        self.raw_length = raw_length
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.mpi_nums = self.comm.Get_size()
        self.mpi_block_nums = int(self.block_num/self.mpi_nums/2) 
        for i in range(150):
            if i**4 > self.G_matrix.shape[0] :
                self.line_index = i
                break
        self.block_index = self.index_length-self.line_index
        self.pr_dict["pt"]=((1-self.pr_dict["pi"])*(1-self.pr_dict["pd"]))
        self.total_msg = int(self.index_length+2*msg_length)
        self.F_dict = {}
        self.fast_mode = fast_mode
        self.output_name = output_name
        if self.rank == 0 :
            get_label()

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
        fin_coding_dict = {}
        check_dict = {}
        for block in range(int(self.block_num/2)):
            for line in range(self.line_nums):
                coding_dict[(block,line)] = []
                check_dict[(block,line)] = 0
        input_name = self.msg_name + "_mpi_coding_" + str(self.rank)
        input_file = open(input_name,"r")
        ##########################################################first row decoding

        total_reads = 0
        error_reads = 0
        for line in tqdm.tqdm(input_file.readlines() ):
            total_reads += 1
            now_code = np.array(eval(line))
            now_code_list = quaternion2binary_matrix(now_code)
            now_code_1 = now_code_list[0]
            now_code_2 = now_code_list[1]

            result_bcjr = log_map_decode(now_code_1,now_code_2,self.total_msg,self.trellis,self.pr_dict,D=2,mode="ex",F_dict=self.F_dict)

            b_code_1,L_array_1,fin_L_array_1 = result_bcjr[0][0],result_bcjr[0][1],result_bcjr[0][2]
            b_code_2,L_array_2,fin_L_array_2 = result_bcjr[1][0],result_bcjr[1][1],result_bcjr[1][2]
            q_code = binary2quaternion_matrix(b_code_1,b_code_2)
            q_block_index = q_code[:self.block_index]
            q_line_index = q_code[self.block_index:self.line_index+self.block_index]
            now_block_index = base4todec(q_block_index)#Changing from binary 4 to binary 10
            now_line_index = base4todec(q_line_index) 

            if (now_block_index,now_line_index) not in coding_dict:
                error_reads +=1

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
 

        msg_block_L_list_1 = [0] + msg_block_L_list_1 + [0]
        msg_block_L_list_2 = [0] + msg_block_L_list_2 + [0]

        #MPI send


        BUFSISE = 1000000000+ MPI.BSEND_OVERHEAD
        buf = bytearray(BUFSISE)


        MPI.Attach_buffer(buf)

        for i in range(self.mpi_nums):
            now_dict = {}
            now_L_list = []

            for key in coding_dict :
                if key[0] >= int(i*self.mpi_block_nums) and key[0] < int((i+1)*self.mpi_block_nums)  :
                    now_dict[key] = coding_dict[key]
            
                
            now_L_list = [msg_block_L_list_1[int(i*self.mpi_block_nums):int((i+1)*self.mpi_block_nums+1)],msg_block_L_list_2[int(i*self.mpi_block_nums):int((i+1)*self.mpi_block_nums+1)]]
            if i == self.rank :
                fin_coding_dict = now_dict
            else:
                send_obj = [now_L_list,now_dict]
                send_reads = self.comm.ibsend(send_obj,dest = i,tag = 0)
                send_reads.wait()
        

        msg_block_L_list_1 = msg_block_L_list_1[int(self.rank*self.mpi_block_nums):int((self.rank+1)*self.mpi_block_nums+1)]
        msg_block_L_list_2 = msg_block_L_list_2[int(self.rank*self.mpi_block_nums):int((self.rank+1)*self.mpi_block_nums+1)]
        

        for i in range(self.mpi_nums):
            if i == self.rank :
                continue
            else:
                recv_req = self.comm.irecv(buf =  1000000000,source = i ,tag = 0)
                recv_obj = recv_req.wait()
                input_dict = recv_obj[1]
                input_msg_list = recv_obj[0]
                for key in input_dict:
                    fin_coding_dict[key] = fin_coding_dict[key] + input_dict[key]
                for i in range(len(msg_block_L_list_1)):
                    if type(msg_block_L_list_1[i]) == int :
                        continue
                    msg_block_L_list_1[i] += input_msg_list[0][i]
                    msg_block_L_list_2[i] += input_msg_list[1][i]

        
        coding_dict = fin_coding_dict
        for key in fin_coding_dict:
            reads_nums = len(coding_dict[key])
            if reads_nums == 0 :
                reads_nums = 1
            now_block_index = key[0]
            now_line_index = key[1]

            if  not  now_block_index  == 0 :
                msg_block_L_list_1[now_block_index-int(self.rank*self.mpi_block_nums)][int(v_column/2+now_line_index),:] = msg_block_L_list_1[now_block_index-int(self.rank*self.mpi_block_nums)][int(v_column/2+now_line_index),:]/reads_nums
                msg_block_L_list_2[now_block_index-int(self.rank*self.mpi_block_nums)][int(v_column/2+now_line_index),:] = msg_block_L_list_2[now_block_index-int(self.rank*self.mpi_block_nums)][int(v_column/2+now_line_index),:]/reads_nums
            if not now_block_index  == msg_block_nums:
                msg_block_L_list_1[now_block_index-int(self.rank*self.mpi_block_nums)+1][now_line_index,:]  = msg_block_L_list_1[now_block_index-int(self.rank*self.mpi_block_nums)+1][now_line_index,:]/reads_nums
                msg_block_L_list_2[now_block_index-int(self.rank*self.mpi_block_nums)+1][now_line_index,:]  = msg_block_L_list_2[now_block_index-int(self.rank*self.mpi_block_nums)+1][now_line_index,:]/reads_nums 

        send_1 = zeros([int(v_column/2),self.msg_length])
        recv_1 = zeros([int(v_column/2),self.msg_length])
        send_2 = zeros([int(v_column/2),self.msg_length])
        recv_2 = zeros([int(v_column/2),self.msg_length])

        #comlun decoding
        for msg_block_index in range(len(msg_block_L_list_1)-1):
            if msg_block_index == 0 and self.rank == 0 :
                continue
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

        if not self.rank == 0 :
            send_1 = msg_block_L_list_1[0][:int(v_column/2),:]
            send_2 = msg_block_L_list_2[0][:int(v_column/2),:]
            send_req_1 = self.comm.Isend(send_1,dest = self.rank - 1,tag = 1)
            send_req_1.Wait()
            send_req_2 = self.comm.Isend(send_2,dest = self.rank - 1,tag = 2)
            send_req_2.Wait()
        if not self.rank == self.mpi_nums -1 :
            recv_req_1 = self.comm.Irecv(recv_1, source=self.rank + 1, tag=1)
            recv_req_1.Wait()
            recv_req_2 = self.comm.Irecv(recv_2, source=self.rank + 1, tag=2)
            recv_req_2.Wait()

            msg_block_L_list_1[-1][:int(v_column/2),:] = recv_1
            msg_block_L_list_2[-1][:int(v_column/2),:] = recv_2

 





        
        output_list = [self.pr_dict]
        for nothing in range(self.iterations-1):

            ex_msg_list_1 = copy.deepcopy(msg_block_L_list_1)
            ex_msg_list_2 = copy.deepcopy(msg_block_L_list_2)
            self.fin_msg_block_list_1 = []
            self.fin_msg_block_list_2 = []
            msg_block_L_list_1 = []
            msg_block_L_list_2 = []  

            for i in range(self.mpi_block_nums+1):
                msg_block_L_list_1.append(zeros([v_column,self.msg_length] ))
                msg_block_L_list_2.append(zeros([v_column,self.msg_length] ))
                self.fin_msg_block_list_1.append(zeros([v_column,self.msg_length] ))
                self.fin_msg_block_list_2.append(zeros([v_column,self.msg_length] ))                
            for key in sorted(coding_dict.keys(),reverse=True) :
                
                block_index = key[0]
                row_index = key[1]
                
                if int(block_index-self.rank*self.mpi_block_nums) == self.mpi_block_nums-2 :
                    if not self.rank == self.mpi_nums-1:
                        send_1 = msg_block_L_list_1[-1][:int(v_column/2),::]
                        send_2 = msg_block_L_list_2[-1][:int(v_column/2),::]
                        send_req_1 = self.comm.Isend(send_1,dest = self.rank +1,tag = 1)
                        send_req_1.Wait()
                        send_req_2 = self.comm.Isend(send_2,dest = self.rank +1,tag = 2)
                        send_req_2.Wait()
                    if not self.rank == 0 :
                        recv_req_1 = self.comm.Irecv(recv_1, source=self.rank - 1, tag=1)
                        recv_req_1.Wait()
                        recv_req_2 = self.comm.Irecv(recv_2, source=self.rank - 1, tag=2)
                        recv_req_2.Wait()
                        msg_block_L_list_1[0][:int(v_column/2),:] = recv_1
                        msg_block_L_list_2[0][:int(v_column/2),:] = recv_2
                if coding_dict[key] == [] :
                    continue
                #get INDEX L-value
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
                    L_int_f_1 = ex_msg_list_1[int(block_index-self.rank*self.mpi_block_nums)][int(v_column/2+row_index),:]
                    L_int_f_2 = ex_msg_list_2[int(block_index-self.rank*self.mpi_block_nums)][int(v_column/2+row_index),:]
                if block_index  == msg_block_nums:
                    L_int_b_1 = zeros([1,self.msg_length])
                    L_int_b_2 = zeros([1,self.msg_length])
                    for i in range(self.msg_length):
                        L_int_b_1[0,i] = 1000
                        L_int_b_2[0,i] = 1000
                else:
                    L_int_b_1 = ex_msg_list_1[int(block_index-self.rank*self.mpi_block_nums)+1][row_index,:]
                    L_int_b_2 = ex_msg_list_2[int(block_index-self.rank*self.mpi_block_nums)+1][row_index,:]

                L_int_1 = np.append(L_int_f_1,L_int_b_1)
                L_int_2 = np.append(L_int_f_2,L_int_b_2)
                L_int_1 = np.append(L_int_1,zeros([1,memory_row]))
                L_int_2 = np.append(L_int_2,zeros([1,memory_row]))
                L_int_1 = np.append(index_L_1,L_int_1)
                L_int_2 = np.append(index_L_2,L_int_2)

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

                        result_bcjr = log_map_decode(codewords_1,codewords_2,self.total_msg,self.trellis,self.pr_dict,D=2,mode="ex",L_int_1=L_int_1,L_int_2=L_int_2)

                        last_read = line 
                        fin_row_list[0] += fin_L_array_1[:-1]
                        fin_row_list[1] += fin_L_array_2[:-1]  
                        now_check_index +=1

                        fin_L_array_1,L_array_1 =result_bcjr[0][2],result_bcjr[0][1]
                        fin_L_array_2,L_array_2 =result_bcjr[1][2],result_bcjr[1][1]

                        if not block_index  == msg_block_nums :
                            msg_block_L_list_1[int(block_index-self.rank*self.mpi_block_nums)+1][row_index,:] += L_array_1[self.msg_length+self.index_length:-memory_row]
                            msg_block_L_list_2[int(block_index-self.rank*self.mpi_block_nums)+1][row_index,:] += L_array_2[self.msg_length+self.index_length:-memory_row]

                        if  not block_index  == 0:
                            msg_block_L_list_1[int(block_index-self.rank*self.mpi_block_nums)][int(v_column/2+row_index),:] += L_array_1[self.index_length:self.msg_length+self.index_length]
                            msg_block_L_list_2[int(block_index-self.rank*self.mpi_block_nums)][int(v_column/2+row_index),:] += L_array_2[self.index_length:self.msg_length+self.index_length]

                    if not block_index  == 0 :
                        msg_block_L_list_1[int(block_index-self.rank*self.mpi_block_nums)][int(v_column/2+row_index),:] = msg_block_L_list_1[int(block_index-self.rank*self.mpi_block_nums)][int(v_column/2+row_index),:]/reads_nums
                        msg_block_L_list_2[int(block_index-self.rank*self.mpi_block_nums)][int(v_column/2+row_index),:] = msg_block_L_list_2[int(block_index-self.rank*self.mpi_block_nums)][int(v_column/2+row_index),:]/reads_nums
                    if not block_index  == msg_block_nums:
                        msg_block_L_list_1[int(block_index-self.rank*self.mpi_block_nums)+1][row_index,:] = msg_block_L_list_1[int(block_index-self.rank*self.mpi_block_nums)+1][row_index,:]/reads_nums
                        msg_block_L_list_2[int(block_index-self.rank*self.mpi_block_nums)+1][row_index,:] = msg_block_L_list_2[int(block_index-self.rank*self.mpi_block_nums)+1][row_index,:]/reads_nums
                    if (get_diff_nums(fin_row_list[0],L_array_1)+get_diff_nums(fin_row_list[1],L_array_2)) == 0:
                        check_dict[key] += 1
                    else:
                        check_dict[key] == 0

  

              
            #column decoding
            for msg_block_index in range(len(msg_block_L_list_1)-1):
                if self.rank == 0 and msg_block_index == 0 :
                    continue
                if msg_block_index == 1:
                    if not self.rank == 0:
                        send_1 = msg_block_L_list_1[0][:int(v_column/2),::]
                        send_2 = msg_block_L_list_2[0][:int(v_column/2),::]
                        send_req_1 = self.comm.Isend(send_1,dest = self.rank -1,tag = 1)
                        send_req_1.Wait()
                        send_req_2 = self.comm.Isend(send_2,dest = self.rank -1,tag = 2)
                        send_req_2.Wait()
                    if not self.rank == self.mpi_nums-1 :
                        recv_req_1 = self.comm.Irecv(recv_1, source=self.rank + 1, tag=1)
                        recv_req_1.Wait()
                        recv_req_2 = self.comm.Irecv(recv_2, source=self.rank + 1, tag=2)
                        recv_req_2.Wait()
                        msg_block_L_list_1[-1][:int(v_column/2),:] = recv_1
                        msg_block_L_list_2[-1][:int(v_column/2),:] = recv_2
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

            self.get_fin_msg(self.fin_msg_block_list_1,self.fin_msg_block_list_2,2)

            if self.rank == 0:

                self.fin_msg_block_list =self.fin_msg_block_list[1:-1]
            else:

                self.fin_msg_block_list = self.fin_msg_block_list[:-1]

            recv_fin = self.comm.gather(self.fin_msg_block_list,root = 0)
            if self.rank == 0 :
                self.fin_msg_block_list = []
                for i in recv_fin :

                    self.fin_msg_block_list += i
                if self.acc_mode == True:
                    now_result = self.Acc()
                    output_list.append(now_result)
        MPI.Detach_buffer()
        if self.rank == 0 :
            name = self.msg_name + "_ACC"
            output_list.append(self.G_matrix)

            if self.output_mode == True :
                output_info = ""
                output_dict = {0:"A",1:"T",2:"G",3:"C"}
                index = 0
                
                for i in range(int(self.block_num/2)):
                    for j in range(self.G_matrix.shape[1]):
                        for k in range(self.msg_length):
                                if index < self.raw_length:
                                    output_info += output_dict[self.fin_msg_block_list[i][j][k]]
                                    index +=1
                np_output_name = self.msg_name+"decode_msg"
                g = open(np_output_name,"w")
                g.write(output_info)
                g.close()
                ref2file(np_output_name,self.output_name)
                    



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
                        error_bit += 1
        print("ACC:",(total_bit-error_bit)/total_bit)
        return (total_bit-error_bit)/total_bit






 