#BCJR algorithm

import numpy as np
from numpy import  zeros
import math
from commpy.utilities import dec2bitarray
import time


def max_star_unit(x,y):
    #Calculating the MAX* function

    output = max(x,y)+math.log(1+math.exp(-abs(x-y)))
    #output = max(x,y)                         
    return output

def max_star_fuc(nums_list):
    #Calculating the MAX* function for the multivariate case

    if len(nums_list) == 1:
        return nums_list[0]
    else:
        x=nums_list[0]
        y=nums_list[1]
        nums_list=nums_list[2:]
        output=max_star_unit(x,y)
        while True:
            if nums_list ==[]:
                break
            unit=nums_list[0]
            nums_list= nums_list[1:]
            output=max_star_unit(output,unit)
        return output

def P_fuc(r_list,x_list,pr_dict,a,b,j):
    #Calculation of the P-function for the pilot calculation of the subsequent gamma function values

    p_c = (pr_dict["pi"]/2)*pr_dict["pd"] + pr_dict["pt"]*(1-pr_dict["ps"])
    p_c_2 = (pr_dict["pi"]/2)*pr_dict["pi"] + pr_dict["pt"]*pr_dict["ps"]
    if j == -1:
        return math.log(pr_dict["pd"])
    elif j >= 0 and r_list[a-1] == x_list[b-1]:
        return math.log(p_c)+j*math.log(pr_dict['pi']/2)
    elif j >=0 :
        return math.log(p_c_2)+j*math.log(pr_dict['pi']/2)

def log_F_fuc(r_list,x_list,pr_dict,a,b,F_dict) :
    #Compute the value of the F function under the log domain

    if a == b ==0 :
        return 0
    elif b == 0 :
        return -10000
    else:
        if (r_list,x_list,a,b) in F_dict:
            return F_dict[(r_list,x_list,a,b)]
        else:
            max_list=[]
            for j in range(-1,a):
                max_list.append(P_fuc(r_list,x_list,pr_dict,a,b,j)+log_F_fuc(r_list,x_list,pr_dict,a-j-1,b-1,F_dict))
            F_dict[(r_list,x_list,a,b)] = max_star_fuc(max_list)
            return F_dict[(r_list,x_list,a,b)]
 
def other_gamma_unit(metrics,tao,current_d,next_d):
    result_gamma = list(metrics[tao,:,:,current_d,next_d].flatten())
    return max_star_fuc(result_gamma)


def _compute_gamma(k,n,tao,d_max,d_min,next_table,output_table,number_states,number_inputs,gamma_metrics,codewords,pr_dict,L_int,gamma_star_metrics,F_dict):
    
    #Calculate the value of gamma and export to the value space of gamma. Also calculate the value of gamma* and export to the value space of gamma*

    next_state_table = next_table
  
    #The following is to do the traversal loop, the purpose of which is to traverse all the positions in the gamma matrix and assign values, where the operation of continue is to remove some of the spatial points that do not meet the conditions
    
    for t in range(1,tao+1):
        for current_state in range(number_states):
            for current_input in range(number_inputs):
                next_state = next_state_table[current_state, current_input]      #当前状态 当前输入下的 下一个状态
                code_symbol = output_table[current_state, current_input]         #当前状态 当前输入下的 输出值v
                x_codeword = dec2bitarray(code_symbol, n)                    #将输出值v转化为二进制变成数组，第二位为数组的大小
                input_codeword = dec2bitarray(current_input,k)               #将输入值u转化为二进制变成数组
                for current_drift in range(d_min,d_max+1):
                    if abs(current_drift) > 2*t  :       #这部分是限制d的绝对值必须小于2τ
                        continue
                    for next_drift_nums in range(-1,2) :
                        next_drift = current_drift+next_drift_nums
                        if next_drift > d_max or next_drift < d_min : #限制next_drift的取值
                            continue
                        r_codeword = codewords[current_drift+(t-1)*n:next_drift+t*n]
                        if r_codeword.size == 0:           
                            continue
                        if next_drift+t*n > codewords.size :
                            continue
                        for i in range(1,k+1) :
                            
                            
                            current_gama_star = 0
                            r_list = tuple(r_codeword.tolist())
                            x_list = tuple(x_codeword.tolist())
                            #以下为最关键的一步：通过F函数计算伽马值
                            current_gama_star += log_F_fuc(r_list,x_list,pr_dict,n+next_drift-current_drift,n,F_dict)

                            for j in range(1,k+1):
                                if j == i :
                                    continue
                                current_gama_star += (1-2*input_codeword[j-1])*L_int[k*(t-1)+j-1]/2

                            gamma_star_metrics[t,current_state,next_state,current_drift-d_min,next_drift-d_min,i] = current_gama_star

                            current_gama = current_gama_star+(1-2*input_codeword[i-1])*L_int[k*(t-1)+i-1]/2

                            gamma_metrics[t,current_state,next_state,current_drift-d_min,next_drift-d_min] = current_gama



def _compute_alpha(n,tao,d_max,d_min,next_table,number_states,number_inputs,alpha_metrics,codewords,gamma_metrics):
    #计算α的值，并导出到α的值空间中

    next_state_table = next_table  #行表示当前状态，列表示输入，值表示输出状态
    for t in range(1,tao+1):
        for current_state in range(number_states):
            for current_input in range(number_inputs):
                next_state = next_state_table[current_state, current_input]
                for current_drift in range(d_min,d_max+1):             
                    if abs(current_drift) > 2*t :       #这部分是限制d的绝对值必须小于2τ
                        continue
                    for next_drift_nums in range(-1,2) :
                        next_drift = current_drift+next_drift_nums
                        if next_drift > d_max or next_drift < d_min :
                            continue
                        r_codeword = codewords[current_drift+(t-1)*n:next_drift+t*n]
                        if r_codeword.size == 0:
                            continue
                        if next_drift > d_max or next_drift < d_min :
                            continue 
                        if alpha_metrics[t-1,current_state,current_drift-d_min] == -10000 :
                            continue
                        if alpha_metrics[t,next_state,next_drift-d_min] == -10000:
                            alpha_metrics[t,next_state,next_drift-d_min] = gamma_metrics[t,current_state,next_state,current_drift-d_min,next_drift-d_min]+alpha_metrics[t-1,current_state,current_drift-d_min]
                        else:
                            alpha_metrics[t,next_state,next_drift-d_min] = max_star_unit(alpha_metrics[t,next_state,next_drift-d_min],gamma_metrics[t,current_state,next_state,current_drift-d_min,next_drift-d_min]+alpha_metrics[t-1,current_state,current_drift-d_min])
               

def _compute_beta(n,tao,d_max,d_min,next_table,number_states,number_inputs,beta_metrics,codewords,gamma_metrics):
    #计算β的值，并导出到β的值空间中
    
    next_state_table = next_table  #行表示当前状态，列表示输入，值表示输出状态

    for t in reversed(range(1,tao)):
        for current_state in range(number_states):
            for current_input in range(number_inputs):
                next_state = next_state_table[current_state, current_input]
                for current_drift in range(d_min,d_max+1):            
                    if abs(current_drift) > 2*t :       #这部分是限制d的绝对值必须小于2τ
                        continue
                    for next_drift_nums in range(-1,2) :
                        next_drift = current_drift+next_drift_nums
                        if next_drift > d_max or next_drift < d_min :
                            continue
                        r_codeword = codewords[current_drift+(t-1)*n:next_drift+t*n]
                        if r_codeword.size == 0:
                            continue
                        if beta_metrics[t+1,next_state,next_drift-d_min] == -10000:
                            continue
                        if beta_metrics[t,current_state,current_drift-d_min] == -10000 :
                            beta_metrics[t,current_state,current_drift-d_min] = gamma_metrics[t+1,current_state,next_state,current_drift-d_min,next_drift-d_min]+beta_metrics[t+1,next_state,next_drift-d_min]
                        else:
                            beta_metrics[t,current_state,current_drift-d_min] = max_star_unit(beta_metrics[t,current_state,current_drift-d_min],(gamma_metrics[t+1,current_state,next_state,current_drift-d_min,next_drift-d_min]+beta_metrics[t+1,next_state,next_drift-d_min]))
def _compute_L(k,n,tao,d_max,d_min,next_table,output_table,number_states,number_inputs,codewords,L_int,alpha_metrics,beta_metrics,gamma_star_metrics,L_metrics,fin_L,output_L):
    #根据α、β、γ的值，来求出LLR（log域）

    next_state_table = next_table  #行表示当前状态，列表示输入，值表示输出状态

    for t in range(1,tao+1):
        
        for current_state in range(number_states):
            for current_input in range(number_inputs):
                next_state = next_state_table[current_state, current_input]      #当前状态 当前输入下的 下一个状态
                code_symbol = output_table[current_state, current_input]         #当前状态 当前输入下的 输出值v
                x_codeword = dec2bitarray(code_symbol, n)                    #将输出值v转化为二进制变成数组，第二位为数组的大小

                input_codeword = dec2bitarray(current_input,k)               #将输入值u转化为二进制变成数组
                for current_drift in range(d_min,d_max+1):
                    if abs(current_drift) > 2*t :       #这部分是限制d的绝对值必须小于2τ
                        continue
                    for next_drift_nums in range(-1,2) :
                        next_drift = current_drift+next_drift_nums   
                        if next_drift > d_max or next_drift < d_min :
                            continue
                        r_codeword = codewords[current_drift+(t-1)*n:next_drift+t*n]
                        if r_codeword.size == 0 :
                            continue   
                        for i in range(1,k+1) :
                            if (
                                alpha_metrics[t-1,current_state,current_drift-d_min]==-10000 or 
                                gamma_star_metrics[t,current_state,next_state,current_drift-d_min,next_drift-d_min,i] ==-10000 or
                                beta_metrics[t,next_state,next_drift-d_min]==-10000
                            )  : 
                                
                                #print(k*(t-1)+i-1,"no")
                                continue
                            if input_codeword[i-1] == 0 :
                                
                                if L_metrics[k*(t-1)+i-1,0] == 0 :
                                    L_metrics[k*(t-1)+i-1,0]=alpha_metrics[t-1,current_state,current_drift-d_min]+gamma_star_metrics[t,current_state,next_state,current_drift-d_min,next_drift-d_min,i]+beta_metrics[t,next_state,next_drift-d_min]
                                else:
                                    L_metrics[k*(t-1)+i-1,0] = max_star_unit(L_metrics[k*(t-1)+i-1,0],
                                (alpha_metrics[t-1,current_state,current_drift-d_min]+gamma_star_metrics[t,current_state,next_state,current_drift-d_min,next_drift-d_min,i]+beta_metrics[t,next_state,next_drift-d_min]))
  
                            if input_codeword[i-1] == 1 :
                                if L_metrics[k*(t-1)+i-1,1] == 0 :
                                    L_metrics[k*(t-1)+i-1,1] = alpha_metrics[t-1,current_state,current_drift-d_min]+gamma_star_metrics[t,current_state,next_state,current_drift-d_min,next_drift-d_min,i]+beta_metrics[t,next_state,next_drift-d_min]
                                else:
                                    L_metrics[k*(t-1)+i-1,1] = max_star_unit(L_metrics[k*(t-1)+i-1,1],
                                (alpha_metrics[t-1,current_state,current_drift-d_min]+gamma_star_metrics[t,current_state,next_state,current_drift-d_min,next_drift-d_min,i]+beta_metrics[t,next_state,next_drift-d_min]))

    for t in range(1,tao+1):
        for i in range(1,k+1) :
            fin_L[k*(t-1)+i-1] = L_metrics[k*(t-1)+i-1,0] - L_metrics[k*(t-1)+i-1,1]  
            output_L[k*(t-1)+i-1] = fin_L[k*(t-1)+i-1] +L_int[k*(t-1)+i-1]+ L_int[k*(t-1)+i-1]





             
def log_map_decode_unit(codewords,trellis,msg_length,next_table,output_table,pr_dict,D=2,mode="nomal",L_int=[],F_dict = {}):
    """
    卷积码编码器

    负责将码字通过log——MAP译码算法译码成原始信息,同时允许输出L值以进行迭代译码。

    --------------------------------------------------------------------------
    输入：

    codewords_list: list [codewords,len_msg.len_U]      码字列表
                    codewords:  经过信道后的码字
                    len_msg :   原始信息位的长度
                    len_U :     结尾清零的信息位长度

    trellis :  class                                    卷积码编码器

    pr_dict : list  {"pt","pr","pi","pd"}               错误率字典
            
    D : int                                             飘移值的冗余阈值

    mode : "nomal","exinfo"                                         暂时没用

    L_int : list                                        先验概率
    ---------------------------------------------------------------------------


    ---------------------------------------------------------------------------
    参数：

    太多了懒得写了，重要参数我一般后边都标了，再不懂就问我吧:D
    ---------------------------------------------------------------------------

    ---------------------------------------------------------------------------
    输出：

    [msg,L_list]                二元列表，其中
                msg :           一维数组，译码后的信息位
                L_list :        二维数组， 译码后的LLR值数组
    ---------------------------------------------------------------------------


    """
    
    n = trellis.n
    k = trellis.k

    st = time.time()
    
    if "pt" not in pr_dict :
        pr_dict["pt"]=((1-pr_dict["pi"])*(1-pr_dict["pd"]))
   
    codewords_list = [codewords,msg_length,msg_length+1]
    rate = float(k)/n   
    L_len = int(codewords_list[2]/rate)    #Length of v
    len_msg=int(codewords_list[2])         #Length of u
    len_msg_true= int(codewords_list[1])
    tao = int(L_len/n)    #Length of τ
    Drift = len(codewords) - L_len #displacement difference
    if  mode== "exinfo" :
        Drift = 0
        D = 0
    if  L_int==[]:
        L_int=zeros(len_msg)

    if Drift >=0 :
        d_max=Drift+D
        d_min=-D
    else:
        d_max=D
        d_min=Drift-D
    
    d_space=d_max-d_min+1  #Size of space for drifting
    #print(d_max,d_min)
    number_states = trellis.number_states #Number of states
    number_inputs = trellis.number_inputs #The number of branches per state in trellis.

    

    alpha_metrics = zeros([tao+1,number_states,d_space])    #Define the value space of alpha, information bit x state bit xd space
    alpha_metrics[:,:,:] = -10000
    alpha_metrics[0,0,00-d_min] = 0
    

    beta_metrics = zeros([tao+1,number_states,d_space])    #Define the value space of beta, information bit x state bit xd space
    beta_metrics[:,:,:] = -10000
    beta_metrics[tao,:,Drift-d_min] = 0


    gamma_metrics = zeros([tao+1,number_states,number_states,d_space,d_space]) #Define the value space of gama, information bit x state bit xd space xk
    gamma_star_metrics = zeros([tao+1,number_states,number_states,d_space,d_space,k+1]) #Define the value space of gama*, information bit x state bit xd space xk
    if mode == "exinfo" :
        gamma_star_metrics = zeros([tao+1,number_states,number_states,d_space,d_space,n+1]) #Define the value space of gama*, information bit x state bit xd space xk      
    gamma_metrics[:,:,:,:,:]=-10000
    gamma_star_metrics[:,:,:,:,:,:]=-10000
    #print((tao+1)*number_states*d_space,)

    if mode == "nomal" :
        L_metrics = zeros([len_msg,2])
        fin_L=zeros(len_msg)
        output_L = zeros(len_msg)
    else:
        L_metrics = zeros([L_len,2])   
        fin_L=zeros(L_len)
        output_L = zeros(L_len)
    #L_metrics[:,:]=-10000


    if mode == "nomal" :
        _compute_gamma(k,n,tao,d_max,d_min,next_table,output_table,
                    number_states,number_inputs,gamma_metrics,codewords,pr_dict,L_int,gamma_star_metrics,F_dict)
        #print("计算γ耗时：",time.time()-st)
        _compute_alpha(n,tao,d_max,d_min,next_table,
                  number_states,number_inputs,alpha_metrics,codewords,gamma_metrics)
        #print("计算α耗时：",time.time()-st)
        _compute_beta(n,tao,d_max,d_min,next_table,
                  number_states,number_inputs,beta_metrics,codewords,
                  gamma_metrics)
        #print("计算β耗时：",time.time()-st)
        _compute_L(k,n,tao,d_max,d_min,next_table,output_table,
                    number_states,number_inputs,codewords,
                   L_int,alpha_metrics,beta_metrics,
                    gamma_star_metrics,L_metrics,fin_L,output_L)



    fin_list=[]
    for i in range(len_msg):
        if output_L[i] > 0 :
            fin_list.append(0)
        else:
            fin_list.append(1)

    return [np.array(fin_list),fin_L,output_L]








