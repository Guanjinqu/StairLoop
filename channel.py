#Channel Simulation Algorithm

import random
import numpy as np

rd = random.Random() 
def channel_model_unit(code, pr_dict):
    """
    Channel simulation function, simulates the effects of a channel on a code sequence.
    
    Parameters:
    code: The original code sequence.
    pr_dict: A dictionary containing the probabilities of different channel errors.
    
    Returns:
    af_code_array: The code sequence after being affected by the channel.
    """
    del_num = 0
    ins_num = 0
    sub_num = 0
    pt_num = 0
    unit_list = [1, 0]

    new_pr_dict = {}
    for key in pr_dict:
        new_pr_dict[key] = 10000*pr_dict[key]

    af_code = []
    if rd.randint(1, 10000) <= new_pr_dict["column"]:
        for i in range(code.size):
            af_code.append(2)
    else:
        for i in range(code.size):
            ins_times = 0
            while ins_times < 1:  
                if rd.randint(1, 10000) <= new_pr_dict["pi"]:
                    af_code.append(random.choice(unit_list))
                    #print("ins:",i,ins_times)
                    ins_num = ins_num + 1
                else:
                    break
            if rd.randint(1, 10000) <= new_pr_dict["pd"]:
                #print("del:",i)
                del_num += 1
                continue
            else:
                pt_num += 1
                if rd.randint(1, 10000) <= new_pr_dict["ps"]:
                    af_code.append(abs(1-code[i]))
                    #print("sub:",i)
                    sub_num += 1
                else:
                    af_code.append(code[i])
                    
    # Convert the list to a numpy array for easier subsequent processing
    af_code_array = np.array(af_code)

    return af_code_array