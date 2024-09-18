
from numpy import array
import commpy.channelcoding.convcode as cc


from load_conf import load_config,get_label
import sys
import copy
import random
random.seed(1)

config_dict = load_config()

if config_dict["high_err_mode"] == True:
    from decode_mpi_h import Decode_mpi_unit
else:
    from decode_mpi import Decode_mpi_unit

H_matrix = config_dict["H_matrix"]
G_matrix = config_dict["G_matrix"]
trellis = config_dict["trellis"]
pr_dict = config_dict["pr_dict"]
name= config_dict["name"]


unit3 = Decode_mpi_unit (
    
    G_matrix ,
    H_matrix,
    trellis,
    msg_name = name,
    pr_dict = pr_dict,
    msg_length = config_dict["msg_length"],
    index_length = config_dict["index_length"] ,
    block_num = config_dict["block_num"],
    iterations = config_dict["iterations"],
    mpi_mode = True,
    output_mode = config_dict["output_mode"],
    acc_mode = config_dict["acc_mode"],
    input_name= config_dict["input_file_name"],
    raw_length = config_dict["raw_length"],
    output_name = config_dict["output_file_name"]
    

)
unit3.fast_decode()


