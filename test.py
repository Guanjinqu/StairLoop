#Test code for StairLoop
from numpy import array
import commpy.channelcoding.convcode as cc
from decode import Decode_unit
from encode import Encode_unit
from LDPC import coding_matrix_systematic,regular_coding_matrix_systematic
from LDPC import log_LDPC_decode
import random
from load_conf import load_config
config_dict = load_config()

H_matrix = config_dict["H_matrix"]
G_matrix = config_dict["G_matrix"]

trellis = config_dict["trellis"]


pr_dict = config_dict["pr_dict"]


name= config_dict["name"]



unit1 = Encode_unit (
    G_matrix ,
    H_matrix,
    trellis,
    unit_name = name,
    msg_length = config_dict["msg_length"],
    index_length = config_dict["index_length"] ,
    block_num = config_dict["block_num"],
    input_file_name = config_dict["input_file_name"],
    input_mode = config_dict["input_mode"]
)

unit1.generate_msg_matrix()
unit1.generate_encode_matrix()
unit1.generate_channel_coding(extra_nums = config_dict["extra_nums"],pr_dict = pr_dict)


unit = Decode_unit (
    
    G_matrix ,
    H_matrix,
    trellis,
    msg_name = name,
    pr_dict = pr_dict,
    msg_length = config_dict["msg_length"],
    index_length = config_dict["index_length"] ,
    block_num = config_dict["block_num"],
    iterations = config_dict["iterations"],
    mpi_mode = False,
    output_mode = False,
    acc_mode = True
)

unit.fast_decode()

