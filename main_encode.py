
from numpy import array
import commpy.channelcoding.convcode as cc
from encode import Encode_unit
from load_conf import load_config,get_label
import numpy as np

from ref2bin import file2ref

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
    block_num =config_dict["block_num"],
    input_file_name = config_dict["input_file_name"],
    input_mode = config_dict["input_mode"]
)

unit1.generate_msg_matrix()
unit1.generate_encode_matrix()
if config_dict["test_mode"] == True:
    unit1.generate_mpi_channel_coding(extra_nums = config_dict["extra_nums"],mpi_nums = config_dict["mpi_nums"],pr_dict = pr_dict)
out_name = name + ".npy"
np.save(out_name,G_matrix)
