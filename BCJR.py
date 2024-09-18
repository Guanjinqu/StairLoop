#BCJR algorithm, note the need to call pyx file contents

from numpy import NaN, array, zeros, exp, log, empty
import commpy.channelcoding.convcode as cc
import random
random.seed(1)
from BCJR_unit import log_map_decode_unit

def log_map_decode(codewords_1,codewords_2,msg_length,trellis,pr_dict,D=2,mode="nomal",L_int_1=[],L_int_2 = [],F_dict = {}):
    '''
    Decoding is performed using the log-MAP algorithm.

    This function performs log MAP decoding based on the provided parameters such as codeword, message length, statechart, and a priori probability dictionary.
    When the mode is set to ‘ex’, the two codewords are decoded independently.

    Parameters: codewords_1
    - codewords_1: first codeword list
    - codewords_2: second codeword list
    - msg_length: length of the message
    - trellis: statechart object containing information about the statechart.
    - pr_dict: dictionary of prior probabilities
    - D: decoding depth
    - mode: decoding mode, ‘nomal’ is normal mode, ‘ex’ is extended mode.
    - L_int_1: internal state list of the first codeword
    - L_int_2: internal state list of the second codeword
    - F_dict: feedback dictionary

    Returns.
    - Returns the decoding result according to the decoding mode. In ‘ex’ mode, it returns a list of two decoding results.
    '''

    if mode == "ex" :
        result1 =  log_map_decode_unit(codewords_1,trellis,msg_length,trellis.next_state_table,trellis.output_table,pr_dict,D,"nomal",L_int_1,F_dict)
        result2 =  log_map_decode_unit(codewords_2,trellis,msg_length,trellis.next_state_table,trellis.output_table,pr_dict,D,"nomal",L_int_2,F_dict)
        return [result1,result2]
    