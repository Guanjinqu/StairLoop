[![PythonVersion](https://img.shields.io/badge/python-3.7-blue)](https://img.shields.io/badge/python-3.7-blue)

# StairLoop
A High Error-Tolerance and High Parallelism DNA Coding Scheme

## Requirements
- scikit-commpy
- mpi4py
- biopython
- pyldpc
- scipy
- cython
- numpy
- tqdm
- pyyaml

## Quick Start
The official version of DNA StairLoop needs to be based on MPI, so you will need to additionally install tools such as MPICH and make sure it adapts to mpi4py.

To make it easier for users to quickly use DNA StairLoop to evaluate its performance. We give the DNA StairLoop test module in single process mode. The test module has similar error correction and iterative decoding capabilities as the official version, but it is not able to parallelise the computation and output of the decoded files, so the decoding speed of the test module will be much slower.

The decoder needs to be compiled before running DNA StairLoop:
> python setup.py build_ext --inplace

Compilation time is usually less than 1 minute.
### DNA StairLoop Test Usage
Run the **test.py** file directly after installing the python package:
> python test.py

StairLoop will simulate encoding 4 blocks and add 3% error to the encoded sequence, and the sequencing depth of 10. And the sequence will be successfully decoded.Code runs usually take less than 10 minutes.

### DNA StairLoop Official Usage

Run **main_encode.py** to encode:
> python main_encode.py 

Execute **mpiexec** for multi-process decoding, where the **-n** argument is the number of processes used:
> mpiexec -n 2 python main_decode.py

After executing the above two commands, DNA StairLoop will encode the "data" file in the folder, pass it through a 3% error rate, 10 sequencing depth simulation channel, and successfully decode and restore the file to "decode_result".

## Customization

### Customize Config
The user can modify the config.yaml file to meet the requirements.


- **Conv_mode [int]** : The coding structure of convolutional codes, where 0 is 1/2 code rate and 1 is 1/3 code rate.
- **LDPC_mode [int]** The coding structure of LDPC codes with code rates of 0, 1, 2, and 3 are 82/100, 64/100, 30/60, and 210/250.
- **acc_mode** : Decode Correct Rate Statistics. When turned on Decode will display the correct rate of decoding.
- **block_num [int]**: The number of blocks. Note that the number of blocks must be a multiple of 4.
- **extra_nums [int]**: Sequencing depth during simulation of channels
- **high_err_mode**: High Error Rate Mode. For some very high error rate cases, turning it on will improve the decoding results.
- **index_length [int]**: Index length. The index length needs to be increased when the number of blocks is large.
- **input_file_name [str]**: File name of the stored file
- **output_file_name [str]**: File name of the decoded output file
- **input_mode**: Input Mode. When turned on StairLoop's encoding will be the specified input file, otherwise the encoded information is randomly generated.
- **iterations [int]**: The number of iterations, increasing the number of iterations for high error rate cases is beneficial to increase the error correction capability.
- **mpi_nums [int]**: Number of processes. Note that this value must be equal to the number of processes specified by mpiexec at the time of decoding.
- **msg_length [int]**: Message length. The length of information in bases for each sequence.
- **name**: Name
- **output_mode**: Output Mode. When turned on StairLoop will output the information to the specified file after decoding.
- **pr_dict [dict]**: Error Rate Dictionary. The user can adjust the error rate of the analogue channel.
- **test_mode**: Test mode. When turned on StairLoop passes through the simulated channel after encoding and generates sequences with an error rate.

### Decode sequencing files
Since DNA StairLoop only accepts sequences in numpy array format. So for fastq format sequencing files, you need to run the **get_mpi_file.py** file first:
> python get_mpi_file.py

Before running, you need to change the contents of lines 56 to 59: 
- **filename**: Sequencing file name
- **output_name**: The name of the output file needs to be the same as the name parameter of config.
- **mpi_nums**: Number of processes
- **unit_nums**: The number of input sequences per process, usually the total number of sequences divided by the number of processes.

**real/get_mpi_file.py** is the processing code for the sequencing file of our in vitro experiment.

## License

DNA StairLoop is licensed under the GNU General Public License, for more information read the LICENSE file or refer to:

http://www.gnu.org/licenses/

## Citation

The paper is under review.
