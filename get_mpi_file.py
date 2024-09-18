
from Bio import SeqIO


base_map = {"A":0,"T":1,"G":2,"C":3}

def get_file(filename,output_name,mpi_nums,unit_nums = 0):
    """
    Processes the sequencing file and splits the content into multiple output files according to the specified parameters, allowing StairLoop to decode it.
    
    Parameters.
    - filename: name of the input file, containing the series data.
    - output_name: base name of the output file to generate the series.
    - mpi_nums: number of output files to be generated.
    - unit_nums: the number of sequences in each output file, default is 0, means it is calculated automatically.
    
    This function reads a fastq format file and divides the sequence data into multiple files according to the specified number of output files.
    If unit_nums is not specified, the approximate number of sequences that each output file should contain will be calculated automatically.
    """
    
    records = SeqIO.parse(filename, "fastq")
    if unit_nums == 0 :
        unit_nums = len(list(records))//mpi_nums
        print(unit_nums)
    
    ex_name = output_name+"_mpi_coding_"
    name_index = 0
    output_name = ex_name+str(name_index)
    output_file = open(output_name,"w")
    i = 1
    for record in records:

        if "N" in record.seq:
            continue

        now_read =  record.seq
        read = [base_map[i] for i in now_read]
        if i < unit_nums:
            output_file.write(str(read))
            output_file.write('\n')
        else:
            if name_index == mpi_nums -1 :
                break
            output_file.close()
            name_index +=1 
            output_name= ex_name+str(name_index)
            output_file = open(output_name,"w")
            output_file.write(str(read))
            output_file.write('\n')
            i = 1

        i +=1
    output_file.close()    


filename = "test.fastq"
output_name ="gen3"
mpi_nums = 10
unit_nums = 200000

get_file(filename,output_name,mpi_nums,unit_nums = unit_nums)
