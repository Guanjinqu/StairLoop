
from Bio import SeqIO


base_map = {"A":0,"T":1,"G":2,"C":3}

def get_file(filename,output_name,mpi_nums,unit_nums = 0):
    
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
        more_file = open("more_ref.txt","r")
        if "N" in record.seq:
            continue
        #print(record.seq)
        if len(str(record.seq)) <152 or len(str(record.seq)) > 158 :
            continue
        now_read = "AACACA" + record.seq[21:-21]+"AC"
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
            
            for line in more_file.readlines():
                output_file.write(str([base_map[i] for i in line[:-2]]))
                output_file.write('\n')
            i = 1

        i +=1
    output_file.close()    


filename = "gen3.assembled.fastq"
output_name ="gen3"
mpi_nums = 10
unit_nums = 200000

get_file(filename,output_name,mpi_nums,unit_nums = unit_nums)
