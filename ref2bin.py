#Interconversion of DNA sequences to original files

from re import S
import struct
import time

df = {"AAAA":0,"AAAT":1,"AAAG":2,"AAAC":3,"AATA":4,"AATT":5,"AATG":6,"AATC":7,"AAGA":8,"AAGT":9,"AAGG":10,"AAGC":11,"AACA":12,"AACT":13,"AACG":14,"AACC":15,"ATAA":16,"ATAT":17,"ATAG":18,"ATAC":19,"ATTA":20,"ATTT":21,"ATTG":22,"ATTC":23,"ATGA":24,"ATGT":25,"ATGG":26,"ATGC":27,"ATCA":28,"ATCT":29,"ATCG":30,"ATCC":31,"AGAA":32,"AGAT":33,"AGAG":34,"AGAC":35,"AGTA":36,"AGTT":37,"AGTG":38,"AGTC":39,"AGGA":40,"AGGT":41,"AGGG":42,"AGGC":43,"AGCA":44,"AGCT":45,"AGCG":46,"AGCC":47,"ACAA":48,"ACAT":49,"ACAG":50,"ACAC":51,"ACTA":52,"ACTT":53,"ACTG":54,"ACTC":55,"ACGA":56,"ACGT":57,"ACGG":58,"ACGC":59,"ACCA":60,"ACCT":61,"ACCG":62,"ACCC":63,"TAAA":64,"TAAT":65,"TAAG":66,"TAAC":67,"TATA":68,"TATT":69,"TATG":70,"TATC":71,"TAGA":72,"TAGT":73,"TAGG":74,"TAGC":75,"TACA":76,"TACT":77,"TACG":78,"TACC":79,"TTAA":80,"TTAT":81,"TTAG":82,"TTAC":83,"TTTA":84,"TTTT":85,"TTTG":86,"TTTC":87,"TTGA":88,"TTGT":89,"TTGG":90,"TTGC":91,"TTCA":92,"TTCT":93,"TTCG":94,"TTCC":95,"TGAA":96,"TGAT":97,"TGAG":98,"TGAC":99,"TGTA":100,"TGTT":101,"TGTG":102,"TGTC":103,"TGGA":104,"TGGT":105,"TGGG":106,"TGGC":107,"TGCA":108,"TGCT":109,"TGCG":110,"TGCC":111,"TCAA":112,"TCAT":113,"TCAG":114,"TCAC":115,"TCTA":116,"TCTT":117,"TCTG":118,"TCTC":119,"TCGA":120,"TCGT":121,"TCGG":122,"TCGC":123,"TCCA":124,"TCCT":125,"TCCG":126,"TCCC":127,"GAAA":128,"GAAT":129,"GAAG":130,"GAAC":131,"GATA":132,"GATT":133,"GATG":134,"GATC":135,"GAGA":136,"GAGT":137,"GAGG":138,"GAGC":139,"GACA":140,"GACT":141,"GACG":142,"GACC":143,"GTAA":144,"GTAT":145,"GTAG":146,"GTAC":147,"GTTA":148,"GTTT":149,"GTTG":150,"GTTC":151,"GTGA":152,"GTGT":153,"GTGG":154,"GTGC":155,"GTCA":156,"GTCT":157,"GTCG":158,"GTCC":159,"GGAA":160,"GGAT":161,"GGAG":162,"GGAC":163,"GGTA":164,"GGTT":165,"GGTG":166,"GGTC":167,"GGGA":168,"GGGT":169,"GGGG":170,"GGGC":171,"GGCA":172,"GGCT":173,"GGCG":174,"GGCC":175,"GCAA":176,"GCAT":177,"GCAG":178,"GCAC":179,"GCTA":180,"GCTT":181,"GCTG":182,"GCTC":183,"GCGA":184,"GCGT":185,"GCGG":186,"GCGC":187,"GCCA":188,"GCCT":189,"GCCG":190,"GCCC":191,"CAAA":192,"CAAT":193,"CAAG":194,"CAAC":195,"CATA":196,"CATT":197,"CATG":198,"CATC":199,"CAGA":200,"CAGT":201,"CAGG":202,"CAGC":203,"CACA":204,"CACT":205,"CACG":206,"CACC":207,"CTAA":208,"CTAT":209,"CTAG":210,"CTAC":211,"CTTA":212,"CTTT":213,"CTTG":214,"CTTC":215,"CTGA":216,"CTGT":217,"CTGG":218,"CTGC":219,"CTCA":220,"CTCT":221,"CTCG":222,"CTCC":223,"CGAA":224,"CGAT":225,"CGAG":226,"CGAC":227,"CGTA":228,"CGTT":229,"CGTG":230,"CGTC":231,"CGGA":232,"CGGT":233,"CGGG":234,"CGGC":235,"CGCA":236,"CGCT":237,"CGCG":238,"CGCC":239,"CCAA":240,"CCAT":241,"CCAG":242,"CCAC":243,"CCTA":244,"CCTT":245,"CCTG":246,"CCTC":247,"CCGA":248,"CCGT":249,"CCGG":250,"CCGC":251,"CCCA":252,"CCCT":253,"CCCG":254,"CCCC":255}
rdf = ('AAAA', 'AAAT', 'AAAG', 'AAAC', 'AATA', 'AATT', 'AATG', 'AATC', 'AAGA', 'AAGT', 'AAGG', 'AAGC', 'AACA', 'AACT', 'AACG', 'AACC', 'ATAA', 'ATAT', 'ATAG', 'ATAC', 'ATTA', 'ATTT', 'ATTG', 'ATTC', 'ATGA', 'ATGT', 'ATGG', 'ATGC', 'ATCA', 'ATCT', 'ATCG', 'ATCC', 'AGAA', 'AGAT', 'AGAG', 'AGAC', 'AGTA', 'AGTT', 'AGTG', 'AGTC', 'AGGA', 'AGGT', 'AGGG', 'AGGC', 'AGCA', 'AGCT', 'AGCG', 'AGCC', 'ACAA', 'ACAT', 'ACAG', 'ACAC', 'ACTA', 'ACTT', 'ACTG', 'ACTC', 'ACGA', 'ACGT', 'ACGG', 'ACGC', 'ACCA', 'ACCT', 'ACCG', 'ACCC', 'TAAA', 'TAAT', 'TAAG', 'TAAC', 'TATA', 'TATT', 'TATG', 'TATC', 'TAGA', 'TAGT', 'TAGG', 'TAGC', 'TACA', 'TACT', 'TACG', 'TACC', 'TTAA', 'TTAT', 'TTAG', 'TTAC', 'TTTA', 'TTTT', 'TTTG', 'TTTC', 'TTGA', 'TTGT', 'TTGG', 'TTGC', 'TTCA', 'TTCT', 'TTCG', 'TTCC', 'TGAA', 'TGAT', 'TGAG', 'TGAC', 'TGTA', 'TGTT', 'TGTG', 'TGTC', 'TGGA', 'TGGT', 'TGGG', 'TGGC', 'TGCA', 'TGCT', 'TGCG', 'TGCC', 'TCAA', 'TCAT', 'TCAG', 'TCAC', 'TCTA', 'TCTT', 'TCTG', 'TCTC', 'TCGA', 'TCGT', 'TCGG', 'TCGC', 'TCCA', 'TCCT', 'TCCG', 'TCCC', 'GAAA', 'GAAT', 'GAAG', 'GAAC', 'GATA', 'GATT', 'GATG', 'GATC', 'GAGA', 'GAGT', 'GAGG', 'GAGC', 'GACA', 'GACT', 'GACG', 'GACC', 'GTAA', 'GTAT', 'GTAG', 'GTAC', 'GTTA', 'GTTT', 'GTTG', 'GTTC', 'GTGA', 'GTGT', 'GTGG', 'GTGC', 'GTCA', 'GTCT', 'GTCG', 'GTCC', 'GGAA', 'GGAT', 'GGAG', 'GGAC', 'GGTA', 'GGTT', 'GGTG', 'GGTC', 'GGGA', 'GGGT', 'GGGG', 'GGGC', 'GGCA', 'GGCT', 'GGCG', 'GGCC', 'GCAA', 'GCAT', 'GCAG', 'GCAC', 'GCTA', 'GCTT', 'GCTG', 'GCTC', 'GCGA', 'GCGT', 'GCGG', 'GCGC', 'GCCA', 'GCCT', 'GCCG', 'GCCC', 'CAAA', 'CAAT', 'CAAG', 'CAAC', 'CATA', 'CATT', 'CATG', 'CATC', 'CAGA', 'CAGT', 'CAGG', 'CAGC', 'CACA', 'CACT', 'CACG', 'CACC', 'CTAA', 'CTAT', 'CTAG', 'CTAC', 'CTTA', 'CTTT', 'CTTG', 'CTTC', 'CTGA', 'CTGT', 'CTGG', 'CTGC', 'CTCA', 'CTCT', 'CTCG', 'CTCC', 'CGAA', 'CGAT', 'CGAG', 'CGAC', 'CGTA', 'CGTT', 'CGTG', 'CGTC', 'CGGA', 'CGGT', 'CGGG', 'CGGC', 'CGCA', 'CGCT', 'CGCG', 'CGCC', 'CCAA', 'CCAT', 'CCAG', 'CCAC', 'CCTA', 'CCTT', 'CCTG', 'CCTC', 'CCGA', 'CCGT', 'CCGG', 'CCGC', 'CCCA', 'CCCT', 'CCCG', 'CCCC')

def pack4f(s):
    return struct.pack("B",df[s])

def unpack4f(data):
    return rdf[data]

def packall(s):
    data = b""
    for i in range(len(s)//4):
        data += pack4f(s[i*4:(i+1)*4])
    return data

def unpackall(data):
    unp_data = ""
    for i in data:
        unp_data += unpack4f(i)
    return unp_data


def file2ref(path_in,path_out):
    """
    Converts a file into quadratic DNA format and stores it.
    path_in: input file address
    path_out:output file address in binary format
    """
    f_in = open(path_in,"rb")
    f_out = open(path_out,"w")
    while True:
        data_chunk = f_in.read(10240)
        if not data_chunk:
            break
        atgc_data = unpackall(data_chunk)
        f_out.write(atgc_data)
    f_in.close()
    f_out.close()



def ref2file(path_in,path_out):
    """
    Decode the DNA sequence file back to the original file.
    path_in:input DNA file address
    path_out:output original file address
    """
    f_in = open(path_in,"r")
    f_out = open(path_out,"wb")

    count = 0
    while True:
        #print(count)
        data_chunk = f_in.read(10240)
        if not data_chunk:
            break
        bin_data = packall(data_chunk)
        f_out.write(bin_data)

        count += 1

    f_in.close()
    f_out.close()


if __name__ == '__main__':
    file2ref("test.txt","test_data")
    #ref2file("gen3decode_msg","data_test.zip")


