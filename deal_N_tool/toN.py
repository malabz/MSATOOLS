import sys, getopt
from numpy import *

def read_fasta(filename, outfile_path):
    with open(outfile_path, 'w') as out:
        with open(filename, 'r') as f:
            temp = ""
            strs = ""
            names = ""
            for line in f:
                if line.startswith('>'):
                    strs=temp
                    if(len(strs)!=0):
                        out.writelines('>' + names + "\n")
                        # for j in range(len(strs)):
                        #     if (strs[j] != 'A' and strs[j] != 'C' and strs[j] != 'G' and strs[j] != 'T'):
                        #         out.write('N')
                        #     else:
                        #         out.write(strs[j])
                        out.write(strs)
                        out.write('\n')
                    names=line[1:-1]
                    temp = ""
                    continue
                if line.startswith('\n'):
                    continue
                temp += line[:-1]
            strs=temp
            out.writelines('>' + names + "\n")
            # for j in range(len(strs)):
            #     if (strs[j] != 'A' and strs[j] != 'C' and strs[j] != 'G' and strs[j] != 'T'):
            #         out.write('N')
            #     else:
            #         out.write(strs[j])
            out.write(strs)
            out.write('\n')


if __name__ == '__main__':
    opts, args = getopt.getopt(sys.argv[1:], "hi:o:")
    infile_path = ""
    outfile_path = ""
    for op, value in opts:
        if op == "-i":
            infile_path = value
            print(infile_path)
        elif op == "-o":
            outfile_path = value
            print(outfile_path)

    # infile_path = "D:\zt\python\code\string\stmsa\data/4/GWHBEBU00000000.genome.id_GWHBEBU00000004.fasta"
    # outfile_path = "D:\zt\python\code\string\stmsa\data/4/4U.fasta"
    read_fasta(infile_path, outfile_path)
