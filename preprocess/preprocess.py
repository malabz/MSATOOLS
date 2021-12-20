"""

Author:zhai1xiao
Date:2021-12-19
Update:2021-12-20

"""

import re
import sys

def readFasta(file_path):
    """
    读取fasta文件

    :param file_path: 路径
    :return: identifiers [] 序列id的列表
             sequences [] 序列的列表
    """
    identifiers = []
    sequences = []
    with open(file_path, 'r') as fin:
        curr = ""
        flag = False
        for line in fin:
            if line[0] == '>':
                if flag:
                    sequences.append(curr.upper())
                curr = ""
                flag = True
                identifiers.append(line[1:].strip())
            else:
                curr += line.strip()
        sequences.append(curr.upper())
    return identifiers, sequences

def preprocess(sequences):
    """
    处理fasta序列中的非法字符（非字母用N代替）

    :param sequences 待处理的序列列表
    :return: processed [] 处理后的序列列表
    """
    processed = []
    for i in sequences:
        tmp = re.sub("[^ACTGNactgn-]", "N", str(i))
        processed.append(tmp)
    return processed

def createFasta(resultName, ids, sequences):
    """
    存储fasta文件

    :param resultName fasta文件名（resultName_processed.fasta）
           ids [] 序列id的列表
           sequences [] 序列的列表

    """
    with open(resultName + "_processed.fasta", "w") as f:
        for i in range(len(sequences)):
            f.write(">" + str(ids[i]) + "\n")
            f.write(str(sequences[i]) + "\n")

if __name__ == '__main__':
    filePath = sys.argv[1]
    resultName = sys.argv[2]
    ids, sequences = readFasta(filePath)
    new = preprocess(sequences)
    createFasta(resultName, ids, sequences)
