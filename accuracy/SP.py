"""

Author:zhai1xiao
Create Date:2021-12-06
Update Date:2021-12-19

"""
import re
import argparse
from itertools import combinations

def read_fasta(file_path):
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
    处理fasta序列中的非法字符（非字母用-代替）

    :param sequences 待处理的序列列表
    :return: processed [] 处理后的序列列表
    """
    processed = []
    for i in sequences:
        tmp = re.sub("[^A-Za-z]", "-", str(i))
        processed.append(tmp)
    return processed

def elementCount(arr):
    """
    利用字典结构统计列表中各个元素的含量

    Args:
        arr 待统计的列表

    Return:
        result
    """
    result = {}
    for i in set(arr):
        result[i] = arr.count(i)
    return result

def evaluate(sequences, matchS, mismatchS, gap1S, gap2S):
    """
    序列打分

    Args:
        sequences  待打分的序列
        matchS     match(碱基相同)的得分
        mismatchS  mismatch(碱基非空，但不相同)的得分
        gap1S      gap比对碱基的得分
        gap2S      gap比对gap的得分

    Return:
        sp score
    """
    match = 0         # nongap==nongap
    mismatch = 0      # nongap!=nongap
    gap_1 = 0         # gap-nongap
    gap_2 = 0         # gap-gap
    for i in range(len(sequences[0])):
        part = []
        for k in range(len(sequences)):
            seed = sequences[k][i]
            part.append(seed)
        countOfSeed = set(part)
        if ("-" not in countOfSeed) and (len(countOfSeed) == 1 or (("N" in countOfSeed) and (len(countOfSeed) == 2))):
            match += (len(sequences)*(len(sequences)-1))/2
        elif ("-" in countOfSeed and len(countOfSeed) == 1):
            gap_2 += (len(sequences)*(len(sequences)-1))/2
        else:
            tmpDict = elementCount(part)
            kind = []
            comGap = 0
            for key in tmpDict.keys():
                if key == '-':
                    comGap = tmpDict[key]
                else:
                    kind.append(tmpDict[key])
            gap_2 += (comGap*(comGap-1)/2)
            gap_1 += (sum(kind) * comGap)
            if len(kind) == 1:
                for i in range(len(kind)):
                    match += kind[i]*(kind[i]-1)/2
            elif len(kind) == 2:
                for i in range(len(kind)):
                    match += kind[i]*(kind[i]-1)/2
                mismatch += kind[0] * kind[1]
            else:
                for i in range(len(kind)):
                    match += kind[i]*(kind[i]-1)/2
                comKind = combinations(kind,2)
                for i in comKind:
                    mismatch += i[0] * i[1]

    score = (match * matchS) + (mismatch * mismatchS) + (gap_1 * gap1S) + (gap_2 * gap2S)
    return score

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, required=True, default='', help='Fasta file path to be scored.')
    parser.add_argument('--match', type=float, default=1, help='match score, default=1.')
    parser.add_argument('--mismatch', type=float, default=-1, help='mismatch score, default=-1.')
    parser.add_argument('--gap1', type=float, default=-2, help='gap-base score, default=-2.')
    parser.add_argument('--gap2', type=float, default=0, help='gap-gap score, default=0.')
    args = parser.parse_args()

    filepath = args.input
    matchScore = args.match
    mismatchScore = args.mismatch
    gap1Score = args.gap1
    gap2Score = args.gap2
    _, sequences = read_fasta(filepath)
    processedSeq = preprocess(sequences)
    sp = evaluate(processedSeq, matchScore, mismatchScore, gap1Score, gap2Score)     # pure sp
    avgSp = sp/(len(sequences)*(len(sequences)-1)/2)                                 # avg sp
    scaledSP = avgSp/len(sequences[0])                                               # scaled sp
    print("SP score: " + str(sp))
    print("Avg SP score: " + str(avgSp))
    print("Scaled SP score: " + str(scaledSP))
    print("Finish!")
