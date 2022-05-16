#encoding: utf-8
import sys
from itertools import combinations
import os
import time
import shutil
import matplotlib.pyplot as plt
from multiprocessing import Pool
import random

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
                    sequences.append(curr)
                curr = ""
                flag = True
                identifiers.append(line[1:].strip())
            else:
                curr += line.strip()
        sequences.append(curr)
    # print("done")
    return identifiers, sequences

def seqLengthDistribution(file_path, dataName):
    """
    1、计算数据集中，序列长度的最小、最大、平均、总和的值和数据集的大小，
    结果保存在dataName.txt中。
    2、绘制序列长度分布图，结果保存在dataName_Distri.png
    :param filepath: 路径
           dataName: 数据集的名字，用于命名结果
    :return:
    """
    _, seqs = read_fasta(file_path)
    length = []
    for i in seqs:
        length.append(len(i))
    minLength = min(length)
    maxLength = max(length)
    sumLength = sum(length)
    avgLength = sumLength/len(length)

    with open(dataName + ".txt", "w") as fout:
        fout.write("Dataset size: " + str(len(length)))
        fout.write("\n")
        fout.write("Min length: " + str(minLength))
        fout.write("\n")
        fout.write("Max length: " + str(maxLength))
        fout.write("\n")
        fout.write("Sum length: " + str(sumLength))
        fout.write("\n")
        fout.write("Avg length: " + str(avgLength))
        fout.write("\n")


    plt.hist(x=length,  # 指定绘图数据
             bins=20,  # 指定直方图中条块的个数
             color='steelblue',  # 指定直方图的填充色
             edgecolor='black'  # 指定直方图的边框色
             )

    # 添加x轴和y轴标签
    plt.xlabel('Length of Sequence')
    plt.ylabel('Count')
    # 添加标题
    plt.title(dataName + ' Sequence length distribution graph')
    plt.savefig(dataName + '_Distri.png')
    plt.close()

def compare(lhs: str, rhs: str) -> float:
    """
    比较两条序列的相似度（匹配字符数量/序列长度）
    :param lhs: 序列1 str
           rhs: 序列2 str
    :return: similarity: 相似度的值 float
    """
    lhs = lhs.lower()
    rhs = rhs.lower()
    n = len(lhs)
    cnt = 0
    for i in range(n):
        if lhs[i] == rhs[i]:
            cnt += 1
    similarity = cnt/n
    return similarity

def create_dir_not_exist(path):
    """
    判断path是否存在，不存在则创建
    :param filepath: 路径
    :return:
    """
    if not os.path.exists(path):
        os.mkdir(path)

def del_file(filepath):
    """
    删除某一目录下的所有文件或文件夹
    :param filepath: 路径
    :return:
    """
    try:
        del_list = os.listdir(filepath)
    except:
        os.remove(filepath)
        return
    for f in del_list:
        file_path = os.path.join(filepath, f)
        if os.path.isfile(file_path):
            os.remove(file_path)
        elif os.path.isdir(file_path):
            shutil.rmtree(file_path)

def dispatch_program(tempdir: str, i: int, j: int, idsi: str, idsj: str, seqi: str, seqj: str) -> float:
    """
    将多序列计算距离的程序采用多线程完成
    :param i: 序列1编号
           j: 序列2编号
           idsi: 序列1标识号
           idsj: 序列2标识号
           seqi: 序列1内容
           seqj: 序列2内容
    :return: ans: 序列相似度
    """
    file_path = tempdir + os.sep + str(i) + "_" + str(j) + ".fasta"
    global records
    with open(file_path, "w") as fout:
        fout.write(">")
        fout.write(idsi)
        fout.write("\n")
        fout.write(seqi)
        fout.write("\n")
        fout.write(">")
        fout.write(idsj)
        fout.write("\n")
        fout.write(seqj)
    outfile_path = file_path + ".aligned"
    os.system(mafft + " " + file_path + " > " + outfile_path + " 2> /dev/null")
    _, current_sequences = read_fasta(outfile_path)
    ans = compare(current_sequences[0], current_sequences[1])
    del_file(file_path)
    del_file(outfile_path)
    return ans


if __name__ == "__main__":
    if len(sys.argv) < 5 or sys.argv[1] == "help":
        print("usage: python", sys.argv[0], "inputfile_name", "mafft_path", 'output_prefix', "mafft_threads")
        print(" ")
        print("positional arguments:")
        print("  inputfile_name    specify the input file")
        print("  mafft_path        the path of executable file of mafft aligner")
        print("  output_prefix     specify the prefix of output files")
        print("  mafft_threads     multi-thread for running mafft alignment")
        print(" ")
        print("print this help information:", sys.argv[0], "help")
        exit(0)
    seqLengthDistribution(sys.argv[1], sys.argv[3])
    start = time.time()
    ids, sqs = read_fasta(sys.argv[1])
    threads = int(sys.argv[4])
    tempdirname = 'temp' + str(random.randint(1000, 32768) * random.randint(197, 2394))
    create_dir_not_exist(tempdirname)
    n = len(ids)
    mafft = sys.argv[2]
    index = []
    for i in range(n):
        index.append(i)
    records_index = list(combinations(index, 2))
    print()
    prei = 0
    sub = []
    records = []
    mypool = Pool(processes = threads)
    mafft_start = time.time()
    for i in records_index:
        sub.append( mypool.apply_async( dispatch_program, ( tempdirname, i[0], i[1], ids[i[0]], ids[i[1]], sqs[i[0]], sqs[i[1]], ) ) )
        if prei != i[0]:
            mypool.close()
            mypool.join()
            for k in sub:
                records.append( k.get() )
            sub = []
            del mypool
            mypool = Pool(processes = threads)
            prei = i[0]
            print("\rSTEP {0:} / {1:}".format(i[0], n - 1))
    if mypool:
        mypool.close()
        mypool.join()
        for k in sub:
            records.append( k.get() )
        del mypool
    mafft_end = time.time()
    plt.hist(x=records,  # 指定绘图数据
             bins=20,  # 指定直方图中条块的个数
             color='steelblue',  # 指定直方图的填充色
             edgecolor='black'  # 指定直方图的边框色
             )
    # 添加x轴和y轴标签
    plt.xlabel('Similarity of Sequence')
    plt.ylabel('Count')
    # 添加标题
    plt.title('Similarity Distribution of ' + sys.argv[3])
    plt.savefig('similarity_' + sys.argv[3] + '.png')
    with open('similarity_' + sys.argv[3] + '.txt', "w") as fout:
        for i in records:
            fout.write(str(i))
            fout.write("\n")
    end = time.time()
    del_file(tempdirname)
    os.removedirs(tempdirname)

    with open(sys.argv[3] + "_time.txt", "w") as fout:
        fout.write("Count of Compare: " + str(len(records_index)))
        fout.write("\n")
        fout.write("MAFFT cost time: " + str(mafft_end - mafft_start))
        fout.write("\n")
        fout.write("Cost time: "+ str(end-start))
        fout.write("\n")
