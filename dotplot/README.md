# 将比对后序列作图的python工具

功能：对.fasta和.maf格式比对后序列作图

前提：序列名称不能带空格，maf不可省略中心序列名称

  包含：
  
       1、python源码                dotplot.py    (python dotplot.py ...)
           
       2、linux可执行程序           dotplot       (./dotplot ...)
           
       3、windows可执行程序         dotplot.exe   (dotplot.exe ...)
           
   

## Usage

```
  1、对fasta格式作图。(有可选参数[指定为中心序列的名字][-number 可容忍的gap数量，默认为0])
     1)  path.fast [-0]                        (默认0号为中心序列，与文件中其他序列分别做图)
     2)  path.fast [center_seq_name] [-0]      (指定中心序列，与文件中其他序列分别做图)
     3)  path.fast seq1_name seq2_name [-0]    (指定文件中该两条作图)
  2、对maf格式作图。
     1)  path.maf center_seq_name              (指定中心序列，与文件中其他序列分别做图)
     2)  path.maf seq1_name seq2_name          (指定文件中该两条作图)
```

## Dependence

python3.7
matplotplib
