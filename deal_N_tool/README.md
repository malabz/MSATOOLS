# cpp为一个处理N的C++脚本程序

（功能：比对前删除N,比对后插回N）
（前提：保证比对前后序列名称不变）
（容许名称带空格，容许交换序列顺序）

## Usage

```
一共有以下3个功能：1+2组合=3
  1、删除N并记录。  
      (试用情景：比对前)
     path1:i-带N的比对前数据.fasta
     path2:o-去N的比对前数据.fasta
     path3:o-记录N的位置信息.tmp
     Win  ：xxx.exe -1 in_original.fasta out_removeN.fasta out_recordN.tmp  
     Linux：./xxx     -1 in_original.fasta out_removeN.fasta out_recordN.tmp  
  
  2、将N插入到结果文件中。 
      (试用情景：比对后)
     path1:i-不带N的比对后结果.fasta
     path2:o-添加N的比对后结果.fasta
     path3:i-记录N的位置信息.tmp 
     Win  ：xxx.exe -2 in_ans_withoutN.fasta out_ans_withN.fasta in_recordN.tmp  
     Linux：./xxx     -2 in_ans_withoutN.fasta out_ans_withN.fasta in_recordN.tmp  
  
  3、从源文件和结果文件中生成带有N的最终结果文件。 
      (试用情景：你自己程序中输入时删除了N，后续没管)
     path1:i-带N的比对前数据.fasta
     path2:i-不带N的比对后结果.fasta
     path3:o-添加N的比对后结果.fasta
     Win  ：xxx.exe -3 in_original.fasta in_ans_withoutN.fasta out_ans_withN.fasta  
     Linux：./xxx     -3 in_original.fasta in_ans_withoutN.fasta out_ans_withN.fasta  

-h 可以打开帮助
```

## toN.py

toN.py 可以把fasta文件中所有简并碱基都转成N，因为带N的序列数据大多软件可以处理，但带其他简并碱基的序列数据可能导致大多软件瘫痪。

usage： python toN.py -i in.fasta -o out.fasta

## Build

g++


