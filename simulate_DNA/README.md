# 根据一条原始序列，生成带有各种结构变异和snp的模拟数据

 * Author: ZhouTong                                                                             
 * Date: 2023/7/20  
## Introduction

```
 * 1函数：生成小插入_1_snp_insert(length):
 * 2函数：生成小删除_2_snp_delete(length):
 * 3函数：生成小snp替换_3_snp_tihuan(length):

* 以下为8种结构变异SV：基础类型：重复，易位(仅染色体内)，倒位，以及两两组合3种 + 大的插入删除
 * 4函数：生成大插入sv_4_big_insert(length):
 * 5函数：生成大删除sv_5_big_delete(length):
 * 6函数：生成串联重复sv_6_duplication(length, times):
 * 61函数：生成散布(易位)重复sv_61_tran_duplication(length):
 * 62函数：生成倒位重复sv_62_inver_duplication(length):
 * 7函数：生成易位sv_7_translocation(length):
 * 8函数：生成倒位sv_8_inversion(length):
 * 9函数：生成倒位易位sv_9_tran_inversion(length):                                                         
```

## Usage

```
 * 两种使用方法：(修改主函数对应字典)
 * 1 通过变异率+长度     generate_maf_fasta
 * 2 通过变异次数+长度   generate_maf_fasta_times
```

## 输入：一个原始DNA的fasta文件
## 输出：一个模拟DNA的fasta文件+一个展示【重复，易位(仅染色体内)，倒位，以及两两组合3种，共6种】的标准答案sv.maf