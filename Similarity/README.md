# 序列相似性的比较
(如您需要可以直接带走.h)
 * Author: ZhouTong                                                                             
 * Date: 2023/4/25  
## Usage

```
 * 1.将序列 seq_center 拆分为长度为 k 的子串，并计算每个子串的哈希值，存储到 Bloom Filter 中；  
 * 2.遍历序列 seq_i，将其拆分为长度为 k 的子串，并计算每个子串的哈希值，                        
 *   然后在 Bloom Filter 中查找是否存在相同的哈希值，若存在，则说明存在相同的子串，统计数量；   
 * 3.比较统计出的相同子串数量。                                                                 
 * 4.这里使用的 MurmurHash 算法计算哈希值，可以快速高效地生成哈希值。而 Bloom Filter            
 *   则可以用较小的内存空间对大规模数据进行快速查询和去重。                                     
 * 5.总体来说，该算法属于基于哈希的近似匹配算法，适用于大规模数据集的相似度查询问题。           

 * 注：返回数值意义不大，只作比较用                                                             
```