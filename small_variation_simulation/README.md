# Small variation simulation of DNA sequence

This script for creating a simulated multi-FASTA alignment file is adopted by [Fenglong Yang](http://lab.malab.cn/~yangfl/) from the script shared by [SNP-sites](https://github.com/sanger-pathogens/snp-sites). The DNA center sequences were simulated randomly with 25% of A, C, T and G. 

Download <a href="http://lab.malab.cn/~tfr/" download="small_variation_simulation_splice.pl">small_variation_simulation_splice.pl</a> 
```
wget http://lab.malab.cn/~tfr/

Usage:

```
perl small_variation_simulation_splice.pl [options]
```

```
Options: 
         -s INT      Number of samples
         -l INT      Number of bases in each genome
	 -sub INT    Number of substitution sites in the alignment
	 -ins INT    Number of insertion sites in the alignment
	 -del INT    Number of deletion sites in the alignment
         -o STR      Output file name
         -h          this help message

```

Example: create an alignment with 20 samples, 1000 bases, 100 substitution sites, 10 insertion sites and 10 deletion sites in each genome：

```
perl small_variation_simulation_splice.pl -o output.aln -s 20 -l 1000 -sub 10 - ins 10 -del 10
```



