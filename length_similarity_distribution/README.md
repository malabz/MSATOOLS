# length and similarity distribution of dataset
## Introduction
Used to calculate the SP value of the result of multiple sequence alignment.

The script automatically replaces the base "U" with "T" in the RNA sequence files and handles illegal characters (characters other than ATCG-) in the file. The solution is to replace illegal characters with "N". **It could only be used to nucleotide sequences.**

By default, pairs with the same nongap characters(match) receive one count, pairs with different nongap characters(mismatch) receive minus one, pairs with one gap(gap1) receive minus two, and pairs with two gaps or consisting of an N and a character(ACTGN) receive zero. **Support user-defined parameters.**

The results:
- SP score: the sum of pair score
- Avg SP score: we divided the SP score by the pair number N * (N – 1) / 2, where N represents the number of sequences.
- Scaled SP score: we divided the Avg SP score by the sequences length L.
## Usage
```bash
usage: python length_similarity_distribution.py inputfile_name mafft_path output_prefix mafft_threads
 
options:
    inputfile_name: specify the input file
        mafft_path: the path of executable file of mafft aligner
     output_prefix: specify the prefix of output files
     mafft_threads: multi-thread for running mafft alignment
```


