# length and similarity distribution of dataset
## Introduction
This python script is used for visualize the distributions of length and similarity given a sequence dataset in FASTA format.
The length is the number of characters of each sequence. The number of length value equals is the number of sequences (n) in dataset.
The similarity between any two sequences is percentage of matched characters to in their pairwise alignment performed by aligner mafft. The number of similarity value equals n(n-1)/2. n is the number of sequences in dataset.
The results:
•	output_prefix_Distri.png: the visualization of length distribution.
•	similarity_output_prefix.png: the visualization of similarity distribution.
•	output_prefix.txt: the statistics of all sequence length in input dataset.
•	similarity_output_prefix.txt: n(n-1)/2 similarity values between any two sequences
•	output_prefix_time.txt: the running time spent by this script.


The results:
- SP score: the sum of pair score
- Avg SP score: we divided the SP score by the pair number N * (N – 1) / 2, where N represents the number of sequences.
- Scaled SP score: we divided the Avg SP score by the sequences length L.
## Usage
```bash
usage: python length_similarity_distribution.py inputfile_name mafft_path output_prefix mafft_threads
 
positional arguments:
  inputfile_name    specify the input file
  mafft_path        the path of executable file of mafft aligner
  output_prefix     specify the prefix of output files
  mafft_threads     multi-thread for running mafft alignment
 
```


