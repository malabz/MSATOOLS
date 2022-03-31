# SP(sum of pair) score
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
python SPscore.py [-h] --input INPUT [--match MATCH] [--mismatch MISMATCH] [--gap1 GAP1] [--gap2 GAP2]

optional arguments:
  -h, --help           show this help message and exit
  --input INPUT        Fasta file path to be scored.
  --match MATCH        match score, default=1.
  --mismatch MISMATCH  mismatch score, default=-1.
  --gap1 GAP1          gap-base score, default=-2.
  --gap2 GAP2          gap-gap score, default=0.
```



