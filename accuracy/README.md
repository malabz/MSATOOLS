# SP(sum of pair) score
## Introduction
Used to calculate the SP value of the result of multiple sequence alignment.

By default, pairs with the same nongap characters(match) receive one count, pairs with different nongap characters(mismatch) receive minus one, pairs with one gap(gap1) receive minus two, and pairs with two gaps(gap2) receive zero. **Support user-defined parameters.**
## Usage
```bash
python SP.py [-h] --input INPUT [--match MATCH] [--mismatch MISMATCH] [--gap1 GAP1] [--gap2 GAP2]

optional arguments:
  -h, --help           show this help message and exit
  --input INPUT        Fasta file path to be scored.
  --match MATCH        match score, default=1.
  --mismatch MISMATCH  mismatch score, default=-1.
  --gap1 GAP1          gap-base score, default=-2.
  --gap2 GAP2          gap-gap score, default=0.
```



