# PreProcess
## Introduction
Used to process fasta files after multiple sequence alignment. 

The script replaces the base "U" with "T" in the RNA sequence files and replaces illegal characters (such as ":") with "N". 

The processed files are stored in resultName_processed.fasta
## Usage
```bash
python preprocess.py filepath resultName
```
The parameters are described as follows:
- filepath: Path of fasta file to be processed
- resultName: Fasta file name after processing


