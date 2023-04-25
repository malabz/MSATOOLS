# 🔧Format conversion tool for bioinformatics

## 1. xmfa2maf
Convert the file from xmfa format to maf format
```
# Compile
  g++ xmfa2maf.cpp -o xmfa2maf
# Run
  ./xmfa2maf -i in.xmfa -o out.maf
```

## 2. aln2fasta
Convert the file from aln format to fasta/fas format
```
# single file
python3 aln2fasta.py PATH/file

# multiple files (based on aln2fasta.py)
sh aln2fasta.sh PATH
```
