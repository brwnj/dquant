#iqseq

De novo identification and quantification of sequence data.

```
usage: dquant.py [-h] [-v] {consensus,quantify} ...

De novo identification and quantification of sequence data.

positional arguments:
  {consensus,quantify}  commands
    quantify
    consensus

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
```

##Quantifying sequences in FASTQ

```
usage: dquant.py quantify [-h] [-m MISMATCH] FASTQ

Find and quantify unique and similar sequences within a FASTQ.

positional arguments:
  FASTQ        reads to process

optional arguments:
  -h, --help   show this help message and exit
  -m MISMATCH  mismatch tolerance when grouping bins [3]
```
