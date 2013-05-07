#iqseq

De novo identification and quantification of sequence data.

```
usage: iqseq.py [-h] [-v] {consensus,quantify} ...

De novo identification and quantification of sequence data.

positional arguments:
  {consensus,quantify}  commands
    quantify            quantify unique and similar sequences
    consensus           build observed sequence library

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
```

##Quantifying sequences in FASTQ

```
usage: iqseq.py quantify [-h] [-c CUTOFF] [-m MISMATCH] FASTQ

Find and quantify unique and similar sequences within a FASTQ.

positional arguments:
  FASTQ        reads to process

optional arguments:
  -h, --help   show this help message and exit
  -c CUTOFF    minimum allowable seq length [18]
  -m MISMATCH  mismatch tolerance when grouping bins [3]
```

##Finding consensus of observed sequences

```
usage: iqseq.py consensus [-h] [-c CUTOFF] [-m MISMATCH] BINS [BINS ...]

Build consensus of sequences across all samples.

positional arguments:
  BINS         results of `quantify`

optional arguments:
  -h, --help   show this help message and exit
  -c CUTOFF    minimum allowable count [10]
  -m MISMATCH  mismatch tolerance when grouping bins [3]
```
