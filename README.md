#Sequence clustering pipeline

De novo identification and quantification of sequence data utilizing sequence
clustering algorithms SEED and CD-HIT (to be implemented soon).

```
usage: iqseq.py [-h] [-v] {consensus,quantify,matrix} ...

De novo identification and quantification of sequence data.

positional arguments:
  {consensus,quantify,matrix}
                        commands
    quantify            quantify unique and similar sequences
    consensus           build observed sequence library
    matrix              generate counts matrix

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
  -c CUTOFF    minimum allowable seq length (default: 18)
  -m MISMATCH  mismatch tolerance when grouping bins (default: 3)
```

##Finding consensus of observed sequences

```
usage: iqseq.py consensus [-h] [-c CUTOFF] [-m MISMATCH] BINS [BINS ...]

Build consensus of sequences across all samples.

positional arguments:
  BINS         results of `quantify`

optional arguments:
  -h, --help   show this help message and exit
  -c CUTOFF    minimum allowable count (default: 100)
  -m MISMATCH  mismatch tolerance when grouping bins (default: 3)
```

##Counts of observed sequences across consensus bins

```
usage: iqseq.py matrix [-h] [-c CUTOFF] [-m MISMATCH] [-n]
                       CONSENSUS COUNTS [COUNTS ...]

Generate counts matrix

positional arguments:
  CONSENSUS    result of `consensus`
  COUNTS       results of `quantify`

optional arguments:
  -h, --help   show this help message and exit
  -c CUTOFF    minimum allowable count for individual sample sequences
               (default: 100)
  -m MISMATCH  mismatch tolerance when grouping bins (default: 3)
  -n           output scaling factor normalized table using method developed
               by Anders and Huber for DESeq:
               http://genomebiology.com/2010/11/10/R106 (default: False)
```