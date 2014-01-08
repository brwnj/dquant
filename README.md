#Sequence clustering pipeline

De novo identification and quantification of sequence data utilizing edit distance
as a metric to group similar sequences. Sequences are paired longest to shortest.

There are many algorithms to accomplish a similar goal and many that better
account for sequence identity matching without the bias using longer sequences
first.

```
usage: sequence_clustering.py [-h] [-v] {consensus,quantify,matrix} ...

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
usage: sequence_clustering.py quantify [-h] [-c CUTOFF] [-m MISMATCH] FASTQ

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
usage: sequence_clustering.py consensus [-h] [-c CUTOFF] [-m MISMATCH]
                                        BINS [BINS ...]

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
usage: sequence_clustering.py matrix [-h] [-c CUTOFF] [-m MISMATCH]
                                     [-n {deseq,totalcount}]
                                     CONSENSUS COUNTS [COUNTS ...]

Generate counts matrix

positional arguments:
  CONSENSUS             result of `consensus`
  COUNTS                results of `quantify`

optional arguments:
  -h, --help            show this help message and exit
  -c CUTOFF             minimum allowable count for individual sample
                        sequences (default: 100)
  -m MISMATCH           mismatch tolerance when grouping bins (default: 3)
  -n {deseq,totalcount}
                        output normalized table using either DESeq or total
                        count method (default: None)
```
