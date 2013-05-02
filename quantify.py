#!/usr/bin/env python
# encoding: utf-8
"""
Finds and quantifies unique and similar sequences among a FASTQ.
"""
import sys
import tools as t
from toolshed import nopen
from collections import Counter
        
def process_exact(fastq):
    """Group identical reads using a Counter. Returns Counter."""
    c = Counter()
    with nopen(fastq) as fh:
        for name, seq, qual in t.read_fastq(fh):
            c.update([seq])
    return c

def main(args):
    print >>sys.stderr, ">> collapsing identical sequences"
    reads = process_exact(args.fastq)
    print >>sys.stderr, ">> collapsing identical subsequences"
    reads = t.process_exact_substring(reads)
    print >>sys.stderr, ">> collapsing similar sequences"
    reads = t.process_similar(reads, args.mismatch)
    for seq, count in reads.iteritems():
        print "%s\t%d" % (seq, count)

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("fastq", metavar="FASTQ", help="reads to process")
    p.add_argument("-m", dest="mismatch", type=int, default=3,
            help="number of mismatches allowed when combining bins [%(default)s]")
    args = p.parse_args()
    main(args)
