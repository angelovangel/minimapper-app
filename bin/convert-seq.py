#!/usr/bin/env python3

# Convert various formats to fasta

import sys
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


in_file = sys.argv[1]
in_format = sys.argv[2]
out_file = sys.argv[3]
out_format = sys.argv[4]


try:
    count = SeqIO.convert(in_file, in_format, out_file, out_format)
    print(f"Successfully converted {count} record(s) from {in_format} to {out_format}. Output: {out_file}")
except Exception as e:
    print(f"Error during conversion: {e}")
    sys.exit(1)
