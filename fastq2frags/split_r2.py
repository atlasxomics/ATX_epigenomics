""" split R2 into new-R2, R3 for downstream processing
"""

import argparse
import logging

from Bio.SeqIO.QualityIO import FastqGeneralIterator
from gzip import open as gzopen

logging.basicConfig(
    format="%(levelname)s - %(asctime)s - %(message)s", level=logging.INFO
)

seq_start = 117
bc2_start, bc2_end = 22, 30
bc1_start, bc1_end = 60, 68

ap = argparse.ArgumentParser()
ap.add_argument("-i", "--input", required=True, help="input file")
ap.add_argument("-o1", "--output_R3", required=True, help="output file R3")
ap.add_argument("-o2", "--output_R2", required=True, help="output file R2")
args = vars(ap.parse_args())
logging.info(f"Arguemnts: {args}")

input_path = args["input"]
output_pathR3 = args["output_R3"]
output_pathR2 = args["output_R2"]

logging.info("R2 parsing initiated...")
with gzopen(input_path, "rt") as in_R2, \
    open(output_pathR3, "w") as out_R3, \
	open(output_pathR2, "w") as out_R2:
    for title, seq, qual in FastqGeneralIterator(in_R2):

        R3_seq = seq[seq_start:]
        R3_qual = qual[seq_start:]
        R2_seq = seq[bc2_start:bc2_end] + seq[bc1_start:bc1_end] 
        R2_qual = qual[bc2_start:bc2_end] + qual[bc1_start:bc1_end]        

        out_R3.write("@%s\n%s\n+\n%s\n" % (title, R3_seq, R3_qual))
        out_R2.write("@%s\n%s\n+\n%s\n" % (title, R2_seq, R2_qual))

logging.info("R2 parsing completed.")
