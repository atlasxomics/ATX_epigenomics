########################################################
# Process R2 for cellranger-atac pipeline
# revised Dec 2022 for spatial ATAC v2
########################################################

import argparse
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from gzip import open as gzopen

seq_start = 117
bc2_start, bc2_end = 22, 30
bc1_start, bc1_end = 60, 68

ap = argparse.ArgumentParser()
ap.add_argument("-i", "--input",
	required=True,
	help="input file")
ap.add_argument("-o1", "--output_R3",
	required=True,
	help="output file R3")
ap.add_argument("-o2",
	"--output_R2",
	required=True,
	help="output file R2")
args = vars(ap.parse_args())

input_file_R2 = args["input"]
output_file_R3 = args["output_R3"]
output_file_R2 = args["output_R2"]

with gzopen(input_file_R2, "rt") as in_handle_R2, \
	open(output_file_R3, "w") as out_handle_R3, \
	open(output_file_R2, "w") as out_handle_R2:
    for title, seq, qual in FastqGeneralIterator(in_handle_R2):
        new_seq_R3 = seq[seq_start:]
        new_qual_R3 = qual[seq_start:]
        barcode = seq[bc2_start:bc2_end] + seq[bc1_start:bc1_end] # !!! BC2 + BC1
        new_qual_R2 = qual[bc2_start:bc2_end] + qual[bc1_start:bc1_end]        
        out_handle_R3.write("@%s\n%s\n+\n%s\n" % (title, new_seq_R3, new_qual_R3))
        out_handle_R2.write("@%s\n%s\n+\n%s\n" % (title, barcode, new_qual_R2))

