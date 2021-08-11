#!/usr/bin/env python3

import sys

CDSFasta = open(sys.argv[1], 'r')
EnsToGene = open(sys.argv[2], 'r')

LineNumber=0

# For each line of the CDS Fasta, go through every line of the EnsToGene file and see if there is a transcript that can be replaced by the gene name

for Line in CDSFasta:
	Line=Line.strip('\n')
	for TransRow in EnsToGene:
		TransRow = TransRow.strip('\n')
		Transcript = TransRow.split('\t')
		Line=Line.replace(Transcript[0], Transcript[1] + "_" + Transcript[0])
	print(Line)
	EnsToGene.seek(0)

CDSFasta.close()
EnsToGene.close()
