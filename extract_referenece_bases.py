#!/usr/bin/env python

import argparse
import ref_base_utils
from datetime import date

parser = argparse.ArgumentParser(description="Generate a reference base information file for Chiasma VCF exports (Illumina manifest and Agena/general format supported)")
parser.add_argument('-f', '--format', required=True, choices=['infinium', 'agena'], help="Format of input manifest, either 'illumina' for Illumina Infinium, or 'agena' for a mart export for Agena/general use.")
parser.add_argument('-i', '--input', required=True, help="Input filename.")
parser.add_argument('-o', '--output', required=True, help= "Output filename.")
parser.add_argument('-g', '--genome', required=True, help= "The source genome for reference base extraction. Fasta format.")
args = parser.parse_args()

current_date = date.today()

header_lines = []
header_lines.append("##fileformat=VCFv4.1")
header_lines.append("fileDate={:%Y%m%d}".format(current_date))
header_lines.append("##source=extract_reference_bases.pl_and_chiasma")
header_lines.append("##reference=file://" + args.genome)
header_lines.append("##phasing=none")

sequences = ref_base_utils.build_seq_dict(args.genome)

