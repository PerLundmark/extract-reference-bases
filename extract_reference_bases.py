#!/usr/bin/env python
"""
A script that extracts reference genome bases for
genotyped positions found in the manifest of a genotyping chip
or similar.

"""
import argparse
import ref_base_utils
from datetime import date

assembly = "b37"

parser = argparse.ArgumentParser(description="Generate a reference base information file for Chiasma VCF exports (Illumina manifest and Agena/general format supported)")
parser.add_argument('-f', '--format', required=True, choices=['infinium', 'agena'], help="Format of input manifest, either 'illumina' for Illumina Infinium, or 'agena' for a mart export for Agena/general use.")
parser.add_argument('-i', '--input', required=True, help="Input filename.")
parser.add_argument('-o', '--output', required=True, help= "Output filename.")
parser.add_argument('-g', '--genome', required=True, help= "The source genome for reference base extraction. Fasta format.")
parser.add_argument('-q', '--quiet', help= "Supress progress printing")

args = parser.parse_args()
current_date = date.today()

#Build VCF header
header_lines = []
header_lines.append("##fileformat=VCFv4.1")
header_lines.append("fileDate={:%Y%m%d}".format(current_date))
header_lines.append("##source=extract_reference_bases.pl_and_chiasma")
header_lines.append("##reference=file://" + args.genome)
header_lines.append("##phasing=none")

(seq_info_list, seq_dict) = ref_base_utils.build_seq_dict(args.genome)

for (current_id, current_length) in seq_info_list:
    header_lines.append('##contig=<ID=' + current_id + ',length=' + str(current_length) + ',assembly=' + assembly + '>')

header_lines.append("##FORMAT=<ID=GT,Number=1,Type=Integer,Description='Genotype'>")

#Read manifest
if(args.format == 'infinium'):
    manifest = ref_base_utils.InfiniumManifest(args.input)
elif(args.format == 'agena'):
    manifest = ref_base_utils.AgenaManifest(args.input)
else:
    raise ValueError("Unknown manifest format")

manifest.read_manifest()
manifest.set_chr_order([row[0] for row in seq_info_list])   #[row[i] for i in manifest_indices]
manifest.sort_markerlist()
markers = manifest.get_markerlist()

#Retrieve reference data
manifest = None


print("Hepp")