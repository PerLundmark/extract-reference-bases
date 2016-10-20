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
header_lines.append("##fileDate={:%Y%m%d}".format(current_date))
header_lines.append("##source=extract_reference_bases.py_and_chiasma")
header_lines.append("##reference=file://" + args.genome)
header_lines.append("##phasing=none")

(seq_info_list, seq_dict) = ref_base_utils.build_seq_dict(args.genome)

for (current_id, current_length) in seq_info_list:
    header_lines.append('##contig=<ID=' + current_id + ',length=' + str(current_length) + ',assembly=' + assembly + '>')

header_lines.append("##FORMAT=<ID=GT,Number=1,Type=Integer,Description='Genotype'>")
header_lines.append("#refbase_data")
header_lines.append("\t".join(['snp_name', 'snp_chr', 'snp_pos', 'snp_var', 'snp_ref_strand', 'ref_name', 'ref_comment', 'ref_base', 'ref_length']))

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

#Retrieve reference data & print output
del(manifest)
try:
    with open(args.output, "w") as output_file:
        for line in header_lines:
            output_file.write(line + '\n')

        marker_counter = 1
        for marker in markers:
            (snp_name, snp_chr, snp_pos, snp_var, snp_ref_strand) = marker
            (clean_snp_chr, clean_snp_pos) = ref_base_utils.clean_chr_and_pos(snp_chr, snp_pos)

            (ref_name, ref_comment, ref_base, ref_length) = ref_base_utils.get_reference_info(seq_dict, clean_snp_chr, clean_snp_pos)

            #output_file.write("\t".join(marker + list(ref_base_utils.get_reference_info(seq_dict, snp_chr, int(snp_pos))))) #fix type of position at build
            output_file.write("\t".join([snp_name, clean_snp_chr, str(clean_snp_pos), snp_var, snp_ref_strand, ref_name, ref_comment, ref_base, str(ref_length)]) + '\n')

            if (marker_counter % 1000 == 0):
                print (marker_counter)

            marker_counter += 1

except IOError:
    print("Unable to write to outputfile", args.output)


print("Hepp")