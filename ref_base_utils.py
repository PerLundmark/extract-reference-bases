"""
Sequence / biopython utility classes and functions
"""

from Bio import SeqIO
from abc import ABC, abstractmethod
import csv

def build_seq_dict(fasta_filename):
    """Builds a dict of Biopython SeqRecord objects from a fasta file (kept in memory)

    Args:
        fasta_filename (str): The file name of the (multi sequence) fasta file to parse.

    Returns:
        A list of tuples with (id, name, length) of sequence records, and a dict of SeqRecord objects
    """
    records = list(SeqIO.parse(fasta_filename, "fasta"))
    records_info = []
    records_dict = {}

    for record in records:
        info_record = (record.id, len(record))
        records_info.append(info_record)
        records_dict[record.id] = record

    return (records_info, records_dict)

def build_indexed_seq_dict(fasta_filename):
    """Builds a dict-like object of Biopython SeqRecord objects from a fasta file (indexed on disk)

    Args:
        fasta_filename (str): The file name of the (multi sequence) fasta file to parse.

    Returns:
        A dict of indexed SeqRecord objects
    """
    return SeqIO.index(fasta_filename, "fasta")

class GenotypingManifest(ABC):
    """Abstract base class for genotyping manifests"""
    def __init__(self, filename):
        self.filename = filename
        self.marker_list = []
        self.chr_order = []

    @abstractmethod
    def read_manifest(self, filename):
        pass

    @abstractmethod
    def get_markerlist(self):
        pass

    @abstractmethod
    def set_chr_order(self, chr_list):
        pass

    @abstractmethod
    def sort_markerlist(self):
        pass

class InfiniumManifest(GenotypingManifest):
    """A Class for Illumina Infinium manifests (csv-format)"""
    def read_manifest(self):
        try:
            with open(self.filename, newline='') as infinium_manifest:
                reader = csv.reader(infinium_manifest, delimiter=',')
                while (next(reader)[0] != "[Assay]"): # Skip header
                    pass

                header = next(reader) #Check positions of fields in manifest (variable order)
                idx_snp_name = header.index('Name')
                idx_snp_chr = header.index('Chr')
                idx_snp_pos = header.index('MapInfo')
                idx_snp_var = header.index('SNP')
                idx_snp_ref_strand = header.index('RefStrand')

                manifest_indices = [idx_snp_name, idx_snp_chr, idx_snp_pos, idx_snp_var, idx_snp_ref_strand]
                self.marker_list = []
                for row in reader:
                    if(row[0] == "[Controls]"): #End of SNP data
                        break
                    self.marker_list.append([row[i] for i in manifest_indices])

        except IOError:
            print("Could not read file: ", self.filename)


    def get_markerlist(self):
        return self.marker_list

    def set_chr_order(self, chr_list):
        self.chr_order = chr_list

    def get_sort_keys(self, item):
        try:
            chr_sort_key = self.chr_order.index(item(1))
        except ValueError:
            chr_sort_key = len(self.chr_order) #Sort 0/unknown/strange chromosome entries behind the rest

        return (chr_sort_key, item[2], item[0])

    def sort_markerlist(self):
#        import functools
#        key = functools.partial(self.get_sort_key, data=something)
#        self.marker_list.sort(key=key)
        self.marker_list.sort(key=self.get_sort_keys)


class AgenaManifest(GenotypingManifest):
    def read_manifest(self, filename):
        pass

    def get_markerlist(self):
        pass
