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
        print ("Reading ", record.id)
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


def get_reference_info(seq_dict, sequence_name, position):
    """ Returns the reference base in a specific sequence and position (in general chromosome and position from a
    reference genome)

    Args:
        seq_dict (dict): The dict of SeqRecords keyed on sequence id.
        sequence_name (str): The name of the sequence/chromosome.
        position (int): The base position to extract.

    Returns:
        A tuple with (name(str), comment(str), reference base(str), reference length(int)) with info from the reference sequence.
    """
    if(sequence_name == "Unknown" or sequence_name == "0"):
        ref_name = "Unknown"
        ref_description = "Chr given as Unknown or 0"
        ref_base = "N"
        ref_length = -1
    elif(position == 0):
        ref_name = seq_dict[sequence_name].name
        ref_description = "No position in manifest"
        ref_base = "N"
        ref_length = len(seq_dict[sequence_name].seq)
    elif(position > len(seq_dict[sequence_name].seq)):
        ref_name = seq_dict[sequence_name].name
        ref_description = "Position in manifest outside reference sequence"
        ref_base = "N"
        ref_length = len(seq_dict[sequence_name].seq)
    else:
        ref_name = seq_dict[sequence_name].name
        ref_description = seq_dict[sequence_name].description
        ref_base = str(seq_dict[sequence_name].seq[position-1:position])
        ref_length = len(seq_dict[sequence_name].seq)

    return (ref_name, ref_description, ref_base, ref_length)

def clean_chr_and_pos(chr, pos):
    """Cleans up the heterogenous Illumina manifest chromsome and position information.
    XY chromosomes are changed to XY since the positions given in the pseudoautosomal region is X chr positions.
    For missing data of different kinds, the function returns 0 for positions and Unknown for chromosomes

    Args:
        chr(str): The SNP chromosome
        pos(int): The chromsome position of the SNP

    Returns:
        A tuple with cleaned up (chromosome, position)
        """

    if (chr == "XY"):
        chr = "X"
    elif (chr == 0):
        chr = "Unknown"

    if (pos == "N/A"):
        pos = 0

    return (chr, pos)


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

                #manifest_indices = [idx_snp_name, idx_snp_chr, idx_snp_pos, idx_snp_var, idx_snp_ref_strand]
                self.marker_list = []
                for row in reader:
                    if(row[0] == "[Controls]"): #End of SNP data
                        break

                    self.marker_list.append([row[idx_snp_name], row[idx_snp_chr], int(row[idx_snp_pos]), row[idx_snp_var], row[idx_snp_ref_strand]])
                    #self.marker_list.append([row[i] for i in manifest_indices])

        except IOError:
            print("Could not read file: ", self.filename)


    def get_markerlist(self):
        return self.marker_list

    def set_chr_order(self, chr_list):
        self.chr_order = chr_list

    def get_sort_keys(self, item):
        try:
            if (item[1] == "XY"): #Sort XY (pseudoautosomal region) as X since that is the base of the coordinates
                item[1] = "X"

            chr_sort_key = self.chr_order.index(item[1])
        except ValueError:
            chr_sort_key = len(self.chr_order) #Sort 0/unknown/strange chromosome entries behind the rest

        return (chr_sort_key, item[2], item[0])

    def sort_markerlist(self):
        self.marker_list.sort(key=self.get_sort_keys)


class AgenaManifest(GenotypingManifest):
    def read_manifest(self):
        try:
            with open(self.filename, newline='') as agena_manifest:
                for line in agena_manifest:
                    if (line[0] == ">"):
                        data = line[1:]
                        header_fields = data.split("|")
                        chr =    header_fields[0]
                        start =  header_fields[1]
                        stop =   header_fields[2]
                        name =   header_fields[3]
                        source = header_fields[4]
                        strand = header_fields[5]
                        var =    header_fields[6]

                        if(strand == "1"):
                            signed_strand = "+"
                        elif(strand == "-1"):
                            signed_strand = "-"
                        else:
                            signed_strand = "NA"

                        if (start != stop):
                            raise ValueError("Start and stop different, variation length > 1 base, start:" + start + " stop:" + stop)
                        self.marker_list.append([name, chr, int(start), var, signed_strand])

        except IOError:
            print("Could not read file: ", self.filename)

        except ValueError as err:
            print(err.message)


    def get_markerlist(self):
        return self.marker_list

    def set_chr_order(self, chr_list):
        self.chr_order = chr_list

    def get_sort_keys(self, item):
        try:
            chr_sort_key = self.chr_order.index(item[1])
        except ValueError:
            chr_sort_key = len(self.chr_order) #Sort 0/unknown/strange chromosome entries behind the rest

        return (chr_sort_key, item[2], item[0])

    def sort_markerlist(self):
        self.marker_list.sort(key=self.get_sort_keys)
