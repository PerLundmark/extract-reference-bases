"""
Sequence / biopython utility functions
"""

from Bio import SeqIO

def build_seq_dict(fasta_filename):
    """Builds a dict of Biopython SeqRecord objects from a fasta file (kept in memory)

    Args:
        fasta_filename (str): The file name of the (multi sequence) fasta file to parse.

    Returns:
        A dict of SeqRecord objects
    """
    return SeqIO.to_dict(SeqIO.parse(fasta_filename, "fasta"))


def build_indexed_seq_dict(fasta_filename):
    """Builds a dict-like object of Biopython SeqRecord objects from a fasta file (indexed on disk)

    Args:
        fasta_filename (str): The file name of the (multi sequence) fasta file to parse.

    Returns:
        A dict of indexed SeqRecord objects
    """
    return SeqIO.index(fasta_filename, "fasta")