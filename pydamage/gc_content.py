from Bio import SeqIO
from Bio.SeqUtils import GC


class fastaFile:
    def __init__(self, fasta):
        """Class representing the fasta file

        Args:
            fasta (str): Path to fasta file
        """
        self.fasta = {}
<<<<<<< HEAD
        for record in SeqIO.parse(fasta_file, "fasta"):
=======
        for record in SeqIO.parse(fasta, "fasta"):
>>>>>>> 204d4fca36b046450dcecc7462719d37557d69ad
            self.fasta[record.id] = record

    def compute_gc(self):
        """Compute GC content

        Args:
            reference (str): name of the reference sequence

        Returns:
            dict: {reference_id : gc_content}
        """
        self.gc = {}
        for ref in self.fasta:
            self.gc[ref] = GC(self.fasta[ref].seq) / 100
        return self.gc
