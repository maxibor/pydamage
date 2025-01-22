import pysam
def parse_fasta(fasta):
    fa_dict = {}
    with pysam.FastxFile(fasta) as fh:
        for entry in fh:
            fa_dict[entry.name] = entry.sequence
    return fa_dict