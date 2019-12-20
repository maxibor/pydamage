import damage
import pysam


bampath = "/Users/borry/Documents/GitHub/pydamage/tests/data/aligned.bam"
bam = pysam.AlignmentFile(bampath, "rb")

references = bam.references
print(references)
wlen = 50
min_al = 10
show_al = False
for ref in references:
    damage.test_ct(bam=bam, ref=ref, wlen=wlen, show_al=show_al, min_al=min_al)
