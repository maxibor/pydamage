# Troubleshooting

## My alignment files don't have a MD tag

You can use [samtools calmd](http://www.htslib.org/doc/samtools-calmd.html) to set the MD tag

Example:

```bash
samtools calmd -b alignment.bam reference.fasta > aln.bam
```
