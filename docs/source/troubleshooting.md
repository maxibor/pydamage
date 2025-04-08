# Troubleshooting

## My alignment files don't have a MD tag

You can use [samtools calmd](http://www.htslib.org/doc/samtools-calmd.html) to set the MD tag

Example:

```bash
samtools calmd -b alignment.bam reference.fasta > aln.bam
```

## I have filtered my BAM file to only retain a few references but pyDamage is still iterating across all references

By default, pyDamage will iterate over every reference that is present in the `SQ` section of the BAM file header and will evaluate whether there are any reads aligned to it. This process is time-consuming particularly when you have pre-filtered your BAM file and many references have no aligned reads.

To speed up the damage analysis, you can use the program [bamAlignCleaner](https://github.com/maxibor/bamAlignCleaner). It will iterate over the BAM file and identify references to which no sequences were aligned. These references will be pruned from the header of the BAM file. By providing a BAM file with only references that contain reads, pyDamage will run faster.
