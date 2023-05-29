This folder is used for the patching of the counts for the DUX4 gene in the raw output of HTSeq-count.
Four parameters are required in total.
1. -f: bed file for the locations of DUX4 genes. Please use the 'Homo_sapiens.GRCh38.V102.withChr.Exon.DUX4.bed' file provided in this folder specifically for the analysis on human samples using GRCh38 as the reference genome.
2. -b: bam file from the alignment.
3. -h: HTSeqcount raw output.
4. -o: the output file name.

Example: 
HTSeqDUX4patch.pl -f Homo_sapiens.GRCh38.V102.withChr.Exon.DUX4.bed -b sample.bam -h sample.HTSeq -o sample.HTSeq.DUX4patch
