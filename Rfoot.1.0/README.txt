Rfoot: Transcriptome-scale identification of RNA-protein complexes from ribosome profiling data

Contact: Zhe Ji (zhe.ji@northwestern.edu)

Materials:
Ribosome profiling datasets in Fastq format;
Genome assembly file in Fasta format;
Ribosomal RNA (rRNA) sequence file in Fasta format;
Transcript definition file in genePred format;
Linux high performance computing cluster;
Perl program installation;
Read mapping software (such as Bowtie (Langmead and Salzberg, 2012) and Tophat (Kim et al., 2013))

1. Download Rfoot package from https://github.com/zhejilab/Rfoot/.

2. Obtain the ribosome profiling dataset with cycloheximide treatment or without drug treatment, trim 3’ adapters of sequencing reads, map trimmed reads to rRNAs, and then align non-rRNA reads to the reference transcriptome and genome.

Example commands:
2a. Download an example ribosome profiling dataset from GEO database using fastq-dump command line as a part of the NIH software sratoolkit, https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/).
fastq-dump -Z SRR1802146 > SRR1802146.fastq
2b. Remove 3’ adapters of sequencing reads. Obtain “removeAdapter.pl” from https://github.com/zhejilab/RibORF.
perl removeAdapter.pl -f SRR1802146.fastq -a CTGTAGGCAC -o adapter.SRR1802146.fastq
2c. Get human rRNA sequences from NCBI database, including 5S rRNA (NR_023363), 5.8S rRNA (NR_003285), 18S rRNA (NR_003286) and 28S rRNA (NR_003287). Put the rRNA sequences in the file “human.ribosomal.rna.fa” with fastq format.
2d. Use Bowtie to index rRNA sequences.
bowtie2-build human.ribosomal.rna.fa hg.ribosome
2e. Align trimmed ribosome profiling reads from step 3b to rRNAs, and obtain non-rRNA reads.
bowtie2 -x hg.ribosome -U adapter.SRR1802146.fastq --un norrna.adapter.SRR1802146.fastq -S ribosome.adapter.SRR1802146.fastq
2f. Align non-rRNA reads to the reference transcriptome and genome, and obtain the alignment file in SAM format. The human reference transcriptome and genome can be obtained from GENCODE database.
tophat --GTF gencode.v28.annotation.gtf --no-convert-bam -o outputDir GRCh38genome.index  norrna.adapter. SRR1802146.fastq

3. Run “Rfoot.pl” to identify transcriptomic non-ribosomal protein-RNA footprints.

Usage: perl Rfoot.pl -i readFile -t transcriptFile -o outputFile [-w windownSize] [-r readnumCutoff] [-f translatedORFFile]
-i readFile: input read mapping file in SAM format;
-t transcriptFile: transcripts of interest in genePred format;
-o outputFile: output file reporting candidate non-ribosomal protein-RNA footprints;
-w windownSize [optional]: scanning window size, default: 60 nt;
-r readnumCutoff [optional]: cutoff of supported read number, default: 10;
-f translatedORFFile [optional]: file with translated ORF regions in genePred format. These regions will be excluded for the consideration of non-ribosomal protein-RNA complexes.

Example commands:
3a. Obtain human reference transcriptome annotation:
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.annotation.gtf.gz
Get gtfToGenePred command line from UCSC Genome Browser, and use the tool to convert GTF file to genePred format:
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred
gtfToGenePred gencode.v28.annotation.gtf gencode.v28.annotation.genePred.txt
3b. Run “Rfoot.pl” to get the prediction:
perl Rfoot.pl -i SRR1802146.mapping.sam -t gencode.v28.annotation.genePred.txt -o candidate.nonribosome.sites.txt -f  translatedORF.genepred.txt

The output file reporting candidate non-ribosomal protein-RNA footprints, with the following columns: “transcript ID”, “chromosome”, “strand”, “footprint start site”, “footprint end site”, “footprint length”, “supported read number”, “position with maximum number of reads”, “maximum number of reads in a location”, “positions with supported reads (location:read number)”.

