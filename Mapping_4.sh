#open a screen (on the chpc, if you don't run a command in screen, it will be killed in 2 minutes)
#start an interactive job on the chpc
#qsub -I -l select=1:ncpus=24:mpiprocs=4 -q smp -P CBBI1195 -l walltime=96:00:00
#load modules (CHPC)
module add chpc/BIOMODULES
module add bwa/0.7.17
module add samtools/1.9
module add picard/2.2.1
module add bcftools/1.6.33

#Tell the script where to find the reference genome by assigning the genome variable to the path to our reference genome:
refgenome=/home/cmeiring/lustre/reference_genome/sis2-181106_HiC.fasta

#Index the reference genome for use by bwa and samtools.
#bwa index $refgenome

#Make directories for all the output files
#mkdir -p sam_files bam_files bcf_files _vcf_files

#The first thing we do is assign the name of the FASTQ file we’re currently working with to a variable called fq1 and tell the script to echo the filename back to us so we can check which file we’re on.
for fq1 in /home/cmeiring/lustre/Wild_dog_samples/OG/WD_20200421/*_forward.fq.gz
  do
  echo 'working with file $fq1'

#We then extract the base name of the file (excluding the path and .fq extension)and assign it to a new variable called base.
  base=$(basename $fq1 _forward.fq.gz)
  echo 'basename is $base'

#We can use the base variable to access both the base_forward.fq and base_reverse.fq input files, and create variables to store the names of our output files. This makes the script easier to read because we don’t need to type out the full name of each of the files: instead, we use the base variable, but add a different extension (e.g. .sam, .bam) for each file produced by our workflow.
  fq1=/home/cmeiring/lustre/Wild_dog_samples/OG/WD_20200421/${base}_forward.fq.gz
  fq2=/home/cmeiring/lustre/Wild_dog_samples/OG/WD_20200421/${base}_reverse.fq.gz
  sam=/home/cmeiring/lustre/Wild_dog_samples/sam_files/${base}.aligned.sam
  bam=/home/cmeiring/lustre/Wild_dog_samples/bam_files/${base}.aligned.bam
  sorted_bam=/home/cmeiring/lustre/Wild_dog_samples/bam_files/${base}.aligned.sorted.bam
  validate=/home/cmeiring/lustre/Wild_dog_samples/bam_files/${base}_validate.txt
  RG_added=/home/cmeiring/lustre/Wild_dog_samples/bam_files/${base}.aligned.bam.rg.bam
  mapping_stats=/home/cmeiring/lustre/Wild_dog_samples/bam_files/${base}.stats.txt
  raw_bcf=/home/cmeiring/lustre/Wild_dog_samples/bcf_files/${base}.raw.bcf
  variants=/home/cmeiring/lustre/Wild_dog_samples/vcf_files/${base}.variants.vcf
  variant_stats=/home/cmeiring/lustre/Wild_dog_samples/vcf_files/${base}.variant.stats???
  final_varaints=/home/cmeiring/lustre/Wild_dog_samples/vcf_files/${base}.final_varaints.vcf

  #1) Align the reads to the reference genome and output a .sam file:
    bwa mem -M -t 32 $refgenome $fq1 $fq2 > $sam
    #-M: mark shorter split hits as secondary (duplicates)
    #-t: number of threads [1]

  #2) Convert the SAM file to BAM format:
    #samtools view -bS $sam > $bam

  #3) Validate the bam files (make sure the .fq files have read groups)
    #picard ValidateSamFile I=$bam O=$validate MODE=SUMMARY

  #4) Add Read Groups if necessary
    #picard AddOrReplaceReadGroups I=$bam O=$RG_added RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20

  #5) Validate the bam files again (make sure the .fq files have read groups)
    #picard ValidateSamFile I=$bam O=$validate MODE=SUMMARY
    #if you have no more errors you can continue

  #6) Sort and index the BAM sam_files
    #samtools sort $RG_added -o $sorted_bam
    #samtools index $sorted_bam

  #7) Generate mapping statistics
    #samtools flagstat $sorted_bam > $mapping_stats

  #8) Use bcftools mpileup to calculate the read coverage of positions in the genome:
    #bcftools mpileup -Ou -E -Q 30 -p -a DP,AD -f $refgenome $sorted_bam
    #-Ou:--output-type TYPE 'b' compressed BCF; 'u' uncompressed BCF;'z' compressed VCF; 'v' uncompressed VCF [v]
    #-E: --redo-BAQ (recalculate BAQ on the fly, ignore existing BQs)
    #-Q: --min-BQ INT (skip bases with baseQ/BAQ smaller than INT [13])
    #-p: --per-sample-mF (apply -m and -F per-sample for increased sensitivity)
    #-a: --annotate LIST (optional tags to output; '?' to list [])
    #-f: --fasta-ref FILE (faidx indexed reference sequence file)

  #9) Call SNPs with bcftools:
    #bcftools call -mv -Ob -o $variants $raw_bcf
    #-m: --multiallelic-caller (alternative model for multiallelic and rare-variant calling (conflicts with -c))
    #-v: --variants-only (output variant sites only)
    #-Ob: --output-type <b|u|z|v> (output type: 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v])

  #10 Get variant statistics
    #bcftools stats $variants -o $variant_stats

  #11 Filter variants

  done
