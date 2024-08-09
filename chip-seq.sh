
###########################################
### Chapter 1: Introduction to ChIP-seq ###
###########################################

# No code


##################################
### Chapter 2: Getting Started ###
##################################

##### Data retrieval from GEO
# Should be adapted for all datasets

# Data from Domcke/Bardet et al. Nature 2015
# ChIP-seq for the TF NRF1 in mouse embryonic stem cells either wild type (WT) or with no DNA methylation (TKO; triple knockout for DNMT1, DNMT3a, DNMT3b)
# Available on GEO GSE67867

# NRF1_CHIP_WT_1 - GSM1891641 - SRR2500883 - single-end
# SRR ID can be found via the SRA link at the bottom of the GEO page

# Create directory for datasets
mkdir data

# Run fastq-dump (~3 h)
fastq-dump --origfmt --outdir data --gzip -A SRR2500883

# Rename sample
mv data/SRR2500883.fastq.gz data/NRF1_CHIP_WT_1.fastq.gz

# Temporary files (e.g. SRR2500883.sra) are stored in the tmp directory defined by TMPDIR so make sure you have enough space and remove them once processed
rm ${TMPDIR}/sra/SRR2500883.sra

# Visualise file
gunzip -c data/NRF1_CHIP_WT_1.fastq.gz | head


##### Running code for multiple samples

# For one sample
sample=NRF1_CHIP_WT_1
ls data/${sample}.fastq.gz

# For several samples
for sample in NRF1_CHIP_WT_1 NRF1_CHIP_WT_2 NRF1_INPUT_WT NRF1_CHIP_TKO_1 NRF1_CHIP_TKO_2 NRF1_INPUT_TKO H3K27AC_CHIP_WT_1 H3K27AC_CHIP_WT_2 H3K27AC_CHIP_TKO_1 H3K27AC_CHIP_TKO_2
do
    ls data/${sample}.fastq.gz
done


##### Unix output redirection

variable="line1\nline2\nline3"
echo -e $variable

# Output redirection as input to the next command using |
echo -e $variable | head -n 1

# Output redirection as input to the current command using <()
head -n 1 <(echo -e $variable)


##### Disk space optimisation

# Extracting ouput of a compressed file without decompressing it
gunzip -c data/NRF1_CHIP_WT_1.fastq.gz | head


##########################################
### Chapter 3: General Quality Control ###
##########################################

##### FastQC
sample=NRF1_CHIP_WT_1

# Run FastQC (~3 min)
fastqc -q -o data data/${sample}.fastq.gz

# Visualise html output file
see data/${sample}_fastqc.html

# Output zipped files containing all information and plot (can be deleted)
ls data/${sample}_fastqc.zip


##### Trim Galore (optional)

# Run Trim Galore (~10 min)
trim_galore -q 20 --stringency 2 -o data data/${sample}.fastq.gz

# List output file
ls data/${sample}_trimmed.fq.gz

# Re-run FastQC to check the trimming performance


####################################
### Chapter 4: Genomic Alignment ###
####################################

##### Genomic alignment with Bowtie 2

# Create directories
mkdir reads
mkdir -p genomes/mm10
mkdir -p indices/mm10
 
# Genome sequence can be downloaded from UCSC
# e.g. for mouse genome (mm10)
wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz -O genomes/mm10/chromFa.tar.gz
tar -zxvf genomes/mm10/chromFa.tar.gz -C genomes/mm10
# Keep only the conventional chromosomes (chr1-19,X,Y,M)
rm genomes/mm10/chromFa.tar.gz genomes/mm10/*random* genomes/mm10/chrUn*
 
# Generate Bowtie 2 index
cd genomes/mm10/
bowtie2-build chr1.fa,chr10.fa,chr11.fa,chr12.fa,chr13.fa,chr14.fa,chr15.fa,chr16.fa,chr17.fa,chr18.fa,chr19.fa,chr2.fa,chr3.fa,chr4.fa,chr5.fa,chr6.fa,chr7.fa,chr8.fa,chr9.fa,chrM.fa,chrX.fa,chrY.fa ../../indices/mm10/mm10
cd ../..

sample=NRF1_CHIP_WT_1

# Run Bowtie 2 (~1 h)
# Input file can be .fastq.gz or _trimmed.fq.gz if trimming was performed
# Option -p n if n cores are available
>bowtie2 -x indices/mm10/mm10 -U data/${sample}.fastq.gz > reads/${sample}.sam
#use cammand
>bowtie2 -x  indices/mm10/mm10 -U data/SRR2500883.fastq.trimmed.gz/SRR2500883_trimmed.fq.gz > reads/83.sam
39473198 reads; of these:
  39473198 (100.00%) were unpaired; of these:
    1672457 (4.24%) aligned 0 times
    27633590 (70.01%) aligned exactly 1 time
    10167151 (25.76%) aligned >1 times
95.76% overall alignment rate

# Filter for -q 10 (~4 min)
samtools view -Sb -q 10 reads/${sample}.sam > reads/${sample}_nonSorted.bam
##use cammand
samtools view -Sb -q 10 reads/83.sam > reads/${sample}_nonSorted.bam

# Sort BAM file by genomic coordinates (~5 min)
samtools sort reads/${sample}_nonSorted.bam > reads/${sample}.bam

# Remove intermediate files
rm reads/${sample}.sam reads/${sample}_nonSorted.bam
 
# One-line alternative (no intermediate files)
bowtie2 -x indices/mm10/mm10 -U data/${sample}.fastq.gz | samtools view -Sb -q 10 | samtools sort > reads/${sample}.bam
 
# Index BAM file to allow fast access
samtools index reads/${sample}.bam
 
# Visualise file in SAM or BED format
samtools view reads/${sample}.bam | head
bamToBed -i reads/${sample}.bam | head


##### Mapping statistics

# Number of raw reads (i.e. every 4th line of FASTQ file)
gunzip -c data/${sample}.fastq.gz | awk 'END{print NR/4}'
gunzip -c data/${sample}.fastq.gz | awk '(NR%4==2)' | wc -l

# Number of aligned reads
bamToBed -i reads/${sample}.bam | wc -l
 
# Number of unique positions chromosome + start
bamToBed -i reads/${sample}.bam | awk '{if($6=="+"){position=$1":"$2}else if($6=="-"){position=$1":"$3};total[position]=1}END{print length(total)}'
 
# One-line command
# sample / raw_reads / mapped_reads / percent_mapped / unique_positions / percent_unique / most_repeated_read / number_repeated / percent_repeated
bamToBed -i reads/${sample}.bam | awk -v OFS="\t" -v sample=$sample -v raw=$(gunzip -c data/${sample}.fastq.gz | awk -v OFS="\t" '(NR%4==2)' | wc -l) 'BEGIN{max=0}{total++;if($6=="+"){position=$1":"$2}else if($6=="-"){position=$1":"$3};count[position]++; if(count[position]>max){max=count[position];maxPos=position}}END{totalPos=length(count); print sample,raw,total,total*100/raw,totalPos,totalPos*100/total,maxPos,count[maxPos],count[maxPos]*100/total}'
USE CAMMAND- 
bamToBed -i reads/${sample}.bam | awk -v OFS="\t" -v sample=$sample -v raw=$(gunzip -c data/SRR2500883.fastq.trimmed.gz/SRR2500883_trimmed.fq.gz | awk -v OFS="\t" '(NR%4==2)' | wc -l) 'BEGIN{max=0}{total++;if($6=="+"){position=$1":"$2}else if($6=="-"){position=$1":"$3};count[position]++; if(count[position]>max){max=count[position];maxPos=position}}END{totalPos=length(count); print sample,raw,total,total*100/raw,totalPos,totalPos*100/total,maxPos,count[maxPos],count[maxPos]*100/total}'
####################################################
### Chapter 5: ChIP-seq-specific Quality Control ###
####################################################

##### ChIPQC

# Open R
R

# Load libraries and set up the environment
library(ChIPQC)
library(BiocParallel)
register(SerialParam())

# Run analysis and create report of results (~10 min)
qc=ChIPQC("NRF1_sample_sheet_without_peaks.csv", "mm10", blacklist="mm10_blacklist.bed", consensus=TRUE, bCount=TRUE)
ChIPQCreport(qc, facet=FALSE, colourBy="Condition")

# Close R (and confirm with y)
q()


###############################
### Chapter 6: Peak Calling ###
###############################

##### Peakzilla for transcription factor data

# Create directory for peaks
mkdir peaks

sample=NRF1_CHIP_WT_1
control=NRF1_INPUT_WT

## Break down of the commands (many intermediate files)

# Run Peakzilla
bamToBed -i reads/${sample}.bam > reads/${sample}.bed
bamToBed -i reads/${control}.bam > reads/${control}.bed
peakzilla.py reads/${sample}.bed reads/${control}.bed -l peaks/${sample}_peakzilla_report.txt > peaks/${sample}_peakzilla.tsv

# Removal of peaks on chrM (density is usually abnormally high along the complete chrM)
# Reorganization of output into BED-like format: chromosome / start (0-based) / end / summit / score / fold_enrichment
awk -v OFS="\t" '(NR>1&&$1!="chrM"){print $1,$2-1,$3,$5,$6,$9}' peaks/${sample}_peakzilla.tsv > peaks/${sample}_peaks_peakzilla_tmp.bed

# Removal of peaks in blacklisted regions
intersectBed -v -a peaks/${sample}_peaks_peakzilla_tmp.bed -b mm10_blacklist.bed.gz > peaks/${sample}_peakzilla.bed

# Remove intermediate files
rm reads/${sample}.bed reads/${control}.bed peaks/${sample}_peakzilla.tsv peaks/${sample}_peaks_peakzilla_tmp.bed

## One-line command (~6 min)
peakzilla.py <(bamToBed -i reads/${sample}.bam) <(bamToBed -i reads/${control}.bam) -l peaks/${sample}_peakzilla_report.txt | awk -v OFS="\t" '(NR>1&&$1!="chrM"){print $1,$2-1,$3,$5,$6,$9}' | intersectBed -v -a stdin -b mm10_blacklist.bed.gz > peaks/${sample}_peaks_peakzilla.bed


##### MACS2 for histone mark data

# Run MACS2 (~7 min)
# q-value threshold can be adapted with the option --broad-cutoff (default 0.1)
macs2 callpeak -t reads/${sample}.bam -c reads/${control}.bam -f BAM -g mm --outdir peaks -n ${sample} --broad 2> peaks/${sample}_macs2_report.txt

# Reorganisation of output into bed-like format: chromosome / start (0-based) / end / summit / score / fold_enrichment
cat peaks/${sample}_peaks.xls | grep -v "#" | awk '(NR>2&&$1!="chrM"){print $1,$2,$3,$4,$6,$7}' | intersectBed -v -a stdin -b mm10_blacklist.bed.gz > peaks/${sample}_peaks_macs.bed

# Removal of intermediate files
rm peaks/${sample}_peaks.xls peaks/${sample}_model.r peaks/${sample}_peaks.broadPeak peaks/${sample}_peaks.gappedPeak


##### Saturation analysis

sample=NRF1_CHIP_WT_1
control=NRF1_INPUT_WT

## Number of peaks found
wc -l peaks/${sample}_peaks_peakzilla.bed

## Saturation analysis
# Does the number of peaks saturate when using all reads available
# When working on MAC OS, use gshuf instead of shuf
# In case seq does not return integers, use seq -f %f
NB=$(samtools idxstats reads/${sample}.bam | awk '{t=t+$3}END{print t}')
L=$sample
for subset in `seq 1000000 1000000 $NB`;do
    L=$L"\t"$(peakzilla.py <(bamToBed -i reads/${sample}.bam | shuf -n $subset --random-source=reads/${sample}.bam) <(bamToBed -i reads/${control}.bam) -l /dev/null | wc -l)
done
echo -e $L > peaks/${sample}_peakzilla_saturation_table.txt

# Open R
R

sample="NRF1_CHIP_WT_1"

# Load saturation table
d=read.table(paste("peaks/",sample,"_peakzilla_saturation_table.txt",sep=""),row.names=1)

# Generate the saturation plot
pdf(paste("peaks/",sample,"_peakzilla_saturation_plot.pdf",sep=""))
par(bg="white")
plot(seq(1,ncol(d),1),d[1,],type="l",xlab="Number of reads in million",ylab="Number of peaks called",main=c("Saturation analysis",rownames(d)[1],paste("Total",d[1,ncol(d)],"peaks")))
dev.off()

# Close R
q()

# Visualise plot
display peaks/${sample}_peakzilla_saturation_plot.pdf &


#####################################
### Chapter 7: Data visualisation ###
#####################################

##### Read densities

# Create a directory for the density tracks
mkdir tracks

sample=NRF1_CHIP_WT_1
 
# Chromosome sizes (to calculate read densities at all positions)
# To be downloaded from UCSC (genome.ucsc.edu/), Downloads, Genome Data
# File is sorted according to chromosome and only conventional chromosomes are kept
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes
sort -k1,1 mm10.chrom.sizes | grep -v chrUn | grep -v random > mm10.chrom.sizes.tmp
mv mm10.chrom.sizes.tmp mm10.chrom.sizes

# Scale to normalise the number of reads at each position to 1,000,000 mapped reads in the library
# Count total number of mapped reads if the index BAI file is available (~1 min)
samtools idxstats reads/${sample}.bam | awk '{if($1!="*"){total=total+$3}}END{print 1000000/total}'
# Otherwise count total number of mapped reads from the BAM file (~ 2 min)
bamToBed -i reads/${sample}.bam | wc -l | awk '{print 1000000/total}'
scale=0.0322351
# One-line command
scale=$(samtools idxstats reads/${sample}.bam | awk '{if($1!="*"){total=total+$3}}END{print 1000000/total}')

# Extend reads checking they don't extend over the limits of chromosomes
extend=200
# Extend reads from the forward strand (since files are sorted according to starts, there is no need to sort the output again)
bamToBed -i reads/${sample}.bam | awk -v OFS="\t" -v E=$extend -v file= mm10.chrom.sizes 'BEGIN{while(getline<file){S[$1]=$2}}($6=="+"){if($2+E>S[$1]){print $1,$2,S[$1]}else{print $1,$2,$2+E}}' > tracks/${sample}_forward.bed
# Extend reads from the reverse strand (since ends are not sorted if reads do not have all the same length e.g. after trimming, we need to sort the output again by chromosome / start)
bamToBed -i reads/${sample}.bam | awk -v OFS="\t" -v E=$extend '($6=="-"){if($3-E<1){print $1,"1",$3}else{print $1,$3-E,$3}}' | sort -k1,1 -k2,2n > tracks/${sample}_reverse.bed
# Merge sorted file (sort -m much faster than just sort since the input files are already sorted)
sort -m -k1,1 -k2,2n tracks/${sample}_forward.bed tracks/${sample}_reverse.bed > tracks/${sample}.bed

# Generate bedGraph file from BED (read coverage at each non-zero position)
genomeCoverageBed -bg -i tracks/${sample}.bed -g mm10.chrom.sizes -scale $scale > tracks/${sample}.bg
# Remove intermediate files
rm tracks/${sample}_forward.bed tracks/${sample}_reverse.bed tracks/${sample}.bed

# One-line command (~6 min)
extend=200
sort -m -k1,1 -k2,2n <(bamToBed -i reads/${sample}.bam | awk -v OFS="\t" -v E=$extend -v file= mm10.chrom.sizes 'BEGIN{while(getline<file){S[$1]=$2}}($6=="+"){if($2+E>S[$1]){print $1,$2,S[$1]}else{print $1,$2,$2+E}}') <(bamToBed -i reads/${sample}.bam | awk -v OFS="\t" -v E=$extend '($6=="-"){if($3-E<1){print $1,"1",$3}else{print $1,$3-E,$3}}' | sort -k1,1 -k2,2n) | genomeCoverageBed -bg -i stdin -g mm10.chrom.sizes -scale $scale > tracks/${sample}.bg

# Generate bigWig file from bedGraph (~2 min)
bedGraphToBigWig tracks/${sample}.bg mm10.chrom.sizes tracks/${sample}.bw

# Remove intermediate file
rm tracks/${sample}.bg


##### Peak regions

# Generate a simple BED file (chromosome / start / end)
cut -f1-3 peaks/${sample}_peaks_peakzilla.bed > peaks/${sample}_peaks_peakzilla_short.bed

# Generate bigBed file from BED
bedToBigBed peaks/${sample}_peaks_peakzilla_short.bed mm10.chrom.sizes tracks/${sample}_peaks_peakzilla.bb
rm peaks/${sample}_peaks_peakzilla_short.bed


##### Genome browser

# Track lines to upload into UCSC

# BAM files (with .bai file in same directory)
# track type=bam name="NRF1_CHIP_WT_1_reads" bigDataUrl=http://chipseq.u-strasbg.fr:8080/reads/NRF1_CHIP_WT_1.bam

# bigWig files
# track type=bigWig name="NRF1_CHIP_WT_1_density" bigDataUrl=http://chipseq.u-strasbg.fr:8080/tracks/NRF1_CHIP_WT_1.bw

# bigBed files
# track type=bigBed name="NRF1_CHIP_WT_1_peaks" bigDataUrl=http://chipseq.u-strasbg.fr:8080/tracks/NRF1_CHIP_WT_1_peaks_peakzilla.bb


#######################################
### Chapter 8: Comparative Analysis ###
#######################################

##### Overlap of peak regions

for sample1 in NRF1_CHIP_WT_1 NRF1_CHIP_WT_2 NRF1_CHIP_TKO_1 NRF1_CHIP_TKO_2
do
    L=$sample1"\t"$(cat peaks/${sample1}_peaks_peakzilla.bed | wc -l)
    for sample2 in NRF1_CHIP_WT_1 NRF1_CHIP_WT_2 NRF1_CHIP_TKO_1 NRF1_CHIP_TKO_2
    do
	    L=$L"\t"$(intersectBed -u -a peaks/${sample1}_peaks_peakzilla.bed -b peaks/${sample2}_peaks_peakzilla.bed | wc -l)
    done
    echo -e $L
done > changes/peak_overlap_table.txt

# Table with: sample / total_number_of_peaks / overlap_for_each_sample
cat changes/peak_overlap_table.txt | column -t

# Transformation into percentages
# Make sure your system is set to decimal number represented as dots and not commas (export LC_ALL=C)
cat changes/peak_overlap_table.txt | awk '{L=$1"\t"$2;for(i=3;i<=NF;i++){L=L"\t"$i*100/$2};print L}' > changes/peak_overlap_percent_table.txt
cat changes/peak_overlap_percent_table.txt | column -t

# Clustering of samples and representation in a heatmap

# Open R and load library
R
library(NMF)

# Load table
x=read.table("changes/peak_overlap_percent_table.txt")
data=x[,3:6]
rownames(data)=colnames(data)=x[,1]

# Generate a heatmap that clusters the samples
pdf("changes/peak_overlap_percent_heatmap.pdf")
aheatmap(data,breaks=50,annRow=data.frame(condition=c("WT","WT","TKO","TKO"),replicate=as.factor(c(1,2,1,2))),annColors=list('condition'=c("blue","red"),'replicate'=c("black","grey")))
dev.off()

# Close R
q()

# Visualise plot
display changes/peak_overlap_percent_heatmap.pdf &

##### Irreproducible discovery rate (IDR)

## Peak calling for IDR

# Call peaks on each replicate with a relaxed threshold (-s 0.1 instead of default 1)
peakzilla.py -s 0.1 <(bamToBed -i reads/NRF1_CHIP_WT_1.bam) <(bamToBed -i reads/NRF1_INPUT_WT.bam) -l peaks/NRF1_CHIP_WT_1_all_peakzilla_report.txt | awk -v OFS="\t" '(NR>1&&$1!="chrM"){print $1,$2-1,$3,$5,$6,"+"}' | intersectBed -v -a stdin -b mm10_blacklist.bed.gz > peaks/NRF1_CHIP_WT_1_all_peaks_peakzilla.bed
peakzilla.py -s 0.1 <(bamToBed -i reads/NRF1_CHIP_WT_2.bam) <(bamToBed -i reads/NRF1_INPUT_WT.bam) -l peaks/NRF1_CHIP_WT_2_all_peakzilla_report.txt | awk -v OFS="\t" '(NR>1&&$1!="chrM"){print $1,$2-1,$3,$5,$6,"+"}' | intersectBed -v -a stdin -b mm10_blacklist.bed.gz > peaks/NRF1_CHIP_WT_2_all_peaks_peakzilla.bed

# Pool BAM files and call peaks on merged file with a relaxed threshold.
samtools merge reads/NRF1_CHIP_WT.bam reads/NRF1_CHIP_WT_1.bam reads/NRF1_CHIP_WT_2.bam

peakzilla.py -s 0.1 <(bamToBed -i reads/NRF1_CHIP_WT.bam) <(bamToBed -i reads/NRF1_INPUT_WT.bam) -l peaks/NRF1_CHIP_WT_1_all_peakzilla_report.txt | awk -v OFS="\t" '(NR>1&&$1!="chrM"){print $1,$2-1,$3,$5,$6,"+"}' | intersectBed -v -a stdin -b mm10_blacklist.bed.gz > peaks/NRF1_CHIP_WT_all_peaks_peakzilla.bed

## Calculating IDR
idr --samples peaks/NRF1_CHIP_WT_1_all_peaks_peakzilla.bed peaks/NRF1_CHIP_WT_2_all_peaks_peakzilla.bed --peak-list peaks/NRF1_CHIP_WT_all_peaks_peakzilla.bed --input-file-type bed --rank score --idr-threshold 0.05 --output-file peaks/NRF1_CHIP_WT_idr.txt --plot

##### Comparison of read densities

## Merging peak regions

# Concatenate peaks from all samples and sort by chr and start
cat peaks/NRF1_CHIP_WT_1_peaks_peakzilla.bed peaks/NRF1_CHIP_WT_2_peaks_peakzilla.bed peaks/NRF1_CHIP_TKO_1_peaks_peakzilla.bed peaks/NRF1_CHIP_TKO_2_peaks_peakzilla.bed | sort -k1,1 -k2,2n > changes/NRF1_all_regions_tmp.txt

# Merge peak regions from all samples
mergeBed -i changes/NRF1_all_regions_tmp.txt > changes/NRF1_all_regions.txt

# Remove intermediate files
rm changes/NRF1_all_regions_tmp.txt

# One-line command
cat peaks/NRF1_CHIP_WT_1_peaks_peakzilla.bed peaks/NRF1_CHIP_WT_2_peaks_peakzilla.bed peaks/NRF1_CHIP_TKO_1_peaks_peakzilla.bed peaks/NRF1_CHIP_TKO_2_peaks_peakzilla.bed | sort -k1,1 -k2,2n | mergeBed -i stdin > changes/NRF1_all_regions.txt

## Counting reads for each sample

# For one sample
intersectBed -c -sorted -a changes/NRF1_all_regions.txt -b reads/NRF1_CHIP_WT_1.bam > changes/NRF1_all_regions_count1.txt

# For all samples (~2 min)
intersectBed -c -sorted -a changes/NRF1_all_regions.txt -b reads/NRF1_CHIP_WT_1.bam | intersectBed -c -sorted -a stdin -b reads/NRF1_CHIP_WT_2.bam | intersectBed -c -sorted -a stdin -b reads/NRF1_CHIP_TKO_1.bam | intersectBed -c -sorted -a stdin -b reads/NRF1_CHIP_TKO_2.bam > changes/NRF1_all_regions_count.txt
cat changes/NRF1_all_regions_count.txt | head | column -t

## Normalising read counts

# Normalised read counts (RPKM table)
# RPKM = reads x (1000000 / total_library) x (1000 / region_size)
# Make sure your system is set to decimal number represented as dots and not commas (export LC_ALL=C)
TOT1=$(samtools idxstats reads/NRF1_CHIP_WT_1.bam | awk '($1!~"*"){t=t+$3}END{print t}')
TOT2=$(samtools idxstats reads/NRF1_CHIP_WT_2.bam | awk '($1!~"*"){t=t+$3}END{print t}')
TOT3=$(samtools idxstats reads/NRF1_CHIP_TKO_1.bam | awk '($1!~"*"){t=t+$3}END{print t}')
TOT4=$(samtools idxstats reads/NRF1_CHIP_TKO_2.bam | awk '($1!~"*"){t=t+$3}END{print t}')
cat changes/NRF1_all_regions_count.txt | awk -v OFS="\t" -vTOT1=$TOT1 -vTOT2=$TOT2 -vTOT3=$TOT3 -vTOT4=$TOT4 '{size=$3-$2;print $1,$2,$3,$4*(1000000/TOT1)*(1000/size),$5*(1000000/TOT2)*(1000/size),$6*(1000000/TOT3)*(1000/size),$7*(1000000/TOT4)*(1000/size)}' > changes/NRF1_all_regions_rpkm.txt
cat changes/NRF1_all_regions_rpkm.txt | head | column -t

# Count reads and generate RPKM tables using for loop (~6 min)
cat changes/NRF1_all_regions.txt > changes/NRF1_all_regions_count.txt
cat changes/NRF1_all_regions.txt > changes/NRF1_all_regions_rpkm.txt
for sample in NRF1_CHIP_WT_1 NRF1_CHIP_WT_2 NRF1_CHIP_TKO_1 NRF1_CHIP_TKO_2
do
    TOT=$(samtools idxstats reads/${sample}.bam | awk '($1!~"*"){t=t+$3}END{print t}')
    intersectBed -c -sorted -a changes/NRF1_all_regions_count.txt -b reads/${sample}.bam > tmp
    mv tmp changes/NRF1_all_regions_count.txt
    paste changes/NRF1_all_regions_rpkm.txt <(cat changes/NRF1_all_regions_count.txt | awk -v OFS="\t" -vTOT=$TOT '{size=$3-$2;print $NF*(1000000/TOT)*(1000/size)}') > tmp
    mv tmp changes/NRF1_all_regions_rpkm.txt
done
cat changes/NRF1_all_regions_count.txt | head | column -t
cat changes/NRF1_all_regions_rpkm.txt | head | column -t

## Comparing read counts

# Open  R
R

# Load RPKM table
data=read.table("changes/NRF1_all_regions_rpkm.txt")
colnames(data)=c("chr","start","end","NRF1_CHIP_WT_1","NRF1_CHIP_WT_2","NRF1_CHIP_TKO_1","NRF1_CHIP_TKO_2")

# Visualise pairwise comparisons in a scatterplot
pdf("changes/NRF1_all_regions_rpkm_cor_scatter.pdf")
par(bg="white",mfrow=c(2,2))
smoothScatter(log2(data$NRF1_CHIP_WT_1),log2(data$NRF1_CHIP_WT_2),xlim=c(-4,10),ylim=c(-4,10),xlab="NRF1_CHIP_WT_1 (log2)",ylab="NRF1_CHIP_WT_2 (log2)",main=paste("PCC =",signif(cor(data$NRF1_CHIP_WT_1,data$NRF1_CHIP_WT_2),2)))
smoothScatter(log2(data$NRF1_CHIP_TKO_1),log2(data$NRF1_CHIP_TKO_2),xlim=c(-4,10),ylim=c(-4,10),xlab="NRF1_CHIP_TKO_1 (log2)",ylab="NRF1_CHIP_TKO_2 (log2)",main=paste("PCC =",signif(cor(data$NRF1_CHIP_TKO_1,data$NRF1_CHIP_TKO_2),2)))
smoothScatter(log2(data$NRF1_CHIP_WT_1),log2(data$NRF1_CHIP_TKO_1),xlim=c(-4,10),ylim=c(-4,10),xlab="NRF1_CHIP_WT_1 (log2)",ylab="NRF1_CHIP_TKO_1 (log2)",main=paste("PCC =",signif(cor(data$NRF1_CHIP_WT_1,data$NRF1_CHIP_TKO_1),2)))
smoothScatter(log2(data$NRF1_CHIP_WT_2),log2(data$NRF1_CHIP_TKO_2),xlim=c(-4,10),ylim=c(-4,10),xlab="NRF1_CHIP_WT_2 (log2)",ylab="NRF1_CHIP_TKO_2 (log2)",main=paste("PCC =",signif(cor(data$NRF1_CHIP_WT_2,data$NRF1_CHIP_TKO_2),2)))
dev.off()

# Compare the RPKM changes between WT and TKO samples across replicated samples
pdf("changes/NRF1_all_regions_rpkm_cor_delta_scatter.pdf")
par(bg="white")
smoothScatter(log2((data$NRF1_CHIP_TKO_1+0.1)/(data$NRF1_CHIP_WT_1+0.1)),log2((data$NRF1_CHIP_TKO_2+0.1)/(data$NRF1_CHIP_WT_2+0.1)),xlim=c(-6,10),ylim=c(-6,10),xlab="NRF1_CHIP_TKO_1 / NRF1_CHIP_WT_1 (log2)",ylab="NRF1_CHIP_TKO_2 / NRF1_CHIP_WT_2 (log2)",main=paste("PCC =",signif(cor((data$NRF1_CHIP_TKO_1+0.1)/(data$NRF1_CHIP_WT_1+0.1),(data$NRF1_CHIP_TKO_2+0.1)/(data$NRF1_CHIP_WT_2+0.1)),2)))
dev.off()

# Close R
q()


# Visualise plots
display changes/NRF1_all_regions_rpkm_cor_scatter.pdf &
display changes/NRF1_all_regions_rpkm_cor_delta_scatter.pdf &


##### DESeq2

# Open R and load library
R
library(DESeq2)

# Load count table (from Section 8.3.2)
x=read.table("changes/NRF1_all_regions_count.txt")
colnames(x)=c("chr","start","end","NRF1_CHIP_WT_1","NRF1_CHIP_WT_2","NRF1_CHIP_TKO_1","NRF1_CHIP_TKO_2")
d=x[,4:7]

# Conditions
colData=data.frame(condition=c("WT","WT","TKO","TKO"))
rownames(colData)=colnames(d)

# DESeq2 matrix
dds=DESeqDataSetFromMatrix(countData=d,colData=colData,design=~condition)
# Filter out regions with no counts (or low counts)
dds=dds[rowSums(counts(dds))>1,]

# Transform data for variance analysis
# Using regularised-logarithm transformation (RLD) when number of samples < 30
rld=rlog(dds,blind=F)
# Using the variance-stabilising transformation (VST) otherwise
vsd=vst(dds,blind=F)

# Principal components analysis (PCA)
pdf("changes/NRF1_all_regions_count_pca.pdf")
plotPCA(rld,intgroup=c("condition"))
dev.off()

# Differential expression (on non-transformed data)
dds=DESeq(dds)

# TKO vs. WT
res=results(dds,contrast=c("condition","TKO","WT"))

# Volcano plot (log2 fold change vs. adjusted p-value)
pdf("changes/NRF1_all_regions_count_TKO_vs_WT_volcano.pdf")
par(bg="white")
plot(res$log2FoldChange,-log10(res$padj),xlim=c(-10,10),ylim=c(0,83),xlab="TKO / WT (log2)",ylab="P-value (-log10)",pch=20)
lines(c(-10,10),c(10,10))
dev.off()

# Table chr / start / end / NRF1_CHIP_WT_1 / NRF1_CHIP_WT_2 / NRF1_CHIP_TKO_1 / NRF1_CHIP_TKO_2 / log2FoldChange / padj
y=cbind(x[,1:3],counts(dds,normalized=T),as.matrix(res)[,c(2,6)])
write.table(y,"changes/NRF1_all_regions_count_TKO_vs_WT_table.txt",quote=F,sep="\t",row.names=F,col.names=T)

# Close R
q()

# Visualise plots
display changes/NRF1_all_regions_count_pca.pdf &
display changes/NRF1_all_regions_count_TKO_vs_WT_volcano.pdf &

# Select differential genes
# Since the adjusted p-value scales with the fold enrichment, we only need a p-value threshold
# e.g. 10^-10
cat changes/NRF1_all_regions_count_TKO_vs_WT_table.txt | head | column -t

# TKO-specific peaks: 2212
# No header, p-value (adjusted) <= 1e-10, log2 fold change >0 (TKO>WT)
cat changes/NRF1_all_regions_count_TKO_vs_WT_table.txt | awk '(NR>1&&$9<=1e-10&&$8>0)' > changes/NRF1_all_regions_count_TKO_spec_table.txt

# WT-specific peaks: 13
# No header, adjusted p-value <= 1e-10, log2 fold change <0 (WT>TKO)
cat changes/NRF1_all_regions_count_TKO_vs_WT_table.txt | awk '(NR>1 && $9<=1e-10 && $8<0)' > changes/NRF1_all_regions_count_WT_spec_table.txt

# Shared peaks: 7068
# No header, absolute log2 fold change < 1 (less than 2 fold change in either direction)
cat changes/NRF1_all_regions_count_TKO_vs_WT_table.txt | awk '(NR>1 && $8<1 && $8>-1)' > changes/NRF1_all_regions_count_shared_table.txt


##### DiffBind

# Open R and load library
R
library(DiffBind)

# Load data (~2 min)
# If peak calling has already been run
NRF1=dba(sampleSheet="NRF1_sample_sheet.csv")
# If peak calling has not yet been run
NRF1=dba(sampleSheet="NRF1_sample_sheet_without_peaks.csv")

plot(NRF1)
NRF1=dba.count(NRF1, minOverlap=2)
NRF1=dba.contrast(NRF1, categories=DBA_CONDITION, minMembers=2)
NRF1=dba.analyze(NRF1, method=c(DBA_DESEQ2))
NRF1.DB=dba.report(NRF1)
dba.plotPCA(NRF1, contrast=1, label=DBA_CONDITION)
dba.plotMA(NRF1)
dba.plotBox(NRF1)
dba.plotHeatmap(NRF1, contrast=1, correlations=FALSE)


######################################
### Chapter 9: Downstream Analyses ###
######################################

##### Genomic location

# Genomic location using the following hierarchy: exon (EXON) > 5'UTR (5UTR) > 3'UTR (3UTR) > intron (INTRON) > 2 kb upstream promoter (P2000) > intergenic (INTER)
head mm10_genomic_features.bed

# Make sure your system is set to decimal number represented as dots and not commas (export LC_ALL=C)

# Define TF peak summit positions
awk -v OFS="\t" '{print $1,$4-1,$4}' peaks/NRF1_CHIP_WT_1_peaks_peakzilla.bed > peaks/NRF1_CHIP_WT_1_peaks_peakzilla_summits.bed
# Genomic location of TF peak summits
intersectBed -wo -a peaks/NRF1_CHIP_WT_1_peaks_peakzilla_summits.bed -b mm10_genomic_features.bed | awk -v OFS="\t" '{F[$7]++;t++}END{for(location in F){print location,F[location]*100/t}}' > peaks/NRF1_CHIP_WT_1_peaks_peakzilla_summits_location.txt
cat peaks/NRF1_CHIP_WT_1_peaks_peakzilla_summits_location.txt

# Repartition of histone peak regions by genomic location
intersectBed -wo -a peaks/H3K27AC_CHIP_WT_1_peaks_macs.bed -b mm10_genomic_features.bed | awk -v OFS="\t" '{F[$10]+=$12;t+=$12}END{for(location in F){print location,F[location]*100/t}}' > peaks/H3K27AC_CHIP_WT_1_peaks_peakzilla_regions_location.txt
cat peaks/H3K27AC_CHIP_WT_1_peaks_peakzilla_regions_location.txt

# Overall repartitioning of genomic location
awk -v OFS="\t" '{F[$4]+=$3-$2;t+=$3-$2}END{for(location in F){print location,F[location]*100/t}}' mm10_genomic_features.bed > peaks/genome_location.txt
cat peaks/genome_location.txt

# Plots
# Open R
R

# Load tables
tf = read.table("peaks/NRF1_CHIP_WT_1_peaks_peakzilla_summits_location.txt")
histone = read.table("peaks/H3K27AC_CHIP_WT_1_peaks_peakzilla_regions_location.txt")
genome = read.table("peaks/genome_location.txt")

# Merge the three tables
d = merge(merge(tf,histone,by=1),genome,by=1)
colnames(d) = c("Location","TF","Histone","Genome")

# Piecharts of the distribution of genomic locations
pdf("peaks/genomic_location_piechart.pdf",height=5,width=15)
par(mfrow=c(1,3),bg="white")
pie(d$TF,labels=d$Location,main="NRF1")
pie(d$Histone,labels=d$Location,main="H3K27AC")
pie(d$Genome,labels=d$Location,main="Genome")
dev.off()

# Barplots of enrichment over genome
pdf("peaks/genomic_location_barplot.pdf")
par(mfrow=c(1,2),bg="white")
barplot(log2(d$TF/d$Genome),names=d$Location,main="NRF1",las=2,ylim=c(-2,6))
barplot(log2(d$Histone/d$Genome),names=d$Location,main="H3K27AC",las=2,ylim=c(-2,6))
dev.off()

# Close R
q()

# Visualise plots
display peaks/genomic_location_piechart.pdf &
display peaks/genomic_location_barplot.pdf &


##### Distance to gene TSS

# TSS positions for all genes
ls mm10_tss.bed

# Distance of TF peak summits to TSS
awk -v OFS="\t" '{print $1,$4-1,$4}' peaks/NRF1_CHIP_WT_1_peaks_peakzilla.bed | closestBed -d -t "first" -a stdin -b mm10_tss.bed | awk '{print $NF}' > peaks/NRF1_CHIP_WT_1_peaks_peakzilla_summits_dist_tss.txt

# Distance of histone peak regions centre to TSS
awk -v OFS="\t" '{c=($2+$3)/2;print $1,c-1,c}' peaks/H3K27AC_CHIP_WT_1_peaks_macs.bed | closestBed -d -t "first" -a stdin -b mm10_tss.bed | awk '{print $NF}' > peaks/H3K27AC_CHIP_WT_1_peaks_macs_center_dist_tss.txt

# Open R
R

# Load the distance to TSS tables
tf=read.table("peaks/NRF1_CHIP_WT_1_peaks_peakzilla_summits_dist_tss.txt")
histone=read.table("peaks/H3K27AC_CHIP_WT_1_peaks_macs_center_dist_tss.txt")

# Histograms of distances of peaks to closest TSS
pdf("peaks/dist_tss_hist.pdf")
par(mfrow=c(2,1),bg="white")
hist(log10(tf[,1]),main="NRF1",xlab="Distance to closet gene TSS (log10)",breaks=seq(0,7,0.2))
hist(log10(histone[,1]),main="H3K27AC",xlab="Distance to closet gene TSS (log10)",breaks=seq(0,7,0.2))
dev.off()

# Close R
q()

# Visualise plot
display peaks/dist_tss_hist.pdf &


##### Assignment to target genes

# TSS positions for all gene transcripts
ls mm10_tss.bed

# Unique list of genes closest to NRF1 peaks (7167 peaks leading to 5595 genes)
closestBed -t "first" -a peaks/NRF1_CHIP_WT_1_peaks_peakzilla.bed -b mm10_tss.bed | awk '{print $10}' | sort -u > peaks/NRF1_CHIP_WT_1_peaks_peakzilla_genes.txt


##### Genomic location annotations using R

# Open R
R

# Extract TSS and promoter regions for the UCSC known genes
require(TxDb.Mmusculus.UCSC.mm10.knownGene)
TxDB=TxDb.Mmusculus.UCSC.mm10.knownGene
library(GenomicRanges)
transcr=transcripts(TxDB, columns=c("gene_id"))
tss=flank(transcr, width=-1)
promoter=promoters(transcr)

# Extract introns of the UCSC known genes
intron=intronsByTranscript(TxDB, use.names=TRUE)

# Extract exons of the UCSC known genes
exons=exons(TxDB, columns=c("gene_id"))

# Extract UTRs of the UCSC known genes
utr5=fiveUTRsByTranscript(TxDB)
utr3=threeUTRsByTranscript(TxDB)

# Distance to nearest TSS using peak summits using peakzilla peaks
peaks=read.table("peaks/NRF1_CHIP_WT_1_peaks_peakzilla.bed", sep="\t", stringsAsFactors=F, header=F)
peakGR=GRanges(seqnames=peaks$V1, ranges=IRanges(start=peaks$V4, end=peaks$V4))
dist=distanceToNearest(peakGR, tss)@elementMetadata$distance

# Histogram of distances of peaks to closest TSS
pdf("peaks/dist_tss_hist_R.pdf")
hist(log10(dist),main="NRF1",xlab="Distance to closet gene TSS (log10)",breaks=seq(0,7,0.2))
dev.off()

# Overlap with features
peakGR$feature="intergenic"
peakGR$feature[queryHits(findOverlaps(peakGR, promoter))]="promoter"
peakGR$feature[queryHits(findOverlaps(peakGR, intron))]="intron"
peakGR$feature[queryHits(findOverlaps(peakGR, utr3))]="3_UTR"
peakGR$feature[queryHits(findOverlaps(peakGR, utr5))]="5_UTR"
peakGR$feature[queryHits(findOverlaps(peakGR, exons))]="exon"
overlaps=table(peakGR$feature)
genome=c(sum(width(ranges(reduce(exons)))),sum(as.numeric(width(ranges(reduce(unlist(intron)))))),sum(width(ranges(reduce(promoter)))), sum(as.numeric(width(ranges(reduce(unlist(utr3)))))), sum(as.numeric(width(ranges(reduce(unlist(utr5)))))))
intergenic=sum(as.numeric(read.table("mm10.chrom.sizes")$V2))-sum(genome)
genome=c(genome[1], intergenic, genome[2:5])
labs=c( "exon","intergenic", "intron", "promoter", "3_UTR", "5_UTR")

pdf("peaks/genomic_location_piechart_R.pdf",height=5,width=15)
par(mfrow=c(1,2),bg="white")
pie(overlaps,labels=names(overlaps),main="NRF1")
pie(genome,labels=labs,main="Genome")
dev.off()

# Alternatively you can use chipseeker to annotate the peaks with overlapping features
library(ChIPseeker)
TxDB=TxDb.Mmusculus.UCSC.mm10.knownGene
require(org.Mm.eg.db)
annot=annotatePeak(peakGR, TxDb=TxDB, annoDb="org.Mm.eg.db")

# Distance to TSS
pdf("peaks/dist_tss_hist_ChIPseeker.pdf")
hist(log10(abs(annot@anno$distanceToTSS)),main="NRF1",xlab="Distance to closet gene TSS (log10)",breaks=seq(0,7,0.2))
dev.off()

# Annotation proportions
annot@annoStat
pdf("peaks/genomic_location_distribution_ChIPseeker.pdf",height=3,width=5)
plotAnnoBar(annot)
dev.off()

# Gene assignment
annot@anno$SYMBOL or annot@anno$ENSEMBL
write.table(annot@anno$SYMBOL)

# Close R
q()


##### De novo motif discovery

sample=NRF1_CHIP_WT_1

# Define regions of 151bp around peak summits
awk -v OFS="\t" '{print $1,$4-75,$4+75}' peaks/${sample}_peaks_peakzilla.bed > peaks/${sample}_peaks_peakzilla_151bp.bed

# HOMER for de novo motif search (~1 h)
findMotifsGenome.pl peaks/${sample}_peaks_peakzilla_151bp.bed genomes/mm10/ motifs -size given

# Open the html page with results
see motifs/homerResults.html &


##### Known motif search

# Download the NRF1 motif from JASPAR
# jaspar.genereg.net/matrix/MA0506.1
wget http://jaspar.genereg.net/api/v1/matrix/MA0506.1.meme -O motifs/NRF1.meme

# Scan genome for motif occurrences using a p-value threshold of 10^-5 and reformat output into BED (~5 min)
mast -hit_list -mt 1e-04 motifs/NRF1.meme genomes/mm10.fa | awk '($1!~/#/){if($2=="+1"){s="+"}else{s="-"};print $1,$3-1,$4,"NRF1",$6,s}' | gzip > motifs/NRF1_mm10.bed.gz

# Number of peak regions: 7167
cat peaks/${sample}_peaks_peakzilla_151bp.bed | wc -l
# Number of peak regions with NRF1 motifs: 5245 (73%)
intersectBed -u -sorted -a peaks/${sample}_peaks_peakzilla_151bp.bed -b motifs/NRF1_mm10.bed.gz | wc -l


##### Integration with additional ChIP-seq datasets

# For each peak selection, for all samples, extract read densities (from the bigwig file) for each position within 5 kb around the peak regions (~2 min)
for peaks in TKO_spec WT_spec shared
do
    for sample in NRF1_CHIP_WT_1 NRF1_CHIP_TKO_1 H3K27AC_CHIP_WT_1 H3K27AC_CHIP_TKO_1
    do
        awk -v OFS="\t" '{center=int(($2+$3)/2);print $1,center-2500,center+2500,$2,$3}' changes/NRF1_all_regions_count_${peaks}_table.txt | bwtool extract -tabs bed stdin tracks/${sample}.bw stdout > changes/NRF1_all_regions_count_${peaks}_density_${sample}.txt
    done
done

# Open R and load library
R
library(gplots)

# Load RPKM table
rpkm=read.table("changes/NRF1_all_regions_rpkm.txt")
colnames(rpkm)=c("chr","start","end","NRF1_CHIP_WT_1","NRF1_CHIP_WT_2","NRF1_CHIP_TKO_1","NRF1_CHIP_TKO_2")

# For each peak selection
for(peaks in c("TKO_spec","shared")){
    png(paste("changes/NRF1_all_regions_count_",peaks,"_density_heatmap.png",sep=""))
    par(bg="white",mfrow=c(1,4))
    
    # For each sample
    for(sample in c("NRF1_CHIP_WT_1","NRF1_CHIP_TKO_1","H3K27AC_CHIP_WT_1","H3K27AC_CHIP_TKO_1")){
    
        # Load the density table
        x=read.table(paste("changes/NRF1_all_regions_count_",peaks,"_density_",sample,".txt",sep=""))
    
        # Add extra column with matching RPKM of region in NRF1_CHIP_WT_1 sample
        y=merge(x,rpkm[,c("chr","start","end","NRF1_CHIP_WT_1")],by.x=c(1,4,5),by.y=c("chr","start","end"))
    
        # Plot image of read densities, ordering regions (row) by RPKM and adjusting colour scale from white to black by 0.1 steps between 0 and 10 and one more step until 100
        image(t(y[order(y$"NRF1_CHIP_WT_1"),7:5006]),axes=F,col=colorpanel(101,"white","black"),breaks=c(seq(0,10,0.1),100),main=sample)
    }
    dev.off()
}

# Plot colour scale
pdf("changes/NRF1_all_regions_count_density_heatmap_scale.pdf",height=3)
par(bg="white")
plot(c(0,110),c(0,1),xlab="Read density",ylab="",pch="",axes=F)
axis(1,at=seq(0,110,10),labels=c(0:10,100))
cols=c(colorpanel(100,"white","black"),rep("black",10))
rect(seq(0,109,1),0,seq(1,110,1),1,border=cols,col=cols)
dev.off()

# Close R
q()

# Visualise plots
display changes/NRF1_all_regions_count_TKO_spec_density_heatmap.png &
display changes/NRF1_all_regions_count_shared_density_heatmap.png &
display changes/NRF1_all_regions_count_density_heatmap_scale.pdf &
