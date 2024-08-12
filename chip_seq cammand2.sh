 Organism:          arabidopsis thaliyana.
SRAtoolkit
fasterq-dump
path files
> export PATH=$HOME/tools/sra-tools/sratoolkit.2.9.6-ubuntu64/bin:$PATH (path installation tools)
fasterq-dump home/scbb/Documents  CHIP_SRR_Acc_List.txt SRR9606582
(only single end data)
1) FastQC
#FastQC is a popular tool to perform quality assessment. As a general rule, read quality decreases towards the 3’ end of reads, and if it becomes too low, bases should be removed to improve mappability.
 > cat list_chip |while read i; do cd $i ; fastqc *.gz;  cd .. ; done
2) Trim_galore
 
 remove=
         .  adapter sequences (adapter trimming)
    • low-quality reads
    • uncalled bases
use cammand-  > cat list_chip| while read i; do cd $i ; trim_galore --stringency 13 -q 30  *1.gz -o trim_galore\_$i ;cd .. ; done
#why use stringency 13 - For adapter trimming, Trim Galore! uses the first 13 bp of Illumina standard adapters ('AGATCGGAAGAGC') by default (suitable for both ends of paired-end libraries and single end), but accepts other adapter sequence..(--illumina Adapter sequence to be trimmed is the first 13bp of the Illumina universal adapter)
#paired-end files, both reads of a read-pair need to be longer than
    <INT> bp to be printed out to validated paired-end files (see option --paired). If only one read #became too short there is the possibility of keeping such unpaired single-end reads (see --retain_unpaired). Default pair-cutoff: 20 bp.
create two files
fastq.gz_trimming_report.txt, val_1.fq.gz 
3) fastqc
 > cat trim_list |while read i; do cd $i; fastqc *.gz  -o ../../trim_qc ; cd .. ; done
#create fastqc graph file and zip file )
#graapph file explain
(For each position, a boxplot is drawn with:

#the median value, represented by the central red line
#the inter-quartile range (25-75%), represented by the yellow box
#the 10% and 90% values in the upper and lower whiskers
#the mean quality, represented by the blue line)

4) bowtie2 
The bowtie2-build indexer
 
 #first download reference genome of Arabidopsis thaliana.
#NCBI - datbase - genome - search txid3702[orgn] (Arabidopsis thaliana = GCF_000001735.4_TAIR10.1_genomic.fna.gz)
bowtie2-build [options]* <reference_in> <bt2_base>  (option= -fThe reference input files (specified as <reference_in>) are FASTA files (usually having extension .fa, .mfa, .fna or similar).
bowtie2-build [options]* <reference_in> <bt2_base>

use cammand =   
> bowtie2-build -f GCF_000001735.4_TAIR10.1_genomic.fna.gz bt2_base

#If your computer has multiple processors/cores, use -p

> cat trim_list |while read i; do cd $i ; bowtie2 -x ../../bowtie2/example/index/bt2_base  *.gz -p 20 -S $i.sam; cd .. ; done (for single end)
82.76% overall alignment rate
10012978 reads; of these:
  10012978 (100.00%) were unpaired; of these:
    1136943 (11.35%) aligned 0 times
    7156785 (71.48%) aligned exactly 1 time
    1719250 (17.17%) aligned >1 times
88.65% overall alignment rate
3500711 reads; of these:
  3500711 (100.00%) were unpaired; of these:
    94490 (2.70%) aligned 0 times
    2386248 (68.16%) aligned exactly 1 time
    1019973 (29.14%) aligned >1 times
97.30% overall alignment rate

create sam files
5) Samtools 
SAM file Convert to BAM

#samtools view – views and converts SAM/BAM/CRAM files 

why use samtools-
* covert sam file to bam format
* remove >1 times aligned seq.
* Remove 0 time aligned seq.
Converting SAM to BAM with samtools “view”

#To do anything meaningful with alignment data from BWA or other aligners (which produce text-based SAM output), we need to first convert the SAM to its binary counterpart, BAM format. The binary format is much easier for computer programs to work with. However, it is consequently very difficult for humans to read. More on that later.

#To convert SAM to BAM, we use the samtools view command. We must specify that our input is in SAM format (by default it expects BAM) using the -S option. We must also say that we want the output to be BAM (by default it produces BAM) with the -b option. Samtools follows the UNIX convention of sending its output to the UNIX STDOUT, so we need to use a redirect operator (“>”) to create a BAM file from the output.

samtools view -S -b  input.sam > output.bam (use cammand) 
> cat list |while read i; do cd $i samtools view -S -b *.sam -o $i.bam; cd .. ; done
samtools merge all bam files

 samtools merge [options] -o out.bam [options] in1.bam ... inN.bam

samtools merge [options] out.bam in1.bam ... inN.bam 
 samtools merge [options] -o out.bam [options] in1.bam ... inN.bam
use cammand =     samtools merge output.bam -b list (list all bam file list name)
sort bam file- 
use cammands- samtools sort -@ 25 sample.bam -o sample.sorted.bam(-@ INT Set number of sorting and compression threads. By default, operation is single-threaded. )
6) sambamba use delete dublicate
 use cammands = samtools rmdup -s output.bam rmdup.bam
 
 ( 6) macs2
 > macs2 callpeak -t sorted.bam -f BAM --name=peak1q --outdir macs2_62(file name) -q 0.01
#create foure file name-
(peak1q_model.pdf  , peak1q_model.r  , peak1q_peaks.narrowPeak  ,  peak1q_peaks.xls  peak1q_summits.bed
create model.pdf file= using macs2 tool , 
Rscript NAME_model.r
Rscript peak1_model.r(create model.r)
Rscript peak1q_model.r( all output files start peak1q_)


 Ok, now let’s do the same peak calling for the rest of our samples:
 
 Command line: 
 > callpeak -t sorted_dup.bam -c SRR9211592_1.sorted_control1.bam -f BAM --name=peak1q --outdir macs2_ct -q 0.01(macs2 control data and treated data)
use cammand - Rscript peak1q_model.r 

  > callpeak -t sorted_dup.bam -c SRR9211592_1.sorted_control1.bam -m 2 50 -n macs2_with_control/SRR9211592_1.sorted_control1.bam
  
gtf to bed =   
> sed '1,4d' GCF_000001735.4_TAIR10.1_genomic.gtf |grep -P "\tgene\t"|awk '{print $1"\t"$4-2000"\t"$5"\t"$10"\t255\t"$7}'|sed 's+"++g' |sed 's+;++'
 > sed '1,4d' GCF_000001735.4_TAIR10.1_genomic.gtf |grep -P "\tgene\t"|awk '{print $1"\t"$4"\t"$5"\t"$10"\t"$7}'|sed 's+"++g' |sed 's+;++' >7.bed
 > awk '{if(!($2>$3)) print $0}' 6.bed  >7.bed
bedtools intersect -a peak1q_summits.bed -b 7.bed 

> cat l1  | while read i ;do cd $i ; bedtools intersect -a /home/scbb/Documents/Dissertation_sarita/Chipseq/samtools_out/macs2_out/7.bed -b $i\_sort.bed >$i\_intesect.bed ;cd .. ;done



1) ######intersect first five(TF) intersect file

(cat l1  | while read i ;do cd $i ; cat *.bed |awk '{print $0"\t+"}' >$i.txt ; cd .. ;done
> cat l1  | while read i ;do cd $i ; bedtools sort -i $i.txt |mergeBed -i -  >$i\_sort.bed ; cd .. ;done

> cat l1  | while read i ;do cd $i ; bedtools intersect -a /home/scbb/Documents/Dissertation_sarita/Chipseq/samtools_out/macs2_out/7.bed -b $i\_sort.bed >$i\_intesect.bed ;cd .. ;done)


2) intersect second five(TF)  

> cat lis | while read i ; do bedtools sort -i $i > $i\_sort.bed && bedtools merge -i $i\_sort.bed > merge.bed ; done(ls *.bed > list.bed)
> cat list.bed | while read i ; do bedtools sort -i $i > $i\_sort.bed && bedtools merge -i $i\_sort.bed ABF1_merge.bed ; done
> cat list.bed | while read i ; do bedtools sort -i $i > $i\_sort.bed && bedtools merge -i $i\_sort.bed ABF3_merge.bed ; done
> cat list.bed | while read i ; do bedtools sort -i $i > $i\_sort.bed && bedtools merge -i $i\_sort.bed ABF4_merge.bed ; done
> cat list.bed | while read i ; do bedtools sort -i $i > $i\_sort.bed && bedtools merge -i $i\_sort.bed  ANAC032_merge.bed ; done
> cat list.bed | while read i ; do bedtools sort -i $i > $i\_sort.bed && bedtools merge -i $i\_sort.bed  ANAC102_merge.bed ; done
intersect second five(TF) *.sort.bed file
use cammand - cat l2  | while read i ;do cd $i ; bedtools intersect -a ../7_mod1.bed -b $i\_sort.bed > $i\_intesect.bed ;cd .. ;done




(gtf - 
> cat 7.bed |sed "s/NC_003070.9/Chr1/g" | sed "s/NC_003071.7/Chr2/g" | sed "s/NC_003074.8/Chr3/g" | sed "s/NC_003075.7/Chr4/g" | sed "s/NC_003076.8/Chr5/g"  | head
> cat 7.bed |sed "s/NC_003070.9/Chr1/g" | sed "s/NC_003071.7/Chr2/g" | sed "s/NC_003074.8/Chr3/g" | sed "s/NC_003075.7/Chr4/g" | sed "s/NC_003076.8/Chr5/g"  > 7_mod.bed
> cat 7_mod.bed |grep -P -i "Chr*" > 7_mod1.bed)

correlation


> cat ABF1_intesect.bed |awk '{print $4}' |sort |uniq > ABF1

> cat ABF3_intesect.bed |awk '{print $4}' |sort |uniq > uniq_ABF3
> cat ABF4_intesect.bed |awk '{print $4}' |sort |uniq > uniq_ABF4


> cat ANAC032_intesect.bed |awk '{print $4}' |sort |uniq > uniq_ANAC032
> cat ANAC102_intesect.bed |awk '{print $4}' |sort |uniq > uniq_ANAC102




 


