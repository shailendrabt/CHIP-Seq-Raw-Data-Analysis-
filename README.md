# CHIP-Seq-Raw-Data-Analysis Workflow 
<img align="right" alt="coding" width ="800" src= "https://crc-pages.pitt.edu/user-manual/_assets/img/advanced-genomics-support/chipseq_02.png">
















##########################################
### Chapter 1: Introduction to ChIP-seq ###
###########################################

Chromatin immunoprecipitation followed by sequencing (ChIP-seq) analysis is a key technology in epigenomic research. This method uses an antibody for a specific DNA-binding protein or a histone modification to identify enriched loci within a genome


##################################
### Chapter 2: Getting Started ###
##################################
##### Trim Galore (optional)

 Run Trim Galore (~10 min)
trim_galore -q 20 --stringency 2 -o data data/${sample}.fastq.gz


##### Data retrieval from GEO
 Should be adapted for all datasets

Data from Domcke/Bardet et al. Nature 2015
ChIP-seq for the TF NRF1 in mouse embryonic stem cells either wild type (WT) or with no DNA methylation (TKO; triple knockout for DNMT1, DNMT3a, DNMT3b)
# Available on GEO GSE67867

# NRF1_CHIP_WT_1 - GSM1891641 - SRR2500883 - single-end
Mapping Cammand LIne -
indexing file -
Download mm10 mouse genome file in terminal 
wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz -O genomes/mm10/chromFa.tar.gz
tar -zxvf genomes/mm10/chromFa.tar.gz -C genomes/mm10
# Keep only the conventional chromosomes (chr1-19,X,Y,M)
rm genomes/mm10/chromFa.tar.gz genomes/mm10/*random* genomes/mm10/chrUn*
 
# Generate Bowtie 2 index
> cd genomes/mm10/
 > bowtie2-build chr1.fa,chr10.fa,chr11.fa,chr12.fa,chr13.fa,chr14.fa,chr15.fa,chr16.fa,chr17.fa,chr18.fa,chr19.fa,chr2.fa,chr3.fa,chr4.fa,chr5.fa,chr6.fa,chr7.fa,chr8.fa,chr9.fa,chrM.fa,chrX.fa,chrY.fa ../../indices/mm10/mm10
> cd ../..
> Mapping- bowtie2 -x  indices/mm10/mm10 -U data/SRR2500883.fastq.trimmed.gz/SRR2500883_trimmed.fq.gz > reads/83.sam











Learning Objectives
Annotating peaks with gene and genomic feature information
Obtaining biological context for identified binding sites using functional enrichment tools
Using the seqqunce data for peaks to identify possible over-represented motifs
Annotation and Functional Analysis
We have identified regions of the genome that are enriched in the number of aligned reads for each of our transcription factors of interest, Nanog and Pou5f1. These enriched regions represent the likely locations of where these proteins bind to the genome. After we obtain a list of peak coordinates, it is important to study the biological implications of the protein–DNA bindings. Certain questions have always been asked: what are the genomic annotations and the functions of these peak regions?
Peak Annotation
Because many cis-regulatory elements are close to TSSs of their targets, it is common to associate each peak to its nearest gene, either upstream or downstream. ChIPseeker is an R Bioconductor package for annotating peaks. Additionally, it has various visualization functions to assess peak coverage over chromosomes and profiles of peaks binding to TSS regions.

Loading data
Peak annotation is generally performed on your high confidence peak calls (after looking at concordance betwee replicates). While we have a confident peak set for our data, this set is rather small and will not result in anything meaningful in our functional analyses. We have generated a set of high confidence peak calls using the full dataset. These were obtained post-IDR analysis, (i.e. concordant peaks between replicates) and are provided in BED format which is optimal input for the ChIPseeker package.
<img align="right" alt="coding" width ="800" src="https://a.storyblok.com/f/196663/640x360/86b5f57ecd/nanopore-sequencing-animation.gif">
Annotation
Many annotation tools use nearest gene methods for assigning a peak to a gene in which the algorithm looks for the nearest TSS to the given genomic coordinates and annotates the peak with that gene. This can be misleading as binding sites might be located between two start sites of different genes.

The annotatePeak function, as part of the ChIPseeker package, uses the nearest gene method described above but also provides parameters to specify a max distance from the TSS. For annotating genomic regions, annotatePeak will not only give the gene information but also reports detail information when genomic region is Exon or Intron. For instance, ‘Exon (uc002sbe.3/9736, exon 69 of 80)’, means that the peak overlaps with the 69th exon of the 80 exons that transcript uc002sbe.3 possess and the corresponding Entrez gene ID is 9736.
Over-representation analysis
There are a plethora of functional enrichment tools that perform some type of over-representation analysis by querying databases containing information about gene function and interactions. Querying these databases for gene function requires the use of a consistent vocabulary to describe gene function. One of the most widely-used vocabularies is the Gene Ontology (GO). This vocabulary was established by the Gene Ontology project, and the words in the vocabulary are referred to as GO terms.

We will be using clusterProfiler to perform over-representation analysis on GO terms associated with our list of significant genes. The tool takes as input a significant gene list and a background gene list and performs statistical enrichment analysis using hypergeometric testing. The basic arguments allow the user to select the appropriate organism and GO ontology (BP, CC, MF) to test.
Functional enrichment: Web-based tools
There are also web-based tool for enrichment analysis on genomic regions, and a popular one is GREAT (Genomic Regions Enrichment of Annotations Tool). GREAT is used to analyze the functional significance of cis-regulatory regions identified by localized measurements of DNA binding events across an entire genome1. It incorporates annotations from 20 different ontologies and is an easy to use tool which generates annotation and downstream functional enrichement results for genomic coordinate files. The utility of GREAT is not limited to ChIP-seq, as it could also be applied to open chromatin, localized epigenomic markers and similar functional data sets, as well as comparative genomics sets.

In the interest of time we will not go into the details of using GREAT, however we have materials linked here if you are interested in testing it out with this dataset. There also demo datasets on the GREAT website that you can use to test out the functionality of the tool.

