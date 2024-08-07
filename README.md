# CHIP-Seq-Raw-Data-Analysis Workflow
<img align="right" alt="coding" width ="800" src="https://hbctraining.github.io/Intro-to-ChIPseq/img/chip_workflow_june2017_step5.png">





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

Annotation
Many annotation tools use nearest gene methods for assigning a peak to a gene in which the algorithm looks for the nearest TSS to the given genomic coordinates and annotates the peak with that gene. This can be misleading as binding sites might be located between two start sites of different genes.

The annotatePeak function, as part of the ChIPseeker package, uses the nearest gene method described above but also provides parameters to specify a max distance from the TSS. For annotating genomic regions, annotatePeak will not only give the gene information but also reports detail information when genomic region is Exon or Intron. For instance, ‘Exon (uc002sbe.3/9736, exon 69 of 80)’, means that the peak overlaps with the 69th exon of the 80 exons that transcript uc002sbe.3 possess and the corresponding Entrez gene ID is 9736.
Over-representation analysis
There are a plethora of functional enrichment tools that perform some type of over-representation analysis by querying databases containing information about gene function and interactions. Querying these databases for gene function requires the use of a consistent vocabulary to describe gene function. One of the most widely-used vocabularies is the Gene Ontology (GO). This vocabulary was established by the Gene Ontology project, and the words in the vocabulary are referred to as GO terms.

We will be using clusterProfiler to perform over-representation analysis on GO terms associated with our list of significant genes. The tool takes as input a significant gene list and a background gene list and performs statistical enrichment analysis using hypergeometric testing. The basic arguments allow the user to select the appropriate organism and GO ontology (BP, CC, MF) to test.
Functional enrichment: Web-based tools
There are also web-based tool for enrichment analysis on genomic regions, and a popular one is GREAT (Genomic Regions Enrichment of Annotations Tool). GREAT is used to analyze the functional significance of cis-regulatory regions identified by localized measurements of DNA binding events across an entire genome1. It incorporates annotations from 20 different ontologies and is an easy to use tool which generates annotation and downstream functional enrichement results for genomic coordinate files. The utility of GREAT is not limited to ChIP-seq, as it could also be applied to open chromatin, localized epigenomic markers and similar functional data sets, as well as comparative genomics sets.

In the interest of time we will not go into the details of using GREAT, however we have materials linked here if you are interested in testing it out with this dataset. There also demo datasets on the GREAT website that you can use to test out the functionality of the tool.

