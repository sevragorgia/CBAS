# CBAS

##Differential expression under contrasting symbiont activity regimes

##Assumptions and requirements

I assume you use a \*nix system. All scripts and analyses ran on a gentoo linux box. No tests or attempts to test the provided scripts on other systems have been done. Support will be provided for linux users only. No support will be provided to MacOS or MS Microsoft users; I have no testing environment for these systems.

I assume you installed git and clone the repository.

You need:

- git
- git-lfs
- R
- Perl

I will try to provide instructions to clone the repo and fetch the large files in this README. Here, I will also provide a brief description of the scripts provided, etc.Note that this repository is intended to be active and to change with time. It is not a dead-end and I plan to upload new data and analyses if I do any in the future. I also intend to modify the manuscript if necessary.

###Provided Scripts

**CBAS_DeSeq2_Analysis.R**

R script providing the steps followed to derived a list of Differentially Expressed Genes (DEGs) from a matrix of counts.

Using this scripts a background set of transcript (necessary for GO term enrichment analysis) can be derived.

**CBAS_TopGO_Enrichment_Analysis.R**

R script used for the GO term enrichment analysis.

**Get_DEG_GO_Annotations.pl**

Perl script with which the GO term annotations for the DEGs and the transcripts in the background set can be extracted from a set GO annotation files (provided in Annotations/GOs)

**Get_GO_Annotaitons.pl**

Perl script to scan swissprot and retrieve the GO annotations associated with a UNIPROT accession. Using this script one can use the results of a blast annotation table against UNIPROT to get the GO terms associated with the UNIPROT match of a transcript and annotate the transcript this way. NOTE: you are annotating against annotations! So, be careful interpreting the results.



***Deriving lists of GO terms for DEGs and their background distributions***

Once lists of DEGs and background transcripts are available the perl scritp *Get_DEG_GO_Annotations.pl* can be used to get the GO annotations.

This is work in progress!
