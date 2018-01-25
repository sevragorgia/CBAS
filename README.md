# CBAS_DE

## Differential expression under contrasting symbiont activity regimes

## Assumptions and requirements

I assume you use a \*nix system. All scripts and analyses ran on a gentoo box. No tests or attempts to test the provided scripts on other operating systems have been done. Support will be provided for linux users **only**. No support will be provided to MacOS or MS Microsoft users; I have no testing environment for these systems. 

I assume you installed git on your computer and want to clone this repository. For this you also need to install **git-lfs**. To run the provided scripts you need:

  - R
  - Perl

If you don't have git-lfs must likely the clone will partially work. All files managed by git-lfs will not be downloaded.

I will try to provide instructions to clone the repo and fetch the large files in this README. 

**Please note** that the scripts provided in this repository are provided without any warranty and may not work on your system. Should that be the case, I'd be happy to help if you use the *Issues* section of this repo to let me know that something is not working for you.

## Description of the folders included.

  - **Annotations**: this folder provides the annotations done for CBAS' reference transcriptome. The following subfolders are provided:
    * All_Bacts: blast results of CBAS' reference transcriptome against **All_Bact**erial genome**s** available at NCBI by the end of 2015. Since sponges harbour a rich microbiome, I blasted against all genomes to see whether a transcript could be of bacterial origin.
    
    * AQU2: blast results of CBAS' reference transcriptome against the *Amphimedon queenslandica* AQU2 protein dataset.
    
    * GOs: GO Terms derived from the best uniprot match to the CBAS' transcripts. Three tables are provided: component, function and process. Each table contains the accesssion of best uniprot match, the CBAS transcript name, and the associated GO terms.
    
    * KEGG: blastKOALA results. Essentially a text file with transcript name and KEGG assignment to an orthology group. **Note** that for this, the transdecoder predicted messenger was used, that is why the transcript names contain |m.XXX (X = {0-9}) at the end.
    
    * Metatable: this folder provide the CBAS' reference transcriptome metatable. This table contains all information about the annotations, transcripts, predicted proteins, etc. for the reference transcriptome. Is a flat csv file that you can use to get (e.g. with grep) all annotations associated with a particular transcript, transcripts containing a given pfam domain, etc. Two versions are provided, the older (in time) version will be removed at some point.
    
    * Pfam: results of pfam_scan.pl
    
    * Transdecoder: transdecoder predicted mRNAs, peptides, CDSs and transdecoder provided track files (not used but provided). The mRNA, peptides and CDSs files are fasta files.
    
    * Uniprot: blast results of CBAS' reference transcriptome against the SWISS-Prot database. These blast results provide the basis for the GO term annotation.
    
  - **Counts**: count matrix used for the DESeq2 analysis. In this folder the RSEM results can be found as well. **Note** that a .metadata text file is provided specifying the details of all the RNA-Seq libraries used. In this file, the *Original_Sample_Name* column specifies the name of the sample as it was used during collection and parts of the lab. work. These names appear in some files and these files **won't** be modified. The reason behind this decision is to keep the data in the repo consistent as far as possible with the Galaxy history associated with the experiments. **Yet, the names given in the counts matrix were changed to make them more readable and easy to handle.** In this matrix, the names used are especified in the column *Name_Used_in_Analysis* in the metadata file. So, you can track everything, I hope. To repeat the DESeq2 analyses the only files needed are: CBAS_Bleaching_RSEM_Expression_Matrix.counts, CBAS_Bleaching_RSEM_Expression_Matrix.info. In this folder you can also find FPKM and TPM matrices.
  
  - **Results**: in this folder a number of results files can be found. In the folder the files *CEGMA_Results*, *Sequencing_QC_Results* and *Trinity_Results.txt* provide an overview of: 1. the completeness of the reference transcriptome as judge by blasting againts the CEGMA gene (**Note** that the results reported in the paper refer to the results of the BUSCO pipeline), 2. the results of the quality control done before assembly and 3. the results provided by trinity (e.g. transcriptome stats). In addition the following subfolders can be found:
  
    * DEGs: provide a number of files containing the list of **D**ifferentially **E**xpressed **G**ene**s** with varied levels of detail. Two different **L**og **F**old **C**hange threshold are provided: |LFC| = 1 and |LFC| = 2. The files starting with OVER_ and UNDER_ provide a copy of the annotations of the best uniprot match for up/down regulated transcripts. I keep this here just to preserve the state of the information I accessed while trying to make sense out of all these data.
    
    * GOs: in this folder, the GO annotation for the DEGs is provided. The background set of transcripts used for the TopGO analysis is also provided here (files starting with Back_). The filenames should be self explanatory. Annotations at two different LFC cutoffs (1 and 2) are provided.
    
    * Pfam: here the pfam annotations for the DEGs and their respective background sets are provided. The filenames should be self explanatory. These files were used for the Pfam enrichment analysis.
  
  - **R_Scripts**: this folder groups a number of R scripts used to conduct some descriptive analyses on the transcriptome, the DEG analysis and the GO/Pfam enrichment analyses. Originally the scripts were independent from each other. Now they are provided as an **Rmarkdown** file that integrates them with a kind of manuscript on the quality of reference transcriptome and the results of the DEG analysis. **Note** that some of the scripts currently provided here will be removed once I am 100% that all the code necessary for all the analyses have been ported to the rmd script. This step only requires a final check, so hopefully these unnecessary, deprecated scripts will be removed soon. Two other files are provided, one is an RData session snapshot and the other includes details of the R environment used, including the versions of the packages loaded. **Note** that different versions of different packages may change the results! Regarding the RData session provided, for this, I ran the rmarkdown file and saved all resulting objects as the RData file provided. It is likely that this RData file will change to a new/final version once all analyses are considered done/closed.
  
  - **Scripts**: here a number of perl scripts (and some old R scripts that will be deleted from this folder) are provided. These were written to do different things. Each script should provided a (moderately) well written help/usage. If any issues/problems are found trying to execute these scripts please use this repository **Issues** tab to ask for feedback/help. I will try to provide a fix or help.
  
  - **Transcripts**: CBAS' reference transcriptome as fasta. A blast database based on this assembly is provided as well.

