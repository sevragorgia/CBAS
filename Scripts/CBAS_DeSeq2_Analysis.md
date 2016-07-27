Changes in gene expression in a common aquarium sponge under different symbiont activity regimes.
================
Sergio Vargas R.
2016-07-26

* * * * *

Preliminary remarks
-------------------

The common blue aquarium sponge is a cyanosponge commonly found in salt water aquaria. Most aquarium hobbiest refer to it as *Collospongia auris*, however its taxonomic affiliation is not clear. I refer to this sponge as the "Common Blue Aquarium Sponge", **CBAS** in short. CBAS is a cyanosponge, it harbours cyanobacteria from the Candidatus *Synechococcus spongiarum* clade.

Culture of the sponge under different light conditions resulted in the observation that when kept under dark conditions the sponge bleaches, turning completely white after ca. 12 weeks under darkness; the "normal" color of the sponge varies between green, blue and violet. Changes in color are accompanied by changes in size and general morphology. For instance, the sponges form elongated projections; the normal shape of the sponge is lobated and plate like.

It is not clear whether the cyanobacteria are still present in the bleached sponge tissues. Bioanalyzer profiles of RNA extracts from bleached sponges clearly showed that the samples lack the rRNA bands belonging the cyanobacteria but qPCR experiments pointed towards the presence of cyanobacteria in DNA extracts from sponges exposed to max. 6 weeeks darkeness [1]. However, this evidence is (in my opinion) not conclusive due to a number of technical problems associated to these experiments; RNA extractions only demonstrate that the symbionts are not active and the qPCR experiments suffered from incosistent amplifying of the samples. Another piece of evidence (weakly) pointing towards this is the capacity of the sponges to regain color after 8-12 weeks re-exposure to light; symbiont uptake from the water column cannot be excluded, however, apparently *S. spongiarum* is not present in the water. All in all, the evidence at hand point to the presence of the symbionts in some kind of inactive (transcriptionl/metabolic) state in the sponge tissues exposed to dark conditions.

Study rational
--------------

Since it is possible to manipulate the trascriptional/metabolic state of the cyanobacterial symbionts through the exposure of the system to darkness, CBAS represents an interesting system to assess how the sponge reacts to symbiont inactivation. I use RNA-Seq data to *de novo* assemble a reference transcriptome for this system and to assess changes in the sponge gene expression profiles associated with different symbiont transcriptional states. I hope these data contributes to our understanding of the molecular mechanisms used by sponges to interact with their microbiomes.

* * * * *

Introduction
============

Materials and Methods
=====================

Sponge culture and experimental setting
---------------------------------------

Sponges are kept in a 200L salt water aquarium under a 12hr light:12hr darkness regime using two T5 fluorescent lights. For all experiments one large sponge was used and a hole-puncher was used to produce explants of 10 mm diameter [2]. The explants were place in a recipient with sand and were allow to recover for 3-5 days under control light conditions. Recovery was visually assesed by inspecting the border of the explants to corroborate the tissue healed at the cutting points. After this, the explants were moved to a section of the aquarium covered with a two sheets of black plastic cloth. The plastic cloth serves to block the light from the T5 lights; the amount of light (in klux) penetrating the plastic sheets was \~0 klux while under normal conditions at least \~15 klux can be measured at the water surface. Twelve weeks after moving the sponges to the dark side of the aquarium, the explants were flash-frozen in liquid nitrogen and kept at -80 °C; at this point in time explants kept under control conditions were also sampled. This experiment was repeated several times[3]; the resulting tissue samples were used to test different RNA extraction protocols and assess whether bleaching can be reproduced in the aquarium.

RNA extraction, libray preparation and sequencing
-------------------------------------------------

Total RNA was extracted using a hybrid protocol that combines a regular CTAB extraction and a spin-column clean-up step. Briefly, using mortar and pestle, samples were pulverized in liquid nitrogen. The tissue powder was lysed in 600 \(\mu\)L warm (56 °C) CTAB-PVP-NaCl buffer containing \(\beta\)-mercaptoethanol for 15-20 minutes with agitation (\~500 rpm). After lysis, one volume acidic Phenol-Chloroform-Isoamyl alcohol (25:24:1) was added to the samples and the extraction was vigorously agitated until the liquid had a milky appearance. The samples were then centrifuged at 14,000 rpm for 15 minutes to separate the phases. 400 \(\mu\)L of the polar phase were then transfered to a new 2mL microcentrifuge tube and the nucleic-acids were precipitated for 10 minutes using one volume isopropanol. The nucleic acids were recovered by centrifugation (14,000rpm for 20 minutes at \~16°C) and the recovered pellets were washed twice with 1 mL cold 80% ethanol. After washing, ethanol droplets were removed with a 10 \(\mu\)L pipette and the pellets were air dried for \~10 minute before resuspending in 30 \(\mu\)L nuclease-free water.

After this initial extraction, the ZR-Duet DNA/RNA MiniPrep [4] was used, following the manufacturers recommendations, to separate DNA and RNA from the sample. The RNA was resuspended in 30 \(\mu\)L nuclease-free water and quality checked initially on 1% agarose gels and finally on a Bioanalyzer 2100 Nano RNA chip. RNA concentration was measured on a Nanodrop 1000. RIN values could be calculated for bleached samples only; control samples have four rRNA peaks corresponding to the 16S, 18S, 23S and 28S rRNA fragments and the Bioanalyzer canno calculate their RIN value. Quality of the control samples was assessed by overlaying bleached samples with RIN and control samples without RIN and assessing their similarity in terms of peak height and baseline level. Five control and four bleached samples were sent on dry-ice to the EMBL Genomics Core Facility [5] where they were used to produce strand-specific libraries with \~110 base pairs (bp). These libraries were multiplexed and pair-end sequenced (50bp reads) in two lanes of a HiSeq 2500 (Illumina).

Transcriptome assembly, annotation
----------------------------------

Reads were quality controled using FastQC [6] and filtered using the BioLite program filter\_illumina.cpp [7]. The surviving read pairs from all libraries were concatenated to produce two fastq files that were used for *de novo* transcriptome assembly in Trinity v2.0.6 (using the --normalize\_reads flag). The resulting contigs (with lenght \>=200bp) were annotated against the Uniprot (SwissProt) [8] and *Amphimedon queenslandica* isoforms (AQU2 proteins) [9] using blastx v2.2.29+ with an expectation cutoff of 0.001. Only the best match per contig was saved and the blast results were saved using blast's XML format converted to a 25 column table [10]. Both tables (i.e. [UNIPROT](https://github.com/sevragorgia/CBAS/tree/master/Annotations/Uniprot) and [AQU2](https://github.com/sevragorgia/CBAS/tree/master/Annotations/AQU2)) are available for download in the [project repository](https://github.com/sevragorgia/CBAS). Gene Ontology [11] annotations for the CBAS transcriptome were obtained by programmatically querying the [QuickGO Webservice](http://www.ebi.ac.uk/QuickGO/) with the transcriptome's UNIPROT annotations [12] and a custom [perl script](https://github.com/sevragorgia/CBAS/blob/master/Scripts/Get_GO_Annotations.pl). For each CBAS transcript the "component", "function" and "process" GO terms associated with its UNIPROT best match were stored in the [project repository](https://github.com/sevragorgia/CBAS/tree/master/Annotations/GOs) as independent tab-separated files that can be easily modified to use as input files in TopGO [13].

In addition, transcripts were translated using the program TransDecoder.LongOrfs [14] and the resulting cds, mRNA, bed, gff3 and pep files stored in the [project repository](https://github.com/sevragorgia/CBAS/tree/master/Annotations/Transdecoder) and used to annotated the transcriptome against the [Pfam](http://pfam.xfam.org/) and [KEEG](http://www.genome.jp/kegg/) databases. For this, the perl script [pfam\_scan.pl](ftp://ftp.ebi.ac.uk/pub/databases/Pfam/Tools/) and the webservice [BlastKOALA](http://www.kegg.jp/blastkoala/) were used; the resulting output files can be found in the [project repository.](https://github.com/sevragorgia/CBAS) Finally, the assembled transcripts were blasted (blastn) against a bacterial genomes database [15].

All annotations were used to build an [annotation meta-table](https://github.com/sevragorgia/CBAS/blob/master/Annotations/Metatable/) using a [custom made perl script](https://github.com/sevragorgia/CBAS/blob/master/Scripts/Create_Transcriptome_Annotation_Table.pl) available in the project repository.

Transcriptome completeness assessment
-------------------------------------

The CBAS transcriptome completeness was assessed by blasting against the CEGMA gene set of [Parra et al. 2007](http://bioinformatics.oxfordjournals.org/content/23/9/1061.abstract)[16] using tblastn with an e-value of \(1e^-19\) as implemented in the scrtipt [find\_cegma\_genes.py](https://bitbucket.org/wrf/galaxy-files/src). For details about the method see [17].

Differential gene expression analysis
-------------------------------------

Using the *de novo* assembled transcriptome, the individual libraries (i.e. control and bleached sponges) were mapped with RSEM and a transcript by sample count matrix was derived from the "gene" counts with the program [abundance\_estimates\_to\_matrix.pl](https://github.com/sevragorgia/CBAS/tree/master/Counts/RSEM) provided as part of the Trinity package [18]. The [count matrix](https://github.com/sevragorgia/CBAS/tree/master/Counts/) was used to find differentially expressed genes in control *vs.* bleached sponges (model = \~Treatment). For this, the package [DeSeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) was used. In brief, transcripts with 0 counts over all samples were removed from the matrix and Principal Component Analysis was used to assess the global expression pattern of control *vs.* bleached sponges and identify potential outliers. Size effects and dispersions were estimated with the DeSe2 methods **estimateSizeFactors** and **estimateDispersions** and differentially expressed transcripts were then infered using a Wald test (DeSeq2 method **nbinomWaldTest**). The resulting p-values were adjusted using the Benjamin-Hochberg correction. Transcripts with **log fold change \(< -2\)** or **log fold change \(> 2\)** and **adjusted p-value \(< 0.01\)** were considered as differentially expressed for further analyses.

A stand-alone R script to replicate the analysis is available in the [project repository](https://github.com/sevragorgia/CBAS/blob/master/Scripts) and can be run using the [count matrix](https://github.com/sevragorgia/CBAS/tree/master/Counts/CBAS_Bleaching_RSEM_Expression_Matrix.counts) and [sample information](https://github.com/sevragorgia/CBAS/tree/master/Counts/CBAS_Bleaching_RSEM_Expression_Matrix.info) provided. And in this version of the manuscript, the code used to calculate something is provided as as R markdown snippets in the Results section.

Gene Ontology term enrichment analysis
--------------------------------------

Pfam domain enrichment analysis
-------------------------------

Results
=======

Discussion
==========

Acknowledgements
================

### Conflict of interest

None declared.

### Author contributions

Sergio Vargas designed the study, conducted the experiments and analyzed the data. Gert Wörheide contributed reagents. Sergio Vargas and Gert Wörheide wrote the manuscript. Sergio Vargas manages the [project repository](https://github.com/sevragorgia/CBAS/).

### Disclaimer

None of the products or companies listed are endorsed or recommended in anyway by the author(s) of this text.

### Scripts and Data availability

All scripts used to analyze the data, as well as some miscellaneous scripts used to, for instance, prepare the annotation table or generate GO-term/Pfam input files for the enrichment analyses can be found in the [project repository](https://github.com/sevragorgia/CBAS/).

References and Notes
====================

[1] I have (repeatedly) failed to prepare libraries for Next Generation Sequencing from bleached tissue, probably due to the presence of secondary metabolites in the extracts.

[2] The 10 mm hole puncher has been replaced by a 6 mm version. [Click here to see the product.](http://www.sigmaaldrich.com/catalog/product/sigma/z708909?lang=de&region=DE)

[3] We, in fact, keep a permanent culture of "bleached" explants that are allowed to bleach for at least 12 weeks, are fixed in liquid nitrogen and stored at -80 °C.

[4] Zymo Research. [Click here to see the product website.](https://www.zymoresearch.de/rna/dna-rna-co-purification/cells-tissue-rna/zr-duet-dna-rna-miniprep)

[5] EMBL GeneCore. [<http://genecore3.genecore.embl.de/genecore3/index.cfm>](http://genecore3.genecore.embl.de/genecore3/index.cfm)

[6] FastQC. [<http://www.bioinformatics.babraham.ac.uk/projects/fastqc/>](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

[7] BioLite: [<https://bitbucket.org/caseywdunn/biolite/overview>](https://bitbucket.org/caseywdunn/biolite/overview)

[8] Downloaded early 2015. [Uniprot download.](http://www.uniprot.org/downloads)

[9] AQU2 isoforms (proteins) [*Amphimedon queenslandica* transcriptome resource.](http://amphimedon.qcloud.qcif.edu.au/downloads.html)

[10] Blast was run from [Galaxy](https://galaxyproject.org/) and the XML file was converted to a 25 column table by the [ncbi\_blast\_plus](http://toolshed.g2.bx.psu.edu/view/devteam/ncbi_blast_plus) tool.

[11] Gene Ontology Consortium. [<http://geneontology.org/>](http://geneontology.org/)

[12] The UNIPROT accession code of the blast best match associated to each transcript in the CBAS transcriptome was used to query the GO annotations associated with the UNIPROT protein. These GO annotations were used to annotated the CBAS transcripts; note that these are annotations based on annotations and should be taken with precaution.

[13] [TopGO](https://bioconductor.org/packages/release/bioc/html/topGO.html) expects GO terms in tab separated {Trascript, GoAnnotation} pairs. The script produces output with a leading column indicating the UNIPROT accession code used for each transcript, which comes handy when controling the correct annotation was used for each transcript and can be easily removed for further analyses.

[14] Transdecoder's Github Repository can be found [here.](https://transdecoder.github.io/)

[15] I honestly do not know what the purpose of this blast run was, I did it for sake of completeness and because the reads were not filtered before assembling the conting. Yet, after seeing the results, it seems as if not too many contigs matched the bacterial database...

[16] Parra et al. 2007. CEGMA: a pipeline to accurately annotate core genes in eukaryotic genomes. [Bioinformatics 23 (9): 1061-1067.](http://dx.doi.org/10.1093/bioinformatics/btm071)

[17] Francis et al. 2013. A comparison across non-model animals suggests an optimal sequencing depth for *de novo* transcriptome assembly. [BMC Genomics 14:167](http://dx.doi.org/10.1186/1471-2164-14-167)

[18] the script is also provided as part of this repository. The Copyright (c) of this script belongs trinityrnaseq (2014), who reserves all rights upon the code. [See the full LICENSE here.](https://github.com/trinityrnaseq/trinityrnaseq/blob/master/LICENSE)