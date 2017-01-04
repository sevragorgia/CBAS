Changes in gene expression in a common aquarium sponge under different symbiont activity regimes.
================
Sergio Vargas R.
2016-07-26

------------------------------------------------------------------------

Preliminary remarks
-------------------

The common blue aquarium sponge is a cyanosponge commonly found in salt water aquaria. Most aquarium hobbiest refer to it as *Collospongia auris*, however its taxonomic affiliation is not clear. I refer to this sponge as the "Common Blue Aquarium Sponge", **CBAS** in short. CBAS is a cyanosponge, it harbours cyanobacteria from the Candidatus *Synechococcus spongiarum* clade.

Culture of the sponge under different light conditions resulted in the observation that when kept under dark conditions the sponge bleaches, turning completely white after ca. 12 weeks under darkness; the "normal" color of the sponge varies between green, blue and violet. Changes in color are accompanied by changes in size and general morphology. For instance, the sponges form elongated projections; the normal shape of the sponge is lobated and plate like.

It is not clear whether the cyanobacteria are still present in the bleached sponge tissues. Bioanalyzer profiles of RNA extracts from bleached sponges clearly showed that the samples lack the rRNA bands belonging the cyanobacteria but qPCR experiments pointed towards the presence of cyanobacteria in DNA extracts from sponges exposed to max. 6 weeeks darkeness [1]. However, this evidence is (in my opinion) not conclusive due to a number of technical problems associated to these experiments; RNA extractions only demonstrate that the symbionts are not active and the qPCR experiments suffered from incosistent amplifying of the samples. Another piece of evidence (weakly) pointing towards this is the capacity of the sponges to regain color after 8-12 weeks re-exposure to light; symbiont uptake from the water column cannot be excluded, however, apparently *S. spongiarum* is not present in the water. All in all, the evidence at hand point to the presence of the symbionts in some kind of inactive (transcriptionl/metabolic) state in the sponge tissues exposed to dark conditions.

Study rational
--------------

Since it is possible to manipulate the trascriptional/metabolic state of the cyanobacterial symbionts through the exposure of the system to darkness, CBAS represents an interesting system to assess how the sponge reacts to symbiont inactivation. I use RNA-Seq data to *de novo* assemble a reference transcriptome for this system and to assess changes in the sponge gene expression profiles associated with different symbiont transcriptional states. I hope these data contributes to our understanding of the molecular mechanisms used by sponges to interact with their microbiomes.

------------------------------------------------------------------------

Introduction
============

Materials and Methods
=====================

Sponge culture and experimental setting
---------------------------------------

Sponges are kept in a 200L salt water aquarium under a 12hr light:12hr darkness regime using two T5 fluorescent lights. For all experiments one large sponge was used and a hole-puncher was used to produce explants of 10 mm diameter [2]. The explants were place in a recipient with sand and were allow to recover for 3-5 days under control light conditions. Recovery was visually assesed by inspecting the border of the explants to corroborate the tissue healed at the cutting points. After this, the explants were moved to a section of the aquarium covered with a two sheets of black plastic cloth. The plastic cloth serves to block the light from the T5 lights; the amount of light (in klux) penetrating the plastic sheets was ~0 klux while under normal conditions at least ~15 klux can be measured at the water surface. Twelve weeks after moving the sponges to the dark side of the aquarium, the explants were flash-frozen in liquid nitrogen and kept at -80 °C; at this point in time explants kept under control conditions were also sampled. This experiment was repeated several times[3]; the resulting tissue samples were used to test different RNA extraction protocols and assess whether bleaching can be reproduced in the aquarium.

RNA extraction, libray preparation and sequencing
-------------------------------------------------

Total RNA was extracted using a hybrid protocol that combines a regular CTAB extraction and a spin-column clean-up step. Briefly, using mortar and pestle, samples were pulverized in liquid nitrogen. The tissue powder was lysed in 600 \(\mu\)L warm (56 °C) CTAB-PVP-NaCl buffer containing \(\beta\)-mercaptoethanol for 15-20 minutes with agitation (~500 rpm). After lysis, one volume acidic Phenol-Chloroform-Isoamyl alcohol (25:24:1) was added to the samples and the extraction was vigorously agitated until the liquid had a milky appearance. The samples were then centrifuged at 14,000 rpm for 15 minutes to separate the phases. 400 \(\mu\)L of the polar phase were then transfered to a new 2mL microcentrifuge tube and the nucleic-acids were precipitated for 10 minutes using one volume isopropanol. The nucleic acids were recovered by centrifugation (14,000rpm for 20 minutes at ~16°C) and the recovered pellets were washed twice with 1 mL cold 80% ethanol. After washing, ethanol droplets were removed with a 10 \(\mu\)L pipette and the pellets were air dried for ~10 minute before resuspending in 30 \(\mu\)L nuclease-free water.

After this initial extraction, the ZR-Duet DNA/RNA MiniPrep [4] was used, following the manufacturers recommendations, to separate DNA and RNA from the sample. The RNA was resuspended in 30 \(\mu\)L nuclease-free water and quality checked initially on 1% agarose gels and finally on a Bioanalyzer 2100 Nano RNA chip. RNA concentration was measured on a Nanodrop 1000. RIN values could be calculated for bleached samples only; control samples have four rRNA peaks corresponding to the 16S, 18S, 23S and 28S rRNA fragments and the Bioanalyzer cannot calculate their RIN value. Quality of the control samples was assessed by overlaying bleached samples with RIN and control samples without RIN and assessing their similarity in terms of peak height and baseline level. Five control and four bleached samples were sent on dry-ice to the EMBL Genomics Core Facility [5] where they were used to produce strand-specific libraries with ~110 base pairs (bp). These libraries were multiplexed and pair-end sequenced (50bp reads) in two lanes of a HiSeq 2500 (Illumina).

Transcriptome assembly, annotation
----------------------------------

Reads were quality controled using FastQC [6] and quality filtered using the BioLite program filter\_illumina.cpp [7]. The surviving read pairs from all libraries were used to produce (using the command cat) two fastq files that were used for *de novo* transcriptome assembly in Trinity v2.0.6 (using the --normalize\_reads flag). The resulting contigs (with lenght &gt;=200bp) were annotated against the Uniprot (SwissProt) [8] and *Amphimedon queenslandica* isoforms (Aqu2 proteins) [9] using blastx v2.2.29+ with an expectation cutoff of 0.001. Only the best match per contig was saved and the blast results were saved using blast's XML format converted to a 25 column table [10]. Both tables (i.e. [UNIPROT](https://github.com/sevragorgia/CBAS/tree/master/Annotations/Uniprot) and [Aqu2](https://github.com/sevragorgia/CBAS/tree/master/Annotations/AQU2)) are available for download in the [project repository](https://github.com/sevragorgia/CBAS). Gene Ontology [11] annotations for the CBAS transcriptome were obtained by programmatically querying the [QuickGO Webservice](http://www.ebi.ac.uk/QuickGO/) with the transcriptome's UNIPROT annotations [12] and a custom [perl script](https://github.com/sevragorgia/CBAS/blob/master/Scripts/Get_GO_Annotations.pl). For each CBAS transcript the "component", "function" and "process" GO terms associated with its UNIPROT best match were stored in the [project repository](https://github.com/sevragorgia/CBAS/tree/master/Annotations/GOs) as independent tab-separated files that can be easily modified to use as input files in TopGO [13].

In addition, transcripts were translated using the program TransDecoder.LongOrfs [14] and the resulting cds, mRNA, bed, gff3 and pep files stored in the [project repository](https://github.com/sevragorgia/CBAS/tree/master/Annotations/Transdecoder) and used to annotate the transcriptome against the [Pfam](http://pfam.xfam.org/) and [KEEG](http://www.genome.jp/kegg/) databases. For this, the perl script [pfam\_scan.pl](ftp://ftp.ebi.ac.uk/pub/databases/Pfam/Tools/) and the webservice [BlastKOALA](http://www.kegg.jp/blastkoala/) were used; the resulting output files can be found in the [project repository.](https://github.com/sevragorgia/CBAS) Finally, the assembled transcripts were blasted (blastn) against a bacterial genomes database [15].

All annotations were used to build an [annotation meta-table](https://github.com/sevragorgia/CBAS/blob/master/Annotations/Metatable/) using a [custom made perl script](https://github.com/sevragorgia/CBAS/blob/master/Scripts/Create_Transcriptome_Annotation_Table.pl) available in the project repository.

Transcriptome completeness assessment
-------------------------------------

The CBAS transcriptome completeness was assessed by blasting against the CEGMA gene set of [Parra et al. 2007](http://bioinformatics.oxfordjournals.org/content/23/9/1061.abstract)[16] using tblastn with an e-value of \(1e^-19\) as implemented in the scrtipt [find\_cegma\_genes.py](https://bitbucket.org/wrf/galaxy-files/src). For details about the method see [17].

Differential gene expression analysis
-------------------------------------

Using the *de novo* assembled transcriptome, the individual libraries (i.e. control and bleached sponges) were mapped with RSEM and a transcript by sample count matrix was derived from the "gene" counts with the program [abundance\_estimates\_to\_matrix.pl](https://github.com/sevragorgia/CBAS/tree/master/Counts/RSEM) provided as part of the Trinity package [18]. The [count matrix](https://github.com/sevragorgia/CBAS/tree/master/Counts/) was used to find differentially expressed genes in control *vs.* bleached sponges (model = ~Treatment). For this, the package [DeSeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) was used. In brief, transcripts with 0 counts over all samples were removed from the matrix and Principal Component Analysis was used to assess the global expression pattern of control *vs.* bleached sponges and identify potential outliers. Size effects and dispersions were estimated with the DeSe2 methods **estimateSizeFactors** and **estimateDispersions** and differentially expressed transcripts were then infered using a Wald test (DeSeq2 method **nbinomWaldTest**). The resulting p-values were adjusted using the Benjamin-Hochberg correction. Transcripts with **log fold change \(< -2\)** or **log fold change \(> 2\)** and **adjusted p-value \(< 0.01\)** were considered as differentially expressed for further analyses.

A stand-alone R script to replicate the analysis is available in the [project repository](https://github.com/sevragorgia/CBAS/blob/master/Scripts) and can be run using the [count matrix](https://github.com/sevragorgia/CBAS/tree/master/Counts/CBAS_Bleaching_RSEM_Expression_Matrix.counts) and [sample information](https://github.com/sevragorgia/CBAS/tree/master/Counts/CBAS_Bleaching_RSEM_Expression_Matrix.info) provided. And in this version of the manuscript, the code used for several calculations is provided embedded as Rmarkdown snippets in the source code of this document.

Gene Ontology Term enrichment analysis
--------------------------------------

In order to assess whether certain Gene Ontology terms are "enriched" among the set of differentially expressed genes a GO-term enrichment analysis was done using the R package [TopGO](https://bioconductor.org/packages/release/bioc/html/topGO.html). For this analysis, "background sets" of not differentially expressed transcripts with similar expression profile to the set of differentially expressed transcripts were created using the method **genefinder** of the package [genefilter](http://bioconductor.org/packages/release/bioc/html/genefilter.html) with the Manhattan distance. Background transcript sets were created for all differentially expressed transcripts (i.e. transcripts with log fold change \(> |2|\) and \(p<0.01\)), for differentially upregulated genes (i.e. transcripts with log fold change \(> 2\) and \(p<0.01\)) and for differentially downregulated genes (i.e. transcripts with log fold change \(< -2\) and \(p<0.01\)). The background sets were then used to test for enrichment of certain GO-terms in the global set of differentially expressed transcripts (log fold change \(> |2|\) and \(p<0.01\)), the set of upregulated genes (log fold change \(> 2\) and \(p<0.01\)) and the set of downregulated genes (log fold change \(< -2\) and \(p<0.01\)). For the enrichment analysis, a Fisher exact test with the *classic* algorithm implemented in the method **runTest** of [TopGo](https://bioconductor.org/packages/release/bioc/html/topGO.html) was used.

Pfam domain enrichment analysis
-------------------------------

In addition to the GO-Term enrichment analysis, and due to the difficulties in assigning a meaning to the GO annotations in our sponge model, a Pfam domain enrichment analysis was also done. Essentially the Pfam enrichment analysis works in a similar way that the GO-term enrichment analysis. First, a background set of genes not differentially expressed but showing similar expression patterns than the DEGs is calculated and used to provide a population of transcript against which the differentially expressed transcripts can be compared. Once the background population of transcripts is available, the frequency with which a given domain is found in the background set of transcripts is compared with the frequency with which the same domain is found in the set of differentially expressed transcripts using an hypergeometric test. The background distribution used for the Pfam enrichment domains was the same used for the analysis of GO-term enrichment.

Results
=======

The transcriptome of the Common Blue Aquarium Sponge, a model cyanosponge
-------------------------------------------------------------------------

Per library, we obtained on average 26,385,834 (\(\pm 3,360,778\)) pairs of reads (50 bp long). After cleaning, each library had an average of 24,942,084 (\(\pm 3,420,083\)) surviving pairs. The concatenated dataset thus had 199,536,668 pairs. Table 1 contains the information on the number of sequenced and surviving reads per library.

**Table 1.** Reads by library before and after quality filtering.

|   Library   | Condition | Sequenced reads | Reads post cleaning | % Kept |
|:-----------:|:---------:|:---------------:|:-------------------:|:------:|
|      1      |  Control  |    33,836,615   |      32,397,726     |  95.75 |
|      2      |  Control  |    26,740,977   |      25,252,825     |  94.43 |
|      3      |  Control  |    22,823,936   |      21,016,976     |  92.08 |
|      4      |  Control  |    24,166,051   |      22,655,601     |  93.75 |
|      5      |  Bleached |    26,579,563   |      25,067,523     |  94.31 |
|      6      |  Bleached |    25,448,522   |      24,229,217     |  95.21 |
|      7      |  Bleached |    22,823,936   |      21,016,976     |  92.08 |
|      8      |  Bleached |    24,166,051   |      22,655,601     |  93.75 |
| Total Reads |           |   211,086,672   |     199,536,668     |  94.53 |

------------------------------------------------------------------------

*De novo* transcriptome assembly using the concatenated dataset resulted in 128,686 transcript \(\gt 200\)bp. The N50 of the transcriptome was 1,281 bp, and the median lenght of the assembled transcripts was 453 bp (MAD=0.2078891). More details about the assembly can be found in Table 2 and the length distribution of the transcripts can be visualized in Figure 1.

**Table 2.** Statistics for the *de novo* assembled CBAS transcriptome

|       Statistic       | Obtained value |
|:---------------------:|:--------------:|
|       GC content      |      40.0      |
|          N50          |      1,281     |
|      Max. length      |     19,309     |
|      Mean length      |       803      |
|     Median length     |       453      |
|      Min. length      |       224      |
| Number of Transcripts |     128,686    |

------------------------------------------------------------------------

<img src="CBAS_DeSeq2_Analysis_files/figure-markdown_github/transcript_length_fig-1.png" style="display: block; margin: auto;" />

**Figure 1.** Transcript length distribution of CBAS reference transcriptome. Maximum length allowed equal 10,000 base-pairs.

------------------------------------------------------------------------

Of the collection of transcripts, 37.75% could be translated into proteins by Transdecoder. The majority of the translated transcripts were *Complete ORFs*, according to Transdecoder. The *Transdecoder ORF Types* distribution of the translated transcripts can be found in Table 3.

**Table 3.** *Transdecoder ORF Type* of the translated transcripts.

|  ORF Type  | Count |   %   |
|:----------:|:-----:|:-----:|
| 3' partial |  6734 | 13.86 |
| 5' partial | 10535 | 21.69 |
|  Internal  | 13473 | 27.74 |
|  Complete  | 17833 | 36.71 |
|    Total   | 48575 |  100  |

------------------------------------------------------------------------

Regarding the annotation of the translated transcripts [19], 29.44% had a matching *Uniprot* annotation. The number of transcripts matching an *Amphimedon queenslandica* (Aqu2) protein was slightly higher, 36.21%. *Gene Ontology Component*, *Function* and *Process* annotations could be retrieved for 26.3%, 26.49% and 26.76% of the transcripts, respectively. Additionally, 23.45% of the transcripts could be annotated with *Pfam* domain information. In contrast, only 5.8% of the assembled transcripts could be annotated against the *KEGG* database.

<img src="CBAS_DeSeq2_Analysis_files/figure-markdown_github/annotation_success_fig-1.png" style="display: block; margin: auto;" /> **Figure 2.** Annotation success by database for transcripts translated by Transdecoder.

------------------------------------------------------------------------

Regarding annotation success, different *Transdecoder ORF Types* had different annotation success rates (Fig. 3) with *complete ORFs* having higher annotation sucess than other *ORF types* and *3' Partial ORFs* having the lowest annotation success.

<img src="CBAS_DeSeq2_Analysis_files/figure-markdown_github/annotation_count_byDB_Fig-1.png" style="display: block; margin: auto;" /> **Figure 3.** Annotation count by database for different types of ORF as defined by Transdecoder.

------------------------------------------------------------------------

For annotation against Uniprot and Aqu2 we used a permissive evalue (0.001) to query these databases. Thus, we risk obtaining annotations that are only partial matches or that are of doubtful homology. Despite our permissive annotation strategy, the mean evalues obtained for the annotations against Uniprot and Aqu2 had clearly lower values than the threshold set (i.e. \(2.519351\times 10^{-7} \pm 1.127118\times 10^{-6}\) and \(2.4146826\times 10^{-7} \pm 1.1139696\times 10^{-6}\), respectively). Interestintly, the mean evalues obtained for both Uniprot and Aqu2 annotations were similar for all *Transdecoder ORF types* (Table 4) and the maximum evalue obtained after querying these two databases was \(1x10^{-05}\) [20].

**Table 4.** Mean evalue by Transdecoder ORF type for annotations obtained from the Uniprot or Aqu2 databases.

| Transdecoder ORF Type | Database |                      Mean evalue (SD)                     |
|:---------------------:|:--------:|:---------------------------------------------------------:|
|       3' partial      |  Uniprot | \(2.1738267\times 10^{-7}\) (\(1.0499955\times 10^{-6}\)) |
|       5' partial      |  Uniprot | \(2.1534162\times 10^{-7}\) (\(1.0780619\times 10^{-6}\)) |
|        Internal       |  Uniprot | \(4.9614047\times 10^{-7}\) (\(1.5576029\times 10^{-6}\)) |
|        Complete       |  Uniprot | \(1.6779664\times 10^{-7}\) (\(8.9792016\times 10^{-7}\)) |
|                       |
|       3' partial      |   Aqu2   | \(2.2338358\times 10^{-7}\) (\(1.0988336\times 10^{-6}\)) |
|       5' partial      |   Aqu2   | \(1.8405011\times 10^{-7}\) (\(9.6965115\times 10^{-7}\)) |
|        Internal       |   Aqu2   | \(4.5642628\times 10^{-7}\) (\(1.5267768\times 10^{-6}\)) |
|        Complete       |   Aqu2   | \(1.6698938\times 10^{-7}\) (\(9.0483509\times 10^{-7}\)) |

------------------------------------------------------------------------

Finally, in terms of completeness, we recovered 243 out of 248 KOGs (i.e. \(97.98\)%) from the CBAS transcriptome. Yet, the number of transcripts classified as *High-confidence, full length matches* was only 175 (i.e. \(70.56\)%). In addition, six transcripts were categorized as *Probable full length matches*, six more were tagged as matching a KOG that was slightly shorter than the query and 4 transcripts matched a KOG that was much longer than the query. Thus, in total \(191\) (i.e. \(77.02\)%) transcripts can be considered high confidence KOG matches. Transcripts matching KOGs that were much shorter and transcripts probably representing missassemblies amount to 22 and 30, respectively (i.e. \(20.97\)%).

Global gene expression patterns change due to symbiont inactivation
-------------------------------------------------------------------

The count matrix used to investigate changes in gene expression in Bleached *vs.* Control sponges had \(100098\) rows with counts. Of these, however, only \(52567\) had counts in all samples (i.e. the transcript was detected in all sequenced libraries). The inspection of the cummulative density distribution and of the density distribution of normalized counts (Fig. 5) for either all rows (not shown) or only the non-zero rows in the matrix showed similar patterns for all libraries indicating that they were sequenced at similar depths and can be compared safely.

<img src="CBAS_DeSeq2_Analysis_files/figure-markdown_github/normalized_counts_cumulativeDist_qc_figure-1.png" style="display: block; margin: auto;" /> **Figure 5.** Cummulative density distribution of normalized counts and density distribution of normalized counts of RNA-Seq libraries prepared samples for *bleached* and *control* **CBAS** samples.

------------------------------------------------------------------------

A pairwise comparison of the individual libraries (Fig 6.) provided further evidence for the technical similarity of the libraries.

<img src="CBAS_DeSeq2_Analysis_files/figure-markdown_github/MDPlots_library_01_qc_figure-1.png" style="display: block; margin: auto;" />

<img src="CBAS_DeSeq2_Analysis_files/figure-markdown_github/MDPlots_library_02_qc_figure-1.png" style="display: block; margin: auto;" />

<img src="CBAS_DeSeq2_Analysis_files/figure-markdown_github/MDPlots_library_03_qc_figure-1.png" style="display: block; margin: auto;" />

<img src="CBAS_DeSeq2_Analysis_files/figure-markdown_github/MDPlots_library_04_qc_figure-1.png" style="display: block; margin: auto;" />

<img src="CBAS_DeSeq2_Analysis_files/figure-markdown_github/MDPlots_library_05_qc_figure-1.png" style="display: block; margin: auto;" />

<img src="CBAS_DeSeq2_Analysis_files/figure-markdown_github/MDPlots_library_06_qc_figure-1.png" style="display: block; margin: auto;" />

<img src="CBAS_DeSeq2_Analysis_files/figure-markdown_github/MDPlots_library_07_qc_figure-1.png" style="display: block; margin: auto;" />

<img src="CBAS_DeSeq2_Analysis_files/figure-markdown_github/MDPlots_library_08_qc_figure-1.png" style="display: block; margin: auto;" /> **Figure 6.** Pairwise library comparisons based on the relation between the *Mean log-fold* and the *Mean (log-transformed) counts*. In general, no trend should be observed in the pairwise comparisons.

------------------------------------------------------------------------

Libraries were also compared in terms of the *FPKM* and *TPM* distribution of the transcripts (Fig. 7). This analysis revealed that two libraries (i.e. CBAS\_Control\_2 and CBAS\_Control\_3) differed markedly from all other sequenced libraries. Thus, we restricted the count matrix to include only *trinity genes* that could be translated using Transdecoder (29,705 transcripts) . After this, the libraries were comparable in terms of the distribution of TPM and FPKM values. The analysis of differential expression was conducted on this reduced data matrix.

<img src="CBAS_DeSeq2_Analysis_files/figure-markdown_github/TPM_figure-1.png" style="display: block; margin: auto;" /> **Figure 7.** TPM and FPKM distribution in datasets including all trinity genes (All transcripts) and only trinity genes that could be translated with trandescoder (Transdecoded transcripts). In general, libraries should have similar TPM/FPKM distributions to be comparable.

------------------------------------------------------------------------

For the analysis of differentially expressed genes, the count matrix was transformed using the *rlog* transformation available in the R package DeSeq2. A cluster analysis based on the between-sample distances calculated using the *rlog* transformed counts grouped samples in two groups matching the *control* and *treatment* groups. Interestingly one treatment sample behaved in an anomalous manner, clustering with neither the treatment nor the control group. This sample (i.e. CBAS\_Bleached\_5) was also found to be different in a Principal Component Analysis (PCA) conducted on the *rlog* transformed count matrix.

------------------------------------------------------------------------

<img src="CBAS_DeSeq2_Analysis_files/figure-markdown_github/sample_heat_pca_plot-1.png" style="display: block; margin: auto;" /><img src="CBAS_DeSeq2_Analysis_files/figure-markdown_github/sample_heat_pca_plot-2.png" style="display: block; margin: auto;" />

**Figure 9.** Cluster and Principal Component Analysis of the *rlog* transformed counts for *bleached vs. control* CBAS specimens.

------------------------------------------------------------------------

If sample *CBAS\_Bleached\_5* is removed from the analysis, *bleached* samples are grouped together while *control* sponges are divided in two groups in both the cluster analysis and the principal component analysis.

<img src="CBAS_DeSeq2_Analysis_files/figure-markdown_github/drop_anomalous_sample_and_replot-1.png" style="display: block; margin: auto;" /><img src="CBAS_DeSeq2_Analysis_files/figure-markdown_github/drop_anomalous_sample_and_replot-2.png" style="display: block; margin: auto;" />

**Figure 10.** Cluster and Principal Component Analysis of the *rlog* transformed counts for *bleached vs. control* CBAS specimens after removing sample *CBAS\_Bleached\_5*.

------------------------------------------------------------------------

In this reduced dataset, a total of 5508 transcripts were differentially expressed between the *control* and *bleached* sponges (Fig. 11). Of these transcripts, 407 were overexpressed (with a log2 fold change higher or equal to 2) in *bleached vs. control* sponges and 1597 were expressed at lower levels (a log2 fold change lower or equal than -2) in *bleached vs. control* specimens.

<img src="CBAS_DeSeq2_Analysis_files/figure-markdown_github/MA_plot-1.png" style="display: block; margin: auto;" />

**Figure 11.** MA-plot showing the relation between log-fold change and mean count for each analyzed CBAS transcript. Transcripts significant at 0.01 are highlighted in red.

------------------------------------------------------------------------

It is worth noting that DeSeq2 was able to adecuately model the variance in the dataset and could correctly calculate the p-values for the transcripts analyzed (Fig 12).

<img src="CBAS_DeSeq2_Analysis_files/figure-markdown_github/probability_plot-1.png" style="display: block; margin: auto;" />

**Figure 12.** Histogram of the p-values obtained for the transcripts analyzed using DeSeq2. A uniform distribution between 0 and 1, with a peak at zero, should be observed if the variance of the dataset could be correctly modelled.

------------------------------------------------------------------------

<img src="CBAS_DeSeq2_Analysis_files/figure-markdown_github/expression_heat_maps-1.png" style="display: block; margin: auto;" /> **Figure 13.** Heatplot showing transcripts with a high log-fold change (LFC \(\geq\) 6) between treatments.

------------------------------------------------------------------------

Inmune response transcripts are highly repressed when symbionts are not active or absent
----------------------------------------------------------------------------------------

<img src="CBAS_DeSeq2_Analysis_files/figure-markdown_github/plot_back_vs_forground_dist-1.png" style="display: block; margin: auto;" />

HERE
====

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

[5] EMBL GeneCore. <http://genecore3.genecore.embl.de/genecore3/index.cfm>

[6] FastQC. <http://www.bioinformatics.babraham.ac.uk/projects/fastqc/>

[7] BioLite: <https://bitbucket.org/caseywdunn/biolite/overview>

[8] Downloaded early 2015. [Uniprot download.](http://www.uniprot.org/downloads)

[9] Aqu2 isoforms (proteins) [*Amphimedon queenslandica* transcriptome resource.](http://amphimedon.qcloud.qcif.edu.au/downloads.html)

[10] Blast was run from [Galaxy](https://galaxyproject.org/) and the XML file was converted to a 25 column table by the [ncbi\_blast\_plus](http://toolshed.g2.bx.psu.edu/view/devteam/ncbi_blast_plus) tool.

[11] Gene Ontology Consortium. <http://geneontology.org/>

[12] The UNIPROT accession code of the blast best match associated to each transcript in the CBAS transcriptome was used to query the GO annotations associated with the UNIPROT protein. These GO annotations were used to annotated the CBAS transcripts; note that these are annotations based on annotations and should be taken with precaution.

[13] TopGO [(https://bioconductor.org/packages/release/bioc/html/topGO.html)](https://bioconductor.org/packages/release/bioc/html/topGO.html) expects GO terms in tab separated {Trascript, GoAnnotation} pairs. The script produces output with a leading column indicating the UNIPROT accession code used for each transcript, which comes handy when controling the correct annotation was used for each transcript and can be easily removed for further analyses.

[14] Transdecoder's Github Repository can be found [here.](https://transdecoder.github.io/)

[15] I honestly do not know what the purpose of this blast run was, I did it for sake of completeness and because the reads were not filtered before assembling the conting. Yet, after seeing the results, it seems as if not too many contigs matched the bacterial database...

[16] Parra et al. 2007. CEGMA: a pipeline to accurately annotate core genes in eukaryotic genomes. [Bioinformatics 23 (9): 1061-1067.](http://dx.doi.org/10.1093/bioinformatics/btm071)

[17] Francis et al. 2013. A comparison across non-model animals suggests an optimal sequencing depth for *de novo* transcriptome assembly. [BMC Genomics 14:167](http://dx.doi.org/10.1186/1471-2164-14-167)

[18] the script is also provided as part of this repository. The Copyright (c) of this script belongs trinityrnaseq (2014), who reserves all rights upon the code. [See the full LICENSE here.](https://github.com/trinityrnaseq/trinityrnaseq/blob/master/LICENSE)

[19] With these graphs were produced for transcripts that could be successfully translated using **Transdecoder**. In the case of the **Uniprot** and **Aqu2** annotations, the transcripts were directly blasted using blastn. Thus could be possible to obtain **Uniprot** and **Aqu2** annotations for transcripts that could not be translated with **Transdecoder**. This is in fact the case, however the number of transcripts that could not *transdecoded* but were annotated is relatively small. In the case of **Uniprot**, the number of transcripts that yielded annotations but were not *transdecoded* was **7535**. For **Aqu2**, only **11915** not *transdecoded* transcripts could be annotated.

[20] Here it is worth noting that the evalue cutoff used for blast against the Uniprot and Aqu2 databases was 0.001. So the annotations obtained had an evalue at least two orders of magnitud smaller than the set cutoff evalue.
