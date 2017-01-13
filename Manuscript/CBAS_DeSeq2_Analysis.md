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

It is worth noting that DeSeq2 was able to adecuately model the variance in the dataset and could correctly calculate the p-values for the analyzed transcripts (Fig 12).

<img src="CBAS_DeSeq2_Analysis_files/figure-markdown_github/probability_plot-1.png" style="display: block; margin: auto;" />

**Figure 12.** Histogram of the p-values obtained for the transcripts analyzed using DeSeq2. A uniform distribution between 0 and 1, with a peak at zero, should be observed if the variance of the dataset could be correctly modelled.

------------------------------------------------------------------------

Inmune response transcripts are highly repressed when symbionts are not active or absent
----------------------------------------------------------------------------------------

For Gene Ontology (GO) enrichment analyses we were able to select a background set of transcripts showing similar expression patterns but no significant change between conditions for all transcripts as well as for overexpressed and underexpressed transcripts (Fig. 13). Enriched Function, Process and Compartment GO terms for over- and underexpressed genes can be found in Tables 5-10.

<img src="CBAS_DeSeq2_Analysis_files/figure-markdown_github/plot_back_vs_forground_dist-1.png" style="display: block; margin: auto;" />

Table 5

| GO.ID        | Term                                        |  Annotated|  Significant|  Expected| classic |
|:-------------|:--------------------------------------------|----------:|------------:|---------:|:--------|
| <GO:0008234> | cysteine-type peptidase activity            |         36|           24|      4.03| 1.3e-15 |
| <GO:0003924> | GTPase activity                             |         36|           17|      4.03| 4.1e-08 |
| <GO:0005200> | structural constituent of cytoskeleton      |         27|           14|      3.02| 1.6e-07 |
| <GO:0004197> | cysteine-type endopeptidase activity        |         17|           11|      1.90| 1.7e-07 |
| <GO:0070011> | peptidase activity, acting on L-amino ac... |        144|           36|     16.10| 6.4e-07 |
| <GO:0008233> | peptidase activity                          |        149|           36|     16.66| 1.6e-06 |
| <GO:0016787> | hydrolase activity                          |        383|           67|     42.83| 1.1e-05 |
| <GO:0005525> | GTP binding                                 |         71|           21|      7.94| 1.2e-05 |
| <GO:0032561> | guanyl ribonucleotide binding               |         72|           21|      8.05| 1.5e-05 |
| <GO:0019001> | guanyl nucleotide binding                   |         74|           21|      8.28| 2.4e-05 |
| <GO:0038024> | cargo receptor activity                     |         37|           13|      4.14| 8.8e-05 |
| <GO:0005044> | scavenger receptor activity                 |         33|           12|      3.69| 0.00011 |
| <GO:0005507> | copper ion binding                          |         14|            7|      1.57| 0.00034 |
| <GO:0004497> | monooxygenase activity                      |         17|            7|      1.90| 0.00142 |
| <GO:0019783> | ubiquitin-like protein-specific protease... |         13|            6|      1.45| 0.00157 |
| <GO:0017111> | nucleoside-triphosphatase activity          |        118|           24|     13.20| 0.00178 |
| <GO:0035375> | zymogen binding                             |         14|            6|      1.57| 0.00250 |
| <GO:0043130> | ubiquitin binding                           |         14|            6|      1.57| 0.00250 |
| <GO:0016462> | pyrophosphatase activity                    |        121|           24|     13.53| 0.00256 |
| <GO:0016817> | hydrolase activity, acting on acid anhyd... |        121|           24|     13.53| 0.00256 |
| <GO:0016818> | hydrolase activity, acting on acid anhyd... |        121|           24|     13.53| 0.00256 |
| <GO:0004843> | ubiquitin-specific protease activity        |         10|            5|      1.12| 0.00261 |
| <GO:0036459> | ubiquitinyl hydrolase activity              |         10|            5|      1.12| 0.00261 |
| <GO:0004175> | endopeptidase activity                      |        122|           24|     13.64| 0.00288 |
| <GO:0032182> | ubiquitin-like protein binding              |         15|            6|      1.68| 0.00379 |
| <GO:0016715> | oxidoreductase activity, acting on paire... |          7|            4|      0.78| 0.00403 |
| <GO:0042393> | histone binding                             |         11|            5|      1.23| 0.00435 |
| <GO:0004435> | phosphatidylinositol phospholipase C act... |          4|            3|      0.45| 0.00505 |
| <GO:0004629> | phospholipase C activity                    |          4|            3|      0.45| 0.00505 |
| <GO:0016682> | oxidoreductase activity, acting on diphe... |          4|            3|      0.45| 0.00505 |
| <GO:0052716> | hydroquinone:oxygen oxidoreductase activ... |          4|            3|      0.45| 0.00505 |
| <GO:0003824> | catalytic activity                          |        814|          107|     91.03| 0.00603 |
| <GO:0060089> | molecular transducer activity               |        137|           25|     15.32| 0.00676 |
| <GO:0005198> | structural molecule activity                |         84|           17|      9.39| 0.00900 |
| <GO:0016705> | oxidoreductase activity, acting on paire... |         23|            7|      2.57| 0.01001 |
| <GO:0016679> | oxidoreductase activity, acting on diphe... |          5|            3|      0.56| 0.01160 |
| <GO:0036094> | small molecule binding                      |        394|           57|     44.06| 0.01191 |
| <GO:0003746> | translation elongation factor activity      |          9|            4|      1.01| 0.01212 |
| <GO:0003979> | UDP-glucose 6-dehydrogenase activity        |          2|            2|      0.22| 0.01244 |
| <GO:0004614> | phosphoglucomutase activity                 |          2|            2|      0.22| 0.01244 |
| <GO:0004796> | thromboxane-A synthase activity             |          2|            2|      0.22| 0.01244 |
| <GO:0016868> | intramolecular transferase activity, pho... |          2|            2|      0.22| 0.01244 |
| <GO:0050811> | GABA receptor binding                       |          2|            2|      0.22| 0.01244 |
| <GO:0016491> | oxidoreductase activity                     |        115|           21|     12.86| 0.01287 |
| <GO:0004872> | receptor activity                           |        126|           22|     14.09| 0.01833 |
| <GO:0001882> | nucleoside binding                          |        329|           48|     36.79| 0.01934 |
| <GO:0001883> | purine nucleoside binding                   |        329|           48|     36.79| 0.01934 |
| <GO:0016853> | isomerase activity                          |         15|            5|      1.68| 0.01947 |
| <GO:0008329> | signaling pattern recognition receptor a... |          6|            3|      0.67| 0.02130 |
| <GO:0038187> | pattern recognition receptor activity       |          6|            3|      0.67| 0.02130 |

Table 6

| GO.ID        | Term                                        |  Annotated|  Significant|  Expected| classic |
|:-------------|:--------------------------------------------|----------:|------------:|---------:|:--------|
| <GO:0006898> | receptor-mediated endocytosis               |         47|           14|      5.14| 0.00028 |
| <GO:0015979> | photosynthesis                              |         19|            8|      2.08| 0.00046 |
| <GO:0019684> | photosynthesis, light reaction              |         15|            7|      1.64| 0.00050 |
| <GO:0019229> | regulation of vasoconstriction              |          3|            3|      0.33| 0.00129 |
| <GO:0042310> | vasoconstriction                            |          3|            3|      0.33| 0.00129 |
| <GO:0045907> | positive regulation of vasoconstriction     |          3|            3|      0.33| 0.00129 |
| <GO:0051289> | protein homotetramerization                 |          9|            5|      0.99| 0.00130 |
| <GO:0006636> | unsaturated fatty acid biosynthetic proc... |          6|            4|      0.66| 0.00174 |
| <GO:0006508> | proteolysis                                 |        240|           40|     26.27| 0.00221 |
| <GO:0010623> | developmental programmed cell death         |         10|            5|      1.09| 0.00236 |
| <GO:0006897> | endocytosis                                 |         97|           20|     10.62| 0.00287 |
| <GO:0021889> | olfactory bulb interneuron differentiati... |          7|            4|      0.77| 0.00372 |
| <GO:0021891> | olfactory bulb interneuron development      |          7|            4|      0.77| 0.00372 |
| <GO:0050920> | regulation of chemotaxis                    |         20|            7|      2.19| 0.00375 |
| <GO:0009698> | phenylpropanoid metabolic process           |          4|            3|      0.44| 0.00475 |
| <GO:0009808> | lignin metabolic process                    |          4|            3|      0.44| 0.00475 |
| <GO:0009814> | defense response, incompatible interacti... |          4|            3|      0.44| 0.00475 |
| <GO:0009816> | defense response to bacterium, incompati... |          4|            3|      0.44| 0.00475 |
| <GO:0046271> | phenylpropanoid catabolic process           |          4|            3|      0.44| 0.00475 |
| <GO:0046274> | lignin catabolic process                    |          4|            3|      0.44| 0.00475 |
| <GO:0030163> | protein catabolic process                   |         81|           17|      8.87| 0.00492 |
| <GO:0050772> | positive regulation of axonogenesis         |         12|            5|      1.31| 0.00619 |
| <GO:0044257> | cellular protein catabolic process          |         70|           15|      7.66| 0.00668 |
| <GO:0016578> | histone deubiquitination                    |          8|            4|      0.88| 0.00682 |
| <GO:0021872> | forebrain generation of neurons             |          8|            4|      0.88| 0.00682 |
| <GO:0021879> | forebrain neuron differentiation            |          8|            4|      0.88| 0.00682 |
| <GO:0021884> | forebrain neuron development                |          8|            4|      0.88| 0.00682 |
| <GO:0070646> | protein modification by small protein re... |         17|            6|      1.86| 0.00696 |
| <GO:0001554> | luteolysis                                  |          5|            3|      0.55| 0.01091 |
| <GO:0009767> | photosynthetic electron transport chain     |          5|            3|      0.55| 0.01091 |
| <GO:0051932> | synaptic transmission, GABAergic            |          5|            3|      0.55| 0.01091 |
| <GO:0061364> | apoptotic process involved in luteolysis    |          5|            3|      0.55| 0.01091 |
| <GO:1903524> | positive regulation of blood circulation    |          5|            3|      0.55| 0.01091 |
| <GO:0021772> | olfactory bulb development                  |          9|            4|      0.99| 0.01123 |
| <GO:0021988> | olfactory lobe development                  |          9|            4|      0.99| 0.01123 |
| <GO:0050770> | regulation of axonogenesis                  |         24|            7|      2.63| 0.01144 |
| <GO:0001516> | prostaglandin biosynthetic process          |          2|            2|      0.22| 0.01192 |
| <GO:0006065> | UDP-glucuronate biosynthetic process        |          2|            2|      0.22| 0.01192 |
| <GO:0009737> | response to abscisic acid                   |          2|            2|      0.22| 0.01192 |
| <GO:0009739> | response to gibberellin                     |          2|            2|      0.22| 0.01192 |
| <GO:0010150> | leaf senescence                             |          2|            2|      0.22| 0.01192 |
| <GO:0010207> | photosystem II assembly                     |          2|            2|      0.22| 0.01192 |
| <GO:0010260> | organ senescence                            |          2|            2|      0.22| 0.01192 |
| <GO:0014909> | smooth muscle cell migration                |          2|            2|      0.22| 0.01192 |
| <GO:0014910> | regulation of smooth muscle cell migrati... |          2|            2|      0.22| 0.01192 |
| <GO:0030002> | cellular anion homeostasis                  |          2|            2|      0.22| 0.01192 |
| <GO:0030320> | cellular monovalent inorganic anion home... |          2|            2|      0.22| 0.01192 |
| <GO:0030644> | cellular chloride ion homeostasis           |          2|            2|      0.22| 0.01192 |
| <GO:0031644> | regulation of neurological system proces... |          2|            2|      0.22| 0.01192 |
| <GO:0031646> | positive regulation of neurological syst... |          2|            2|      0.22| 0.01192 |

Table 7

| GO.ID        | Term                               |  Annotated|  Significant|  Expected| classic |
|:-------------|:-----------------------------------|----------:|------------:|---------:|:--------|
| <GO:0005618> | cell wall                          |         17|            9|      1.97| 3.3e-05 |
| <GO:0030312> | external encapsulating structure   |         19|            9|      2.20| 0.00010 |
| <GO:0005576> | extracellular region               |        344|           60|     39.92| 0.00016 |
| <GO:0005615> | extracellular space                |         70|           17|      8.12| 0.00177 |
| <GO:0042588> | zymogen granule                    |         14|            6|      1.62| 0.00304 |
| <GO:0042589> | zymogen granule membrane           |         14|            6|      1.62| 0.00304 |
| <GO:0005874> | microtubule                        |         45|           12|      5.22| 0.00373 |
| <GO:0005773> | vacuole                            |         69|           16|      8.01| 0.00401 |
| <GO:0045335> | phagocytic vesicle                 |         16|            6|      1.86| 0.00663 |
| <GO:0030667> | secretory granule membrane         |         21|            7|      2.44| 0.00708 |
| <GO:0043083> | synaptic cleft                     |          8|            4|      0.93| 0.00843 |
| <GO:0044434> | chloroplast part                   |         23|            7|      2.67| 0.01221 |
| <GO:0044435> | plastid part                       |         23|            7|      2.67| 0.01221 |
| <GO:0000325> | plant-type vacuole                 |          2|            2|      0.23| 0.01340 |
| <GO:0010282> | senescence-associated vacuole      |          2|            2|      0.23| 0.01340 |
| <GO:0005578> | proteinaceous extracellular matrix |         65|           14|      7.54| 0.01377 |
| <GO:0031012> | extracellular matrix               |         66|           14|      7.66| 0.01573 |
| <GO:0005581> | collagen trimer                    |         14|            5|      1.62| 0.01659 |
| <GO:0009507> | chloroplast                        |         30|            8|      3.48| 0.01728 |
| <GO:0009536> | plastid                            |         30|            8|      3.48| 0.01728 |
| <GO:0030424> | axon                               |         61|           13|      7.08| 0.01898 |
| <GO:0030673> | axolemma                           |         10|            4|      1.16| 0.02096 |
| <GO:0032589> | neuron projection membrane         |         10|            4|      1.16| 0.02096 |
| <GO:0019898> | extrinsic component of membrane    |         37|            9|      4.29| 0.02152 |
| <GO:0009535> | chloroplast thylakoid membrane     |         15|            5|      1.74| 0.02259 |
| <GO:0055035> | plastid thylakoid membrane         |         15|            5|      1.74| 0.02259 |
| <GO:0009534> | chloroplast thylakoid              |         16|            5|      1.86| 0.02985 |
| <GO:0031256> | leading edge membrane              |         16|            5|      1.86| 0.02985 |
| <GO:0031976> | plastid thylakoid                  |         16|            5|      1.86| 0.02985 |
| <GO:0005764> | lysosome                           |         52|           11|      6.03| 0.03160 |
| <GO:0033267> | axon part                          |         34|            8|      3.95| 0.03574 |
| <GO:0034357> | photosynthetic membrane            |         17|            5|      1.97| 0.03842 |
| <GO:0042651> | thylakoid membrane                 |         17|            5|      1.97| 0.03842 |
| <GO:0044436> | thylakoid part                     |         17|            5|      1.97| 0.03842 |
| <GO:0000323> | lytic vacuole                      |         54|           11|      6.27| 0.04072 |
| <GO:0030670> | phagocytic vesicle membrane        |         12|            4|      1.39| 0.04104 |
| <GO:0030141> | secretory granule                  |         35|            8|      4.06| 0.04188 |
| <GO:0031410> | cytoplasmic vesicle                |        111|           19|     12.88| 0.04694 |
| <GO:0009579> | thylakoid                          |         18|            5|      2.09| 0.04835 |
| <GO:0030139> | endocytic vesicle                  |         25|            6|      2.90| 0.06045 |
| <GO:0000786> | nucleosome                         |          4|            2|      0.46| 0.06860 |
| <GO:0044815> | DNA packaging complex              |          4|            2|      0.46| 0.06860 |
| <GO:0005788> | endoplasmic reticulum lumen        |         14|            4|      1.62| 0.06905 |
| <GO:0043195> | terminal bouton                    |         14|            4|      1.62| 0.06905 |
| <GO:0009279> | cell outer membrane                |          9|            3|      1.04| 0.07614 |
| <GO:0009523> | photosystem II                     |          9|            3|      1.04| 0.07614 |
| <GO:0044304> | main axon                          |         15|            4|      1.74| 0.08594 |
| <GO:0005737> | cytoplasm                          |        987|          123|    114.53| 0.09142 |
| <GO:0000790> | nuclear chromatin                  |         10|            3|      1.16| 0.09979 |
| <GO:0009521> | photosystem                        |         10|            3|      1.16| 0.09979 |

Table 8

| GO.ID        | Term                                        |  Annotated|  Significant|  Expected| classic |
|:-------------|:--------------------------------------------|----------:|------------:|---------:|:--------|
| <GO:0003774> | motor activity                              |        113|           62|     15.28| 8.5e-26 |
| <GO:0003777> | microtubule motor activity                  |         80|           51|     10.81| 1.3e-25 |
| <GO:0005201> | extracellular matrix structural constitu... |         46|           22|      6.22| 1.8e-08 |
| <GO:0004550> | nucleoside diphosphate kinase activity      |         11|           10|      1.49| 1.9e-08 |
| <GO:0004867> | serine-type endopeptidase inhibitor acti... |         24|           15|      3.24| 3.2e-08 |
| <GO:0004866> | endopeptidase inhibitor activity            |         36|           18|      4.87| 1.6e-07 |
| <GO:0030414> | peptidase inhibitor activity                |         40|           19|      5.41| 2.0e-07 |
| <GO:0061135> | endopeptidase regulator activity            |         38|           18|      5.14| 4.5e-07 |
| <GO:0046906> | tetrapyrrole binding                        |         91|           29|     12.30| 4.5e-06 |
| <GO:0015085> | calcium ion transmembrane transporter ac... |         36|           15|      4.87| 2.8e-05 |
| <GO:0016887> | ATPase activity                             |        302|           66|     40.83| 3.1e-05 |
| <GO:0061134> | peptidase regulator activity                |         53|           19|      7.16| 3.2e-05 |
| <GO:0005432> | calcium:sodium antiporter activity          |         12|            8|      1.62| 3.2e-05 |
| <GO:0015368> | calcium:cation antiporter activity          |         12|            8|      1.62| 3.2e-05 |
| <GO:0017111> | nucleoside-triphosphatase activity          |        471|           94|     63.67| 3.3e-05 |
| <GO:0020037> | heme binding                                |         72|           23|      9.73| 4.2e-05 |
| <GO:0008191> | metalloendopeptidase inhibitor activity     |          5|            5|      0.68| 4.5e-05 |
| <GO:0010576> | metalloenzyme regulator activity            |          5|            5|      0.68| 4.5e-05 |
| <GO:0048551> | metalloenzyme inhibitor activity            |          5|            5|      0.68| 4.5e-05 |
| <GO:0015081> | sodium ion transmembrane transporter act... |         23|           11|      3.11| 7.4e-05 |
| <GO:0016776> | phosphotransferase activity, phosphate g... |         20|           10|      2.70| 9.9e-05 |
| <GO:0016462> | pyrophosphatase activity                    |        486|           94|     65.70| 0.00011 |
| <GO:0016817> | hydrolase activity, acting on acid anhyd... |        487|           94|     65.83| 0.00012 |
| <GO:0016818> | hydrolase activity, acting on acid anhyd... |        487|           94|     65.83| 0.00012 |
| <GO:0015491> | cation:cation antiporter activity           |         14|            8|      1.89| 0.00015 |
| <GO:0004497> | monooxygenase activity                      |         55|           18|      7.44| 0.00020 |
| <GO:0015298> | solute:cation antiporter activity           |         18|            9|      2.43| 0.00022 |
| <GO:0016682> | oxidoreductase activity, acting on diphe... |          6|            5|      0.81| 0.00024 |
| <GO:0052716> | hydroquinone:oxygen oxidoreductase activ... |          6|            5|      0.81| 0.00024 |
| <GO:0016712> | oxidoreductase activity, acting on paire... |         12|            7|      1.62| 0.00034 |
| <GO:0070330> | aromatase activity                          |         12|            7|      1.62| 0.00034 |
| <GO:0004796> | thromboxane-A synthase activity             |          7|            5|      0.95| 0.00074 |
| <GO:0016679> | oxidoreductase activity, acting on diphe... |          7|            5|      0.95| 0.00074 |
| <GO:0005518> | collagen binding                            |         25|           10|      3.38| 0.00092 |
| <GO:0019205> | nucleobase-containing compound kinase ac... |         25|           10|      3.38| 0.00092 |
| <GO:0004857> | enzyme inhibitor activity                   |         87|           23|     11.76| 0.00093 |
| <GO:0072509> | divalent inorganic cation transmembrane ... |         53|           16|      7.16| 0.00121 |
| <GO:0005516> | calmodulin binding                          |         94|           24|     12.71| 0.00123 |
| <GO:0004161> | dimethylallyltranstransferase activity      |          5|            4|      0.68| 0.00148 |
| <GO:0004337> | geranyltranstransferase activity            |          5|            4|      0.68| 0.00148 |
| <GO:0004735> | pyrroline-5-carboxylate reductase activi... |          5|            4|      0.68| 0.00148 |
| <GO:0030332> | cyclin binding                              |         11|            6|      1.49| 0.00151 |
| <GO:0004017> | adenylate kinase activity                   |          8|            5|      1.08| 0.00175 |
| <GO:0005220> | inositol 1,4,5-trisphosphate-sensitive c... |          8|            5|      1.08| 0.00175 |
| <GO:0043394> | proteoglycan binding                        |         15|            7|      2.03| 0.00193 |
| <GO:0030246> | carbohydrate binding                        |        114|           27|     15.41| 0.00211 |
| <GO:0004872> | receptor activity                           |        401|           74|     54.21| 0.00245 |
| <GO:0015643> | toxic substance binding                     |          3|            3|      0.41| 0.00246 |
| <GO:0018585> | fluorene oxygenase activity                 |          3|            3|      0.41| 0.00246 |
| <GO:0048407> | platelet-derived growth factor binding      |          3|            3|      0.41| 0.00246 |

Table 9

| GO.ID        | Term                                        |  Annotated|  Significant|  Expected| classic    |
|:-------------|:--------------------------------------------|----------:|------------:|---------:|:-----------|
| <GO:0007018> | microtubule-based movement                  |        158|           93|     21.54| &lt; 1e-30 |
| <GO:0007017> | microtubule-based process                   |        332|          133|     45.26| &lt; 1e-30 |
| <GO:0060271> | cilium morphogenesis                        |        180|           92|     24.54| &lt; 1e-30 |
| <GO:0044782> | cilium organization                         |        174|           89|     23.72| &lt; 1e-30 |
| <GO:0042384> | cilium assembly                             |        164|           84|     22.36| &lt; 1e-30 |
| <GO:0035082> | axoneme assembly                            |         59|           45|      8.04| 7.8e-28    |
| <GO:0001578> | microtubule bundle formation                |         66|           46|      9.00| 1.4e-25    |
| <GO:0003341> | cilium movement                             |         57|           42|      7.77| 5.1e-25    |
| <GO:0010927> | cellular component assembly involved in ... |        210|           88|     28.63| 6.7e-25    |
| <GO:0030031> | cell projection assembly                    |        211|           86|     28.76| 2.4e-23    |
| <GO:0001539> | cilium or flagellum-dependent cell motil... |         32|           29|      4.36| 1.7e-22    |
| <GO:0006928> | movement of cell or subcellular componen... |        659|          172|     89.83| 7.8e-20    |
| <GO:0048646> | anatomical structure formation involved ... |        482|          135|     65.70| 2.7e-18    |
| <GO:0048858> | cell projection morphogenesis               |        449|          126|     61.20| 3.8e-17    |
| <GO:0032990> | cell part morphogenesis                     |        452|          126|     61.61| 6.8e-17    |
| <GO:0070925> | organelle assembly                          |        326|          100|     44.44| 1.6e-16    |
| <GO:0000226> | microtubule cytoskeleton organization       |        210|           73|     28.63| 2.0e-15    |
| <GO:0030030> | cell projection organization                |        630|          155|     85.88| 2.9e-15    |
| <GO:0000902> | cell morphogenesis                          |        560|          140|     76.33| 2.2e-14    |
| <GO:0070286> | axonemal dynein complex assembly            |         30|           22|      4.09| 1.4e-13    |
| <GO:0009653> | anatomical structure morphogenesis          |       1009|          213|    137.54| 3.8e-13    |
| <GO:0032989> | cellular component morphogenesis            |        597|          143|     81.38| 4.5e-13    |
| <GO:0060285> | cilium-dependent cell motility              |         20|           17|      2.73| 1.3e-12    |
| <GO:0060294> | cilium movement involved in cell motilit... |         14|           12|      1.91| 2.7e-09    |
| <GO:0036158> | outer dynein arm assembly                   |         17|           13|      2.32| 7.2e-09    |
| <GO:0003351> | epithelial cilium movement                  |         13|           11|      1.77| 1.7e-08    |
| <GO:0048869> | cellular developmental process              |       1327|          243|    180.89| 3.0e-08    |
| <GO:0042073> | intraciliary transport                      |         21|           14|      2.86| 3.2e-08    |
| <GO:0048856> | anatomical structure development            |       1740|          304|    237.18| 3.8e-08    |
| <GO:0036159> | inner dynein arm assembly                   |         12|           10|      1.64| 1.1e-07    |
| <GO:0007286> | spermatid development                       |         39|           19|      5.32| 1.4e-07    |
| <GO:0048515> | spermatid differentiation                   |         40|           19|      5.45| 2.3e-07    |
| <GO:0044767> | single-organism developmental process       |       1963|          330|    267.58| 5.3e-07    |
| <GO:0032502> | developmental process                       |       1979|          332|    269.76| 6.0e-07    |
| <GO:0007130> | synaptonemal complex assembly               |          7|            7|      0.95| 8.6e-07    |
| <GO:0007131> | reciprocal meiotic recombination            |         23|           13|      3.14| 1.6e-06    |
| <GO:0035825> | reciprocal DNA recombination                |         23|           13|      3.14| 1.6e-06    |
| <GO:0045143> | homologous chromosome segregation           |         23|           13|      3.14| 1.6e-06    |
| <GO:0070192> | chromosome organization involved in meio... |         23|           13|      3.14| 1.6e-06    |
| <GO:0007129> | synapsis                                    |         17|           11|      2.32| 1.6e-06    |
| <GO:0021591> | ventricular system development              |         30|           15|      4.09| 1.9e-06    |
| <GO:0007010> | cytoskeleton organization                   |        409|           89|     55.75| 2.1e-06    |
| <GO:0030198> | extracellular matrix organization           |        112|           34|     15.27| 2.9e-06    |
| <GO:0043062> | extracellular structure organization        |        118|           35|     16.08| 3.7e-06    |
| <GO:0007283> | spermatogenesis                             |        148|           41|     20.17| 3.8e-06    |
| <GO:0007127> | meiosis I                                   |         39|           17|      5.32| 4.4e-06    |
| <GO:0022607> | cellular component assembly                 |        901|          167|    122.82| 4.5e-06    |
| <GO:0048232> | male gamete generation                      |        150|           41|     20.45| 5.5e-06    |
| <GO:0070193> | synaptonemal complex organization           |          8|            7|      1.09| 6.0e-06    |
| <GO:0007140> | male meiosis                                |         22|           12|      3.00| 6.6e-06    |

Table 10

| GO.ID        | Term                                        |  Annotated|  Significant|  Expected| classic    |
|:-------------|:--------------------------------------------|----------:|------------:|---------:|:-----------|
| <GO:0005929> | cilium                                      |        296|          156|     41.66| &lt; 1e-30 |
| <GO:0005930> | axoneme                                     |         95|           73|     13.37| &lt; 1e-30 |
| <GO:0097014> | ciliary cytoplasm                           |         95|           73|     13.37| &lt; 1e-30 |
| <GO:0044441> | ciliary part                                |        179|          101|     25.19| &lt; 1e-30 |
| <GO:0032838> | cell projection cytoplasm                   |        112|           74|     15.76| &lt; 1e-30 |
| <GO:0044447> | axoneme part                                |         48|           44|      6.75| &lt; 1e-30 |
| <GO:0030286> | dynein complex                              |         68|           50|      9.57| 7.6e-29    |
| <GO:0042995> | cell projection                             |        801|          223|    112.72| 1.6e-28    |
| <GO:0031514> | motile cilium                               |         60|           45|      8.44| 1.1e-26    |
| <GO:0044463> | cell projection part                        |        415|          132|     58.40| 5.3e-22    |
| <GO:0044430> | cytoskeletal part                           |        598|          166|     84.15| 1.3e-20    |
| <GO:0005858> | axonemal dynein complex                     |         32|           28|      4.50| 2.0e-20    |
| <GO:0005856> | cytoskeleton                                |        855|          213|    120.32| 3.2e-20    |
| <GO:0015630> | microtubule cytoskeleton                    |        507|          146|     71.35| 1.1e-19    |
| <GO:0005874> | microtubule                                 |        211|           79|     29.69| 5.8e-18    |
| <GO:0005868> | cytoplasmic dynein complex                  |         38|           29|      5.35| 6.1e-18    |
| <GO:0005875> | microtubule associated complex              |        118|           55|     16.61| 8.5e-18    |
| <GO:0005578> | proteinaceous extracellular matrix          |        259|           74|     36.45| 4.1e-10    |
| <GO:0031012> | extracellular matrix                        |        264|           74|     37.15| 1.1e-09    |
| <GO:0036126> | sperm flagellum                             |         21|           15|      2.96| 3.6e-09    |
| <GO:0005814> | centriole                                   |         60|           27|      8.44| 5.9e-09    |
| <GO:0097223> | sperm part                                  |         40|           21|      5.63| 9.6e-09    |
| <GO:1990716> | axonemal central apparatus                  |          9|            9|      1.27| 2.1e-08    |
| <GO:0044450> | microtubule organizing center part          |         76|           30|     10.70| 3.5e-08    |
| <GO:0005581> | collagen trimer                             |         52|           23|      7.32| 1.2e-07    |
| <GO:0036156> | inner dynein arm                            |         10|            9|      1.41| 1.8e-07    |
| <GO:0005813> | centrosome                                  |        212|           57|     29.83| 4.5e-07    |
| <GO:0030990> | intraciliary transport particle             |         17|           11|      2.39| 2.2e-06    |
| <GO:0005815> | microtubule organizing center               |        276|           67|     38.84| 2.5e-06    |
| <GO:0072372> | primary cilium                              |         96|           31|     13.51| 3.7e-06    |
| <GO:0036064> | ciliary basal body                          |         54|           21|      7.60| 5.3e-06    |
| <GO:0005583> | fibrillar collagen trimer                   |          6|            6|      0.84| 7.7e-06    |
| <GO:0098643> | banded collagen fibril                      |          6|            6|      0.84| 7.7e-06    |
| <GO:0043228> | non-membrane-bounded organelle              |       1512|          264|    212.78| 1.1e-05    |
| <GO:0043232> | intracellular non-membrane-bounded organ... |       1512|          264|    212.78| 1.1e-05    |
| <GO:0098644> | complex of collagen trimers                 |          9|            7|      1.27| 3.0e-05    |
| <GO:0097228> | sperm principal piece                       |          5|            5|      0.70| 5.5e-05    |
| <GO:1990718> | axonemal central pair projection            |          5|            5|      0.70| 5.5e-05    |
| <GO:0005615> | extracellular space                         |        286|           64|     40.25| 6.7e-05    |
| <GO:0030992> | intraciliary transport particle B           |         10|            7|      1.41| 8.7e-05    |
| <GO:0036157> | outer dynein arm                            |         10|            7|      1.41| 8.7e-05    |
| <GO:0030141> | secretory granule                           |        114|           31|     16.04| 0.00016    |
| <GO:0000794> | condensed nuclear chromosome                |         36|           14|      5.07| 0.00020    |
| <GO:0000795> | synaptonemal complex                        |         14|            8|      1.97| 0.00020    |
| <GO:0030089> | phycobilisome                               |          4|            4|      0.56| 0.00039    |
| <GO:0005576> | extracellular region                        |       1190|          204|    167.47| 0.00051    |
| <GO:0042383> | sarcolemma                                  |         53|           17|      7.46| 0.00063    |
| <GO:0034357> | photosynthetic membrane                     |         40|           14|      5.63| 0.00071    |
| <GO:0044436> | thylakoid part                              |         41|           14|      5.77| 0.00094    |
| <GO:0031512> | motile primary cilium                       |         10|            6|      1.41| 0.00097    |

PFAM
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
