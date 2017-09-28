------------------------------------------------------------------------

Preliminary remarks
-------------------

The common blue aquarium sponge is a cyanosponge commonly found in salt
water aquaria. Most aquarium hobbiest refer to it as *Collospongia
auris*, however its taxonomic affiliation is not clear. I refer to this
sponge as the "Common Blue Aquarium Sponge", **CBAS** in short. CBAS is
a cyanosponge, it harbours cyanobacteria from the Candidatus
*Synechococcus spongiarum* clade.

Culture of the sponge under different light conditions resulted in the
observation that when kept under dark conditions the sponge bleaches,
turning completely white after ca. 12 weeks under darkness; the "normal"
color of the sponge varies between green, blue and violet. Changes in
color are accompanied by changes in size and general morphology. For
instance, the sponges form elongated projections; the normal shape of
the sponge is lobated and plate like.

It is not clear whether the cyanobacteria are still present in the
bleached sponge tissues. Bioanalyzer profiles of RNA extracts from
bleached sponges clearly showed that the samples lack the rRNA bands
belonging the cyanobacteria but qPCR experiments pointed towards the
presence of cyanobacteria in DNA extracts from sponges exposed to max. 6
weeeks darkeness [1]. However, this evidence is (in my opinion) not
conclusive due to a number of technical problems associated to these
experiments; RNA extractions only demonstrate that the symbionts are not
active and the qPCR experiments suffered from incosistent amplifying of
the samples. Another piece of evidence (weakly) pointing towards this is
the capacity of the sponges to regain color after 8-12 weeks re-exposure
to light; symbiont uptake from the water column cannot be excluded,
however, apparently *S. spongiarum* is not present in the water. All in
all, the evidence at hand point to the presence of the symbionts in some
kind of inactive (transcriptionl/metabolic) state in the sponge tissues
exposed to dark conditions.

Study rational
--------------

Since it is possible to manipulate the trascriptional/metabolic state of
the cyanobacterial symbionts through the exposure of the system to
darkness, CBAS represents an interesting system to assess how the sponge
reacts to symbiont inactivation. I use RNA-Seq data to *de novo*
assemble a reference transcriptome for this system and to assess changes
in the sponge gene expression profiles associated with different
symbiont transcriptional states. I hope these data contributes to our
understanding of the molecular mechanisms used by sponges to interact
with their microbiomes.

------------------------------------------------------------------------

Introduction
============

Materials and Methods
=====================

Sponge culture and experimental setting
---------------------------------------

Sponges are kept in a 200L salt water aquarium under a 12hr light:12hr
darkness regime using two T5 fluorescent lights. For all experiments one
large sponge was used and a hole-puncher was used to produce explants of
10 mm diameter [2]. The explants were place in a recipient with sand and
were allow to recover for 3-5 days under control light conditions.
Recovery was visually assesed by inspecting the border of the explants
to corroborate the tissue healed at the cutting points. After this, the
explants were moved to a section of the aquarium covered with a two
sheets of black plastic cloth. The plastic cloth serves to block the
light from the T5 lights; the amount of light (in klux) penetrating the
plastic sheets was ~0 klux while under normal conditions at least ~15
klux can be measured at the water surface. Twelve weeks after moving the
sponges to the dark side of the aquarium, the explants were flash-frozen
in liquid nitrogen and kept at -80 °C; at this point in time explants
kept under control conditions were also sampled. This experiment was
repeated several times[3]; the resulting tissue samples were used to
test different RNA extraction protocols and assess whether bleaching can
be reproduced in the aquarium.

RNA extraction, libray preparation and sequencing
-------------------------------------------------

Total RNA was extracted using a hybrid protocol that combines a regular
CTAB extraction and a spin-column clean-up step. Briefly, using mortar
and pestle, samples were pulverized in liquid nitrogen. The tissue
powder was lysed in 600 *μ*L warm (56 °C) CTAB-PVP-NaCl buffer
containing *β*-mercaptoethanol for 15-20 minutes with agitation (~500
rpm). After lysis, one volume acidic Phenol-Chloroform-Isoamyl alcohol
(25:24:1) was added to the samples and the extraction was vigorously
agitated until the liquid had a milky appearance. The samples were then
centrifuged at 14,000 rpm for 15 minutes to separate the phases. 400
*μ*L of the polar phase were then transfered to a new 2mL
microcentrifuge tube and the nucleic-acids were precipitated for 10
minutes using one volume isopropanol. The nucleic acids were recovered
by centrifugation (14,000rpm for 20 minutes at ~16°C) and the recovered
pellets were washed twice with 1 mL cold 80% ethanol. After washing,
ethanol droplets were removed with a 10 *μ*L pipette and the pellets
were air dried for ~10 minute before resuspending in 30 *μ*L
nuclease-free water.

After this initial extraction, the ZR-Duet DNA/RNA MiniPrep [4] was
used, following the manufacturers recommendations, to separate DNA and
RNA from the sample. The RNA was resuspended in 30 *μ*L nuclease-free
water and quality checked initially on 1% agarose gels and finally on a
Bioanalyzer 2100 Nano RNA chip. RNA concentration was measured on a
Nanodrop 1000. RIN values could be calculated for bleached samples only;
control samples have four rRNA peaks corresponding to the 16S, 18S, 23S
and 28S rRNA fragments and the Bioanalyzer cannot calculate their RIN
value. Quality of the control samples was assessed by overlaying
bleached samples with RIN and control samples without RIN and assessing
their similarity in terms of peak height and baseline level. Five
control and four bleached samples were sent on dry-ice to the EMBL
Genomics Core Facility [5] where they were used to produce
strand-specific libraries with ~110 base pairs (bp). These libraries
were multiplexed and pair-end sequenced (50bp reads) in two lanes of a
HiSeq 2500 (Illumina).

Transcriptome assembly, annotation
----------------------------------

Reads were quality controled using FastQC [6] and quality filtered using
the BioLite program filter\_illumina.cpp [7]. The surviving read pairs
from all libraries were used to produce (using the command cat) two
fastq files that were used for *de novo* transcriptome assembly in
Trinity v2.0.6 (using the --normalize\_reads flag). The resulting
contigs (with lenght &gt;=200bp) were annotated against the Uniprot
(SwissProt) [8] and *Amphimedon queenslandica* isoforms (Aqu2 proteins)
[9] using blastx v2.2.29+ with an expectation cutoff of 0.001. Only the
best match per contig was saved and the blast results were saved using
blast's XML format converted to a 25 column table [10]. Both tables
(i.e.
[UNIPROT](https://github.com/sevragorgia/CBAS/tree/master/Annotations/Uniprot)
and
[Aqu2](https://github.com/sevragorgia/CBAS/tree/master/Annotations/AQU2))
are available for download in the [project
repository](https://github.com/sevragorgia/CBAS). Gene Ontology [11]
annotations for the CBAS transcriptome were obtained by programmatically
querying the [QuickGO Webservice](http://www.ebi.ac.uk/QuickGO/) with
the transcriptome's UNIPROT annotations [12] and a custom [perl
script](https://github.com/sevragorgia/CBAS/blob/master/Scripts/Get_GO_Annotations.pl).
For each CBAS transcript the "component", "function" and "process" GO
terms associated with its UNIPROT best match were stored in the [project
repository](https://github.com/sevragorgia/CBAS/tree/master/Annotations/GOs)
as independent tab-separated files that can be easily modified to use as
input files in TopGO [13].

In addition, transcripts were translated using the program
TransDecoder.LongOrfs [14] and the resulting cds, mRNA, bed, gff3 and
pep files stored in the [project
repository](https://github.com/sevragorgia/CBAS/tree/master/Annotations/Transdecoder)
and used to annotate the transcriptome against the
[Pfam](http://pfam.xfam.org/) and [KEEG](http://www.genome.jp/kegg/)
databases. For this, the perl script
[pfam\_scan.pl](ftp://ftp.ebi.ac.uk/pub/databases/Pfam/Tools/) and the
webservice [BlastKOALA](http://www.kegg.jp/blastkoala/) were used; the
resulting output files can be found in the [project
repository.](https://github.com/sevragorgia/CBAS) Finally, the assembled
transcripts were blasted (blastn) against a bacterial genomes database
[15].

All annotations were used to build an [annotation
meta-table](https://github.com/sevragorgia/CBAS/blob/master/Annotations/Metatable/)
using a [custom made perl
script](https://github.com/sevragorgia/CBAS/blob/master/Scripts/Create_Transcriptome_Annotation_Table.pl)
available in the project repository.

Transcriptome completeness assessment
-------------------------------------

The CBAS transcriptome completeness was assessed by blasting against the
CEGMA gene set of [Parra et al.
2007](http://bioinformatics.oxfordjournals.org/content/23/9/1061.abstract)[16]
using tblastn with an e-value of 1*e*<sup>−</sup>19 as implemented in
the scrtipt
[find\_cegma\_genes.py](https://bitbucket.org/wrf/galaxy-files/src). For
details about the method see [17].

Differential gene expression analysis
-------------------------------------

Using the *de novo* assembled transcriptome, the individual libraries
(i.e. control and bleached sponges) were mapped with RSEM and a
transcript by sample count matrix was derived from the "gene" counts
with the program
[abundance\_estimates\_to\_matrix.pl](https://github.com/sevragorgia/CBAS/tree/master/Counts/RSEM)
provided as part of the Trinity package [18]. The [count
matrix](https://github.com/sevragorgia/CBAS/tree/master/Counts/) was
used to find differentially expressed genes in control *vs.* bleached
sponges (model = ~Treatment). For this, the package
[DeSeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
was used. In brief, transcripts with 0 counts over all samples were
removed from the matrix and Principal Component Analysis was used to
assess the global expression pattern of control *vs.* bleached sponges
and identify potential outliers. Size effects and dispersions were
estimated with the DeSe2 methods **estimateSizeFactors** and
**estimateDispersions** and differentially expressed transcripts were
then infered using a Wald test (DeSeq2 method **nbinomWaldTest**). The
resulting p-values were adjusted using the Benjamin-Hochberg correction.
Transcripts with **log fold change &lt; − 2** or **log fold change
&gt;2** and **adjusted p-value &lt;0.01** were considered as
differentially expressed for further analyses.

A stand-alone R script to replicate the analysis is available in the
[project
repository](https://github.com/sevragorgia/CBAS/blob/master/Scripts) and
can be run using the [count
matrix](https://github.com/sevragorgia/CBAS/tree/master/Counts/CBAS_Bleaching_RSEM_Expression_Matrix.counts)
and [sample
information](https://github.com/sevragorgia/CBAS/tree/master/Counts/CBAS_Bleaching_RSEM_Expression_Matrix.info)
provided. And in this version of the manuscript, the code used for
several calculations is provided embedded as Rmarkdown snippets in the
source code of this document.

Gene Ontology Term enrichment analysis
--------------------------------------

In order to assess whether certain Gene Ontology terms are "enriched"
among the set of differentially expressed genes a GO-term enrichment
analysis was done using the R package
[TopGO](https://bioconductor.org/packages/release/bioc/html/topGO.html).
For this analysis, "background sets" of not differentially expressed
transcripts with similar expression profile to the set of differentially
expressed transcripts were created using the method **genefinder** of
the package
[genefilter](http://bioconductor.org/packages/release/bioc/html/genefilter.html)
with the Manhattan distance. Background transcript sets were created for
all differentially expressed transcripts (i.e. transcripts with log fold
change &gt;|2| and *p* &lt; 0.01), for differentially upregulated genes
(i.e. transcripts with log fold change &gt;2 and *p* &lt; 0.01) and for
differentially downregulated genes (i.e. transcripts with log fold
change &lt; − 2 and *p* &lt; 0.01). The background sets were then used
to test for enrichment of certain GO-terms in the global set of
differentially expressed transcripts (log fold change &gt;|2| and
*p* &lt; 0.01), the set of upregulated genes (log fold change &gt;2 and
*p* &lt; 0.01) and the set of downregulated genes (log fold change
&lt; − 2 and *p* &lt; 0.01). For the enrichment analysis, a Fisher exact
test with the *classic* algorithm implemented in the method **runTest**
of
[TopGo](https://bioconductor.org/packages/release/bioc/html/topGO.html)
was used.

Pfam domain enrichment analysis
-------------------------------

In addition to the GO-Term enrichment analysis, and due to the
difficulties in assigning a meaning to the GO annotations in our sponge
model, a Pfam domain enrichment analysis was also done. Essentially the
Pfam enrichment analysis works in a similar way that the GO-term
enrichment analysis. First, a background set of genes not differentially
expressed but showing similar expression patterns than the DEGs is
calculated and used to provide a population of transcripts against which
the differentially expressed transcripts can be compared. Once the
background population of transcripts is available, the frequency with
which a given domain is found in the background set of transcripts is
compared with the frequency with which the same domain is found in the
set of differentially expressed transcripts using an hypergeometric
test. To account for multiple comparisons, the p-values are adjusted
using the Benjamini-Hochberg correctino. The background distribution
used for the Pfam enrichment domains was the same used for the analysis
of GO-term enrichment.

Results
=======

The transcriptome of the Common Blue Aquarium Sponge, a model cyanosponge
-------------------------------------------------------------------------

Per library, we obtained on average 26,385,834 (±3, 360, 778) pairs of
reads (50 bp long). After cleaning, each library had an average of
24,942,084 (±3, 420, 083) surviving pairs. The concatenated dataset thus
had 199,536,668 pairs. Table 1 contains the information on the number of
sequenced and surviving reads per library.

**Table 1.** Reads by library before and after quality filtering.

<table>
<thead>
<tr class="header">
<th align="center">Library</th>
<th align="center">Condition</th>
<th align="center">Sequenced reads</th>
<th align="center">Reads post cleaning</th>
<th align="center">% Kept</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="center">1</td>
<td align="center">Control</td>
<td align="center">33,836,615</td>
<td align="center">32,397,726</td>
<td align="center">95.75</td>
</tr>
<tr class="even">
<td align="center">2</td>
<td align="center">Control</td>
<td align="center">26,740,977</td>
<td align="center">25,252,825</td>
<td align="center">94.43</td>
</tr>
<tr class="odd">
<td align="center">3</td>
<td align="center">Control</td>
<td align="center">22,823,936</td>
<td align="center">21,016,976</td>
<td align="center">92.08</td>
</tr>
<tr class="even">
<td align="center">4</td>
<td align="center">Control</td>
<td align="center">24,166,051</td>
<td align="center">22,655,601</td>
<td align="center">93.75</td>
</tr>
<tr class="odd">
<td align="center">5</td>
<td align="center">Bleached</td>
<td align="center">26,579,563</td>
<td align="center">25,067,523</td>
<td align="center">94.31</td>
</tr>
<tr class="even">
<td align="center">6</td>
<td align="center">Bleached</td>
<td align="center">25,448,522</td>
<td align="center">24,229,217</td>
<td align="center">95.21</td>
</tr>
<tr class="odd">
<td align="center">7</td>
<td align="center">Bleached</td>
<td align="center">22,823,936</td>
<td align="center">21,016,976</td>
<td align="center">92.08</td>
</tr>
<tr class="even">
<td align="center">8</td>
<td align="center">Bleached</td>
<td align="center">24,166,051</td>
<td align="center">22,655,601</td>
<td align="center">93.75</td>
</tr>
<tr class="odd">
<td align="center">Total Reads</td>
<td align="center"></td>
<td align="center">211,086,672</td>
<td align="center">199,536,668</td>
<td align="center">94.53</td>
</tr>
</tbody>
</table>

------------------------------------------------------------------------

*De novo* transcriptome assembly using the concatenated dataset resulted
in 128,686 transcript &gt;200bp. The N50 of the transcriptome was 1,281
bp, and the median lenght of the assembled transcripts was 453 bp
(MAD=0.2078891). More details about the assembly can be found in Table 2
and the length distribution of the transcripts can be visualized in
Figure 1.

**Table 2.** Statistics for the *de novo* assembled CBAS transcriptome

<table>
<thead>
<tr class="header">
<th align="center">Statistic</th>
<th align="center">Obtained value</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="center">GC content</td>
<td align="center">40.0</td>
</tr>
<tr class="even">
<td align="center">N50</td>
<td align="center">1,281</td>
</tr>
<tr class="odd">
<td align="center">Max. length</td>
<td align="center">19,309</td>
</tr>
<tr class="even">
<td align="center">Mean length</td>
<td align="center">803</td>
</tr>
<tr class="odd">
<td align="center">Median length</td>
<td align="center">453</td>
</tr>
<tr class="even">
<td align="center">Min. length</td>
<td align="center">224</td>
</tr>
<tr class="odd">
<td align="center">Number of Transcripts</td>
<td align="center">128,686</td>
</tr>
</tbody>
</table>

------------------------------------------------------------------------

<img src="CBAS_DeSeq2_Analysis_files/figure-markdown_strict/transcript_length_fig-1.png" style="display: block; margin: auto;" />

**Figure 1.** Transcript length distribution of CBAS reference
transcriptome. Maximum length allowed equal 10,000 base-pairs.

------------------------------------------------------------------------

Of the collection of transcripts, 37.75% could be translated into
proteins by Transdecoder. The majority of the translated transcripts
were *Complete ORFs*, according to Transdecoder. The *Transdecoder ORF
Types* distribution of the translated transcripts can be found in Table
3.

**Table 3.** *Transdecoder ORF Type* of the translated transcripts.

<table>
<thead>
<tr class="header">
<th align="center">ORF Type</th>
<th align="center">Count</th>
<th align="center">%</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="center">3' partial</td>
<td align="center">6734</td>
<td align="center">13.86</td>
</tr>
<tr class="even">
<td align="center">5' partial</td>
<td align="center">10535</td>
<td align="center">21.69</td>
</tr>
<tr class="odd">
<td align="center">Internal</td>
<td align="center">13473</td>
<td align="center">27.74</td>
</tr>
<tr class="even">
<td align="center">Complete</td>
<td align="center">17833</td>
<td align="center">36.71</td>
</tr>
<tr class="odd">
<td align="center">Total</td>
<td align="center">48575</td>
<td align="center">100</td>
</tr>
</tbody>
</table>

------------------------------------------------------------------------

Regarding the annotation of the translated transcripts [19], 29.44% had
a matching *Uniprot* annotation. The number of transcripts matching an
*Amphimedon queenslandica* (Aqu2) protein was slightly higher, 36.21%.
*Gene Ontology Component*, *Function* and *Process* annotations could be
retrieved for 26.3%, 26.49% and 26.76% of the transcripts, respectively.
Additionally, 23.45% of the transcripts could be annotated with *Pfam*
domain information. In contrast, only 5.8% of the assembled transcripts
could be annotated against the *KEGG* database.

<img src="CBAS_DeSeq2_Analysis_files/figure-markdown_strict/annotation_success_fig-1.png" style="display: block; margin: auto;" />
**Figure 2.** Annotation success by database for transcripts translated
by Transdecoder.

------------------------------------------------------------------------

Regarding annotation success, different *Transdecoder ORF Types* had
different annotation success rates (Fig. 3) with *complete ORFs* having
higher annotation sucess than other *ORF types* and *3' Partial ORFs*
having the lowest annotation success.

<img src="CBAS_DeSeq2_Analysis_files/figure-markdown_strict/annotation_count_byDB_Fig-1.png" style="display: block; margin: auto;" />
**Figure 3.** Annotation count by database for different types of ORF as
defined by Transdecoder.

------------------------------------------------------------------------

For annotation against Uniprot and Aqu2 we used a permissive evalue
(0.001) to query these databases. Thus, we risk obtaining annotations
that are only partial matches or that are of doubtful homology. Despite
our permissive annotation strategy, the mean evalues obtained for the
annotations against Uniprot and Aqu2 had clearly lower values than the
threshold set (i.e.
2.519351 × 10<sup>−7</sup> ± 1.127118 × 10<sup>−6</sup> and
2.4146826 × 10<sup>−7</sup> ± 1.1139696 × 10<sup>−6</sup>,
respectively). Interestintly, the mean evalues obtained for both Uniprot
and Aqu2 annotations were similar for all *Transdecoder ORF types*
(Table 4) and the maximum evalue obtained after querying these two
databases was 1*x*10<sup>−05</sup> [20].

**Table 4.** Mean evalue by Transdecoder ORF type for annotations
obtained from the Uniprot or Aqu2 databases.

<table>
<thead>
<tr class="header">
<th align="center">Transdecoder ORF Type</th>
<th align="center">Database</th>
<th align="center">Mean evalue (SD)</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="center">3' partial</td>
<td align="center">Uniprot</td>
<td align="center"><span class="math inline">2.1738267 × 10<sup>−7</sup></span> (<span class="math inline">1.0499955 × 10<sup>−6</sup></span>)</td>
</tr>
<tr class="even">
<td align="center">5' partial</td>
<td align="center">Uniprot</td>
<td align="center"><span class="math inline">2.1534162 × 10<sup>−7</sup></span> (<span class="math inline">1.0780619 × 10<sup>−6</sup></span>)</td>
</tr>
<tr class="odd">
<td align="center">Internal</td>
<td align="center">Uniprot</td>
<td align="center"><span class="math inline">4.9614047 × 10<sup>−7</sup></span> (<span class="math inline">1.5576029 × 10<sup>−6</sup></span>)</td>
</tr>
<tr class="even">
<td align="center">Complete</td>
<td align="center">Uniprot</td>
<td align="center"><span class="math inline">1.6779664 × 10<sup>−7</sup></span> (<span class="math inline">8.9792016 × 10<sup>−7</sup></span>)</td>
</tr>
<tr class="odd">
<td align="center"></td>
<td align="center"></td>
</tr>
<tr class="even">
<td align="center">3' partial</td>
<td align="center">Aqu2</td>
<td align="center"><span class="math inline">2.2338358 × 10<sup>−7</sup></span> (<span class="math inline">1.0988336 × 10<sup>−6</sup></span>)</td>
</tr>
<tr class="odd">
<td align="center">5' partial</td>
<td align="center">Aqu2</td>
<td align="center"><span class="math inline">1.8405011 × 10<sup>−7</sup></span> (<span class="math inline">9.6965115 × 10<sup>−7</sup></span>)</td>
</tr>
<tr class="even">
<td align="center">Internal</td>
<td align="center">Aqu2</td>
<td align="center"><span class="math inline">4.5642628 × 10<sup>−7</sup></span> (<span class="math inline">1.5267768 × 10<sup>−6</sup></span>)</td>
</tr>
<tr class="odd">
<td align="center">Complete</td>
<td align="center">Aqu2</td>
<td align="center"><span class="math inline">1.6698938 × 10<sup>−7</sup></span> (<span class="math inline">9.0483509 × 10<sup>−7</sup></span>)</td>
</tr>
</tbody>
</table>

------------------------------------------------------------------------

Finally, in terms of completeness, we recovered 243 out of 248 KOGs
(i.e. 97.98%) from the CBAS transcriptome. Yet, the number of
transcripts classified as *High-confidence, full length matches* was
only 175 (i.e. 70.56%). In addition, six transcripts were categorized as
*Probable full length matches*, six more were tagged as matching a KOG
that was slightly shorter than the query and 4 transcripts matched a KOG
that was much longer than the query. Thus, in total 191 (i.e. 77.02%)
transcripts can be considered high confidence KOG matches. Transcripts
matching KOGs that were much shorter and transcripts probably
representing missassemblies amount to 22 and 30, respectively (i.e.
20.97%).

Global gene expression patterns change due to symbiont inactivation
-------------------------------------------------------------------

The count matrix used to investigate changes in gene expression in
Bleached *vs.* Control sponges had 100098 rows with counts. Of these,
however, only 52567 had counts in all samples (i.e. the transcript was
detected in all sequenced libraries). The inspection of the cummulative
density distribution and of the density distribution of normalized
counts (Fig. 5) for either all rows (not shown) or only the non-zero
rows in the matrix showed similar patterns for all libraries indicating
that they were sequenced at similar depths and can be compared safely.

<img src="CBAS_DeSeq2_Analysis_files/figure-markdown_strict/normalized_counts_cumulativeDist_qc_figure-1.png" style="display: block; margin: auto;" />
**Figure 5.** Cummulative density distribution of normalized counts and
density distribution of normalized counts of RNA-Seq libraries prepared
samples for *bleached* and *control* **CBAS** samples.

------------------------------------------------------------------------

A pairwise comparison of the individual libraries (Fig 6.) provided
further evidence for the technical similarity of the libraries.

<img src="CBAS_DeSeq2_Analysis_files/figure-markdown_strict/MDPlots_library_01_qc_figure-1.png" style="display: block; margin: auto;" />

<img src="CBAS_DeSeq2_Analysis_files/figure-markdown_strict/MDPlots_library_02_qc_figure-1.png" style="display: block; margin: auto;" />

<img src="CBAS_DeSeq2_Analysis_files/figure-markdown_strict/MDPlots_library_03_qc_figure-1.png" style="display: block; margin: auto;" />

<img src="CBAS_DeSeq2_Analysis_files/figure-markdown_strict/MDPlots_library_04_qc_figure-1.png" style="display: block; margin: auto;" />

<img src="CBAS_DeSeq2_Analysis_files/figure-markdown_strict/MDPlots_library_05_qc_figure-1.png" style="display: block; margin: auto;" />

<img src="CBAS_DeSeq2_Analysis_files/figure-markdown_strict/MDPlots_library_06_qc_figure-1.png" style="display: block; margin: auto;" />

<img src="CBAS_DeSeq2_Analysis_files/figure-markdown_strict/MDPlots_library_07_qc_figure-1.png" style="display: block; margin: auto;" />

<img src="CBAS_DeSeq2_Analysis_files/figure-markdown_strict/MDPlots_library_08_qc_figure-1.png" style="display: block; margin: auto;" />
**Figure 6.** Pairwise library comparisons based on the relation between
the *Mean log-fold* and the *Mean (log-transformed) counts*. In general,
no trend should be observed in the pairwise comparisons.

------------------------------------------------------------------------

Libraries were also compared in terms of the *FPKM* and *TPM*
distribution of the transcripts (Fig. 7). This analysis revealed that
two libraries (i.e. CBAS\_Control\_2 and CBAS\_Control\_3) differed
markedly from all other sequenced libraries. Thus, we restricted the
count matrix to include only *trinity genes* that could be translated
using Transdecoder (29,705 transcripts) . After this, the libraries were
comparable in terms of the distribution of TPM and FPKM values. The
analysis of differential expression was conducted on this reduced data
matrix.

<img src="CBAS_DeSeq2_Analysis_files/figure-markdown_strict/TPM_figure-1.png" style="display: block; margin: auto;" />
**Figure 7.** TPM and FPKM distribution in datasets including all
trinity genes (All transcripts) and only trinity genes that could be
translated with trandescoder (Transdecoded transcripts). In general,
libraries should have similar TPM/FPKM distributions to be comparable.

------------------------------------------------------------------------

For the analysis of differentially expressed genes, the count matrix was
transformed using the *rlog* transformation available in the R package
DeSeq2. A cluster analysis based on the between-sample distances
calculated using the *rlog* transformed counts grouped samples in two
groups matching the *control* and *treatment* groups. Interestingly one
treatment sample behaved in an anomalous manner, clustering with neither
the treatment nor the control group. This sample (i.e.
CBAS\_Bleached\_5) was also found to be different in a Principal
Component Analysis (PCA) conducted on the *rlog* transformed count
matrix.

------------------------------------------------------------------------

<img src="CBAS_DeSeq2_Analysis_files/figure-markdown_strict/sample_heat_pca_plot-1.png" style="display: block; margin: auto;" /><img src="CBAS_DeSeq2_Analysis_files/figure-markdown_strict/sample_heat_pca_plot-2.png" style="display: block; margin: auto;" />

**Figure 9.** Cluster and Principal Component Analysis of the *rlog*
transformed counts for *bleached vs. control* CBAS specimens.

------------------------------------------------------------------------

If sample *CBAS\_Bleached\_5* is removed from the analysis, *bleached*
samples are grouped together while *control* sponges are divided in two
groups in both the cluster analysis and the principal component
analysis.

<img src="CBAS_DeSeq2_Analysis_files/figure-markdown_strict/drop_anomalous_sample_and_replot-1.png" style="display: block; margin: auto;" /><img src="CBAS_DeSeq2_Analysis_files/figure-markdown_strict/drop_anomalous_sample_and_replot-2.png" style="display: block; margin: auto;" />

**Figure 10.** Cluster and Principal Component Analysis of the *rlog*
transformed counts for *bleached vs. control* CBAS specimens after
removing sample *CBAS\_Bleached\_5*.

------------------------------------------------------------------------

We detected significant differences in the global expression patters
bewteen treated and control sponges (Adonis PseudoF = 5.1164, df = 1,
p-value = 0.039). In this reduced dataset, a total of 5508 transcripts
were differentially expressed between the *control* and *bleached*
sponges (Fig. 11). Of these transcripts, 1166 were overexpressed (with a
log2 fold change higher or equal to 1) in *bleached vs. control* sponges
and 2563 were expressed at lower levels (a log2 fold change lower or
equal than -1) in *bleached vs. control* specimens.

<img src="CBAS_DeSeq2_Analysis_files/figure-markdown_strict/MA_plot-1.png" style="display: block; margin: auto;" /><img src="CBAS_DeSeq2_Analysis_files/figure-markdown_strict/MA_plot-2.png" style="display: block; margin: auto;" />

**Figure 11.** MA-plot showing the relation between log-fold change and
mean count for each analyzed CBAS transcript. Transcripts significant at
0.01 are highlighted in red.

------------------------------------------------------------------------

It is worth noting that DeSeq2 was able to adecuately model the variance
in the dataset and could correctly calculate the p-values for the
analyzed transcripts (Fig 12).

<img src="CBAS_DeSeq2_Analysis_files/figure-markdown_strict/probability_plot-1.png" style="display: block; margin: auto;" />

**Figure 12.** Histogram of the p-values obtained for the transcripts
analyzed using DeSeq2. A uniform distribution between 0 and 1, with a
peak at zero, should be observed if the variance of the dataset could be
correctly modelled.

------------------------------------------------------------------------

The expression of immune response related transcripts changes between symbiontic states
---------------------------------------------------------------------------------------

For Gene Ontology (GO) enrichment analyses we were able to select a
background set of transcripts showing similar expression patterns but no
significant change between conditions for all transcripts as well as for
overexpressed and underexpressed transcripts (Fig. 13). Enriched
Function, Process and Compartment GO terms for over- and underexpressed
genes can be found in Tables 5-10.

<img src="CBAS_DeSeq2_Analysis_files/figure-markdown_strict/plot_back_vs_forground_dist-1.png" style="display: block; margin: auto;" />

Table 5: GO-Term functions overrepresented in the set of overexpressed
DEGs.

<table>
<thead>
<tr class="header">
<th></th>
<th align="left">GO.ID</th>
<th align="left">Term</th>
<th align="right">Annotated</th>
<th align="right">Significant</th>
<th align="right">Expected</th>
<th align="left">classic</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>12</td>
<td align="left"><a href="GO:0004197" class="uri">GO:0004197</a></td>
<td align="left">cysteine-type endopeptidase activity</td>
<td align="right">20</td>
<td align="right">14</td>
<td align="right">5.59</td>
<td align="left">0.00010</td>
</tr>
<tr class="even">
<td>13</td>
<td align="left"><a href="GO:0016787" class="uri">GO:0016787</a></td>
<td align="left">hydrolase activity</td>
<td align="right">473</td>
<td align="right">162</td>
<td align="right">132.26</td>
<td align="left">0.00031</td>
</tr>
<tr class="odd">
<td>14</td>
<td align="left"><a href="GO:0004843" class="uri">GO:0004843</a></td>
<td align="left">ubiquitin-specific protease activity</td>
<td align="right">15</td>
<td align="right">11</td>
<td align="right">4.19</td>
<td align="left">0.00032</td>
</tr>
<tr class="even">
<td>15</td>
<td align="left"><a href="GO:0036459" class="uri">GO:0036459</a></td>
<td align="left">ubiquitinyl hydrolase activity</td>
<td align="right">15</td>
<td align="right">11</td>
<td align="right">4.19</td>
<td align="left">0.00032</td>
</tr>
<tr class="odd">
<td>16</td>
<td align="left"><a href="GO:0004872" class="uri">GO:0004872</a></td>
<td align="left">receptor activity</td>
<td align="right">166</td>
<td align="right">66</td>
<td align="right">46.42</td>
<td align="left">0.00039</td>
</tr>
<tr class="even">
<td>17</td>
<td align="left"><a href="GO:0019783" class="uri">GO:0019783</a></td>
<td align="left">ubiquitin-like protein-specific protease...</td>
<td align="right">18</td>
<td align="right">12</td>
<td align="right">5.03</td>
<td align="left">0.00067</td>
</tr>
<tr class="odd">
<td>18</td>
<td align="left"><a href="GO:0035375" class="uri">GO:0035375</a></td>
<td align="left">zymogen binding</td>
<td align="right">19</td>
<td align="right">12</td>
<td align="right">5.31</td>
<td align="left">0.00137</td>
</tr>
<tr class="even">
<td>19</td>
<td align="left"><a href="GO:0005200" class="uri">GO:0005200</a></td>
<td align="left">structural constituent of cytoskeleton</td>
<td align="right">29</td>
<td align="right">16</td>
<td align="right">8.11</td>
<td align="left">0.00172</td>
</tr>
<tr class="odd">
<td>20</td>
<td align="left"><a href="GO:0004871" class="uri">GO:0004871</a></td>
<td align="left">signal transducer activity</td>
<td align="right">128</td>
<td align="right">51</td>
<td align="right">35.79</td>
<td align="left">0.00176</td>
</tr>
<tr class="even">
<td>21</td>
<td align="left"><a href="GO:0003735" class="uri">GO:0003735</a></td>
<td align="left">structural constituent of ribosome</td>
<td align="right">52</td>
<td align="right">24</td>
<td align="right">14.54</td>
<td align="left">0.00346</td>
</tr>
<tr class="odd">
<td>22</td>
<td align="left"><a href="GO:0004175" class="uri">GO:0004175</a></td>
<td align="left">endopeptidase activity</td>
<td align="right">154</td>
<td align="right">58</td>
<td align="right">43.06</td>
<td align="left">0.00407</td>
</tr>
<tr class="even">
<td>23</td>
<td align="left"><a href="GO:0004715" class="uri">GO:0004715</a></td>
<td align="left">non-membrane spanning protein tyrosine k...</td>
<td align="right">19</td>
<td align="right">11</td>
<td align="right">5.31</td>
<td align="left">0.00569</td>
</tr>
<tr class="odd">
<td>24</td>
<td align="left"><a href="GO:0005507" class="uri">GO:0005507</a></td>
<td align="left">copper ion binding</td>
<td align="right">17</td>
<td align="right">10</td>
<td align="right">4.75</td>
<td align="left">0.00721</td>
</tr>
<tr class="even">
<td>25</td>
<td align="left"><a href="GO:0042393" class="uri">GO:0042393</a></td>
<td align="left">histone binding</td>
<td align="right">15</td>
<td align="right">9</td>
<td align="right">4.19</td>
<td align="left">0.00910</td>
</tr>
<tr class="odd">
<td>26</td>
<td align="left"><a href="GO:0043167" class="uri">GO:0043167</a></td>
<td align="left">ion binding</td>
<td align="right">986</td>
<td align="right">299</td>
<td align="right">275.70</td>
<td align="left">0.00931</td>
</tr>
<tr class="even">
<td>27</td>
<td align="left"><a href="GO:0038023" class="uri">GO:0038023</a></td>
<td align="left">signaling receptor activity</td>
<td align="right">110</td>
<td align="right">42</td>
<td align="right">30.76</td>
<td align="left">0.01074</td>
</tr>
<tr class="odd">
<td>28</td>
<td align="left"><a href="GO:0004497" class="uri">GO:0004497</a></td>
<td align="left">monooxygenase activity</td>
<td align="right">21</td>
<td align="right">11</td>
<td align="right">5.87</td>
<td align="left">0.01505</td>
</tr>
<tr class="even">
<td>29</td>
<td align="left"><a href="GO:0004888" class="uri">GO:0004888</a></td>
<td align="left">transmembrane signaling receptor activit...</td>
<td align="right">100</td>
<td align="right">38</td>
<td align="right">27.96</td>
<td align="left">0.01626</td>
</tr>
<tr class="odd">
<td>30</td>
<td align="left"><a href="GO:0005201" class="uri">GO:0005201</a></td>
<td align="left">extracellular matrix structural constitu...</td>
<td align="right">14</td>
<td align="right">8</td>
<td align="right">3.91</td>
<td align="left">0.02025</td>
</tr>
<tr class="even">
<td>31</td>
<td align="left"><a href="GO:0005007" class="uri">GO:0005007</a></td>
<td align="left">fibroblast growth factor-activated recep...</td>
<td align="right">7</td>
<td align="right">5</td>
<td align="right">1.96</td>
<td align="left">0.02097</td>
</tr>
<tr class="odd">
<td>32</td>
<td align="left"><a href="GO:0005218" class="uri">GO:0005218</a></td>
<td align="left">intracellular ligand-gated calcium chann...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.84</td>
<td align="left">0.02177</td>
</tr>
<tr class="even">
<td>33</td>
<td align="left"><a href="GO:0005219" class="uri">GO:0005219</a></td>
<td align="left">ryanodine-sensitive calcium-release chan...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.84</td>
<td align="left">0.02177</td>
</tr>
<tr class="odd">
<td>34</td>
<td align="left"><a href="GO:0016716" class="uri">GO:0016716</a></td>
<td align="left">oxidoreductase activity, acting on paire...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.84</td>
<td align="left">0.02177</td>
</tr>
<tr class="even">
<td>35</td>
<td align="left"><a href="GO:0017147" class="uri">GO:0017147</a></td>
<td align="left">Wnt-protein binding</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.84</td>
<td align="left">0.02177</td>
</tr>
<tr class="odd">
<td>36</td>
<td align="left"><a href="GO:0032500" class="uri">GO:0032500</a></td>
<td align="left">muramyl dipeptide binding</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.84</td>
<td align="left">0.02177</td>
</tr>
<tr class="even">
<td>37</td>
<td align="left"><a href="GO:0070064" class="uri">GO:0070064</a></td>
<td align="left">proline-rich region binding</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.84</td>
<td align="left">0.02177</td>
</tr>
<tr class="odd">
<td>38</td>
<td align="left"><a href="GO:0001882" class="uri">GO:0001882</a></td>
<td align="left">nucleoside binding</td>
<td align="right">409</td>
<td align="right">131</td>
<td align="right">114.36</td>
<td align="left">0.02294</td>
</tr>
<tr class="even">
<td>39</td>
<td align="left"><a href="GO:0001883" class="uri">GO:0001883</a></td>
<td align="left">purine nucleoside binding</td>
<td align="right">409</td>
<td align="right">131</td>
<td align="right">114.36</td>
<td align="left">0.02294</td>
</tr>
<tr class="odd">
<td>40</td>
<td align="left"><a href="GO:0015278" class="uri">GO:0015278</a></td>
<td align="left">calcium-release channel activity</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.40</td>
<td align="left">0.02357</td>
</tr>
<tr class="even">
<td>41</td>
<td align="left"><a href="GO:0016682" class="uri">GO:0016682</a></td>
<td align="left">oxidoreductase activity, acting on diphe...</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.40</td>
<td align="left">0.02357</td>
</tr>
<tr class="odd">
<td>42</td>
<td align="left"><a href="GO:0016860" class="uri">GO:0016860</a></td>
<td align="left">intramolecular oxidoreductase activity</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.40</td>
<td align="left">0.02357</td>
</tr>
<tr class="even">
<td>43</td>
<td align="left"><a href="GO:0043130" class="uri">GO:0043130</a></td>
<td align="left">ubiquitin binding</td>
<td align="right">17</td>
<td align="right">9</td>
<td align="right">4.75</td>
<td align="left">0.02532</td>
</tr>
<tr class="odd">
<td>44</td>
<td align="left"><a href="GO:0004713" class="uri">GO:0004713</a></td>
<td align="left">protein tyrosine kinase activity</td>
<td align="right">54</td>
<td align="right">22</td>
<td align="right">15.10</td>
<td align="left">0.02734</td>
</tr>
<tr class="even">
<td>45</td>
<td align="left"><a href="GO:0032549" class="uri">GO:0032549</a></td>
<td align="left">ribonucleoside binding</td>
<td align="right">408</td>
<td align="right">130</td>
<td align="right">114.08</td>
<td align="left">0.02807</td>
</tr>
<tr class="odd">
<td>46</td>
<td align="left"><a href="GO:0032550" class="uri">GO:0032550</a></td>
<td align="left">purine ribonucleoside binding</td>
<td align="right">408</td>
<td align="right">130</td>
<td align="right">114.08</td>
<td align="left">0.02807</td>
</tr>
<tr class="even">
<td>47</td>
<td align="left"><a href="GO:0035639" class="uri">GO:0035639</a></td>
<td align="left">purine ribonucleoside triphosphate bindi...</td>
<td align="right">408</td>
<td align="right">130</td>
<td align="right">114.08</td>
<td align="left">0.02807</td>
</tr>
<tr class="odd">
<td>48</td>
<td align="left"><a href="GO:0017111" class="uri">GO:0017111</a></td>
<td align="left">nucleoside-triphosphatase activity</td>
<td align="right">145</td>
<td align="right">51</td>
<td align="right">40.54</td>
<td align="left">0.02936</td>
</tr>
<tr class="even">
<td>49</td>
<td align="left"><a href="GO:0016462" class="uri">GO:0016462</a></td>
<td align="left">pyrophosphatase activity</td>
<td align="right">149</td>
<td align="right">52</td>
<td align="right">41.66</td>
<td align="left">0.03244</td>
</tr>
<tr class="odd">
<td>50</td>
<td align="left"><a href="GO:0016817" class="uri">GO:0016817</a></td>
<td align="left">hydrolase activity, acting on acid anhyd...</td>
<td align="right">149</td>
<td align="right">52</td>
<td align="right">41.66</td>
<td align="left">0.03244</td>
</tr>
<tr class="even">
<td>51</td>
<td align="left"><a href="GO:0016818" class="uri">GO:0016818</a></td>
<td align="left">hydrolase activity, acting on acid anhyd...</td>
<td align="right">149</td>
<td align="right">52</td>
<td align="right">41.66</td>
<td align="left">0.03244</td>
</tr>
<tr class="odd">
<td>52</td>
<td align="left"><a href="GO:0032553" class="uri">GO:0032553</a></td>
<td align="left">ribonucleotide binding</td>
<td align="right">411</td>
<td align="right">130</td>
<td align="right">114.92</td>
<td align="left">0.03573</td>
</tr>
<tr class="even">
<td>53</td>
<td align="left"><a href="GO:0032555" class="uri">GO:0032555</a></td>
<td align="left">purine ribonucleotide binding</td>
<td align="right">411</td>
<td align="right">130</td>
<td align="right">114.92</td>
<td align="left">0.03573</td>
</tr>
<tr class="odd">
<td>54</td>
<td align="left"><a href="GO:0001948" class="uri">GO:0001948</a></td>
<td align="left">glycoprotein binding</td>
<td align="right">18</td>
<td align="right">9</td>
<td align="right">5.03</td>
<td align="left">0.03841</td>
</tr>
<tr class="even">
<td>55</td>
<td align="left"><a href="GO:0032182" class="uri">GO:0032182</a></td>
<td align="left">ubiquitin-like protein binding</td>
<td align="right">18</td>
<td align="right">9</td>
<td align="right">5.03</td>
<td align="left">0.03841</td>
</tr>
<tr class="odd">
<td>56</td>
<td align="left"><a href="GO:0017076" class="uri">GO:0017076</a></td>
<td align="left">purine nucleotide binding</td>
<td align="right">413</td>
<td align="right">130</td>
<td align="right">115.48</td>
<td align="left">0.04171</td>
</tr>
<tr class="even">
<td>57</td>
<td align="left"><a href="GO:0003887" class="uri">GO:0003887</a></td>
<td align="left">DNA-directed DNA polymerase activity</td>
<td align="right">47</td>
<td align="right">19</td>
<td align="right">13.14</td>
<td align="left">0.04220</td>
</tr>
<tr class="odd">
<td>58</td>
<td align="left"><a href="GO:0008046" class="uri">GO:0008046</a></td>
<td align="left">axon guidance receptor activity</td>
<td align="right">13</td>
<td align="right">7</td>
<td align="right">3.64</td>
<td align="left">0.04317</td>
</tr>
<tr class="even">
<td>59</td>
<td align="left"><a href="GO:0015085" class="uri">GO:0015085</a></td>
<td align="left">calcium ion transmembrane transporter ac...</td>
<td align="right">8</td>
<td align="right">5</td>
<td align="right">2.24</td>
<td align="left">0.04319</td>
</tr>
<tr class="odd">
<td>60</td>
<td align="left"><a href="GO:0051018" class="uri">GO:0051018</a></td>
<td align="left">protein kinase A binding</td>
<td align="right">8</td>
<td align="right">5</td>
<td align="right">2.24</td>
<td align="left">0.04319</td>
</tr>
</tbody>
</table>

Table 6: GO-Term processes overrepresented in the set of overexpressed
DEGs.

<table>
<thead>
<tr class="header">
<th></th>
<th align="left">GO.ID</th>
<th align="left">Term</th>
<th align="right">Annotated</th>
<th align="right">Significant</th>
<th align="right">Expected</th>
<th align="left">classic</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>7</td>
<td align="left"><a href="GO:0031344" class="uri">GO:0031344</a></td>
<td align="left">regulation of cell projection organizati...</td>
<td align="right">69</td>
<td align="right">33</td>
<td align="right">19.15</td>
<td align="left">0.00024</td>
</tr>
<tr class="even">
<td>8</td>
<td align="left"><a href="GO:0030595" class="uri">GO:0030595</a></td>
<td align="left">leukocyte chemotaxis</td>
<td align="right">15</td>
<td align="right">11</td>
<td align="right">4.16</td>
<td align="left">0.00030</td>
</tr>
<tr class="odd">
<td>9</td>
<td align="left"><a href="GO:0050731" class="uri">GO:0050731</a></td>
<td align="left">positive regulation of peptidyl-tyrosine...</td>
<td align="right">13</td>
<td align="right">10</td>
<td align="right">3.61</td>
<td align="left">0.00031</td>
</tr>
<tr class="even">
<td>10</td>
<td align="left"><a href="GO:0050920" class="uri">GO:0050920</a></td>
<td align="left">regulation of chemotaxis</td>
<td align="right">31</td>
<td align="right">18</td>
<td align="right">8.60</td>
<td align="left">0.00034</td>
</tr>
<tr class="odd">
<td>11</td>
<td align="left"><a href="GO:0006898" class="uri">GO:0006898</a></td>
<td align="left">receptor-mediated endocytosis</td>
<td align="right">60</td>
<td align="right">29</td>
<td align="right">16.65</td>
<td align="left">0.00045</td>
</tr>
<tr class="even">
<td>12</td>
<td align="left"><a href="GO:0050900" class="uri">GO:0050900</a></td>
<td align="left">leukocyte migration</td>
<td align="right">39</td>
<td align="right">21</td>
<td align="right">10.82</td>
<td align="left">0.00046</td>
</tr>
<tr class="odd">
<td>13</td>
<td align="left"><a href="GO:0050767" class="uri">GO:0050767</a></td>
<td align="left">regulation of neurogenesis</td>
<td align="right">77</td>
<td align="right">35</td>
<td align="right">21.37</td>
<td align="left">0.00052</td>
</tr>
<tr class="even">
<td>14</td>
<td align="left"><a href="GO:0009719" class="uri">GO:0009719</a></td>
<td align="left">response to endogenous stimulus</td>
<td align="right">130</td>
<td align="right">53</td>
<td align="right">36.08</td>
<td align="left">0.00062</td>
</tr>
<tr class="odd">
<td>15</td>
<td align="left"><a href="GO:0051259" class="uri">GO:0051259</a></td>
<td align="left">protein oligomerization</td>
<td align="right">61</td>
<td align="right">29</td>
<td align="right">16.93</td>
<td align="left">0.00064</td>
</tr>
<tr class="even">
<td>16</td>
<td align="left"><a href="GO:0021772" class="uri">GO:0021772</a></td>
<td align="left">olfactory bulb development</td>
<td align="right">16</td>
<td align="right">11</td>
<td align="right">4.44</td>
<td align="left">0.00072</td>
</tr>
<tr class="odd">
<td>17</td>
<td align="left"><a href="GO:0021988" class="uri">GO:0021988</a></td>
<td align="left">olfactory lobe development</td>
<td align="right">16</td>
<td align="right">11</td>
<td align="right">4.44</td>
<td align="left">0.00072</td>
</tr>
<tr class="even">
<td>18</td>
<td align="left"><a href="GO:0050730" class="uri">GO:0050730</a></td>
<td align="left">regulation of peptidyl-tyrosine phosphor...</td>
<td align="right">14</td>
<td align="right">10</td>
<td align="right">3.89</td>
<td align="left">0.00081</td>
</tr>
<tr class="odd">
<td>19</td>
<td align="left"><a href="GO:0006887" class="uri">GO:0006887</a></td>
<td align="left">exocytosis</td>
<td align="right">51</td>
<td align="right">25</td>
<td align="right">14.15</td>
<td align="left">0.00087</td>
</tr>
<tr class="even">
<td>20</td>
<td align="left"><a href="GO:0002688" class="uri">GO:0002688</a></td>
<td align="left">regulation of leukocyte chemotaxis</td>
<td align="right">10</td>
<td align="right">8</td>
<td align="right">2.78</td>
<td align="left">0.00087</td>
</tr>
<tr class="odd">
<td>21</td>
<td align="left"><a href="GO:0043552" class="uri">GO:0043552</a></td>
<td align="left">positive regulation of phosphatidylinosi...</td>
<td align="right">10</td>
<td align="right">8</td>
<td align="right">2.78</td>
<td align="left">0.00087</td>
</tr>
<tr class="even">
<td>22</td>
<td align="left"><a href="GO:0090218" class="uri">GO:0090218</a></td>
<td align="left">positive regulation of lipid kinase acti...</td>
<td align="right">10</td>
<td align="right">8</td>
<td align="right">2.78</td>
<td align="left">0.00087</td>
</tr>
<tr class="odd">
<td>23</td>
<td align="left"><a href="GO:1903727" class="uri">GO:1903727</a></td>
<td align="left">positive regulation of phospholipid meta...</td>
<td align="right">10</td>
<td align="right">8</td>
<td align="right">2.78</td>
<td align="left">0.00087</td>
</tr>
<tr class="even">
<td>24</td>
<td align="left"><a href="GO:0051289" class="uri">GO:0051289</a></td>
<td align="left">protein homotetramerization</td>
<td align="right">12</td>
<td align="right">9</td>
<td align="right">3.33</td>
<td align="left">0.00088</td>
</tr>
<tr class="odd">
<td>25</td>
<td align="left"><a href="GO:0030111" class="uri">GO:0030111</a></td>
<td align="left">regulation of Wnt signaling pathway</td>
<td align="right">38</td>
<td align="right">20</td>
<td align="right">10.55</td>
<td align="left">0.00092</td>
</tr>
<tr class="even">
<td>26</td>
<td align="left"><a href="GO:0051240" class="uri">GO:0051240</a></td>
<td align="left">positive regulation of multicellular org...</td>
<td align="right">111</td>
<td align="right">46</td>
<td align="right">30.80</td>
<td align="left">0.00094</td>
</tr>
<tr class="odd">
<td>27</td>
<td align="left"><a href="GO:0006935" class="uri">GO:0006935</a></td>
<td align="left">chemotaxis</td>
<td align="right">100</td>
<td align="right">42</td>
<td align="right">27.75</td>
<td align="left">0.00114</td>
</tr>
<tr class="even">
<td>28</td>
<td align="left"><a href="GO:0045664" class="uri">GO:0045664</a></td>
<td align="left">regulation of neuron differentiation</td>
<td align="right">66</td>
<td align="right">30</td>
<td align="right">18.32</td>
<td align="left">0.00132</td>
</tr>
<tr class="odd">
<td>29</td>
<td align="left"><a href="GO:0010975" class="uri">GO:0010975</a></td>
<td align="left">regulation of neuron projection developm...</td>
<td align="right">55</td>
<td align="right">26</td>
<td align="right">15.26</td>
<td align="left">0.00136</td>
</tr>
<tr class="even">
<td>30</td>
<td align="left"><a href="GO:1901652" class="uri">GO:1901652</a></td>
<td align="left">response to peptide</td>
<td align="right">34</td>
<td align="right">18</td>
<td align="right">9.44</td>
<td align="left">0.00152</td>
</tr>
<tr class="odd">
<td>31</td>
<td align="left"><a href="GO:1901653" class="uri">GO:1901653</a></td>
<td align="left">cellular response to peptide</td>
<td align="right">29</td>
<td align="right">16</td>
<td align="right">8.05</td>
<td align="left">0.00157</td>
</tr>
<tr class="even">
<td>32</td>
<td align="left"><a href="GO:0034135" class="uri">GO:0034135</a></td>
<td align="left">regulation of toll-like receptor 2 signa...</td>
<td align="right">5</td>
<td align="right">5</td>
<td align="right">1.39</td>
<td align="left">0.00162</td>
</tr>
<tr class="odd">
<td>33</td>
<td align="left"><a href="GO:0034136" class="uri">GO:0034136</a></td>
<td align="left">negative regulation of toll-like recepto...</td>
<td align="right">5</td>
<td align="right">5</td>
<td align="right">1.39</td>
<td align="left">0.00162</td>
</tr>
<tr class="even">
<td>34</td>
<td align="left"><a href="GO:0050766" class="uri">GO:0050766</a></td>
<td align="left">positive regulation of phagocytosis</td>
<td align="right">5</td>
<td align="right">5</td>
<td align="right">1.39</td>
<td align="left">0.00162</td>
</tr>
<tr class="odd">
<td>35</td>
<td align="left"><a href="GO:1902106" class="uri">GO:1902106</a></td>
<td align="left">negative regulation of leukocyte differe...</td>
<td align="right">5</td>
<td align="right">5</td>
<td align="right">1.39</td>
<td align="left">0.00162</td>
</tr>
<tr class="even">
<td>36</td>
<td align="left"><a href="GO:0060284" class="uri">GO:0060284</a></td>
<td align="left">regulation of cell development</td>
<td align="right">87</td>
<td align="right">37</td>
<td align="right">24.14</td>
<td align="left">0.00171</td>
</tr>
<tr class="odd">
<td>37</td>
<td align="left"><a href="GO:0035567" class="uri">GO:0035567</a></td>
<td align="left">non-canonical Wnt signaling pathway</td>
<td align="right">22</td>
<td align="right">13</td>
<td align="right">6.11</td>
<td align="left">0.00190</td>
</tr>
<tr class="even">
<td>38</td>
<td align="left"><a href="GO:0031346" class="uri">GO:0031346</a></td>
<td align="left">positive regulation of cell projection o...</td>
<td align="right">32</td>
<td align="right">17</td>
<td align="right">8.88</td>
<td align="left">0.00196</td>
</tr>
<tr class="odd">
<td>39</td>
<td align="left"><a href="GO:0060326" class="uri">GO:0060326</a></td>
<td align="left">cell chemotaxis</td>
<td align="right">27</td>
<td align="right">15</td>
<td align="right">7.49</td>
<td align="left">0.00201</td>
</tr>
<tr class="even">
<td>40</td>
<td align="left"><a href="GO:0032755" class="uri">GO:0032755</a></td>
<td align="left">positive regulation of interleukin-6 pro...</td>
<td align="right">13</td>
<td align="right">9</td>
<td align="right">3.61</td>
<td align="left">0.00215</td>
</tr>
<tr class="odd">
<td>41</td>
<td align="left"><a href="GO:0060401" class="uri">GO:0060401</a></td>
<td align="left">cytosolic calcium ion transport</td>
<td align="right">13</td>
<td align="right">9</td>
<td align="right">3.61</td>
<td align="left">0.00215</td>
</tr>
<tr class="even">
<td>42</td>
<td align="left"><a href="GO:0060402" class="uri">GO:0060402</a></td>
<td align="left">calcium ion transport into cytosol</td>
<td align="right">13</td>
<td align="right">9</td>
<td align="right">3.61</td>
<td align="left">0.00215</td>
</tr>
<tr class="odd">
<td>43</td>
<td align="left"><a href="GO:1901700" class="uri">GO:1901700</a></td>
<td align="left">response to oxygen-containing compound</td>
<td align="right">130</td>
<td align="right">51</td>
<td align="right">36.08</td>
<td align="left">0.00219</td>
</tr>
<tr class="even">
<td>44</td>
<td align="left"><a href="GO:0051960" class="uri">GO:0051960</a></td>
<td align="left">regulation of nervous system development</td>
<td align="right">94</td>
<td align="right">39</td>
<td align="right">26.09</td>
<td align="left">0.00225</td>
</tr>
<tr class="odd">
<td>45</td>
<td align="left"><a href="GO:0042330" class="uri">GO:0042330</a></td>
<td align="left">taxis</td>
<td align="right">103</td>
<td align="right">42</td>
<td align="right">28.58</td>
<td align="left">0.00229</td>
</tr>
<tr class="even">
<td>46</td>
<td align="left"><a href="GO:0034121" class="uri">GO:0034121</a></td>
<td align="left">regulation of toll-like receptor signali...</td>
<td align="right">7</td>
<td align="right">6</td>
<td align="right">1.94</td>
<td align="left">0.00239</td>
</tr>
<tr class="odd">
<td>47</td>
<td align="left"><a href="GO:0002685" class="uri">GO:0002685</a></td>
<td align="left">regulation of leukocyte migration</td>
<td align="right">11</td>
<td align="right">8</td>
<td align="right">3.05</td>
<td align="left">0.00242</td>
</tr>
<tr class="even">
<td>48</td>
<td align="left"><a href="GO:0007405" class="uri">GO:0007405</a></td>
<td align="left">neuroblast proliferation</td>
<td align="right">11</td>
<td align="right">8</td>
<td align="right">3.05</td>
<td align="left">0.00242</td>
</tr>
<tr class="odd">
<td>49</td>
<td align="left"><a href="GO:0021889" class="uri">GO:0021889</a></td>
<td align="left">olfactory bulb interneuron differentiati...</td>
<td align="right">11</td>
<td align="right">8</td>
<td align="right">3.05</td>
<td align="left">0.00242</td>
</tr>
<tr class="even">
<td>50</td>
<td align="left"><a href="GO:0021891" class="uri">GO:0021891</a></td>
<td align="left">olfactory bulb interneuron development</td>
<td align="right">11</td>
<td align="right">8</td>
<td align="right">3.05</td>
<td align="left">0.00242</td>
</tr>
<tr class="odd">
<td>51</td>
<td align="left"><a href="GO:0032677" class="uri">GO:0032677</a></td>
<td align="left">regulation of interleukin-8 production</td>
<td align="right">11</td>
<td align="right">8</td>
<td align="right">3.05</td>
<td align="left">0.00242</td>
</tr>
<tr class="even">
<td>52</td>
<td align="left"><a href="GO:0032757" class="uri">GO:0032757</a></td>
<td align="left">positive regulation of interleukin-8 pro...</td>
<td align="right">11</td>
<td align="right">8</td>
<td align="right">3.05</td>
<td align="left">0.00242</td>
</tr>
<tr class="odd">
<td>53</td>
<td align="left"><a href="GO:0043550" class="uri">GO:0043550</a></td>
<td align="left">regulation of lipid kinase activity</td>
<td align="right">11</td>
<td align="right">8</td>
<td align="right">3.05</td>
<td align="left">0.00242</td>
</tr>
<tr class="even">
<td>54</td>
<td align="left"><a href="GO:0043551" class="uri">GO:0043551</a></td>
<td align="left">regulation of phosphatidylinositol 3-kin...</td>
<td align="right">11</td>
<td align="right">8</td>
<td align="right">3.05</td>
<td align="left">0.00242</td>
</tr>
<tr class="odd">
<td>55</td>
<td align="left"><a href="GO:1903725" class="uri">GO:1903725</a></td>
<td align="left">regulation of phospholipid metabolic pro...</td>
<td align="right">11</td>
<td align="right">8</td>
<td align="right">3.05</td>
<td align="left">0.00242</td>
</tr>
<tr class="even">
<td>56</td>
<td align="left"><a href="GO:0006953" class="uri">GO:0006953</a></td>
<td align="left">acute-phase response</td>
<td align="right">9</td>
<td align="right">7</td>
<td align="right">2.50</td>
<td align="left">0.00256</td>
</tr>
<tr class="odd">
<td>57</td>
<td align="left"><a href="GO:0007204" class="uri">GO:0007204</a></td>
<td align="left">positive regulation of cytosolic calcium...</td>
<td align="right">18</td>
<td align="right">11</td>
<td align="right">5.00</td>
<td align="left">0.00299</td>
</tr>
<tr class="even">
<td>58</td>
<td align="left"><a href="GO:0016055" class="uri">GO:0016055</a></td>
<td align="left">Wnt signaling pathway</td>
<td align="right">52</td>
<td align="right">24</td>
<td align="right">14.43</td>
<td align="left">0.00310</td>
</tr>
<tr class="odd">
<td>59</td>
<td align="left"><a href="GO:0007411" class="uri">GO:0007411</a></td>
<td align="left">axon guidance</td>
<td align="right">72</td>
<td align="right">31</td>
<td align="right">19.98</td>
<td align="left">0.00320</td>
</tr>
<tr class="even">
<td>60</td>
<td align="left"><a href="GO:0097485" class="uri">GO:0097485</a></td>
<td align="left">neuron projection guidance</td>
<td align="right">72</td>
<td align="right">31</td>
<td align="right">19.98</td>
<td align="left">0.00320</td>
</tr>
<tr class="odd">
<td>61</td>
<td align="left"><a href="GO:0040012" class="uri">GO:0040012</a></td>
<td align="left">regulation of locomotion</td>
<td align="right">58</td>
<td align="right">26</td>
<td align="right">16.10</td>
<td align="left">0.00349</td>
</tr>
<tr class="even">
<td>62</td>
<td align="left"><a href="GO:0043405" class="uri">GO:0043405</a></td>
<td align="left">regulation of MAP kinase activity</td>
<td align="right">36</td>
<td align="right">18</td>
<td align="right">9.99</td>
<td align="left">0.00350</td>
</tr>
<tr class="odd">
<td>63</td>
<td align="left"><a href="GO:0050795" class="uri">GO:0050795</a></td>
<td align="left">regulation of behavior</td>
<td align="right">36</td>
<td align="right">18</td>
<td align="right">9.99</td>
<td align="left">0.00350</td>
</tr>
<tr class="even">
<td>64</td>
<td align="left"><a href="GO:0043043" class="uri">GO:0043043</a></td>
<td align="left">peptide biosynthetic process</td>
<td align="right">108</td>
<td align="right">43</td>
<td align="right">29.97</td>
<td align="left">0.00351</td>
</tr>
<tr class="odd">
<td>65</td>
<td align="left"><a href="GO:0048839" class="uri">GO:0048839</a></td>
<td align="left">inner ear development</td>
<td align="right">39</td>
<td align="right">19</td>
<td align="right">10.82</td>
<td align="left">0.00395</td>
</tr>
<tr class="even">
<td>66</td>
<td align="left"><a href="GO:0010976" class="uri">GO:0010976</a></td>
<td align="left">positive regulation of neuron projection...</td>
<td align="right">26</td>
<td align="right">14</td>
<td align="right">7.22</td>
<td align="left">0.00418</td>
</tr>
<tr class="odd">
<td>67</td>
<td align="left"><a href="GO:0030163" class="uri">GO:0030163</a></td>
<td align="left">protein catabolic process</td>
<td align="right">103</td>
<td align="right">41</td>
<td align="right">28.58</td>
<td align="left">0.00436</td>
</tr>
<tr class="even">
<td>68</td>
<td align="left"><a href="GO:0014074" class="uri">GO:0014074</a></td>
<td align="left">response to purine-containing compound</td>
<td align="right">14</td>
<td align="right">9</td>
<td align="right">3.89</td>
<td align="left">0.00456</td>
</tr>
<tr class="odd">
<td>69</td>
<td align="left"><a href="GO:0097529" class="uri">GO:0097529</a></td>
<td align="left">myeloid leukocyte migration</td>
<td align="right">14</td>
<td align="right">9</td>
<td align="right">3.89</td>
<td align="left">0.00456</td>
</tr>
<tr class="even">
<td>70</td>
<td align="left"><a href="GO:0009725" class="uri">GO:0009725</a></td>
<td align="left">response to hormone</td>
<td align="right">59</td>
<td align="right">26</td>
<td align="right">16.37</td>
<td align="left">0.00465</td>
</tr>
<tr class="odd">
<td>71</td>
<td align="left"><a href="GO:0043434" class="uri">GO:0043434</a></td>
<td align="left">response to peptide hormone</td>
<td align="right">29</td>
<td align="right">15</td>
<td align="right">8.05</td>
<td align="left">0.00503</td>
</tr>
<tr class="even">
<td>72</td>
<td align="left"><a href="GO:0016579" class="uri">GO:0016579</a></td>
<td align="left">protein deubiquitination</td>
<td align="right">19</td>
<td align="right">11</td>
<td align="right">5.27</td>
<td align="left">0.00534</td>
</tr>
<tr class="odd">
<td>73</td>
<td align="left"><a href="GO:0051239" class="uri">GO:0051239</a></td>
<td align="left">regulation of multicellular organismal p...</td>
<td align="right">260</td>
<td align="right">90</td>
<td align="right">72.15</td>
<td align="left">0.00541</td>
</tr>
<tr class="even">
<td>74</td>
<td align="left"><a href="GO:0048667" class="uri">GO:0048667</a></td>
<td align="left">cell morphogenesis involved in neuron di...</td>
<td align="right">98</td>
<td align="right">39</td>
<td align="right">27.20</td>
<td align="left">0.00542</td>
</tr>
<tr class="odd">
<td>75</td>
<td align="left"><a href="GO:0051094" class="uri">GO:0051094</a></td>
<td align="left">positive regulation of developmental pro...</td>
<td align="right">98</td>
<td align="right">39</td>
<td align="right">27.20</td>
<td align="left">0.00542</td>
</tr>
<tr class="even">
<td>76</td>
<td align="left"><a href="GO:0000186" class="uri">GO:0000186</a></td>
<td align="left">activation of MAPKK activity</td>
<td align="right">12</td>
<td align="right">8</td>
<td align="right">3.33</td>
<td align="left">0.00551</td>
</tr>
<tr class="odd">
<td>77</td>
<td align="left"><a href="GO:0021872" class="uri">GO:0021872</a></td>
<td align="left">forebrain generation of neurons</td>
<td align="right">12</td>
<td align="right">8</td>
<td align="right">3.33</td>
<td align="left">0.00551</td>
</tr>
<tr class="even">
<td>78</td>
<td align="left"><a href="GO:0021879" class="uri">GO:0021879</a></td>
<td align="left">forebrain neuron differentiation</td>
<td align="right">12</td>
<td align="right">8</td>
<td align="right">3.33</td>
<td align="left">0.00551</td>
</tr>
<tr class="odd">
<td>79</td>
<td align="left"><a href="GO:0021884" class="uri">GO:0021884</a></td>
<td align="left">forebrain neuron development</td>
<td align="right">12</td>
<td align="right">8</td>
<td align="right">3.33</td>
<td align="left">0.00551</td>
</tr>
<tr class="even">
<td>80</td>
<td align="left"><a href="GO:0032637" class="uri">GO:0032637</a></td>
<td align="left">interleukin-8 production</td>
<td align="right">12</td>
<td align="right">8</td>
<td align="right">3.33</td>
<td align="left">0.00551</td>
</tr>
<tr class="odd">
<td>81</td>
<td align="left"><a href="GO:0051208" class="uri">GO:0051208</a></td>
<td align="left">sequestering of calcium ion</td>
<td align="right">12</td>
<td align="right">8</td>
<td align="right">3.33</td>
<td align="left">0.00551</td>
</tr>
<tr class="even">
<td>82</td>
<td align="left"><a href="GO:0051209" class="uri">GO:0051209</a></td>
<td align="left">release of sequestered calcium ion into ...</td>
<td align="right">12</td>
<td align="right">8</td>
<td align="right">3.33</td>
<td align="left">0.00551</td>
</tr>
<tr class="odd">
<td>83</td>
<td align="left"><a href="GO:0051282" class="uri">GO:0051282</a></td>
<td align="left">regulation of sequestering of calcium io...</td>
<td align="right">12</td>
<td align="right">8</td>
<td align="right">3.33</td>
<td align="left">0.00551</td>
</tr>
<tr class="even">
<td>84</td>
<td align="left"><a href="GO:0051283" class="uri">GO:0051283</a></td>
<td align="left">negative regulation of sequestering of c...</td>
<td align="right">12</td>
<td align="right">8</td>
<td align="right">3.33</td>
<td align="left">0.00551</td>
</tr>
<tr class="odd">
<td>85</td>
<td align="left"><a href="GO:0031644" class="uri">GO:0031644</a></td>
<td align="left">regulation of neurological system proces...</td>
<td align="right">4</td>
<td align="right">4</td>
<td align="right">1.11</td>
<td align="left">0.00588</td>
</tr>
<tr class="even">
<td>86</td>
<td align="left"><a href="GO:0071313" class="uri">GO:0071313</a></td>
<td align="left">cellular response to caffeine</td>
<td align="right">4</td>
<td align="right">4</td>
<td align="right">1.11</td>
<td align="left">0.00588</td>
</tr>
<tr class="odd">
<td>87</td>
<td align="left"><a href="GO:0071415" class="uri">GO:0071415</a></td>
<td align="left">cellular response to purine-containing c...</td>
<td align="right">4</td>
<td align="right">4</td>
<td align="right">1.11</td>
<td align="left">0.00588</td>
</tr>
<tr class="even">
<td>88</td>
<td align="left"><a href="GO:0090022" class="uri">GO:0090022</a></td>
<td align="left">regulation of neutrophil chemotaxis</td>
<td align="right">4</td>
<td align="right">4</td>
<td align="right">1.11</td>
<td align="left">0.00588</td>
</tr>
<tr class="odd">
<td>89</td>
<td align="left"><a href="GO:1902622" class="uri">GO:1902622</a></td>
<td align="left">regulation of neutrophil migration</td>
<td align="right">4</td>
<td align="right">4</td>
<td align="right">1.11</td>
<td align="left">0.00588</td>
</tr>
<tr class="even">
<td>90</td>
<td align="left"><a href="GO:0071900" class="uri">GO:0071900</a></td>
<td align="left">regulation of protein serine/threonine k...</td>
<td align="right">43</td>
<td align="right">20</td>
<td align="right">11.93</td>
<td align="left">0.00609</td>
</tr>
<tr class="odd">
<td>91</td>
<td align="left"><a href="GO:0018108" class="uri">GO:0018108</a></td>
<td align="left">peptidyl-tyrosine phosphorylation</td>
<td align="right">63</td>
<td align="right">27</td>
<td align="right">17.48</td>
<td align="left">0.00628</td>
</tr>
<tr class="even">
<td>92</td>
<td align="left"><a href="GO:0030178" class="uri">GO:0030178</a></td>
<td align="left">negative regulation of Wnt signaling pat...</td>
<td align="right">27</td>
<td align="right">14</td>
<td align="right">7.49</td>
<td align="left">0.00651</td>
</tr>
<tr class="odd">
<td>93</td>
<td align="left"><a href="GO:0043406" class="uri">GO:0043406</a></td>
<td align="left">positive regulation of MAP kinase activi...</td>
<td align="right">27</td>
<td align="right">14</td>
<td align="right">7.49</td>
<td align="left">0.00651</td>
</tr>
<tr class="even">
<td>94</td>
<td align="left"><a href="GO:0071495" class="uri">GO:0071495</a></td>
<td align="left">cellular response to endogenous stimulus</td>
<td align="right">105</td>
<td align="right">41</td>
<td align="right">29.14</td>
<td align="left">0.00651</td>
</tr>
<tr class="odd">
<td>95</td>
<td align="left"><a href="GO:0006412" class="uri">GO:0006412</a></td>
<td align="left">translation</td>
<td align="right">102</td>
<td align="right">40</td>
<td align="right">28.31</td>
<td align="left">0.00658</td>
</tr>
<tr class="even">
<td>96</td>
<td align="left"><a href="GO:0006508" class="uri">GO:0006508</a></td>
<td align="left">proteolysis</td>
<td align="right">298</td>
<td align="right">101</td>
<td align="right">82.70</td>
<td align="left">0.00662</td>
</tr>
<tr class="odd">
<td>97</td>
<td align="left"><a href="GO:0044257" class="uri">GO:0044257</a></td>
<td align="left">cellular protein catabolic process</td>
<td align="right">90</td>
<td align="right">36</td>
<td align="right">24.98</td>
<td align="left">0.00677</td>
</tr>
<tr class="even">
<td>98</td>
<td align="left"><a href="GO:0048814" class="uri">GO:0048814</a></td>
<td align="left">regulation of dendrite morphogenesis</td>
<td align="right">17</td>
<td align="right">10</td>
<td align="right">4.72</td>
<td align="left">0.00680</td>
</tr>
<tr class="odd">
<td>99</td>
<td align="left"><a href="GO:0032868" class="uri">GO:0032868</a></td>
<td align="left">response to insulin</td>
<td align="right">22</td>
<td align="right">12</td>
<td align="right">6.11</td>
<td align="left">0.00694</td>
</tr>
<tr class="even">
<td>100</td>
<td align="left"><a href="GO:0060828" class="uri">GO:0060828</a></td>
<td align="left">regulation of canonical Wnt signaling pa...</td>
<td align="right">22</td>
<td align="right">12</td>
<td align="right">6.11</td>
<td align="left">0.00694</td>
</tr>
<tr class="odd">
<td>101</td>
<td align="left"><a href="GO:0070646" class="uri">GO:0070646</a></td>
<td align="left">protein modification by small protein re...</td>
<td align="right">22</td>
<td align="right">12</td>
<td align="right">6.11</td>
<td align="left">0.00694</td>
</tr>
<tr class="even">
<td>102</td>
<td align="left"><a href="GO:0010769" class="uri">GO:0010769</a></td>
<td align="left">regulation of cell morphogenesis involve...</td>
<td align="right">52</td>
<td align="right">23</td>
<td align="right">14.43</td>
<td align="left">0.00721</td>
</tr>
<tr class="odd">
<td>103</td>
<td align="left"><a href="GO:0030593" class="uri">GO:0030593</a></td>
<td align="left">neutrophil chemotaxis</td>
<td align="right">8</td>
<td align="right">6</td>
<td align="right">2.22</td>
<td align="left">0.00733</td>
</tr>
<tr class="even">
<td>104</td>
<td align="left"><a href="GO:0070536" class="uri">GO:0070536</a></td>
<td align="left">protein K63-linked deubiquitination</td>
<td align="right">8</td>
<td align="right">6</td>
<td align="right">2.22</td>
<td align="left">0.00733</td>
</tr>
<tr class="odd">
<td>105</td>
<td align="left"><a href="GO:1990266" class="uri">GO:1990266</a></td>
<td align="left">neutrophil migration</td>
<td align="right">8</td>
<td align="right">6</td>
<td align="right">2.22</td>
<td align="left">0.00733</td>
</tr>
<tr class="even">
<td>106</td>
<td align="left"><a href="GO:0002532" class="uri">GO:0002532</a></td>
<td align="left">production of molecular mediator involve...</td>
<td align="right">6</td>
<td align="right">5</td>
<td align="right">1.67</td>
<td align="left">0.00750</td>
</tr>
<tr class="odd">
<td>107</td>
<td align="left"><a href="GO:0002825" class="uri">GO:0002825</a></td>
<td align="left">regulation of T-helper 1 type immune res...</td>
<td align="right">6</td>
<td align="right">5</td>
<td align="right">1.67</td>
<td align="left">0.00750</td>
</tr>
<tr class="even">
<td>108</td>
<td align="left"><a href="GO:0002827" class="uri">GO:0002827</a></td>
<td align="left">positive regulation of T-helper 1 type i...</td>
<td align="right">6</td>
<td align="right">5</td>
<td align="right">1.67</td>
<td align="left">0.00750</td>
</tr>
<tr class="odd">
<td>109</td>
<td align="left"><a href="GO:0007260" class="uri">GO:0007260</a></td>
<td align="left">tyrosine phosphorylation of STAT protein</td>
<td align="right">6</td>
<td align="right">5</td>
<td align="right">1.67</td>
<td align="left">0.00750</td>
</tr>
<tr class="even">
<td>110</td>
<td align="left"><a href="GO:0034122" class="uri">GO:0034122</a></td>
<td align="left">negative regulation of toll-like recepto...</td>
<td align="right">6</td>
<td align="right">5</td>
<td align="right">1.67</td>
<td align="left">0.00750</td>
</tr>
<tr class="odd">
<td>111</td>
<td align="left"><a href="GO:0035397" class="uri">GO:0035397</a></td>
<td align="left">helper T cell enhancement of adaptive im...</td>
<td align="right">6</td>
<td align="right">5</td>
<td align="right">1.67</td>
<td align="left">0.00750</td>
</tr>
<tr class="even">
<td>112</td>
<td align="left"><a href="GO:0044243" class="uri">GO:0044243</a></td>
<td align="left">multicellular organismal catabolic proce...</td>
<td align="right">6</td>
<td align="right">5</td>
<td align="right">1.67</td>
<td align="left">0.00750</td>
</tr>
<tr class="odd">
<td>113</td>
<td align="left"><a href="GO:0050764" class="uri">GO:0050764</a></td>
<td align="left">regulation of phagocytosis</td>
<td align="right">6</td>
<td align="right">5</td>
<td align="right">1.67</td>
<td align="left">0.00750</td>
</tr>
<tr class="even">
<td>114</td>
<td align="left"><a href="GO:0072604" class="uri">GO:0072604</a></td>
<td align="left">interleukin-6 secretion</td>
<td align="right">6</td>
<td align="right">5</td>
<td align="right">1.67</td>
<td align="left">0.00750</td>
</tr>
<tr class="odd">
<td>115</td>
<td align="left"><a href="GO:0072606" class="uri">GO:0072606</a></td>
<td align="left">interleukin-8 secretion</td>
<td align="right">6</td>
<td align="right">5</td>
<td align="right">1.67</td>
<td align="left">0.00750</td>
</tr>
<tr class="even">
<td>116</td>
<td align="left"><a href="GO:2000482" class="uri">GO:2000482</a></td>
<td align="left">regulation of interleukin-8 secretion</td>
<td align="right">6</td>
<td align="right">5</td>
<td align="right">1.67</td>
<td align="left">0.00750</td>
</tr>
<tr class="odd">
<td>117</td>
<td align="left"><a href="GO:2000484" class="uri">GO:2000484</a></td>
<td align="left">positive regulation of interleukin-8 sec...</td>
<td align="right">6</td>
<td align="right">5</td>
<td align="right">1.67</td>
<td align="left">0.00750</td>
</tr>
<tr class="even">
<td>118</td>
<td align="left"><a href="GO:2000778" class="uri">GO:2000778</a></td>
<td align="left">positive regulation of interleukin-6 sec...</td>
<td align="right">6</td>
<td align="right">5</td>
<td align="right">1.67</td>
<td align="left">0.00750</td>
</tr>
<tr class="odd">
<td>119</td>
<td align="left"><a href="GO:0043604" class="uri">GO:0043604</a></td>
<td align="left">amide biosynthetic process</td>
<td align="right">115</td>
<td align="right">44</td>
<td align="right">31.91</td>
<td align="left">0.00753</td>
</tr>
<tr class="even">
<td>120</td>
<td align="left"><a href="GO:0016192" class="uri">GO:0016192</a></td>
<td align="left">vesicle-mediated transport</td>
<td align="right">204</td>
<td align="right">72</td>
<td align="right">56.61</td>
<td align="left">0.00769</td>
</tr>
<tr class="odd">
<td>121</td>
<td align="left"><a href="GO:0043583" class="uri">GO:0043583</a></td>
<td align="left">ear development</td>
<td align="right">41</td>
<td align="right">19</td>
<td align="right">11.38</td>
<td align="left">0.00782</td>
</tr>
<tr class="even">
<td>122</td>
<td align="left"><a href="GO:0032940" class="uri">GO:0032940</a></td>
<td align="left">secretion by cell</td>
<td align="right">103</td>
<td align="right">40</td>
<td align="right">28.58</td>
<td align="left">0.00799</td>
</tr>
<tr class="odd">
<td>123</td>
<td align="left"><a href="GO:0018212" class="uri">GO:0018212</a></td>
<td align="left">peptidyl-tyrosine modification</td>
<td align="right">64</td>
<td align="right">27</td>
<td align="right">17.76</td>
<td align="left">0.00809</td>
</tr>
<tr class="even">
<td>124</td>
<td align="left"><a href="GO:0022604" class="uri">GO:0022604</a></td>
<td align="left">regulation of cell morphogenesis</td>
<td align="right">67</td>
<td align="right">28</td>
<td align="right">18.59</td>
<td align="left">0.00821</td>
</tr>
<tr class="odd">
<td>125</td>
<td align="left"><a href="GO:0032870" class="uri">GO:0032870</a></td>
<td align="left">cellular response to hormone stimulus</td>
<td align="right">44</td>
<td align="right">20</td>
<td align="right">12.21</td>
<td align="left">0.00834</td>
</tr>
<tr class="even">
<td>126</td>
<td align="left"><a href="GO:0051603" class="uri">GO:0051603</a></td>
<td align="left">proteolysis involved in cellular protein...</td>
<td align="right">88</td>
<td align="right">35</td>
<td align="right">24.42</td>
<td align="left">0.00838</td>
</tr>
<tr class="odd">
<td>127</td>
<td align="left"><a href="GO:0060070" class="uri">GO:0060070</a></td>
<td align="left">canonical Wnt signaling pathway</td>
<td align="right">25</td>
<td align="right">13</td>
<td align="right">6.94</td>
<td align="left">0.00843</td>
</tr>
<tr class="even">
<td>128</td>
<td align="left"><a href="GO:0071375" class="uri">GO:0071375</a></td>
<td align="left">cellular response to peptide hormone sti...</td>
<td align="right">25</td>
<td align="right">13</td>
<td align="right">6.94</td>
<td align="left">0.00843</td>
</tr>
<tr class="odd">
<td>129</td>
<td align="left"><a href="GO:0045597" class="uri">GO:0045597</a></td>
<td align="left">positive regulation of cell differentiat...</td>
<td align="right">79</td>
<td align="right">32</td>
<td align="right">21.92</td>
<td align="left">0.00844</td>
</tr>
<tr class="even">
<td>130</td>
<td align="left"><a href="GO:0036336" class="uri">GO:0036336</a></td>
<td align="left">dendritic cell migration</td>
<td align="right">15</td>
<td align="right">9</td>
<td align="right">4.16</td>
<td align="left">0.00862</td>
</tr>
<tr class="odd">
<td>131</td>
<td align="left"><a href="GO:0050773" class="uri">GO:0050773</a></td>
<td align="left">regulation of dendrite development</td>
<td align="right">20</td>
<td align="right">11</td>
<td align="right">5.55</td>
<td align="left">0.00895</td>
</tr>
<tr class="even">
<td>132</td>
<td align="left"><a href="GO:0051480" class="uri">GO:0051480</a></td>
<td align="left">cytosolic calcium ion homeostasis</td>
<td align="right">20</td>
<td align="right">11</td>
<td align="right">5.55</td>
<td align="left">0.00895</td>
</tr>
<tr class="odd">
<td>133</td>
<td align="left"><a href="GO:0045596" class="uri">GO:0045596</a></td>
<td align="left">negative regulation of cell differentiat...</td>
<td align="right">53</td>
<td align="right">23</td>
<td align="right">14.71</td>
<td align="left">0.00951</td>
</tr>
<tr class="even">
<td>134</td>
<td align="left"><a href="GO:0010770" class="uri">GO:0010770</a></td>
<td align="left">positive regulation of cell morphogenesi...</td>
<td align="right">28</td>
<td align="right">14</td>
<td align="right">7.77</td>
<td align="left">0.00978</td>
</tr>
<tr class="odd">
<td>135</td>
<td align="left"><a href="GO:0045666" class="uri">GO:0045666</a></td>
<td align="left">positive regulation of neuron differenti...</td>
<td align="right">28</td>
<td align="right">14</td>
<td align="right">7.77</td>
<td align="left">0.00978</td>
</tr>
<tr class="even">
<td>136</td>
<td align="left"><a href="GO:0022603" class="uri">GO:0022603</a></td>
<td align="left">regulation of anatomical structure morph...</td>
<td align="right">101</td>
<td align="right">39</td>
<td align="right">28.03</td>
<td align="left">0.00979</td>
</tr>
<tr class="odd">
<td>137</td>
<td align="left"><a href="GO:0016358" class="uri">GO:0016358</a></td>
<td align="left">dendrite development</td>
<td align="right">39</td>
<td align="right">18</td>
<td align="right">10.82</td>
<td align="left">0.01005</td>
</tr>
<tr class="even">
<td>138</td>
<td align="left"><a href="GO:0042742" class="uri">GO:0042742</a></td>
<td align="left">defense response to bacterium</td>
<td align="right">39</td>
<td align="right">18</td>
<td align="right">10.82</td>
<td align="left">0.01005</td>
</tr>
<tr class="odd">
<td>139</td>
<td align="left"><a href="GO:0051962" class="uri">GO:0051962</a></td>
<td align="left">positive regulation of nervous system de...</td>
<td align="right">42</td>
<td align="right">19</td>
<td align="right">11.66</td>
<td align="left">0.01068</td>
</tr>
<tr class="even">
<td>140</td>
<td align="left"><a href="GO:0045595" class="uri">GO:0045595</a></td>
<td align="left">regulation of cell differentiation</td>
<td align="right">155</td>
<td align="right">56</td>
<td align="right">43.01</td>
<td align="left">0.01089</td>
</tr>
<tr class="odd">
<td>141</td>
<td align="left"><a href="GO:0051238" class="uri">GO:0051238</a></td>
<td align="left">sequestering of metal ion</td>
<td align="right">13</td>
<td align="right">8</td>
<td align="right">3.61</td>
<td align="left">0.01089</td>
</tr>
<tr class="even">
<td>142</td>
<td align="left"><a href="GO:2000050" class="uri">GO:2000050</a></td>
<td align="left">regulation of non-canonical Wnt signalin...</td>
<td align="right">13</td>
<td align="right">8</td>
<td align="right">3.61</td>
<td align="left">0.01089</td>
</tr>
<tr class="odd">
<td>143</td>
<td align="left"><a href="GO:2000051" class="uri">GO:2000051</a></td>
<td align="left">negative regulation of non-canonical Wnt...</td>
<td align="right">13</td>
<td align="right">8</td>
<td align="right">3.61</td>
<td align="left">0.01089</td>
</tr>
<tr class="even">
<td>144</td>
<td align="left"><a href="GO:0032635" class="uri">GO:0032635</a></td>
<td align="left">interleukin-6 production</td>
<td align="right">31</td>
<td align="right">15</td>
<td align="right">8.60</td>
<td align="left">0.01097</td>
</tr>
<tr class="odd">
<td>145</td>
<td align="left"><a href="GO:0032675" class="uri">GO:0032675</a></td>
<td align="left">regulation of interleukin-6 production</td>
<td align="right">31</td>
<td align="right">15</td>
<td align="right">8.60</td>
<td align="left">0.01097</td>
</tr>
<tr class="even">
<td>146</td>
<td align="left"><a href="GO:0071902" class="uri">GO:0071902</a></td>
<td align="left">positive regulation of protein serine/th...</td>
<td align="right">31</td>
<td align="right">15</td>
<td align="right">8.60</td>
<td align="left">0.01097</td>
</tr>
<tr class="odd">
<td>147</td>
<td align="left"><a href="GO:0051241" class="uri">GO:0051241</a></td>
<td align="left">negative regulation of multicellular org...</td>
<td align="right">111</td>
<td align="right">42</td>
<td align="right">30.80</td>
<td align="left">0.01116</td>
</tr>
<tr class="even">
<td>148</td>
<td align="left"><a href="GO:0050808" class="uri">GO:0050808</a></td>
<td align="left">synapse organization</td>
<td align="right">45</td>
<td align="right">20</td>
<td align="right">12.49</td>
<td align="left">0.01122</td>
</tr>
<tr class="odd">
<td>149</td>
<td align="left"><a href="GO:0045087" class="uri">GO:0045087</a></td>
<td align="left">innate immune response</td>
<td align="right">99</td>
<td align="right">38</td>
<td align="right">27.47</td>
<td align="left">0.01197</td>
</tr>
<tr class="even">
<td>150</td>
<td align="left"><a href="GO:0007416" class="uri">GO:0007416</a></td>
<td align="left">synapse assembly</td>
<td align="right">34</td>
<td align="right">16</td>
<td align="right">9.44</td>
<td align="left">0.01202</td>
</tr>
<tr class="odd">
<td>151</td>
<td align="left"><a href="GO:0014070" class="uri">GO:0014070</a></td>
<td align="left">response to organic cyclic compound</td>
<td align="right">54</td>
<td align="right">23</td>
<td align="right">14.99</td>
<td align="left">0.01235</td>
</tr>
<tr class="even">
<td>152</td>
<td align="left"><a href="GO:0044057" class="uri">GO:0044057</a></td>
<td align="left">regulation of system process</td>
<td align="right">26</td>
<td align="right">13</td>
<td align="right">7.22</td>
<td align="left">0.01268</td>
</tr>
<tr class="odd">
<td>153</td>
<td align="left"><a href="GO:0050707" class="uri">GO:0050707</a></td>
<td align="left">regulation of cytokine secretion</td>
<td align="right">26</td>
<td align="right">13</td>
<td align="right">7.22</td>
<td align="left">0.01268</td>
</tr>
<tr class="even">
<td>154</td>
<td align="left"><a href="GO:0002684" class="uri">GO:0002684</a></td>
<td align="left">positive regulation of immune system pro...</td>
<td align="right">84</td>
<td align="right">33</td>
<td align="right">23.31</td>
<td align="left">0.01275</td>
</tr>
<tr class="odd">
<td>155</td>
<td align="left"><a href="GO:0006909" class="uri">GO:0006909</a></td>
<td align="left">phagocytosis</td>
<td align="right">37</td>
<td align="right">17</td>
<td align="right">10.27</td>
<td align="left">0.01291</td>
</tr>
<tr class="even">
<td>156</td>
<td align="left"><a href="GO:0051093" class="uri">GO:0051093</a></td>
<td align="left">negative regulation of developmental pro...</td>
<td align="right">75</td>
<td align="right">30</td>
<td align="right">20.81</td>
<td align="left">0.01300</td>
</tr>
<tr class="odd">
<td>157</td>
<td align="left"><a href="GO:0007267" class="uri">GO:0007267</a></td>
<td align="left">cell-cell signaling</td>
<td align="right">69</td>
<td align="right">28</td>
<td align="right">19.15</td>
<td align="left">0.01303</td>
</tr>
<tr class="even">
<td>158</td>
<td align="left"><a href="GO:0002682" class="uri">GO:0002682</a></td>
<td align="left">regulation of immune system process</td>
<td align="right">112</td>
<td align="right">42</td>
<td align="right">31.08</td>
<td align="left">0.01326</td>
</tr>
<tr class="odd">
<td>159</td>
<td align="left"><a href="GO:0006461" class="uri">GO:0006461</a></td>
<td align="left">protein complex assembly</td>
<td align="right">131</td>
<td align="right">48</td>
<td align="right">36.35</td>
<td align="left">0.01353</td>
</tr>
<tr class="even">
<td>160</td>
<td align="left"><a href="GO:0070271" class="uri">GO:0070271</a></td>
<td align="left">protein complex biogenesis</td>
<td align="right">131</td>
<td align="right">48</td>
<td align="right">36.35</td>
<td align="left">0.01353</td>
</tr>
<tr class="odd">
<td>161</td>
<td align="left"><a href="GO:0002695" class="uri">GO:0002695</a></td>
<td align="left">negative regulation of leukocyte activat...</td>
<td align="right">11</td>
<td align="right">7</td>
<td align="right">3.05</td>
<td align="left">0.01364</td>
</tr>
<tr class="even">
<td>162</td>
<td align="left"><a href="GO:0035385" class="uri">GO:0035385</a></td>
<td align="left">Roundabout signaling pathway</td>
<td align="right">11</td>
<td align="right">7</td>
<td align="right">3.05</td>
<td align="left">0.01364</td>
</tr>
<tr class="odd">
<td>163</td>
<td align="left"><a href="GO:0046683" class="uri">GO:0046683</a></td>
<td align="left">response to organophosphorus</td>
<td align="right">11</td>
<td align="right">7</td>
<td align="right">3.05</td>
<td align="left">0.01364</td>
</tr>
<tr class="even">
<td>164</td>
<td align="left"><a href="GO:0050866" class="uri">GO:0050866</a></td>
<td align="left">negative regulation of cell activation</td>
<td align="right">11</td>
<td align="right">7</td>
<td align="right">3.05</td>
<td align="left">0.01364</td>
</tr>
<tr class="odd">
<td>165</td>
<td align="left"><a href="GO:0090263" class="uri">GO:0090263</a></td>
<td align="left">positive regulation of canonical Wnt sig...</td>
<td align="right">11</td>
<td align="right">7</td>
<td align="right">3.05</td>
<td align="left">0.01364</td>
</tr>
<tr class="even">
<td>166</td>
<td align="left"><a href="GO:1901224" class="uri">GO:1901224</a></td>
<td align="left">positive regulation of NIK/NF-kappaB sig...</td>
<td align="right">11</td>
<td align="right">7</td>
<td align="right">3.05</td>
<td align="left">0.01364</td>
</tr>
<tr class="odd">
<td>167</td>
<td align="left"><a href="GO:0010720" class="uri">GO:0010720</a></td>
<td align="left">positive regulation of cell development</td>
<td align="right">40</td>
<td align="right">18</td>
<td align="right">11.10</td>
<td align="left">0.01368</td>
</tr>
<tr class="even">
<td>168</td>
<td align="left"><a href="GO:0006518" class="uri">GO:0006518</a></td>
<td align="left">peptide metabolic process</td>
<td align="right">119</td>
<td align="right">44</td>
<td align="right">33.02</td>
<td align="left">0.01495</td>
</tr>
<tr class="odd">
<td>169</td>
<td align="left"><a href="GO:0098742" class="uri">GO:0098742</a></td>
<td align="left">cell-cell adhesion via plasma-membrane a...</td>
<td align="right">52</td>
<td align="right">22</td>
<td align="right">14.43</td>
<td align="left">0.01560</td>
</tr>
<tr class="even">
<td>170</td>
<td align="left"><a href="GO:0006952" class="uri">GO:0006952</a></td>
<td align="left">defense response</td>
<td align="right">148</td>
<td align="right">53</td>
<td align="right">41.07</td>
<td align="left">0.01588</td>
</tr>
<tr class="odd">
<td>171</td>
<td align="left"><a href="GO:0019538" class="uri">GO:0019538</a></td>
<td align="left">protein metabolic process</td>
<td align="right">764</td>
<td align="right">233</td>
<td align="right">212.02</td>
<td align="left">0.01606</td>
</tr>
<tr class="even">
<td>172</td>
<td align="left"><a href="GO:0007423" class="uri">GO:0007423</a></td>
<td align="left">sensory organ development</td>
<td align="right">70</td>
<td align="right">28</td>
<td align="right">19.43</td>
<td align="left">0.01617</td>
</tr>
<tr class="odd">
<td>173</td>
<td align="left"><a href="GO:0048738" class="uri">GO:0048738</a></td>
<td align="left">cardiac muscle tissue development</td>
<td align="right">24</td>
<td align="right">12</td>
<td align="right">6.66</td>
<td align="left">0.01647</td>
</tr>
<tr class="even">
<td>174</td>
<td align="left"><a href="GO:0071407" class="uri">GO:0071407</a></td>
<td align="left">cellular response to organic cyclic comp...</td>
<td align="right">24</td>
<td align="right">12</td>
<td align="right">6.66</td>
<td align="left">0.01647</td>
</tr>
<tr class="odd">
<td>175</td>
<td align="left"><a href="GO:0048812" class="uri">GO:0048812</a></td>
<td align="left">neuron projection morphogenesis</td>
<td align="right">104</td>
<td align="right">39</td>
<td align="right">28.86</td>
<td align="left">0.01675</td>
</tr>
<tr class="even">
<td>176</td>
<td align="left"><a href="GO:0007214" class="uri">GO:0007214</a></td>
<td align="left">gamma-aminobutyric acid signaling pathwa...</td>
<td align="right">9</td>
<td align="right">6</td>
<td align="right">2.50</td>
<td align="left">0.01686</td>
</tr>
<tr class="odd">
<td>177</td>
<td align="left"><a href="GO:0019748" class="uri">GO:0019748</a></td>
<td align="left">secondary metabolic process</td>
<td align="right">9</td>
<td align="right">6</td>
<td align="right">2.50</td>
<td align="left">0.01686</td>
</tr>
<tr class="even">
<td>178</td>
<td align="left"><a href="GO:0022028" class="uri">GO:0022028</a></td>
<td align="left">tangential migration from the subventric...</td>
<td align="right">9</td>
<td align="right">6</td>
<td align="right">2.50</td>
<td align="left">0.01686</td>
</tr>
<tr class="odd">
<td>179</td>
<td align="left"><a href="GO:0032963" class="uri">GO:0032963</a></td>
<td align="left">collagen metabolic process</td>
<td align="right">9</td>
<td align="right">6</td>
<td align="right">2.50</td>
<td align="left">0.01686</td>
</tr>
<tr class="even">
<td>180</td>
<td align="left"><a href="GO:0044259" class="uri">GO:0044259</a></td>
<td align="left">multicellular organismal macromolecule m...</td>
<td align="right">9</td>
<td align="right">6</td>
<td align="right">2.50</td>
<td align="left">0.01686</td>
</tr>
<tr class="odd">
<td>181</td>
<td align="left"><a href="GO:0051964" class="uri">GO:0051964</a></td>
<td align="left">negative regulation of synapse assembly</td>
<td align="right">9</td>
<td align="right">6</td>
<td align="right">2.50</td>
<td align="left">0.01686</td>
</tr>
<tr class="even">
<td>182</td>
<td align="left"><a href="GO:0071621" class="uri">GO:0071621</a></td>
<td align="left">granulocyte chemotaxis</td>
<td align="right">9</td>
<td align="right">6</td>
<td align="right">2.50</td>
<td align="left">0.01686</td>
</tr>
<tr class="odd">
<td>183</td>
<td align="left"><a href="GO:1902115" class="uri">GO:1902115</a></td>
<td align="left">regulation of organelle assembly</td>
<td align="right">9</td>
<td align="right">6</td>
<td align="right">2.50</td>
<td align="left">0.01686</td>
</tr>
<tr class="even">
<td>184</td>
<td align="left"><a href="GO:0030334" class="uri">GO:0030334</a></td>
<td align="left">regulation of cell migration</td>
<td align="right">38</td>
<td align="right">17</td>
<td align="right">10.55</td>
<td align="left">0.01751</td>
</tr>
<tr class="odd">
<td>185</td>
<td align="left"><a href="GO:0007610" class="uri">GO:0007610</a></td>
<td align="left">behavior</td>
<td align="right">83</td>
<td align="right">32</td>
<td align="right">23.03</td>
<td align="left">0.01901</td>
</tr>
<tr class="even">
<td>186</td>
<td align="left"><a href="GO:1901701" class="uri">GO:1901701</a></td>
<td align="left">cellular response to oxygen-containing c...</td>
<td align="right">80</td>
<td align="right">31</td>
<td align="right">22.20</td>
<td align="left">0.01927</td>
</tr>
<tr class="odd">
<td>187</td>
<td align="left"><a href="GO:0007268" class="uri">GO:0007268</a></td>
<td align="left">synaptic transmission</td>
<td align="right">47</td>
<td align="right">20</td>
<td align="right">13.04</td>
<td align="left">0.01932</td>
</tr>
<tr class="even">
<td>188</td>
<td align="left"><a href="GO:0002764" class="uri">GO:0002764</a></td>
<td align="left">immune response-regulating signaling pat...</td>
<td align="right">50</td>
<td align="right">21</td>
<td align="right">13.88</td>
<td align="left">0.01967</td>
</tr>
<tr class="odd">
<td>189</td>
<td align="left"><a href="GO:0000904" class="uri">GO:0000904</a></td>
<td align="left">cell morphogenesis involved in different...</td>
<td align="right">124</td>
<td align="right">45</td>
<td align="right">34.41</td>
<td align="left">0.01992</td>
</tr>
<tr class="even">
<td>190</td>
<td align="left"><a href="GO:0002221" class="uri">GO:0002221</a></td>
<td align="left">pattern recognition receptor signaling p...</td>
<td align="right">30</td>
<td align="right">14</td>
<td align="right">8.33</td>
<td align="left">0.02002</td>
</tr>
<tr class="odd">
<td>191</td>
<td align="left"><a href="GO:0042088" class="uri">GO:0042088</a></td>
<td align="left">T-helper 1 type immune response</td>
<td align="right">7</td>
<td align="right">5</td>
<td align="right">1.94</td>
<td align="left">0.02028</td>
</tr>
<tr class="even">
<td>192</td>
<td align="left"><a href="GO:0050830" class="uri">GO:0050830</a></td>
<td align="left">defense response to Gram-positive bacter...</td>
<td align="right">7</td>
<td align="right">5</td>
<td align="right">1.94</td>
<td align="left">0.02028</td>
</tr>
<tr class="odd">
<td>193</td>
<td align="left"><a href="GO:0090102" class="uri">GO:0090102</a></td>
<td align="left">cochlea development</td>
<td align="right">7</td>
<td align="right">5</td>
<td align="right">1.94</td>
<td align="left">0.02028</td>
</tr>
<tr class="even">
<td>194</td>
<td align="left"><a href="GO:1903524" class="uri">GO:1903524</a></td>
<td align="left">positive regulation of blood circulation</td>
<td align="right">7</td>
<td align="right">5</td>
<td align="right">1.94</td>
<td align="left">0.02028</td>
</tr>
<tr class="odd">
<td>195</td>
<td align="left"><a href="GO:0000727" class="uri">GO:0000727</a></td>
<td align="left">double-strand break repair via break-ind...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="even">
<td>196</td>
<td align="left"><a href="GO:0002053" class="uri">GO:0002053</a></td>
<td align="left">positive regulation of mesenchymal cell ...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="odd">
<td>197</td>
<td align="left"><a href="GO:0002227" class="uri">GO:0002227</a></td>
<td align="left">innate immune response in mucosa</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="even">
<td>198</td>
<td align="left"><a href="GO:0002251" class="uri">GO:0002251</a></td>
<td align="left">organ or tissue specific immune response</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="odd">
<td>199</td>
<td align="left"><a href="GO:0002374" class="uri">GO:0002374</a></td>
<td align="left">cytokine secretion involved in immune re...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="even">
<td>200</td>
<td align="left"><a href="GO:0002385" class="uri">GO:0002385</a></td>
<td align="left">mucosal immune response</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="odd">
<td>201</td>
<td align="left"><a href="GO:0002534" class="uri">GO:0002534</a></td>
<td align="left">cytokine production involved in inflamma...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="even">
<td>202</td>
<td align="left"><a href="GO:0002548" class="uri">GO:0002548</a></td>
<td align="left">monocyte chemotaxis</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="odd">
<td>203</td>
<td align="left"><a href="GO:0002710" class="uri">GO:0002710</a></td>
<td align="left">negative regulation of T cell mediated i...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="even">
<td>204</td>
<td align="left"><a href="GO:0002714" class="uri">GO:0002714</a></td>
<td align="left">positive regulation of B cell mediated i...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="odd">
<td>205</td>
<td align="left"><a href="GO:0002815" class="uri">GO:0002815</a></td>
<td align="left">biosynthetic process of antibacterial pe...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="even">
<td>206</td>
<td align="left"><a href="GO:0002816" class="uri">GO:0002816</a></td>
<td align="left">regulation of biosynthetic process of an...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="odd">
<td>207</td>
<td align="left"><a href="GO:0002861" class="uri">GO:0002861</a></td>
<td align="left">regulation of inflammatory response to a...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="even">
<td>208</td>
<td align="left"><a href="GO:0002862" class="uri">GO:0002862</a></td>
<td align="left">negative regulation of inflammatory resp...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="odd">
<td>209</td>
<td align="left"><a href="GO:0002891" class="uri">GO:0002891</a></td>
<td align="left">positive regulation of immunoglobulin me...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="even">
<td>210</td>
<td align="left"><a href="GO:0002925" class="uri">GO:0002925</a></td>
<td align="left">positive regulation of humoral immune re...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="odd">
<td>211</td>
<td align="left"><a href="GO:0006582" class="uri">GO:0006582</a></td>
<td align="left">melanin metabolic process</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="even">
<td>212</td>
<td align="left"><a href="GO:0006965" class="uri">GO:0006965</a></td>
<td align="left">positive regulation of biosynthetic proc...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="odd">
<td>213</td>
<td align="left"><a href="GO:0010464" class="uri">GO:0010464</a></td>
<td align="left">regulation of mesenchymal cell prolifera...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="even">
<td>214</td>
<td align="left"><a href="GO:0010591" class="uri">GO:0010591</a></td>
<td align="left">regulation of lamellipodium assembly</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="odd">
<td>215</td>
<td align="left"><a href="GO:0010669" class="uri">GO:0010669</a></td>
<td align="left">epithelial structure maintenance</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="even">
<td>216</td>
<td align="left"><a href="GO:0010863" class="uri">GO:0010863</a></td>
<td align="left">positive regulation of phospholipase C a...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="odd">
<td>217</td>
<td align="left"><a href="GO:0019229" class="uri">GO:0019229</a></td>
<td align="left">regulation of vasoconstriction</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="even">
<td>218</td>
<td align="left"><a href="GO:0030277" class="uri">GO:0030277</a></td>
<td align="left">maintenance of gastrointestinal epitheli...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="odd">
<td>219</td>
<td align="left"><a href="GO:0032074" class="uri">GO:0032074</a></td>
<td align="left">negative regulation of nuclease activity</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="even">
<td>220</td>
<td align="left"><a href="GO:0032490" class="uri">GO:0032490</a></td>
<td align="left">detection of molecule of bacterial origi...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="odd">
<td>221</td>
<td align="left"><a href="GO:0032494" class="uri">GO:0032494</a></td>
<td align="left">response to peptidoglycan</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="even">
<td>222</td>
<td align="left"><a href="GO:0032498" class="uri">GO:0032498</a></td>
<td align="left">detection of muramyl dipeptide</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="odd">
<td>223</td>
<td align="left"><a href="GO:0032499" class="uri">GO:0032499</a></td>
<td align="left">detection of peptidoglycan</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="even">
<td>224</td>
<td align="left"><a href="GO:0032602" class="uri">GO:0032602</a></td>
<td align="left">chemokine production</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="odd">
<td>225</td>
<td align="left"><a href="GO:0032620" class="uri">GO:0032620</a></td>
<td align="left">interleukin-17 production</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="even">
<td>226</td>
<td align="left"><a href="GO:0032660" class="uri">GO:0032660</a></td>
<td align="left">regulation of interleukin-17 production</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="odd">
<td>227</td>
<td align="left"><a href="GO:0032701" class="uri">GO:0032701</a></td>
<td align="left">negative regulation of interleukin-18 pr...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="even">
<td>228</td>
<td align="left"><a href="GO:0032703" class="uri">GO:0032703</a></td>
<td align="left">negative regulation of interleukin-2 pro...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="odd">
<td>229</td>
<td align="left"><a href="GO:0032740" class="uri">GO:0032740</a></td>
<td align="left">positive regulation of interleukin-17 pr...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="even">
<td>230</td>
<td align="left"><a href="GO:0033198" class="uri">GO:0033198</a></td>
<td align="left">response to ATP</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="odd">
<td>231</td>
<td align="left"><a href="GO:0042310" class="uri">GO:0042310</a></td>
<td align="left">vasoconstriction</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="even">
<td>232</td>
<td align="left"><a href="GO:0042438" class="uri">GO:0042438</a></td>
<td align="left">melanin biosynthetic process</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="odd">
<td>233</td>
<td align="left"><a href="GO:0045168" class="uri">GO:0045168</a></td>
<td align="left">cell-cell signaling involved in cell fat...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="even">
<td>234</td>
<td align="left"><a href="GO:0045581" class="uri">GO:0045581</a></td>
<td align="left">negative regulation of T cell differenti...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="odd">
<td>235</td>
<td align="left"><a href="GO:0045620" class="uri">GO:0045620</a></td>
<td align="left">negative regulation of lymphocyte differ...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="even">
<td>236</td>
<td align="left"><a href="GO:0045907" class="uri">GO:0045907</a></td>
<td align="left">positive regulation of vasoconstriction</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="odd">
<td>237</td>
<td align="left"><a href="GO:0046629" class="uri">GO:0046629</a></td>
<td align="left">gamma-delta T cell activation</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="even">
<td>238</td>
<td align="left"><a href="GO:0046643" class="uri">GO:0046643</a></td>
<td align="left">regulation of gamma-delta T cell activat...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="odd">
<td>239</td>
<td align="left"><a href="GO:0046645" class="uri">GO:0046645</a></td>
<td align="left">positive regulation of gamma-delta T cel...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="even">
<td>240</td>
<td align="left"><a href="GO:0048339" class="uri">GO:0048339</a></td>
<td align="left">paraxial mesoderm development</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="odd">
<td>241</td>
<td align="left"><a href="GO:0051930" class="uri">GO:0051930</a></td>
<td align="left">regulation of sensory perception of pain</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="even">
<td>242</td>
<td align="left"><a href="GO:0051931" class="uri">GO:0051931</a></td>
<td align="left">regulation of sensory perception</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="odd">
<td>243</td>
<td align="left"><a href="GO:0061430" class="uri">GO:0061430</a></td>
<td align="left">bone trabecula morphogenesis</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="even">
<td>244</td>
<td align="left"><a href="GO:0070431" class="uri">GO:0070431</a></td>
<td align="left">nucleotide-binding oligomerization domai...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="odd">
<td>245</td>
<td align="left"><a href="GO:0071224" class="uri">GO:0071224</a></td>
<td align="left">cellular response to peptidoglycan</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="even">
<td>246</td>
<td align="left"><a href="GO:0071608" class="uri">GO:0071608</a></td>
<td align="left">macrophage inflammatory protein-1 alpha ...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="odd">
<td>247</td>
<td align="left"><a href="GO:0090103" class="uri">GO:0090103</a></td>
<td align="left">cochlea morphogenesis</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="even">
<td>248</td>
<td align="left"><a href="GO:0098792" class="uri">GO:0098792</a></td>
<td align="left">xenophagy</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="odd">
<td>249</td>
<td align="left"><a href="GO:1900015" class="uri">GO:1900015</a></td>
<td align="left">regulation of cytokine production involv...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="even">
<td>250</td>
<td align="left"><a href="GO:1900017" class="uri">GO:1900017</a></td>
<td align="left">positive regulation of cytokine producti...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="odd">
<td>251</td>
<td align="left"><a href="GO:1900274" class="uri">GO:1900274</a></td>
<td align="left">regulation of phospholipase C activity</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="even">
<td>252</td>
<td align="left"><a href="GO:1902017" class="uri">GO:1902017</a></td>
<td align="left">regulation of cilium assembly</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="odd">
<td>253</td>
<td align="left"><a href="GO:1902743" class="uri">GO:1902743</a></td>
<td align="left">regulation of lamellipodium organization</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.83</td>
<td align="left">0.02128</td>
</tr>
<tr class="even">
<td>254</td>
<td align="left"><a href="GO:0051963" class="uri">GO:0051963</a></td>
<td align="left">regulation of synapse assembly</td>
<td align="right">22</td>
<td align="right">11</td>
<td align="right">6.11</td>
<td align="left">0.02145</td>
</tr>
<tr class="odd">
<td>255</td>
<td align="left"><a href="GO:2000021" class="uri">GO:2000021</a></td>
<td align="left">regulation of ion homeostasis</td>
<td align="right">22</td>
<td align="right">11</td>
<td align="right">6.11</td>
<td align="left">0.02145</td>
</tr>
<tr class="even">
<td>256</td>
<td align="left"><a href="GO:0050896" class="uri">GO:0050896</a></td>
<td align="left">response to stimulus</td>
<td align="right">729</td>
<td align="right">222</td>
<td align="right">202.31</td>
<td align="left">0.02149</td>
</tr>
<tr class="odd">
<td>257</td>
<td align="left"><a href="GO:0046903" class="uri">GO:0046903</a></td>
<td align="left">secretion</td>
<td align="right">115</td>
<td align="right">42</td>
<td align="right">31.91</td>
<td align="left">0.02157</td>
</tr>
<tr class="even">
<td>258</td>
<td align="left"><a href="GO:0043603" class="uri">GO:0043603</a></td>
<td align="left">cellular amide metabolic process</td>
<td align="right">131</td>
<td align="right">47</td>
<td align="right">36.35</td>
<td align="left">0.02183</td>
</tr>
<tr class="odd">
<td>259</td>
<td align="left"><a href="GO:0001933" class="uri">GO:0001933</a></td>
<td align="left">negative regulation of protein phosphory...</td>
<td align="right">36</td>
<td align="right">16</td>
<td align="right">9.99</td>
<td align="left">0.02243</td>
</tr>
<tr class="even">
<td>260</td>
<td align="left"><a href="GO:0002275" class="uri">GO:0002275</a></td>
<td align="left">myeloid cell activation involved in immu...</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.39</td>
<td align="left">0.02292</td>
</tr>
<tr class="odd">
<td>261</td>
<td align="left"><a href="GO:0002279" class="uri">GO:0002279</a></td>
<td align="left">mast cell activation involved in immune ...</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.39</td>
<td align="left">0.02292</td>
</tr>
<tr class="even">
<td>262</td>
<td align="left"><a href="GO:0002444" class="uri">GO:0002444</a></td>
<td align="left">myeloid leukocyte mediated immunity</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.39</td>
<td align="left">0.02292</td>
</tr>
<tr class="odd">
<td>263</td>
<td align="left"><a href="GO:0002448" class="uri">GO:0002448</a></td>
<td align="left">mast cell mediated immunity</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.39</td>
<td align="left">0.02292</td>
</tr>
<tr class="even">
<td>264</td>
<td align="left"><a href="GO:0002579" class="uri">GO:0002579</a></td>
<td align="left">positive regulation of antigen processin...</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.39</td>
<td align="left">0.02292</td>
</tr>
<tr class="odd">
<td>265</td>
<td align="left"><a href="GO:0002698" class="uri">GO:0002698</a></td>
<td align="left">negative regulation of immune effector p...</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.39</td>
<td align="left">0.02292</td>
</tr>
<tr class="even">
<td>266</td>
<td align="left"><a href="GO:0002886" class="uri">GO:0002886</a></td>
<td align="left">regulation of myeloid leukocyte mediated...</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.39</td>
<td align="left">0.02292</td>
</tr>
<tr class="odd">
<td>267</td>
<td align="left"><a href="GO:0007586" class="uri">GO:0007586</a></td>
<td align="left">digestion</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.39</td>
<td align="left">0.02292</td>
</tr>
<tr class="even">
<td>268</td>
<td align="left"><a href="GO:0010042" class="uri">GO:0010042</a></td>
<td align="left">response to manganese ion</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.39</td>
<td align="left">0.02292</td>
</tr>
<tr class="odd">
<td>269</td>
<td align="left"><a href="GO:0010517" class="uri">GO:0010517</a></td>
<td align="left">regulation of phospholipase activity</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.39</td>
<td align="left">0.02292</td>
</tr>
<tr class="even">
<td>270</td>
<td align="left"><a href="GO:0010518" class="uri">GO:0010518</a></td>
<td align="left">positive regulation of phospholipase act...</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.39</td>
<td align="left">0.02292</td>
</tr>
<tr class="odd">
<td>271</td>
<td align="left"><a href="GO:0014003" class="uri">GO:0014003</a></td>
<td align="left">oligodendrocyte development</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.39</td>
<td align="left">0.02292</td>
</tr>
<tr class="even">
<td>272</td>
<td align="left"><a href="GO:0014041" class="uri">GO:0014041</a></td>
<td align="left">regulation of neuron maturation</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.39</td>
<td align="left">0.02292</td>
</tr>
<tr class="odd">
<td>273</td>
<td align="left"><a href="GO:0030574" class="uri">GO:0030574</a></td>
<td align="left">collagen catabolic process</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.39</td>
<td align="left">0.02292</td>
</tr>
<tr class="even">
<td>274</td>
<td align="left"><a href="GO:0031000" class="uri">GO:0031000</a></td>
<td align="left">response to caffeine</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.39</td>
<td align="left">0.02292</td>
</tr>
<tr class="odd">
<td>275</td>
<td align="left"><a href="GO:0032418" class="uri">GO:0032418</a></td>
<td align="left">lysosome localization</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.39</td>
<td align="left">0.02292</td>
</tr>
<tr class="even">
<td>276</td>
<td align="left"><a href="GO:0032735" class="uri">GO:0032735</a></td>
<td align="left">positive regulation of interleukin-12 pr...</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.39</td>
<td align="left">0.02292</td>
</tr>
<tr class="odd">
<td>277</td>
<td align="left"><a href="GO:0033003" class="uri">GO:0033003</a></td>
<td align="left">regulation of mast cell activation</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.39</td>
<td align="left">0.02292</td>
</tr>
<tr class="even">
<td>278</td>
<td align="left"><a href="GO:0033006" class="uri">GO:0033006</a></td>
<td align="left">regulation of mast cell activation invol...</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.39</td>
<td align="left">0.02292</td>
</tr>
<tr class="odd">
<td>279</td>
<td align="left"><a href="GO:0035640" class="uri">GO:0035640</a></td>
<td align="left">exploration behavior</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.39</td>
<td align="left">0.02292</td>
</tr>
<tr class="even">
<td>280</td>
<td align="left"><a href="GO:0035641" class="uri">GO:0035641</a></td>
<td align="left">locomotory exploration behavior</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.39</td>
<td align="left">0.02292</td>
</tr>
<tr class="odd">
<td>281</td>
<td align="left"><a href="GO:0042136" class="uri">GO:0042136</a></td>
<td align="left">neurotransmitter biosynthetic process</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.39</td>
<td align="left">0.02292</td>
</tr>
<tr class="even">
<td>282</td>
<td align="left"><a href="GO:0042509" class="uri">GO:0042509</a></td>
<td align="left">regulation of tyrosine phosphorylation o...</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.39</td>
<td align="left">0.02292</td>
</tr>
<tr class="odd">
<td>283</td>
<td align="left"><a href="GO:0042531" class="uri">GO:0042531</a></td>
<td align="left">positive regulation of tyrosine phosphor...</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.39</td>
<td align="left">0.02292</td>
</tr>
<tr class="even">
<td>284</td>
<td align="left"><a href="GO:0043299" class="uri">GO:0043299</a></td>
<td align="left">leukocyte degranulation</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.39</td>
<td align="left">0.02292</td>
</tr>
<tr class="odd">
<td>285</td>
<td align="left"><a href="GO:0043300" class="uri">GO:0043300</a></td>
<td align="left">regulation of leukocyte degranulation</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.39</td>
<td align="left">0.02292</td>
</tr>
<tr class="even">
<td>286</td>
<td align="left"><a href="GO:0043303" class="uri">GO:0043303</a></td>
<td align="left">mast cell degranulation</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.39</td>
<td align="left">0.02292</td>
</tr>
<tr class="odd">
<td>287</td>
<td align="left"><a href="GO:0043304" class="uri">GO:0043304</a></td>
<td align="left">regulation of mast cell degranulation</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.39</td>
<td align="left">0.02292</td>
</tr>
<tr class="even">
<td>288</td>
<td align="left"><a href="GO:0043331" class="uri">GO:0043331</a></td>
<td align="left">response to dsRNA</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.39</td>
<td align="left">0.02292</td>
</tr>
<tr class="odd">
<td>289</td>
<td align="left"><a href="GO:0045055" class="uri">GO:0045055</a></td>
<td align="left">regulated secretory pathway</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.39</td>
<td align="left">0.02292</td>
</tr>
<tr class="even">
<td>290</td>
<td align="left"><a href="GO:0045576" class="uri">GO:0045576</a></td>
<td align="left">mast cell activation</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.39</td>
<td align="left">0.02292</td>
</tr>
<tr class="odd">
<td>291</td>
<td align="left"><a href="GO:0046189" class="uri">GO:0046189</a></td>
<td align="left">phenol-containing compound biosynthetic ...</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.39</td>
<td align="left">0.02292</td>
</tr>
<tr class="even">
<td>292</td>
<td align="left"><a href="GO:0051882" class="uri">GO:0051882</a></td>
<td align="left">mitochondrial depolarization</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.39</td>
<td align="left">0.02292</td>
</tr>
<tr class="odd">
<td>293</td>
<td align="left"><a href="GO:0051900" class="uri">GO:0051900</a></td>
<td align="left">regulation of mitochondrial depolarizati...</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.39</td>
<td align="left">0.02292</td>
</tr>
<tr class="even">
<td>294</td>
<td align="left"><a href="GO:0060191" class="uri">GO:0060191</a></td>
<td align="left">regulation of lipase activity</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.39</td>
<td align="left">0.02292</td>
</tr>
<tr class="odd">
<td>295</td>
<td align="left"><a href="GO:0060193" class="uri">GO:0060193</a></td>
<td align="left">positive regulation of lipase activity</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.39</td>
<td align="left">0.02292</td>
</tr>
<tr class="even">
<td>296</td>
<td align="left"><a href="GO:0071287" class="uri">GO:0071287</a></td>
<td align="left">cellular response to manganese ion</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.39</td>
<td align="left">0.02292</td>
</tr>
<tr class="odd">
<td>297</td>
<td align="left"><a href="GO:0071312" class="uri">GO:0071312</a></td>
<td align="left">cellular response to alkaloid</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.39</td>
<td align="left">0.02292</td>
</tr>
<tr class="even">
<td>298</td>
<td align="left"><a href="GO:0071622" class="uri">GO:0071622</a></td>
<td align="left">regulation of granulocyte chemotaxis</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.39</td>
<td align="left">0.02292</td>
</tr>
<tr class="odd">
<td>299</td>
<td align="left"><a href="GO:0071867" class="uri">GO:0071867</a></td>
<td align="left">response to monoamine</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.39</td>
<td align="left">0.02292</td>
</tr>
<tr class="even">
<td>300</td>
<td align="left"><a href="GO:0071868" class="uri">GO:0071868</a></td>
<td align="left">cellular response to monoamine stimulus</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.39</td>
<td align="left">0.02292</td>
</tr>
<tr class="odd">
<td>301</td>
<td align="left"><a href="GO:0071869" class="uri">GO:0071869</a></td>
<td align="left">response to catecholamine</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.39</td>
<td align="left">0.02292</td>
</tr>
<tr class="even">
<td>302</td>
<td align="left"><a href="GO:0071870" class="uri">GO:0071870</a></td>
<td align="left">cellular response to catecholamine stimu...</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.39</td>
<td align="left">0.02292</td>
</tr>
<tr class="odd">
<td>303</td>
<td align="left"><a href="GO:1903305" class="uri">GO:1903305</a></td>
<td align="left">regulation of regulated secretory pathwa...</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.39</td>
<td align="left">0.02292</td>
</tr>
<tr class="even">
<td>304</td>
<td align="left"><a href="GO:1903429" class="uri">GO:1903429</a></td>
<td align="left">regulation of cell maturation</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.39</td>
<td align="left">0.02292</td>
</tr>
<tr class="odd">
<td>305</td>
<td align="left"><a href="GO:0040011" class="uri">GO:0040011</a></td>
<td align="left">locomotion</td>
<td align="right">203</td>
<td align="right">69</td>
<td align="right">56.33</td>
<td align="left">0.02314</td>
</tr>
<tr class="even">
<td>306</td>
<td align="left"><a href="GO:0009605" class="uri">GO:0009605</a></td>
<td align="left">response to external stimulus</td>
<td align="right">283</td>
<td align="right">93</td>
<td align="right">78.54</td>
<td align="left">0.02325</td>
</tr>
<tr class="odd">
<td>307</td>
<td align="left"><a href="GO:2000145" class="uri">GO:2000145</a></td>
<td align="left">regulation of cell motility</td>
<td align="right">39</td>
<td align="right">17</td>
<td align="right">10.82</td>
<td align="left">0.02328</td>
</tr>
<tr class="even">
<td>308</td>
<td align="left"><a href="GO:0061351" class="uri">GO:0061351</a></td>
<td align="left">neural precursor cell proliferation</td>
<td align="right">25</td>
<td align="right">12</td>
<td align="right">6.94</td>
<td align="left">0.02390</td>
</tr>
<tr class="odd">
<td>309</td>
<td align="left"><a href="GO:0019882" class="uri">GO:0019882</a></td>
<td align="left">antigen processing and presentation</td>
<td align="right">17</td>
<td align="right">9</td>
<td align="right">4.72</td>
<td align="left">0.02411</td>
</tr>
<tr class="even">
<td>310</td>
<td align="left"><a href="GO:1901222" class="uri">GO:1901222</a></td>
<td align="left">regulation of NIK/NF-kappaB signaling</td>
<td align="right">17</td>
<td align="right">9</td>
<td align="right">4.72</td>
<td align="left">0.02411</td>
</tr>
<tr class="odd">
<td>311</td>
<td align="left"><a href="GO:0051270" class="uri">GO:0051270</a></td>
<td align="left">regulation of cellular component movemen...</td>
<td align="right">48</td>
<td align="right">20</td>
<td align="right">13.32</td>
<td align="left">0.02478</td>
</tr>
<tr class="even">
<td>312</td>
<td align="left"><a href="GO:0007040" class="uri">GO:0007040</a></td>
<td align="left">lysosome organization</td>
<td align="right">12</td>
<td align="right">7</td>
<td align="right">3.33</td>
<td align="left">0.02503</td>
</tr>
<tr class="odd">
<td>313</td>
<td align="left"><a href="GO:0016199" class="uri">GO:0016199</a></td>
<td align="left">axon midline choice point recognition</td>
<td align="right">12</td>
<td align="right">7</td>
<td align="right">3.33</td>
<td align="left">0.02503</td>
</tr>
<tr class="even">
<td>314</td>
<td align="left"><a href="GO:0030177" class="uri">GO:0030177</a></td>
<td align="left">positive regulation of Wnt signaling pat...</td>
<td align="right">12</td>
<td align="right">7</td>
<td align="right">3.33</td>
<td align="left">0.02503</td>
</tr>
<tr class="odd">
<td>315</td>
<td align="left"><a href="GO:0080171" class="uri">GO:0080171</a></td>
<td align="left">lytic vacuole organization</td>
<td align="right">12</td>
<td align="right">7</td>
<td align="right">3.33</td>
<td align="left">0.02503</td>
</tr>
<tr class="even">
<td>316</td>
<td align="left"><a href="GO:0050663" class="uri">GO:0050663</a></td>
<td align="left">cytokine secretion</td>
<td align="right">28</td>
<td align="right">13</td>
<td align="right">7.77</td>
<td align="left">0.02589</td>
</tr>
<tr class="odd">
<td>317</td>
<td align="left"><a href="GO:2000026" class="uri">GO:2000026</a></td>
<td align="left">regulation of multicellular organismal d...</td>
<td align="right">181</td>
<td align="right">62</td>
<td align="right">50.23</td>
<td align="left">0.02609</td>
</tr>
<tr class="even">
<td>318</td>
<td align="left"><a href="GO:0050769" class="uri">GO:0050769</a></td>
<td align="left">positive regulation of neurogenesis</td>
<td align="right">31</td>
<td align="right">14</td>
<td align="right">8.60</td>
<td align="left">0.02747</td>
</tr>
<tr class="odd">
<td>319</td>
<td align="left"><a href="GO:0051865" class="uri">GO:0051865</a></td>
<td align="left">protein autoubiquitination</td>
<td align="right">31</td>
<td align="right">14</td>
<td align="right">8.60</td>
<td align="left">0.02747</td>
</tr>
<tr class="even">
<td>320</td>
<td align="left"><a href="GO:0001894" class="uri">GO:0001894</a></td>
<td align="left">tissue homeostasis</td>
<td align="right">20</td>
<td align="right">10</td>
<td align="right">5.55</td>
<td align="left">0.02799</td>
</tr>
<tr class="odd">
<td>321</td>
<td align="left"><a href="GO:0032869" class="uri">GO:0032869</a></td>
<td align="left">cellular response to insulin stimulus</td>
<td align="right">20</td>
<td align="right">10</td>
<td align="right">5.55</td>
<td align="left">0.02799</td>
</tr>
<tr class="even">
<td>322</td>
<td align="left"><a href="GO:0002250" class="uri">GO:0002250</a></td>
<td align="left">adaptive immune response</td>
<td align="right">34</td>
<td align="right">15</td>
<td align="right">9.44</td>
<td align="left">0.02872</td>
</tr>
<tr class="odd">
<td>323</td>
<td align="left"><a href="GO:0002376" class="uri">GO:0002376</a></td>
<td align="left">immune system process</td>
<td align="right">195</td>
<td align="right">66</td>
<td align="right">54.11</td>
<td align="left">0.02874</td>
</tr>
<tr class="even">
<td>324</td>
<td align="left"><a href="GO:0060627" class="uri">GO:0060627</a></td>
<td align="left">regulation of vesicle-mediated transport</td>
<td align="right">37</td>
<td align="right">16</td>
<td align="right">10.27</td>
<td align="left">0.02967</td>
</tr>
<tr class="odd">
<td>325</td>
<td align="left"><a href="GO:0022008" class="uri">GO:0022008</a></td>
<td align="left">neurogenesis</td>
<td align="right">192</td>
<td align="right">65</td>
<td align="right">53.28</td>
<td align="left">0.02974</td>
</tr>
<tr class="even">
<td>326</td>
<td align="left"><a href="GO:0048585" class="uri">GO:0048585</a></td>
<td align="left">negative regulation of response to stimu...</td>
<td align="right">143</td>
<td align="right">50</td>
<td align="right">39.68</td>
<td align="left">0.03010</td>
</tr>
<tr class="odd">
<td>327</td>
<td align="left"><a href="GO:0021537" class="uri">GO:0021537</a></td>
<td align="left">telencephalon development</td>
<td align="right">40</td>
<td align="right">17</td>
<td align="right">11.10</td>
<td align="left">0.03038</td>
</tr>
<tr class="even">
<td>328</td>
<td align="left"><a href="GO:0048871" class="uri">GO:0048871</a></td>
<td align="left">multicellular organismal homeostasis</td>
<td align="right">23</td>
<td align="right">11</td>
<td align="right">6.38</td>
<td align="left">0.03109</td>
</tr>
<tr class="odd">
<td>329</td>
<td align="left"><a href="GO:0046777" class="uri">GO:0046777</a></td>
<td align="left">protein autophosphorylation</td>
<td align="right">46</td>
<td align="right">19</td>
<td align="right">12.77</td>
<td align="left">0.03118</td>
</tr>
<tr class="even">
<td>330</td>
<td align="left"><a href="GO:0033674" class="uri">GO:0033674</a></td>
<td align="left">positive regulation of kinase activity</td>
<td align="right">52</td>
<td align="right">21</td>
<td align="right">14.43</td>
<td align="left">0.03138</td>
</tr>
<tr class="odd">
<td>331</td>
<td align="left"><a href="GO:0001649" class="uri">GO:0001649</a></td>
<td align="left">osteoblast differentiation</td>
<td align="right">15</td>
<td align="right">8</td>
<td align="right">4.16</td>
<td align="left">0.03157</td>
</tr>
<tr class="even">
<td>332</td>
<td align="left"><a href="GO:0038095" class="uri">GO:0038095</a></td>
<td align="left">Fc-epsilon receptor signaling pathway</td>
<td align="right">15</td>
<td align="right">8</td>
<td align="right">4.16</td>
<td align="left">0.03157</td>
</tr>
<tr class="odd">
<td>333</td>
<td align="left"><a href="GO:0050921" class="uri">GO:0050921</a></td>
<td align="left">positive regulation of chemotaxis</td>
<td align="right">15</td>
<td align="right">8</td>
<td align="right">4.16</td>
<td align="left">0.03157</td>
</tr>
<tr class="even">
<td>334</td>
<td align="left"><a href="GO:0016578" class="uri">GO:0016578</a></td>
<td align="left">histone deubiquitination</td>
<td align="right">10</td>
<td align="right">6</td>
<td align="right">2.78</td>
<td align="left">0.03240</td>
</tr>
<tr class="odd">
<td>335</td>
<td align="left"><a href="GO:0044236" class="uri">GO:0044236</a></td>
<td align="left">multicellular organismal metabolic proce...</td>
<td align="right">10</td>
<td align="right">6</td>
<td align="right">2.78</td>
<td align="left">0.03240</td>
</tr>
<tr class="even">
<td>336</td>
<td align="left"><a href="GO:0050777" class="uri">GO:0050777</a></td>
<td align="left">negative regulation of immune response</td>
<td align="right">10</td>
<td align="right">6</td>
<td align="right">2.78</td>
<td align="left">0.03240</td>
</tr>
<tr class="odd">
<td>337</td>
<td align="left"><a href="GO:0050829" class="uri">GO:0050829</a></td>
<td align="left">defense response to Gram-negative bacter...</td>
<td align="right">10</td>
<td align="right">6</td>
<td align="right">2.78</td>
<td align="left">0.03240</td>
</tr>
<tr class="even">
<td>338</td>
<td align="left"><a href="GO:0051250" class="uri">GO:0051250</a></td>
<td align="left">negative regulation of lymphocyte activa...</td>
<td align="right">10</td>
<td align="right">6</td>
<td align="right">2.78</td>
<td align="left">0.03240</td>
</tr>
<tr class="odd">
<td>339</td>
<td align="left"><a href="GO:0060688" class="uri">GO:0060688</a></td>
<td align="left">regulation of morphogenesis of a branchi...</td>
<td align="right">10</td>
<td align="right">6</td>
<td align="right">2.78</td>
<td align="left">0.03240</td>
</tr>
<tr class="even">
<td>340</td>
<td align="left"><a href="GO:0097530" class="uri">GO:0097530</a></td>
<td align="left">granulocyte migration</td>
<td align="right">10</td>
<td align="right">6</td>
<td align="right">2.78</td>
<td align="left">0.03240</td>
</tr>
<tr class="odd">
<td>341</td>
<td align="left"><a href="GO:0007409" class="uri">GO:0007409</a></td>
<td align="left">axonogenesis</td>
<td align="right">83</td>
<td align="right">31</td>
<td align="right">23.03</td>
<td align="left">0.03310</td>
</tr>
<tr class="even">
<td>342</td>
<td align="left"><a href="GO:0000187" class="uri">GO:0000187</a></td>
<td align="left">activation of MAPK activity</td>
<td align="right">18</td>
<td align="right">9</td>
<td align="right">5.00</td>
<td align="left">0.03666</td>
</tr>
<tr class="odd">
<td>343</td>
<td align="left"><a href="GO:0008286" class="uri">GO:0008286</a></td>
<td align="left">insulin receptor signaling pathway</td>
<td align="right">18</td>
<td align="right">9</td>
<td align="right">5.00</td>
<td align="left">0.03666</td>
</tr>
<tr class="even">
<td>344</td>
<td align="left"><a href="GO:0002218" class="uri">GO:0002218</a></td>
<td align="left">activation of innate immune response</td>
<td align="right">32</td>
<td align="right">14</td>
<td align="right">8.88</td>
<td align="left">0.03679</td>
</tr>
<tr class="odd">
<td>345</td>
<td align="left"><a href="GO:0002758" class="uri">GO:0002758</a></td>
<td align="left">innate immune response-activating signal...</td>
<td align="right">32</td>
<td align="right">14</td>
<td align="right">8.88</td>
<td align="left">0.03679</td>
</tr>
<tr class="even">
<td>346</td>
<td align="left"><a href="GO:0048813" class="uri">GO:0048813</a></td>
<td align="left">dendrite morphogenesis</td>
<td align="right">32</td>
<td align="right">14</td>
<td align="right">8.88</td>
<td align="left">0.03679</td>
</tr>
<tr class="odd">
<td>347</td>
<td align="left"><a href="GO:0051235" class="uri">GO:0051235</a></td>
<td align="left">maintenance of location</td>
<td align="right">32</td>
<td align="right">14</td>
<td align="right">8.88</td>
<td align="left">0.03679</td>
</tr>
<tr class="even">
<td>348</td>
<td align="left"><a href="GO:0051272" class="uri">GO:0051272</a></td>
<td align="left">positive regulation of cellular componen...</td>
<td align="right">32</td>
<td align="right">14</td>
<td align="right">8.88</td>
<td align="left">0.03679</td>
</tr>
<tr class="odd">
<td>349</td>
<td align="left"><a href="GO:0010033" class="uri">GO:0010033</a></td>
<td align="left">response to organic substance</td>
<td align="right">224</td>
<td align="right">74</td>
<td align="right">62.16</td>
<td align="left">0.03718</td>
</tr>
<tr class="even">
<td>350</td>
<td align="left"><a href="GO:0031175" class="uri">GO:0031175</a></td>
<td align="left">neuron projection development</td>
<td align="right">122</td>
<td align="right">43</td>
<td align="right">33.86</td>
<td align="left">0.03745</td>
</tr>
<tr class="odd">
<td>351</td>
<td align="left"><a href="GO:0044707" class="uri">GO:0044707</a></td>
<td align="left">single-multicellular organism process</td>
<td align="right">642</td>
<td align="right">195</td>
<td align="right">178.16</td>
<td align="left">0.03843</td>
</tr>
<tr class="even">
<td>352</td>
<td align="left"><a href="GO:0042326" class="uri">GO:0042326</a></td>
<td align="left">negative regulation of phosphorylation</td>
<td align="right">38</td>
<td align="right">16</td>
<td align="right">10.55</td>
<td align="left">0.03852</td>
</tr>
<tr class="odd">
<td>353</td>
<td align="left"><a href="GO:0051961" class="uri">GO:0051961</a></td>
<td align="left">negative regulation of nervous system de...</td>
<td align="right">38</td>
<td align="right">16</td>
<td align="right">10.55</td>
<td align="left">0.03852</td>
</tr>
<tr class="even">
<td>354</td>
<td align="left"><a href="GO:0006954" class="uri">GO:0006954</a></td>
<td align="left">inflammatory response</td>
<td align="right">56</td>
<td align="right">22</td>
<td align="right">15.54</td>
<td align="left">0.03857</td>
</tr>
<tr class="odd">
<td>355</td>
<td align="left"><a href="GO:0044265" class="uri">GO:0044265</a></td>
<td align="left">cellular macromolecule catabolic process</td>
<td align="right">103</td>
<td align="right">37</td>
<td align="right">28.58</td>
<td align="left">0.03888</td>
</tr>
<tr class="even">
<td>356</td>
<td align="left"><a href="GO:0002757" class="uri">GO:0002757</a></td>
<td align="left">immune response-activating signal transd...</td>
<td align="right">41</td>
<td align="right">17</td>
<td align="right">11.38</td>
<td align="left">0.03896</td>
</tr>
<tr class="odd">
<td>357</td>
<td align="left"><a href="GO:0001525" class="uri">GO:0001525</a></td>
<td align="left">angiogenesis</td>
<td align="right">50</td>
<td align="right">20</td>
<td align="right">13.88</td>
<td align="left">0.03913</td>
</tr>
<tr class="even">
<td>358</td>
<td align="left"><a href="GO:0007156" class="uri">GO:0007156</a></td>
<td align="left">homophilic cell adhesion via plasma memb...</td>
<td align="right">47</td>
<td align="right">19</td>
<td align="right">13.04</td>
<td align="left">0.03923</td>
</tr>
<tr class="odd">
<td>359</td>
<td align="left"><a href="GO:0001505" class="uri">GO:0001505</a></td>
<td align="left">regulation of neurotransmitter levels</td>
<td align="right">21</td>
<td align="right">10</td>
<td align="right">5.83</td>
<td align="left">0.04054</td>
</tr>
<tr class="even">
<td>360</td>
<td align="left"><a href="GO:0051262" class="uri">GO:0051262</a></td>
<td align="left">protein tetramerization</td>
<td align="right">21</td>
<td align="right">10</td>
<td align="right">5.83</td>
<td align="left">0.04054</td>
</tr>
<tr class="odd">
<td>361</td>
<td align="left"><a href="GO:0048699" class="uri">GO:0048699</a></td>
<td align="left">generation of neurons</td>
<td align="right">185</td>
<td align="right">62</td>
<td align="right">51.34</td>
<td align="left">0.04111</td>
</tr>
<tr class="even">
<td>362</td>
<td align="left"><a href="GO:0050793" class="uri">GO:0050793</a></td>
<td align="left">regulation of developmental process</td>
<td align="right">215</td>
<td align="right">71</td>
<td align="right">59.67</td>
<td align="left">0.04125</td>
</tr>
<tr class="odd">
<td>363</td>
<td align="left"><a href="GO:0016198" class="uri">GO:0016198</a></td>
<td align="left">axon choice point recognition</td>
<td align="right">13</td>
<td align="right">7</td>
<td align="right">3.61</td>
<td align="left">0.04152</td>
</tr>
<tr class="even">
<td>364</td>
<td align="left"><a href="GO:0019730" class="uri">GO:0019730</a></td>
<td align="left">antimicrobial humoral response</td>
<td align="right">13</td>
<td align="right">7</td>
<td align="right">3.61</td>
<td align="left">0.04152</td>
</tr>
<tr class="odd">
<td>365</td>
<td align="left"><a href="GO:0019731" class="uri">GO:0019731</a></td>
<td align="left">antibacterial humoral response</td>
<td align="right">13</td>
<td align="right">7</td>
<td align="right">3.61</td>
<td align="left">0.04152</td>
</tr>
<tr class="even">
<td>366</td>
<td align="left"><a href="GO:0022029" class="uri">GO:0022029</a></td>
<td align="left">telencephalon cell migration</td>
<td align="right">13</td>
<td align="right">7</td>
<td align="right">3.61</td>
<td align="left">0.04152</td>
</tr>
<tr class="odd">
<td>367</td>
<td align="left"><a href="GO:0032615" class="uri">GO:0032615</a></td>
<td align="left">interleukin-12 production</td>
<td align="right">13</td>
<td align="right">7</td>
<td align="right">3.61</td>
<td align="left">0.04152</td>
</tr>
<tr class="even">
<td>368</td>
<td align="left"><a href="GO:0032655" class="uri">GO:0032655</a></td>
<td align="left">regulation of interleukin-12 production</td>
<td align="right">13</td>
<td align="right">7</td>
<td align="right">3.61</td>
<td align="left">0.04152</td>
</tr>
<tr class="odd">
<td>369</td>
<td align="left"><a href="GO:0022408" class="uri">GO:0022408</a></td>
<td align="left">negative regulation of cell-cell adhesio...</td>
<td align="right">8</td>
<td align="right">5</td>
<td align="right">2.22</td>
<td align="left">0.04186</td>
</tr>
<tr class="even">
<td>370</td>
<td align="left"><a href="GO:0048709" class="uri">GO:0048709</a></td>
<td align="left">oligodendrocyte differentiation</td>
<td align="right">8</td>
<td align="right">5</td>
<td align="right">2.22</td>
<td align="left">0.04186</td>
</tr>
<tr class="odd">
<td>371</td>
<td align="left"><a href="GO:0050864" class="uri">GO:0050864</a></td>
<td align="left">regulation of B cell activation</td>
<td align="right">8</td>
<td align="right">5</td>
<td align="right">2.22</td>
<td align="left">0.04186</td>
</tr>
<tr class="even">
<td>372</td>
<td align="left"><a href="GO:0016477" class="uri">GO:0016477</a></td>
<td align="left">cell migration</td>
<td align="right">123</td>
<td align="right">43</td>
<td align="right">34.13</td>
<td align="left">0.04280</td>
</tr>
<tr class="odd">
<td>373</td>
<td align="left"><a href="GO:0001816" class="uri">GO:0001816</a></td>
<td align="left">cytokine production</td>
<td align="right">72</td>
<td align="right">27</td>
<td align="right">19.98</td>
<td align="left">0.04285</td>
</tr>
<tr class="even">
<td>374</td>
<td align="left"><a href="GO:0048705" class="uri">GO:0048705</a></td>
<td align="left">skeletal system morphogenesis</td>
<td align="right">24</td>
<td align="right">11</td>
<td align="right">6.66</td>
<td align="left">0.04345</td>
</tr>
<tr class="odd">
<td>375</td>
<td align="left"><a href="GO:0009968" class="uri">GO:0009968</a></td>
<td align="left">negative regulation of signal transducti...</td>
<td align="right">120</td>
<td align="right">42</td>
<td align="right">33.30</td>
<td align="left">0.04417</td>
</tr>
<tr class="even">
<td>376</td>
<td align="left"><a href="GO:0001764" class="uri">GO:0001764</a></td>
<td align="left">neuron migration</td>
<td align="right">27</td>
<td align="right">12</td>
<td align="right">7.49</td>
<td align="left">0.04560</td>
</tr>
<tr class="odd">
<td>377</td>
<td align="left"><a href="GO:0032103" class="uri">GO:0032103</a></td>
<td align="left">positive regulation of response to exter...</td>
<td align="right">27</td>
<td align="right">12</td>
<td align="right">7.49</td>
<td align="left">0.04560</td>
</tr>
<tr class="even">
<td>378</td>
<td align="left"><a href="GO:0071241" class="uri">GO:0071241</a></td>
<td align="left">cellular response to inorganic substance</td>
<td align="right">27</td>
<td align="right">12</td>
<td align="right">7.49</td>
<td align="left">0.04560</td>
</tr>
<tr class="odd">
<td>379</td>
<td align="left"><a href="GO:0072089" class="uri">GO:0072089</a></td>
<td align="left">stem cell proliferation</td>
<td align="right">27</td>
<td align="right">12</td>
<td align="right">7.49</td>
<td align="left">0.04560</td>
</tr>
<tr class="even">
<td>380</td>
<td align="left"><a href="GO:0050776" class="uri">GO:0050776</a></td>
<td align="left">regulation of immune response</td>
<td align="right">82</td>
<td align="right">30</td>
<td align="right">22.76</td>
<td align="left">0.04708</td>
</tr>
<tr class="odd">
<td>381</td>
<td align="left"><a href="GO:0014706" class="uri">GO:0014706</a></td>
<td align="left">striated muscle tissue development</td>
<td align="right">30</td>
<td align="right">13</td>
<td align="right">8.33</td>
<td align="left">0.04712</td>
</tr>
<tr class="even">
<td>382</td>
<td align="left"><a href="GO:0040017" class="uri">GO:0040017</a></td>
<td align="left">positive regulation of locomotion</td>
<td align="right">33</td>
<td align="right">14</td>
<td align="right">9.16</td>
<td align="left">0.04816</td>
</tr>
<tr class="odd">
<td>383</td>
<td align="left"><a href="GO:0060537" class="uri">GO:0060537</a></td>
<td align="left">muscle tissue development</td>
<td align="right">33</td>
<td align="right">14</td>
<td align="right">9.16</td>
<td align="left">0.04816</td>
</tr>
<tr class="even">
<td>384</td>
<td align="left"><a href="GO:0002526" class="uri">GO:0002526</a></td>
<td align="left">acute inflammatory response</td>
<td align="right">16</td>
<td align="right">8</td>
<td align="right">4.44</td>
<td align="left">0.04820</td>
</tr>
<tr class="odd">
<td>385</td>
<td align="left"><a href="GO:0017157" class="uri">GO:0017157</a></td>
<td align="left">regulation of exocytosis</td>
<td align="right">16</td>
<td align="right">8</td>
<td align="right">4.44</td>
<td align="left">0.04820</td>
</tr>
<tr class="even">
<td>386</td>
<td align="left"><a href="GO:0021885" class="uri">GO:0021885</a></td>
<td align="left">forebrain cell migration</td>
<td align="right">16</td>
<td align="right">8</td>
<td align="right">4.44</td>
<td align="left">0.04820</td>
</tr>
<tr class="odd">
<td>387</td>
<td align="left"><a href="GO:0032845" class="uri">GO:0032845</a></td>
<td align="left">negative regulation of homeostatic proce...</td>
<td align="right">16</td>
<td align="right">8</td>
<td align="right">4.44</td>
<td align="left">0.04820</td>
</tr>
<tr class="even">
<td>388</td>
<td align="left"><a href="GO:0033559" class="uri">GO:0033559</a></td>
<td align="left">unsaturated fatty acid metabolic process</td>
<td align="right">16</td>
<td align="right">8</td>
<td align="right">4.44</td>
<td align="left">0.04820</td>
</tr>
<tr class="odd">
<td>389</td>
<td align="left"><a href="GO:0045834" class="uri">GO:0045834</a></td>
<td align="left">positive regulation of lipid metabolic p...</td>
<td align="right">16</td>
<td align="right">8</td>
<td align="right">4.44</td>
<td align="left">0.04820</td>
</tr>
<tr class="even">
<td>390</td>
<td align="left"><a href="GO:0008284" class="uri">GO:0008284</a></td>
<td align="left">positive regulation of cell proliferatio...</td>
<td align="right">48</td>
<td align="right">19</td>
<td align="right">13.32</td>
<td align="left">0.04870</td>
</tr>
<tr class="odd">
<td>391</td>
<td align="left"><a href="GO:0032147" class="uri">GO:0032147</a></td>
<td align="left">activation of protein kinase activity</td>
<td align="right">36</td>
<td align="right">15</td>
<td align="right">9.99</td>
<td align="left">0.04880</td>
</tr>
<tr class="even">
<td>392</td>
<td align="left"><a href="GO:0007067" class="uri">GO:0007067</a></td>
<td align="left">mitotic nuclear division</td>
<td align="right">45</td>
<td align="right">18</td>
<td align="right">12.49</td>
<td align="left">0.04902</td>
</tr>
<tr class="odd">
<td>393</td>
<td align="left"><a href="GO:0043410" class="uri">GO:0043410</a></td>
<td align="left">positive regulation of MAPK cascade</td>
<td align="right">45</td>
<td align="right">18</td>
<td align="right">12.49</td>
<td align="left">0.04902</td>
</tr>
<tr class="even">
<td>394</td>
<td align="left"><a href="GO:0002253" class="uri">GO:0002253</a></td>
<td align="left">activation of immune response</td>
<td align="right">42</td>
<td align="right">17</td>
<td align="right">11.66</td>
<td align="left">0.04917</td>
</tr>
</tbody>
</table>

Table 7: GO-Term compartments overrepresented in the set of
overexpressed DEGs.

<table>
<thead>
<tr class="header">
<th></th>
<th align="left">GO.ID</th>
<th align="left">Term</th>
<th align="right">Annotated</th>
<th align="right">Significant</th>
<th align="right">Expected</th>
<th align="left">classic</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>4</td>
<td align="left"><a href="GO:0005578" class="uri">GO:0005578</a></td>
<td align="left">proteinaceous extracellular matrix</td>
<td align="right">92</td>
<td align="right">42</td>
<td align="right">25.80</td>
<td align="left">0.00016</td>
</tr>
<tr class="even">
<td>5</td>
<td align="left"><a href="GO:0005604" class="uri">GO:0005604</a></td>
<td align="left">basement membrane</td>
<td align="right">41</td>
<td align="right">22</td>
<td align="right">11.50</td>
<td align="left">0.00042</td>
</tr>
<tr class="odd">
<td>6</td>
<td align="left"><a href="GO:0005840" class="uri">GO:0005840</a></td>
<td align="left">ribosome</td>
<td align="right">54</td>
<td align="right">27</td>
<td align="right">15.14</td>
<td align="left">0.00043</td>
</tr>
<tr class="even">
<td>7</td>
<td align="left"><a href="GO:0044420" class="uri">GO:0044420</a></td>
<td align="left">extracellular matrix component</td>
<td align="right">47</td>
<td align="right">24</td>
<td align="right">13.18</td>
<td align="left">0.00061</td>
</tr>
<tr class="odd">
<td>8</td>
<td align="left"><a href="GO:0005615" class="uri">GO:0005615</a></td>
<td align="left">extracellular space</td>
<td align="right">92</td>
<td align="right">40</td>
<td align="right">25.80</td>
<td align="left">0.00082</td>
</tr>
<tr class="even">
<td>9</td>
<td align="left"><a href="GO:0098552" class="uri">GO:0098552</a></td>
<td align="left">side of membrane</td>
<td align="right">51</td>
<td align="right">25</td>
<td align="right">14.30</td>
<td align="left">0.00102</td>
</tr>
<tr class="odd">
<td>10</td>
<td align="left"><a href="GO:0019897" class="uri">GO:0019897</a></td>
<td align="left">extrinsic component of plasma membrane</td>
<td align="right">28</td>
<td align="right">16</td>
<td align="right">7.85</td>
<td align="left">0.00107</td>
</tr>
<tr class="even">
<td>11</td>
<td align="left"><a href="GO:0005618" class="uri">GO:0005618</a></td>
<td align="left">cell wall</td>
<td align="right">21</td>
<td align="right">13</td>
<td align="right">5.89</td>
<td align="left">0.00115</td>
</tr>
<tr class="odd">
<td>12</td>
<td align="left"><a href="GO:0042588" class="uri">GO:0042588</a></td>
<td align="left">zymogen granule</td>
<td align="right">19</td>
<td align="right">12</td>
<td align="right">5.33</td>
<td align="left">0.00141</td>
</tr>
<tr class="even">
<td>13</td>
<td align="left"><a href="GO:0042589" class="uri">GO:0042589</a></td>
<td align="left">zymogen granule membrane</td>
<td align="right">19</td>
<td align="right">12</td>
<td align="right">5.33</td>
<td align="left">0.00141</td>
</tr>
<tr class="odd">
<td>14</td>
<td align="left"><a href="GO:0030312" class="uri">GO:0030312</a></td>
<td align="left">external encapsulating structure</td>
<td align="right">23</td>
<td align="right">13</td>
<td align="right">6.45</td>
<td align="left">0.00363</td>
</tr>
<tr class="even">
<td>15</td>
<td align="left"><a href="GO:0005581" class="uri">GO:0005581</a></td>
<td align="left">collagen trimer</td>
<td align="right">21</td>
<td align="right">12</td>
<td align="right">5.89</td>
<td align="left">0.00460</td>
</tr>
<tr class="odd">
<td>16</td>
<td align="left"><a href="GO:0071944" class="uri">GO:0071944</a></td>
<td align="left">cell periphery</td>
<td align="right">496</td>
<td align="right">162</td>
<td align="right">139.07</td>
<td align="left">0.00467</td>
</tr>
<tr class="even">
<td>17</td>
<td align="left"><a href="GO:0022627" class="uri">GO:0022627</a></td>
<td align="left">cytosolic small ribosomal subunit</td>
<td align="right">12</td>
<td align="right">8</td>
<td align="right">3.36</td>
<td align="left">0.00591</td>
</tr>
<tr class="odd">
<td>18</td>
<td align="left"><a href="GO:0043083" class="uri">GO:0043083</a></td>
<td align="left">synaptic cleft</td>
<td align="right">12</td>
<td align="right">8</td>
<td align="right">3.36</td>
<td align="left">0.00591</td>
</tr>
<tr class="even">
<td>19</td>
<td align="left"><a href="GO:0019867" class="uri">GO:0019867</a></td>
<td align="left">outer membrane</td>
<td align="right">24</td>
<td align="right">13</td>
<td align="right">6.73</td>
<td align="left">0.00593</td>
</tr>
<tr class="odd">
<td>20</td>
<td align="left"><a href="GO:0043195" class="uri">GO:0043195</a></td>
<td align="left">terminal bouton</td>
<td align="right">22</td>
<td align="right">12</td>
<td align="right">6.17</td>
<td align="left">0.00759</td>
</tr>
<tr class="even">
<td>21</td>
<td align="left"><a href="GO:0045335" class="uri">GO:0045335</a></td>
<td align="left">phagocytic vesicle</td>
<td align="right">22</td>
<td align="right">12</td>
<td align="right">6.17</td>
<td align="left">0.00759</td>
</tr>
<tr class="odd">
<td>22</td>
<td align="left"><a href="GO:0044391" class="uri">GO:0044391</a></td>
<td align="left">ribosomal subunit</td>
<td align="right">38</td>
<td align="right">18</td>
<td align="right">10.65</td>
<td align="left">0.00812</td>
</tr>
<tr class="even">
<td>23</td>
<td align="left"><a href="GO:0005886" class="uri">GO:0005886</a></td>
<td align="left">plasma membrane</td>
<td align="right">461</td>
<td align="right">150</td>
<td align="right">129.26</td>
<td align="left">0.00814</td>
</tr>
<tr class="odd">
<td>24</td>
<td align="left"><a href="GO:0015935" class="uri">GO:0015935</a></td>
<td align="left">small ribosomal subunit</td>
<td align="right">15</td>
<td align="right">9</td>
<td align="right">4.21</td>
<td align="left">0.00928</td>
</tr>
<tr class="even">
<td>25</td>
<td align="left"><a href="GO:0098562" class="uri">GO:0098562</a></td>
<td align="left">cytoplasmic side of membrane</td>
<td align="right">39</td>
<td align="right">18</td>
<td align="right">10.94</td>
<td align="left">0.01124</td>
</tr>
<tr class="odd">
<td>26</td>
<td align="left"><a href="GO:0031982" class="uri">GO:0031982</a></td>
<td align="left">vesicle</td>
<td align="right">325</td>
<td align="right">108</td>
<td align="right">91.13</td>
<td align="left">0.01376</td>
</tr>
<tr class="even">
<td>27</td>
<td align="left"><a href="GO:0030667" class="uri">GO:0030667</a></td>
<td align="left">secretory granule membrane</td>
<td align="right">26</td>
<td align="right">13</td>
<td align="right">7.29</td>
<td align="left">0.01387</td>
</tr>
<tr class="odd">
<td>28</td>
<td align="left"><a href="GO:0012506" class="uri">GO:0012506</a></td>
<td align="left">vesicle membrane</td>
<td align="right">66</td>
<td align="right">27</td>
<td align="right">18.51</td>
<td align="left">0.01493</td>
</tr>
<tr class="even">
<td>29</td>
<td align="left"><a href="GO:0030424" class="uri">GO:0030424</a></td>
<td align="left">axon</td>
<td align="right">78</td>
<td align="right">31</td>
<td align="right">21.87</td>
<td align="left">0.01506</td>
</tr>
<tr class="odd">
<td>30</td>
<td align="left"><a href="GO:0022626" class="uri">GO:0022626</a></td>
<td align="left">cytosolic ribosome</td>
<td align="right">29</td>
<td align="right">14</td>
<td align="right">8.13</td>
<td align="left">0.01558</td>
</tr>
<tr class="even">
<td>31</td>
<td align="left"><a href="GO:0044421" class="uri">GO:0044421</a></td>
<td align="left">extracellular region part</td>
<td align="right">310</td>
<td align="right">103</td>
<td align="right">86.92</td>
<td align="left">0.01630</td>
</tr>
<tr class="odd">
<td>32</td>
<td align="left"><a href="GO:0033267" class="uri">GO:0033267</a></td>
<td align="left">axon part</td>
<td align="right">46</td>
<td align="right">20</td>
<td align="right">12.90</td>
<td align="left">0.01665</td>
</tr>
<tr class="even">
<td>33</td>
<td align="left"><a href="GO:0043679" class="uri">GO:0043679</a></td>
<td align="left">axon terminus</td>
<td align="right">24</td>
<td align="right">12</td>
<td align="right">6.73</td>
<td align="left">0.01790</td>
</tr>
<tr class="odd">
<td>34</td>
<td align="left"><a href="GO:0044306" class="uri">GO:0044306</a></td>
<td align="left">neuron projection terminus</td>
<td align="right">24</td>
<td align="right">12</td>
<td align="right">6.73</td>
<td align="left">0.01790</td>
</tr>
<tr class="even">
<td>35</td>
<td align="left"><a href="GO:0030659" class="uri">GO:0030659</a></td>
<td align="left">cytoplasmic vesicle membrane</td>
<td align="right">64</td>
<td align="right">26</td>
<td align="right">17.94</td>
<td align="left">0.01851</td>
</tr>
<tr class="odd">
<td>36</td>
<td align="left"><a href="GO:0044433" class="uri">GO:0044433</a></td>
<td align="left">cytoplasmic vesicle part</td>
<td align="right">73</td>
<td align="right">29</td>
<td align="right">20.47</td>
<td align="left">0.01856</td>
</tr>
<tr class="even">
<td>37</td>
<td align="left"><a href="GO:0031410" class="uri">GO:0031410</a></td>
<td align="left">cytoplasmic vesicle</td>
<td align="right">142</td>
<td align="right">51</td>
<td align="right">39.82</td>
<td align="left">0.02051</td>
</tr>
<tr class="odd">
<td>38</td>
<td align="left"><a href="GO:0000811" class="uri">GO:0000811</a></td>
<td align="left">GINS complex</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.84</td>
<td align="left">0.02195</td>
</tr>
<tr class="even">
<td>39</td>
<td align="left"><a href="GO:0031261" class="uri">GO:0031261</a></td>
<td align="left">DNA replication preinitiation complex</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.84</td>
<td align="left">0.02195</td>
</tr>
<tr class="odd">
<td>40</td>
<td align="left"><a href="GO:0031298" class="uri">GO:0031298</a></td>
<td align="left">replication fork protection complex</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.84</td>
<td align="left">0.02195</td>
</tr>
<tr class="even">
<td>41</td>
<td align="left"><a href="GO:0033162" class="uri">GO:0033162</a></td>
<td align="left">melanosome membrane</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.84</td>
<td align="left">0.02195</td>
</tr>
<tr class="odd">
<td>42</td>
<td align="left"><a href="GO:0045009" class="uri">GO:0045009</a></td>
<td align="left">chitosome</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.84</td>
<td align="left">0.02195</td>
</tr>
<tr class="even">
<td>43</td>
<td align="left"><a href="GO:0032839" class="uri">GO:0032839</a></td>
<td align="left">dendrite cytoplasm</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.40</td>
<td align="left">0.02381</td>
</tr>
<tr class="odd">
<td>44</td>
<td align="left"><a href="GO:0009897" class="uri">GO:0009897</a></td>
<td align="left">external side of plasma membrane</td>
<td align="right">12</td>
<td align="right">7</td>
<td align="right">3.36</td>
<td align="left">0.02647</td>
</tr>
<tr class="even">
<td>45</td>
<td align="left"><a href="GO:0005788" class="uri">GO:0005788</a></td>
<td align="left">endoplasmic reticulum lumen</td>
<td align="right">20</td>
<td align="right">10</td>
<td align="right">5.61</td>
<td align="left">0.03002</td>
</tr>
<tr class="odd">
<td>46</td>
<td align="left"><a href="GO:0030139" class="uri">GO:0030139</a></td>
<td align="left">endocytic vesicle</td>
<td align="right">34</td>
<td align="right">15</td>
<td align="right">9.53</td>
<td align="left">0.03141</td>
</tr>
<tr class="even">
<td>47</td>
<td align="left"><a href="GO:0009986" class="uri">GO:0009986</a></td>
<td align="left">cell surface</td>
<td align="right">67</td>
<td align="right">26</td>
<td align="right">18.79</td>
<td align="left">0.03418</td>
</tr>
<tr class="odd">
<td>48</td>
<td align="left"><a href="GO:0044456" class="uri">GO:0044456</a></td>
<td align="left">synapse part</td>
<td align="right">64</td>
<td align="right">25</td>
<td align="right">17.94</td>
<td align="left">0.03449</td>
</tr>
<tr class="even">
<td>49</td>
<td align="left"><a href="GO:0016023" class="uri">GO:0016023</a></td>
<td align="left">cytoplasmic membrane-bounded vesicle</td>
<td align="right">127</td>
<td align="right">45</td>
<td align="right">35.61</td>
<td align="left">0.03639</td>
</tr>
<tr class="odd">
<td>50</td>
<td align="left"><a href="GO:0031988" class="uri">GO:0031988</a></td>
<td align="left">membrane-bounded vesicle</td>
<td align="right">307</td>
<td align="right">99</td>
<td align="right">86.08</td>
<td align="left">0.04311</td>
</tr>
<tr class="even">
<td>51</td>
<td align="left"><a href="GO:0030141" class="uri">GO:0030141</a></td>
<td align="left">secretory granule</td>
<td align="right">44</td>
<td align="right">18</td>
<td align="right">12.34</td>
<td align="left">0.04312</td>
</tr>
<tr class="odd">
<td>52</td>
<td align="left"><a href="GO:0000922" class="uri">GO:0000922</a></td>
<td align="left">spindle pole</td>
<td align="right">8</td>
<td align="right">5</td>
<td align="right">2.24</td>
<td align="left">0.04368</td>
</tr>
<tr class="even">
<td>53</td>
<td align="left"><a href="GO:0009279" class="uri">GO:0009279</a></td>
<td align="left">cell outer membrane</td>
<td align="right">13</td>
<td align="right">7</td>
<td align="right">3.65</td>
<td align="left">0.04379</td>
</tr>
</tbody>
</table>

Table 8: GO-Term functions overrepresented in the set of underexpressed
DEGs.

<table>
<thead>
<tr class="header">
<th></th>
<th align="left">GO.ID</th>
<th align="left">Term</th>
<th align="right">Annotated</th>
<th align="right">Significant</th>
<th align="right">Expected</th>
<th align="left">classic</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>24</td>
<td align="left"><a href="GO:0030246" class="uri">GO:0030246</a></td>
<td align="left">carbohydrate binding</td>
<td align="right">120</td>
<td align="right">45</td>
<td align="right">26.76</td>
<td align="left">0.00010</td>
</tr>
<tr class="even">
<td>25</td>
<td align="left"><a href="GO:0016462" class="uri">GO:0016462</a></td>
<td align="left">pyrophosphatase activity</td>
<td align="right">508</td>
<td align="right">148</td>
<td align="right">113.26</td>
<td align="left">0.00011</td>
</tr>
<tr class="odd">
<td>26</td>
<td align="left"><a href="GO:0016817" class="uri">GO:0016817</a></td>
<td align="left">hydrolase activity, acting on acid anhyd...</td>
<td align="right">509</td>
<td align="right">148</td>
<td align="right">113.49</td>
<td align="left">0.00012</td>
</tr>
<tr class="even">
<td>27</td>
<td align="left"><a href="GO:0016818" class="uri">GO:0016818</a></td>
<td align="left">hydrolase activity, acting on acid anhyd...</td>
<td align="right">509</td>
<td align="right">148</td>
<td align="right">113.49</td>
<td align="left">0.00012</td>
</tr>
<tr class="odd">
<td>28</td>
<td align="left"><a href="GO:0015238" class="uri">GO:0015238</a></td>
<td align="left">drug transmembrane transporter activity</td>
<td align="right">27</td>
<td align="right">15</td>
<td align="right">6.02</td>
<td align="left">0.00017</td>
</tr>
<tr class="even">
<td>29</td>
<td align="left"><a href="GO:0004497" class="uri">GO:0004497</a></td>
<td align="left">monooxygenase activity</td>
<td align="right">59</td>
<td align="right">25</td>
<td align="right">13.15</td>
<td align="left">0.00043</td>
</tr>
<tr class="odd">
<td>30</td>
<td align="left"><a href="GO:0090484" class="uri">GO:0090484</a></td>
<td align="left">drug transporter activity</td>
<td align="right">29</td>
<td align="right">15</td>
<td align="right">6.47</td>
<td align="left">0.00048</td>
</tr>
<tr class="even">
<td>31</td>
<td align="left"><a href="GO:0005220" class="uri">GO:0005220</a></td>
<td align="left">inositol 1,4,5-trisphosphate-sensitive c...</td>
<td align="right">11</td>
<td align="right">8</td>
<td align="right">2.45</td>
<td align="left">0.00051</td>
</tr>
<tr class="odd">
<td>32</td>
<td align="left"><a href="GO:0008191" class="uri">GO:0008191</a></td>
<td align="left">metalloendopeptidase inhibitor activity</td>
<td align="right">5</td>
<td align="right">5</td>
<td align="right">1.11</td>
<td align="left">0.00055</td>
</tr>
<tr class="even">
<td>33</td>
<td align="left"><a href="GO:0010576" class="uri">GO:0010576</a></td>
<td align="left">metalloenzyme regulator activity</td>
<td align="right">5</td>
<td align="right">5</td>
<td align="right">1.11</td>
<td align="left">0.00055</td>
</tr>
<tr class="odd">
<td>34</td>
<td align="left"><a href="GO:0048551" class="uri">GO:0048551</a></td>
<td align="left">metalloenzyme inhibitor activity</td>
<td align="right">5</td>
<td align="right">5</td>
<td align="right">1.11</td>
<td align="left">0.00055</td>
</tr>
<tr class="even">
<td>35</td>
<td align="left"><a href="GO:0004888" class="uri">GO:0004888</a></td>
<td align="left">transmembrane signaling receptor activit...</td>
<td align="right">292</td>
<td align="right">89</td>
<td align="right">65.10</td>
<td align="left">0.00055</td>
</tr>
<tr class="odd">
<td>36</td>
<td align="left"><a href="GO:0016524" class="uri">GO:0016524</a></td>
<td align="left">latrotoxin receptor activity</td>
<td align="right">9</td>
<td align="right">7</td>
<td align="right">2.01</td>
<td align="left">0.00063</td>
</tr>
<tr class="even">
<td>37</td>
<td align="left"><a href="GO:0015085" class="uri">GO:0015085</a></td>
<td align="left">calcium ion transmembrane transporter ac...</td>
<td align="right">39</td>
<td align="right">18</td>
<td align="right">8.70</td>
<td align="left">0.00080</td>
</tr>
<tr class="odd">
<td>38</td>
<td align="left"><a href="GO:0015491" class="uri">GO:0015491</a></td>
<td align="left">cation:cation antiporter activity</td>
<td align="right">14</td>
<td align="right">9</td>
<td align="right">3.12</td>
<td align="left">0.00088</td>
</tr>
<tr class="even">
<td>39</td>
<td align="left"><a href="GO:0016646" class="uri">GO:0016646</a></td>
<td align="left">oxidoreductase activity, acting on the C...</td>
<td align="right">14</td>
<td align="right">9</td>
<td align="right">3.12</td>
<td align="left">0.00088</td>
</tr>
<tr class="odd">
<td>40</td>
<td align="left"><a href="GO:0061134" class="uri">GO:0061134</a></td>
<td align="left">peptidase regulator activity</td>
<td align="right">56</td>
<td align="right">23</td>
<td align="right">12.49</td>
<td align="left">0.00120</td>
</tr>
<tr class="even">
<td>41</td>
<td align="left"><a href="GO:0005432" class="uri">GO:0005432</a></td>
<td align="left">calcium:sodium antiporter activity</td>
<td align="right">12</td>
<td align="right">8</td>
<td align="right">2.68</td>
<td align="left">0.00124</td>
</tr>
<tr class="odd">
<td>42</td>
<td align="left"><a href="GO:0015368" class="uri">GO:0015368</a></td>
<td align="left">calcium:cation antiporter activity</td>
<td align="right">12</td>
<td align="right">8</td>
<td align="right">2.68</td>
<td align="left">0.00124</td>
</tr>
<tr class="even">
<td>43</td>
<td align="left"><a href="GO:0016705" class="uri">GO:0016705</a></td>
<td align="left">oxidoreductase activity, acting on paire...</td>
<td align="right">83</td>
<td align="right">31</td>
<td align="right">18.51</td>
<td align="left">0.00126</td>
</tr>
<tr class="odd">
<td>44</td>
<td align="left"><a href="GO:0015081" class="uri">GO:0015081</a></td>
<td align="left">sodium ion transmembrane transporter act...</td>
<td align="right">23</td>
<td align="right">12</td>
<td align="right">5.13</td>
<td align="left">0.00162</td>
</tr>
<tr class="even">
<td>45</td>
<td align="left"><a href="GO:0008574" class="uri">GO:0008574</a></td>
<td align="left">ATP-dependent microtubule motor activity...</td>
<td align="right">10</td>
<td align="right">7</td>
<td align="right">2.23</td>
<td align="left">0.00170</td>
</tr>
<tr class="odd">
<td>46</td>
<td align="left"><a href="GO:0016491" class="uri">GO:0016491</a></td>
<td align="left">oxidoreductase activity</td>
<td align="right">425</td>
<td align="right">120</td>
<td align="right">94.76</td>
<td align="left">0.00176</td>
</tr>
<tr class="even">
<td>47</td>
<td align="left"><a href="GO:0004796" class="uri">GO:0004796</a></td>
<td align="left">thromboxane-A synthase activity</td>
<td align="right">8</td>
<td align="right">6</td>
<td align="right">1.78</td>
<td align="left">0.00224</td>
</tr>
<tr class="odd">
<td>48</td>
<td align="left"><a href="GO:0003989" class="uri">GO:0003989</a></td>
<td align="left">acetyl-CoA carboxylase activity</td>
<td align="right">4</td>
<td align="right">4</td>
<td align="right">0.89</td>
<td align="left">0.00246</td>
</tr>
<tr class="even">
<td>49</td>
<td align="left"><a href="GO:0004161" class="uri">GO:0004161</a></td>
<td align="left">dimethylallyltranstransferase activity</td>
<td align="right">6</td>
<td align="right">5</td>
<td align="right">1.34</td>
<td align="left">0.00268</td>
</tr>
<tr class="odd">
<td>50</td>
<td align="left"><a href="GO:0004337" class="uri">GO:0004337</a></td>
<td align="left">geranyltranstransferase activity</td>
<td align="right">6</td>
<td align="right">5</td>
<td align="right">1.34</td>
<td align="left">0.00268</td>
</tr>
<tr class="even">
<td>51</td>
<td align="left"><a href="GO:0016682" class="uri">GO:0016682</a></td>
<td align="left">oxidoreductase activity, acting on diphe...</td>
<td align="right">6</td>
<td align="right">5</td>
<td align="right">1.34</td>
<td align="left">0.00268</td>
</tr>
<tr class="odd">
<td>52</td>
<td align="left"><a href="GO:0052716" class="uri">GO:0052716</a></td>
<td align="left">hydroquinone:oxygen oxidoreductase activ...</td>
<td align="right">6</td>
<td align="right">5</td>
<td align="right">1.34</td>
<td align="left">0.00268</td>
</tr>
<tr class="even">
<td>53</td>
<td align="left"><a href="GO:0004857" class="uri">GO:0004857</a></td>
<td align="left">enzyme inhibitor activity</td>
<td align="right">90</td>
<td align="right">32</td>
<td align="right">20.07</td>
<td align="left">0.00268</td>
</tr>
<tr class="odd">
<td>54</td>
<td align="left"><a href="GO:0016820" class="uri">GO:0016820</a></td>
<td align="left">hydrolase activity, acting on acid anhyd...</td>
<td align="right">80</td>
<td align="right">29</td>
<td align="right">17.84</td>
<td align="left">0.00301</td>
</tr>
<tr class="even">
<td>55</td>
<td align="left"><a href="GO:0042626" class="uri">GO:0042626</a></td>
<td align="left">ATPase activity, coupled to transmembran...</td>
<td align="right">80</td>
<td align="right">29</td>
<td align="right">17.84</td>
<td align="left">0.00301</td>
</tr>
<tr class="odd">
<td>56</td>
<td align="left"><a href="GO:0005516" class="uri">GO:0005516</a></td>
<td align="left">calmodulin binding</td>
<td align="right">98</td>
<td align="right">34</td>
<td align="right">21.85</td>
<td align="left">0.00319</td>
</tr>
<tr class="even">
<td>57</td>
<td align="left"><a href="GO:0043394" class="uri">GO:0043394</a></td>
<td align="left">proteoglycan binding</td>
<td align="right">16</td>
<td align="right">9</td>
<td align="right">3.57</td>
<td align="left">0.00325</td>
</tr>
<tr class="odd">
<td>58</td>
<td align="left"><a href="GO:0016879" class="uri">GO:0016879</a></td>
<td align="left">ligase activity, forming carbon-nitrogen...</td>
<td align="right">28</td>
<td align="right">13</td>
<td align="right">6.24</td>
<td align="left">0.00397</td>
</tr>
<tr class="even">
<td>59</td>
<td align="left"><a href="GO:0015399" class="uri">GO:0015399</a></td>
<td align="left">primary active transmembrane transporter...</td>
<td align="right">82</td>
<td align="right">29</td>
<td align="right">18.28</td>
<td align="left">0.00454</td>
</tr>
<tr class="odd">
<td>60</td>
<td align="left"><a href="GO:0015405" class="uri">GO:0015405</a></td>
<td align="left">P-P-bond-hydrolysis-driven transmembrane...</td>
<td align="right">82</td>
<td align="right">29</td>
<td align="right">18.28</td>
<td align="left">0.00454</td>
</tr>
<tr class="even">
<td>61</td>
<td align="left"><a href="GO:0043492" class="uri">GO:0043492</a></td>
<td align="left">ATPase activity, coupled to movement of ...</td>
<td align="right">82</td>
<td align="right">29</td>
<td align="right">18.28</td>
<td align="left">0.00454</td>
</tr>
<tr class="odd">
<td>62</td>
<td align="left"><a href="GO:0008017" class="uri">GO:0008017</a></td>
<td align="left">microtubule binding</td>
<td align="right">93</td>
<td align="right">32</td>
<td align="right">20.74</td>
<td align="left">0.00479</td>
</tr>
<tr class="even">
<td>63</td>
<td align="left"><a href="GO:0015631" class="uri">GO:0015631</a></td>
<td align="left">tubulin binding</td>
<td align="right">126</td>
<td align="right">41</td>
<td align="right">28.09</td>
<td align="left">0.00485</td>
</tr>
<tr class="odd">
<td>64</td>
<td align="left"><a href="GO:0032403" class="uri">GO:0032403</a></td>
<td align="left">protein complex binding</td>
<td align="right">364</td>
<td align="right">102</td>
<td align="right">81.16</td>
<td align="left">0.00490</td>
</tr>
<tr class="even">
<td>65</td>
<td align="left"><a href="GO:0005506" class="uri">GO:0005506</a></td>
<td align="left">iron ion binding</td>
<td align="right">76</td>
<td align="right">27</td>
<td align="right">16.94</td>
<td align="left">0.00565</td>
</tr>
<tr class="odd">
<td>66</td>
<td align="left"><a href="GO:0016776" class="uri">GO:0016776</a></td>
<td align="left">phosphotransferase activity, phosphate g...</td>
<td align="right">20</td>
<td align="right">10</td>
<td align="right">4.46</td>
<td align="left">0.00590</td>
</tr>
<tr class="even">
<td>67</td>
<td align="left"><a href="GO:0004725" class="uri">GO:0004725</a></td>
<td align="left">protein tyrosine phosphatase activity</td>
<td align="right">70</td>
<td align="right">25</td>
<td align="right">15.61</td>
<td align="left">0.00702</td>
</tr>
<tr class="odd">
<td>68</td>
<td align="left"><a href="GO:0044877" class="uri">GO:0044877</a></td>
<td align="left">macromolecular complex binding</td>
<td align="right">516</td>
<td align="right">138</td>
<td align="right">115.05</td>
<td align="left">0.00738</td>
</tr>
<tr class="even">
<td>69</td>
<td align="left"><a href="GO:0016679" class="uri">GO:0016679</a></td>
<td align="left">oxidoreductase activity, acting on diphe...</td>
<td align="right">7</td>
<td align="right">5</td>
<td align="right">1.56</td>
<td align="left">0.00765</td>
</tr>
<tr class="odd">
<td>70</td>
<td align="left"><a href="GO:0004930" class="uri">GO:0004930</a></td>
<td align="left">G-protein coupled receptor activity</td>
<td align="right">123</td>
<td align="right">39</td>
<td align="right">27.42</td>
<td align="left">0.00951</td>
</tr>
<tr class="even">
<td>71</td>
<td align="left"><a href="GO:0051015" class="uri">GO:0051015</a></td>
<td align="left">actin filament binding</td>
<td align="right">68</td>
<td align="right">24</td>
<td align="right">15.16</td>
<td align="left">0.00961</td>
</tr>
<tr class="odd">
<td>72</td>
<td align="left"><a href="GO:0000150" class="uri">GO:0000150</a></td>
<td align="left">recombinase activity</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.11</td>
<td align="left">0.01012</td>
</tr>
<tr class="even">
<td>73</td>
<td align="left"><a href="GO:0004735" class="uri">GO:0004735</a></td>
<td align="left">pyrroline-5-carboxylate reductase activi...</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.11</td>
<td align="left">0.01012</td>
</tr>
<tr class="odd">
<td>74</td>
<td align="left"><a href="GO:0019238" class="uri">GO:0019238</a></td>
<td align="left">cyclohydrolase activity</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.11</td>
<td align="left">0.01012</td>
</tr>
<tr class="even">
<td>75</td>
<td align="left"><a href="GO:0043236" class="uri">GO:0043236</a></td>
<td align="left">laminin binding</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.11</td>
<td align="left">0.01012</td>
</tr>
<tr class="odd">
<td>76</td>
<td align="left"><a href="GO:0046912" class="uri">GO:0046912</a></td>
<td align="left">transferase activity, transferring acyl ...</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.11</td>
<td align="left">0.01012</td>
</tr>
<tr class="even">
<td>77</td>
<td align="left"><a href="GO:0000981" class="uri">GO:0000981</a></td>
<td align="left">RNA polymerase II transcription factor a...</td>
<td align="right">51</td>
<td align="right">19</td>
<td align="right">11.37</td>
<td align="left">0.01078</td>
</tr>
<tr class="odd">
<td>78</td>
<td align="left"><a href="GO:0004089" class="uri">GO:0004089</a></td>
<td align="left">carbonate dehydratase activity</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01107</td>
</tr>
<tr class="even">
<td>79</td>
<td align="left"><a href="GO:0004329" class="uri">GO:0004329</a></td>
<td align="left">formate-tetrahydrofolate ligase activity</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01107</td>
</tr>
<tr class="odd">
<td>80</td>
<td align="left"><a href="GO:0004477" class="uri">GO:0004477</a></td>
<td align="left">methenyltetrahydrofolate cyclohydrolase ...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01107</td>
</tr>
<tr class="even">
<td>81</td>
<td align="left"><a href="GO:0004488" class="uri">GO:0004488</a></td>
<td align="left">methylenetetrahydrofolate dehydrogenase ...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01107</td>
</tr>
<tr class="odd">
<td>82</td>
<td align="left"><a href="GO:0010521" class="uri">GO:0010521</a></td>
<td align="left">telomerase inhibitor activity</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01107</td>
</tr>
<tr class="even">
<td>83</td>
<td align="left"><a href="GO:0015643" class="uri">GO:0015643</a></td>
<td align="left">toxic substance binding</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01107</td>
</tr>
<tr class="odd">
<td>84</td>
<td align="left"><a href="GO:0018585" class="uri">GO:0018585</a></td>
<td align="left">fluorene oxygenase activity</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01107</td>
</tr>
<tr class="even">
<td>85</td>
<td align="left"><a href="GO:0019136" class="uri">GO:0019136</a></td>
<td align="left">deoxynucleoside kinase activity</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01107</td>
</tr>
<tr class="odd">
<td>86</td>
<td align="left"><a href="GO:0048407" class="uri">GO:0048407</a></td>
<td align="left">platelet-derived growth factor binding</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01107</td>
</tr>
<tr class="even">
<td>87</td>
<td align="left"><a href="GO:0070739" class="uri">GO:0070739</a></td>
<td align="left">protein-glutamic acid ligase activity</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01107</td>
</tr>
<tr class="odd">
<td>88</td>
<td align="left"><a href="GO:0070740" class="uri">GO:0070740</a></td>
<td align="left">tubulin-glutamic acid ligase activity</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01107</td>
</tr>
<tr class="even">
<td>89</td>
<td align="left"><a href="GO:0022804" class="uri">GO:0022804</a></td>
<td align="left">active transmembrane transporter activit...</td>
<td align="right">155</td>
<td align="right">47</td>
<td align="right">34.56</td>
<td align="left">0.01164</td>
</tr>
<tr class="odd">
<td>90</td>
<td align="left"><a href="GO:0003697" class="uri">GO:0003697</a></td>
<td align="left">single-stranded DNA binding</td>
<td align="right">28</td>
<td align="right">12</td>
<td align="right">6.24</td>
<td align="left">0.01199</td>
</tr>
<tr class="even">
<td>91</td>
<td align="left"><a href="GO:0004017" class="uri">GO:0004017</a></td>
<td align="left">adenylate kinase activity</td>
<td align="right">8</td>
<td align="right">5</td>
<td align="right">1.78</td>
<td align="left">0.01667</td>
</tr>
<tr class="odd">
<td>92</td>
<td align="left"><a href="GO:0043015" class="uri">GO:0043015</a></td>
<td align="left">gamma-tubulin binding</td>
<td align="right">8</td>
<td align="right">5</td>
<td align="right">1.78</td>
<td align="left">0.01667</td>
</tr>
<tr class="even">
<td>93</td>
<td align="left"><a href="GO:0009055" class="uri">GO:0009055</a></td>
<td align="left">electron carrier activity</td>
<td align="right">43</td>
<td align="right">16</td>
<td align="right">9.59</td>
<td align="left">0.01873</td>
</tr>
<tr class="odd">
<td>94</td>
<td align="left"><a href="GO:0016645" class="uri">GO:0016645</a></td>
<td align="left">oxidoreductase activity, acting on the C...</td>
<td align="right">23</td>
<td align="right">10</td>
<td align="right">5.13</td>
<td align="left">0.01898</td>
</tr>
<tr class="even">
<td>95</td>
<td align="left"><a href="GO:0008559" class="uri">GO:0008559</a></td>
<td align="left">xenobiotic-transporting ATPase activity</td>
<td align="right">11</td>
<td align="right">6</td>
<td align="right">2.45</td>
<td align="left">0.01980</td>
</tr>
<tr class="odd">
<td>96</td>
<td align="left"><a href="GO:0030332" class="uri">GO:0030332</a></td>
<td align="left">cyclin binding</td>
<td align="right">11</td>
<td align="right">6</td>
<td align="right">2.45</td>
<td align="left">0.01980</td>
</tr>
<tr class="even">
<td>97</td>
<td align="left"><a href="GO:0042910" class="uri">GO:0042910</a></td>
<td align="left">xenobiotic transporter activity</td>
<td align="right">11</td>
<td align="right">6</td>
<td align="right">2.45</td>
<td align="left">0.01980</td>
</tr>
<tr class="odd">
<td>98</td>
<td align="left"><a href="GO:0004713" class="uri">GO:0004713</a></td>
<td align="left">protein tyrosine kinase activity</td>
<td align="right">144</td>
<td align="right">43</td>
<td align="right">32.11</td>
<td align="left">0.02003</td>
</tr>
<tr class="even">
<td>99</td>
<td align="left"><a href="GO:0008138" class="uri">GO:0008138</a></td>
<td align="left">protein tyrosine/serine/threonine phosph...</td>
<td align="right">20</td>
<td align="right">9</td>
<td align="right">4.46</td>
<td align="left">0.02007</td>
</tr>
<tr class="odd">
<td>100</td>
<td align="left"><a href="GO:0015278" class="uri">GO:0015278</a></td>
<td align="left">calcium-release channel activity</td>
<td align="right">20</td>
<td align="right">9</td>
<td align="right">4.46</td>
<td align="left">0.02007</td>
</tr>
<tr class="even">
<td>101</td>
<td align="left"><a href="GO:0005217" class="uri">GO:0005217</a></td>
<td align="left">intracellular ligand-gated ion channel a...</td>
<td align="right">17</td>
<td align="right">8</td>
<td align="right">3.79</td>
<td align="left">0.02081</td>
</tr>
<tr class="odd">
<td>102</td>
<td align="left"><a href="GO:0015299" class="uri">GO:0015299</a></td>
<td align="left">solute:proton antiporter activity</td>
<td align="right">6</td>
<td align="right">4</td>
<td align="right">1.34</td>
<td align="left">0.02501</td>
</tr>
<tr class="even">
<td>103</td>
<td align="left"><a href="GO:0050840" class="uri">GO:0050840</a></td>
<td align="left">extracellular matrix binding</td>
<td align="right">6</td>
<td align="right">4</td>
<td align="right">1.34</td>
<td align="left">0.02501</td>
</tr>
<tr class="odd">
<td>104</td>
<td align="left"><a href="GO:0051920" class="uri">GO:0051920</a></td>
<td align="left">peroxiredoxin activity</td>
<td align="right">6</td>
<td align="right">4</td>
<td align="right">1.34</td>
<td align="left">0.02501</td>
</tr>
<tr class="even">
<td>105</td>
<td align="left"><a href="GO:0097367" class="uri">GO:0097367</a></td>
<td align="left">carbohydrate derivative binding</td>
<td align="right">1467</td>
<td align="right">355</td>
<td align="right">327.08</td>
<td align="left">0.02544</td>
</tr>
<tr class="odd">
<td>106</td>
<td align="left"><a href="GO:0043167" class="uri">GO:0043167</a></td>
<td align="left">ion binding</td>
<td align="right">3120</td>
<td align="right">728</td>
<td align="right">695.63</td>
<td align="left">0.02716</td>
</tr>
<tr class="even">
<td>107</td>
<td align="left"><a href="GO:0030554" class="uri">GO:0030554</a></td>
<td align="left">adenyl nucleotide binding</td>
<td align="right">1157</td>
<td align="right">283</td>
<td align="right">257.96</td>
<td align="left">0.02836</td>
</tr>
<tr class="odd">
<td>108</td>
<td align="left"><a href="GO:0030674" class="uri">GO:0030674</a></td>
<td align="left">protein binding, bridging</td>
<td align="right">56</td>
<td align="right">19</td>
<td align="right">12.49</td>
<td align="left">0.03032</td>
</tr>
<tr class="even">
<td>109</td>
<td align="left"><a href="GO:0072509" class="uri">GO:0072509</a></td>
<td align="left">divalent inorganic cation transmembrane ...</td>
<td align="right">56</td>
<td align="right">19</td>
<td align="right">12.49</td>
<td align="left">0.03032</td>
</tr>
<tr class="odd">
<td>110</td>
<td align="left"><a href="GO:0016863" class="uri">GO:0016863</a></td>
<td align="left">intramolecular oxidoreductase activity, ...</td>
<td align="right">9</td>
<td align="right">5</td>
<td align="right">2.01</td>
<td align="left">0.03069</td>
</tr>
<tr class="even">
<td>111</td>
<td align="left"><a href="GO:0031683" class="uri">GO:0031683</a></td>
<td align="left">G-protein beta/gamma-subunit complex bin...</td>
<td align="right">9</td>
<td align="right">5</td>
<td align="right">2.01</td>
<td align="left">0.03069</td>
</tr>
<tr class="odd">
<td>112</td>
<td align="left"><a href="GO:0034987" class="uri">GO:0034987</a></td>
<td align="left">immunoglobulin receptor binding</td>
<td align="right">9</td>
<td align="right">5</td>
<td align="right">2.01</td>
<td align="left">0.03069</td>
</tr>
<tr class="even">
<td>113</td>
<td align="left"><a href="GO:0034988" class="uri">GO:0034988</a></td>
<td align="left">Fc-gamma receptor I complex binding</td>
<td align="right">9</td>
<td align="right">5</td>
<td align="right">2.01</td>
<td align="left">0.03069</td>
</tr>
<tr class="odd">
<td>114</td>
<td align="left"><a href="GO:0032559" class="uri">GO:0032559</a></td>
<td align="left">adenyl ribonucleotide binding</td>
<td align="right">1155</td>
<td align="right">282</td>
<td align="right">257.52</td>
<td align="left">0.03116</td>
</tr>
<tr class="even">
<td>115</td>
<td align="left"><a href="GO:0016209" class="uri">GO:0016209</a></td>
<td align="left">antioxidant activity</td>
<td align="right">49</td>
<td align="right">17</td>
<td align="right">10.92</td>
<td align="left">0.03174</td>
</tr>
<tr class="odd">
<td>116</td>
<td align="left"><a href="GO:0019201" class="uri">GO:0019201</a></td>
<td align="left">nucleotide kinase activity</td>
<td align="right">15</td>
<td align="right">7</td>
<td align="right">3.34</td>
<td align="left">0.03177</td>
</tr>
<tr class="even">
<td>117</td>
<td align="left"><a href="GO:0016860" class="uri">GO:0016860</a></td>
<td align="left">intramolecular oxidoreductase activity</td>
<td align="right">28</td>
<td align="right">11</td>
<td align="right">6.24</td>
<td align="left">0.03180</td>
</tr>
<tr class="odd">
<td>118</td>
<td align="left"><a href="GO:0035375" class="uri">GO:0035375</a></td>
<td align="left">zymogen binding</td>
<td align="right">35</td>
<td align="right">13</td>
<td align="right">7.80</td>
<td align="left">0.03294</td>
</tr>
<tr class="even">
<td>119</td>
<td align="left"><a href="GO:0004095" class="uri">GO:0004095</a></td>
<td align="left">carnitine O-palmitoyltransferase activit...</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03687</td>
</tr>
<tr class="odd">
<td>120</td>
<td align="left"><a href="GO:0016416" class="uri">GO:0016416</a></td>
<td align="left">O-palmitoyltransferase activity</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03687</td>
</tr>
<tr class="even">
<td>121</td>
<td align="left"><a href="GO:0016775" class="uri">GO:0016775</a></td>
<td align="left">phosphotransferase activity, nitrogenous...</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03687</td>
</tr>
<tr class="odd">
<td>122</td>
<td align="left"><a href="GO:0019206" class="uri">GO:0019206</a></td>
<td align="left">nucleoside kinase activity</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03687</td>
</tr>
<tr class="even">
<td>123</td>
<td align="left"><a href="GO:1901681" class="uri">GO:1901681</a></td>
<td align="left">sulfur compound binding</td>
<td align="right">65</td>
<td align="right">21</td>
<td align="right">14.49</td>
<td align="left">0.03999</td>
</tr>
<tr class="odd">
<td>124</td>
<td align="left"><a href="GO:0004601" class="uri">GO:0004601</a></td>
<td align="left">peroxidase activity</td>
<td align="right">36</td>
<td align="right">13</td>
<td align="right">8.03</td>
<td align="left">0.04136</td>
</tr>
<tr class="even">
<td>125</td>
<td align="left"><a href="GO:0008201" class="uri">GO:0008201</a></td>
<td align="left">heparin binding</td>
<td align="right">36</td>
<td align="right">13</td>
<td align="right">8.03</td>
<td align="left">0.04136</td>
</tr>
<tr class="odd">
<td>126</td>
<td align="left"><a href="GO:0016684" class="uri">GO:0016684</a></td>
<td align="left">oxidoreductase activity, acting on perox...</td>
<td align="right">36</td>
<td align="right">13</td>
<td align="right">8.03</td>
<td align="left">0.04136</td>
</tr>
<tr class="even">
<td>127</td>
<td align="left"><a href="GO:0038024" class="uri">GO:0038024</a></td>
<td align="left">cargo receptor activity</td>
<td align="right">85</td>
<td align="right">26</td>
<td align="right">18.95</td>
<td align="left">0.04648</td>
</tr>
<tr class="odd">
<td>128</td>
<td align="left"><a href="GO:0004075" class="uri">GO:0004075</a></td>
<td align="left">biotin carboxylase activity</td>
<td align="right">7</td>
<td align="right">4</td>
<td align="right">1.56</td>
<td align="left">0.04816</td>
</tr>
<tr class="even">
<td>129</td>
<td align="left"><a href="GO:0030552" class="uri">GO:0030552</a></td>
<td align="left">cAMP binding</td>
<td align="right">7</td>
<td align="right">4</td>
<td align="right">1.56</td>
<td align="left">0.04816</td>
</tr>
<tr class="odd">
<td>130</td>
<td align="left"><a href="GO:0000155" class="uri">GO:0000155</a></td>
<td align="left">phosphorelay sensor kinase activity</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04968</td>
</tr>
<tr class="even">
<td>131</td>
<td align="left"><a href="GO:0003878" class="uri">GO:0003878</a></td>
<td align="left">ATP citrate synthase activity</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04968</td>
</tr>
<tr class="odd">
<td>132</td>
<td align="left"><a href="GO:0003904" class="uri">GO:0003904</a></td>
<td align="left">deoxyribodipyrimidine photo-lyase activi...</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04968</td>
</tr>
<tr class="even">
<td>133</td>
<td align="left"><a href="GO:0003913" class="uri">GO:0003913</a></td>
<td align="left">DNA photolyase activity</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04968</td>
</tr>
<tr class="odd">
<td>134</td>
<td align="left"><a href="GO:0004316" class="uri">GO:0004316</a></td>
<td align="left">3-oxoacyl-[acyl-carrier-protein] reducta...</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04968</td>
</tr>
<tr class="even">
<td>135</td>
<td align="left"><a href="GO:0004421" class="uri">GO:0004421</a></td>
<td align="left">hydroxymethylglutaryl-CoA synthase activ...</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04968</td>
</tr>
<tr class="odd">
<td>136</td>
<td align="left"><a href="GO:0004466" class="uri">GO:0004466</a></td>
<td align="left">long-chain-acyl-CoA dehydrogenase activi...</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04968</td>
</tr>
<tr class="even">
<td>137</td>
<td align="left"><a href="GO:0004478" class="uri">GO:0004478</a></td>
<td align="left">methionine adenosyltransferase activity</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04968</td>
</tr>
<tr class="odd">
<td>138</td>
<td align="left"><a href="GO:0004487" class="uri">GO:0004487</a></td>
<td align="left">methylenetetrahydrofolate dehydrogenase ...</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04968</td>
</tr>
<tr class="even">
<td>139</td>
<td align="left"><a href="GO:0004499" class="uri">GO:0004499</a></td>
<td align="left">N,N-dimethylaniline monooxygenase activi...</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04968</td>
</tr>
<tr class="odd">
<td>140</td>
<td align="left"><a href="GO:0004512" class="uri">GO:0004512</a></td>
<td align="left">inositol-3-phosphate synthase activity</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04968</td>
</tr>
<tr class="even">
<td>141</td>
<td align="left"><a href="GO:0004603" class="uri">GO:0004603</a></td>
<td align="left">phenylethanolamine N-methyltransferase a...</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04968</td>
</tr>
<tr class="odd">
<td>142</td>
<td align="left"><a href="GO:0004673" class="uri">GO:0004673</a></td>
<td align="left">protein histidine kinase activity</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04968</td>
</tr>
<tr class="even">
<td>143</td>
<td align="left"><a href="GO:0008395" class="uri">GO:0008395</a></td>
<td align="left">steroid hydroxylase activity</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04968</td>
</tr>
<tr class="odd">
<td>144</td>
<td align="left"><a href="GO:0016215" class="uri">GO:0016215</a></td>
<td align="left">acyl-CoA desaturase activity</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04968</td>
</tr>
<tr class="even">
<td>145</td>
<td align="left"><a href="GO:0016936" class="uri">GO:0016936</a></td>
<td align="left">galactoside binding</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04968</td>
</tr>
<tr class="odd">
<td>146</td>
<td align="left"><a href="GO:0017099" class="uri">GO:0017099</a></td>
<td align="left">very-long-chain-acyl-CoA dehydrogenase a...</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04968</td>
</tr>
<tr class="even">
<td>147</td>
<td align="left"><a href="GO:0030160" class="uri">GO:0030160</a></td>
<td align="left">GKAP/Homer scaffold activity</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04968</td>
</tr>
<tr class="odd">
<td>148</td>
<td align="left"><a href="GO:0033925" class="uri">GO:0033925</a></td>
<td align="left">mannosyl-glycoprotein endo-beta-N-acetyl...</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04968</td>
</tr>
<tr class="even">
<td>149</td>
<td align="left"><a href="GO:0043515" class="uri">GO:0043515</a></td>
<td align="left">kinetochore binding</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04968</td>
</tr>
<tr class="odd">
<td>150</td>
<td align="left"><a href="GO:0045156" class="uri">GO:0045156</a></td>
<td align="left">electron transporter, transferring elect...</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04968</td>
</tr>
<tr class="even">
<td>151</td>
<td align="left"><a href="GO:0047237" class="uri">GO:0047237</a></td>
<td align="left">glucuronylgalactosylproteoglycan 4-beta-...</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04968</td>
</tr>
<tr class="odd">
<td>152</td>
<td align="left"><a href="GO:0055100" class="uri">GO:0055100</a></td>
<td align="left">adiponectin binding</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04968</td>
</tr>
<tr class="even">
<td>153</td>
<td align="left"><a href="GO:0061676" class="uri">GO:0061676</a></td>
<td align="left">importin-alpha family protein binding</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04968</td>
</tr>
<tr class="odd">
<td>154</td>
<td align="left"><a href="GO:0070735" class="uri">GO:0070735</a></td>
<td align="left">protein-glycine ligase activity</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04968</td>
</tr>
<tr class="even">
<td>155</td>
<td align="left"><a href="GO:0097003" class="uri">GO:0097003</a></td>
<td align="left">adipokinetic hormone receptor activity</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04968</td>
</tr>
<tr class="odd">
<td>156</td>
<td align="left"><a href="GO:0097109" class="uri">GO:0097109</a></td>
<td align="left">neuroligin family protein binding</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04968</td>
</tr>
<tr class="even">
<td>157</td>
<td align="left"><a href="GO:0098634" class="uri">GO:0098634</a></td>
<td align="left">protein binding involved in cell-matrix ...</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04968</td>
</tr>
<tr class="odd">
<td>158</td>
<td align="left"><a href="GO:0098639" class="uri">GO:0098639</a></td>
<td align="left">collagen binding involved in cell-matrix...</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04968</td>
</tr>
<tr class="even">
<td>159</td>
<td align="left"><a href="GO:2001065" class="uri">GO:2001065</a></td>
<td align="left">mannan binding</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04968</td>
</tr>
<tr class="odd">
<td>160</td>
<td align="left"><a href="GO:0015077" class="uri">GO:0015077</a></td>
<td align="left">monovalent inorganic cation transmembran...</td>
<td align="right">74</td>
<td align="right">23</td>
<td align="right">16.50</td>
<td align="left">0.04971</td>
</tr>
</tbody>
</table>

Table 9: GO-Term processes overrepresented in the set of underexpressed
DEGs.

<table>
<thead>
<tr class="header">
<th></th>
<th align="left">GO.ID</th>
<th align="left">Term</th>
<th align="right">Annotated</th>
<th align="right">Significant</th>
<th align="right">Expected</th>
<th align="left">classic</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>1</td>
<td align="left"><a href="GO:0007017" class="uri">GO:0007017</a></td>
<td align="left">microtubule-based process</td>
<td align="right">356</td>
<td align="right">181</td>
<td align="right">79.47</td>
<td align="left">&lt; 1e-30</td>
</tr>
<tr class="even">
<td>2</td>
<td align="left"><a href="GO:0007018" class="uri">GO:0007018</a></td>
<td align="left">microtubule-based movement</td>
<td align="right">167</td>
<td align="right">109</td>
<td align="right">37.28</td>
<td align="left">&lt; 1e-30</td>
</tr>
<tr class="odd">
<td>3</td>
<td align="left"><a href="GO:0060271" class="uri">GO:0060271</a></td>
<td align="left">cilium morphogenesis</td>
<td align="right">200</td>
<td align="right">120</td>
<td align="right">44.65</td>
<td align="left">&lt; 1e-30</td>
</tr>
<tr class="even">
<td>120</td>
<td align="left"><a href="GO:0006633" class="uri">GO:0006633</a></td>
<td align="left">fatty acid biosynthetic process</td>
<td align="right">58</td>
<td align="right">26</td>
<td align="right">12.95</td>
<td align="left">0.00011</td>
</tr>
<tr class="odd">
<td>121</td>
<td align="left"><a href="GO:0035108" class="uri">GO:0035108</a></td>
<td align="left">limb morphogenesis</td>
<td align="right">52</td>
<td align="right">24</td>
<td align="right">11.61</td>
<td align="left">0.00012</td>
</tr>
<tr class="even">
<td>122</td>
<td align="left"><a href="GO:0009163" class="uri">GO:0009163</a></td>
<td align="left">nucleoside biosynthetic process</td>
<td align="right">49</td>
<td align="right">23</td>
<td align="right">10.94</td>
<td align="left">0.00012</td>
</tr>
<tr class="odd">
<td>123</td>
<td align="left"><a href="GO:0060122" class="uri">GO:0060122</a></td>
<td align="left">inner ear receptor stereocilium organiza...</td>
<td align="right">29</td>
<td align="right">16</td>
<td align="right">6.47</td>
<td align="left">0.00012</td>
</tr>
<tr class="even">
<td>124</td>
<td align="left"><a href="GO:1902850" class="uri">GO:1902850</a></td>
<td align="left">microtubule cytoskeleton organization in...</td>
<td align="right">29</td>
<td align="right">16</td>
<td align="right">6.47</td>
<td align="left">0.00012</td>
</tr>
<tr class="odd">
<td>125</td>
<td align="left"><a href="GO:0048598" class="uri">GO:0048598</a></td>
<td align="left">embryonic morphogenesis</td>
<td align="right">228</td>
<td align="right">75</td>
<td align="right">50.90</td>
<td align="left">0.00012</td>
</tr>
<tr class="even">
<td>126</td>
<td align="left"><a href="GO:0019367" class="uri">GO:0019367</a></td>
<td align="left">fatty acid elongation, saturated fatty a...</td>
<td align="right">6</td>
<td align="right">6</td>
<td align="right">1.34</td>
<td align="left">0.00012</td>
</tr>
<tr class="odd">
<td>127</td>
<td align="left"><a href="GO:0046949" class="uri">GO:0046949</a></td>
<td align="left">fatty-acyl-CoA biosynthetic process</td>
<td align="right">6</td>
<td align="right">6</td>
<td align="right">1.34</td>
<td align="left">0.00012</td>
</tr>
<tr class="even">
<td>128</td>
<td align="left"><a href="GO:0007155" class="uri">GO:0007155</a></td>
<td align="left">cell adhesion</td>
<td align="right">621</td>
<td align="right">176</td>
<td align="right">138.62</td>
<td align="left">0.00013</td>
</tr>
<tr class="odd">
<td>129</td>
<td align="left"><a href="GO:0098813" class="uri">GO:0098813</a></td>
<td align="left">nuclear chromosome segregation</td>
<td align="right">104</td>
<td align="right">40</td>
<td align="right">23.22</td>
<td align="left">0.00013</td>
</tr>
<tr class="even">
<td>130</td>
<td align="left"><a href="GO:0022610" class="uri">GO:0022610</a></td>
<td align="left">biological adhesion</td>
<td align="right">622</td>
<td align="right">176</td>
<td align="right">138.85</td>
<td align="left">0.00014</td>
</tr>
<tr class="odd">
<td>131</td>
<td align="left"><a href="GO:0007140" class="uri">GO:0007140</a></td>
<td align="left">male meiosis</td>
<td align="right">24</td>
<td align="right">14</td>
<td align="right">5.36</td>
<td align="left">0.00014</td>
</tr>
<tr class="even">
<td>132</td>
<td align="left"><a href="GO:0009799" class="uri">GO:0009799</a></td>
<td align="left">specification of symmetry</td>
<td align="right">65</td>
<td align="right">28</td>
<td align="right">14.51</td>
<td align="left">0.00014</td>
</tr>
<tr class="odd">
<td>133</td>
<td align="left"><a href="GO:0009855" class="uri">GO:0009855</a></td>
<td align="left">determination of bilateral symmetry</td>
<td align="right">65</td>
<td align="right">28</td>
<td align="right">14.51</td>
<td align="left">0.00014</td>
</tr>
<tr class="even">
<td>134</td>
<td align="left"><a href="GO:0030855" class="uri">GO:0030855</a></td>
<td align="left">epithelial cell differentiation</td>
<td align="right">218</td>
<td align="right">72</td>
<td align="right">48.66</td>
<td align="left">0.00014</td>
</tr>
<tr class="odd">
<td>135</td>
<td align="left"><a href="GO:0051094" class="uri">GO:0051094</a></td>
<td align="left">positive regulation of developmental pro...</td>
<td align="right">331</td>
<td align="right">102</td>
<td align="right">73.89</td>
<td align="left">0.00015</td>
</tr>
<tr class="even">
<td>136</td>
<td align="left"><a href="GO:0055114" class="uri">GO:0055114</a></td>
<td align="left">oxidation-reduction process</td>
<td align="right">511</td>
<td align="right">148</td>
<td align="right">114.07</td>
<td align="left">0.00015</td>
</tr>
<tr class="odd">
<td>137</td>
<td align="left"><a href="GO:0051321" class="uri">GO:0051321</a></td>
<td align="left">meiotic cell cycle</td>
<td align="right">98</td>
<td align="right">38</td>
<td align="right">21.88</td>
<td align="left">0.00016</td>
</tr>
<tr class="even">
<td>138</td>
<td align="left"><a href="GO:0009790" class="uri">GO:0009790</a></td>
<td align="left">embryo development</td>
<td align="right">417</td>
<td align="right">124</td>
<td align="right">93.09</td>
<td align="left">0.00016</td>
</tr>
<tr class="odd">
<td>139</td>
<td align="left"><a href="GO:0051017" class="uri">GO:0051017</a></td>
<td align="left">actin filament bundle assembly</td>
<td align="right">41</td>
<td align="right">20</td>
<td align="right">9.15</td>
<td align="left">0.00017</td>
</tr>
<tr class="even">
<td>140</td>
<td align="left"><a href="GO:0061326" class="uri">GO:0061326</a></td>
<td align="left">renal tubule development</td>
<td align="right">41</td>
<td align="right">20</td>
<td align="right">9.15</td>
<td align="left">0.00017</td>
</tr>
<tr class="odd">
<td>141</td>
<td align="left"><a href="GO:0001656" class="uri">GO:0001656</a></td>
<td align="left">metanephros development</td>
<td align="right">50</td>
<td align="right">23</td>
<td align="right">11.16</td>
<td align="left">0.00017</td>
</tr>
<tr class="even">
<td>142</td>
<td align="left"><a href="GO:0048736" class="uri">GO:0048736</a></td>
<td align="left">appendage development</td>
<td align="right">85</td>
<td align="right">34</td>
<td align="right">18.97</td>
<td align="left">0.00017</td>
</tr>
<tr class="odd">
<td>143</td>
<td align="left"><a href="GO:0001764" class="uri">GO:0001764</a></td>
<td align="left">neuron migration</td>
<td align="right">72</td>
<td align="right">30</td>
<td align="right">16.07</td>
<td align="left">0.00017</td>
</tr>
<tr class="even">
<td>144</td>
<td align="left"><a href="GO:0001658" class="uri">GO:0001658</a></td>
<td align="left">branching involved in ureteric bud morph...</td>
<td align="right">27</td>
<td align="right">15</td>
<td align="right">6.03</td>
<td align="left">0.00017</td>
</tr>
<tr class="odd">
<td>145</td>
<td align="left"><a href="GO:0006855" class="uri">GO:0006855</a></td>
<td align="left">drug transmembrane transport</td>
<td align="right">27</td>
<td align="right">15</td>
<td align="right">6.03</td>
<td align="left">0.00017</td>
</tr>
<tr class="even">
<td>146</td>
<td align="left"><a href="GO:0090307" class="uri">GO:0090307</a></td>
<td align="left">mitotic spindle assembly</td>
<td align="right">27</td>
<td align="right">15</td>
<td align="right">6.03</td>
<td align="left">0.00017</td>
</tr>
<tr class="odd">
<td>147</td>
<td align="left"><a href="GO:0006560" class="uri">GO:0006560</a></td>
<td align="left">proline metabolic process</td>
<td align="right">8</td>
<td align="right">7</td>
<td align="right">1.79</td>
<td align="left">0.00018</td>
</tr>
<tr class="even">
<td>148</td>
<td align="left"><a href="GO:0035161" class="uri">GO:0035161</a></td>
<td align="left">imaginal disc lineage restriction</td>
<td align="right">8</td>
<td align="right">7</td>
<td align="right">1.79</td>
<td align="left">0.00018</td>
</tr>
<tr class="odd">
<td>149</td>
<td align="left"><a href="GO:0035337" class="uri">GO:0035337</a></td>
<td align="left">fatty-acyl-CoA metabolic process</td>
<td align="right">8</td>
<td align="right">7</td>
<td align="right">1.79</td>
<td align="left">0.00018</td>
</tr>
<tr class="even">
<td>150</td>
<td align="left"><a href="GO:0022402" class="uri">GO:0022402</a></td>
<td align="left">cell cycle process</td>
<td align="right">529</td>
<td align="right">152</td>
<td align="right">118.09</td>
<td align="left">0.00019</td>
</tr>
<tr class="odd">
<td>151</td>
<td align="left"><a href="GO:0072028" class="uri">GO:0072028</a></td>
<td align="left">nephron morphogenesis</td>
<td align="right">30</td>
<td align="right">16</td>
<td align="right">6.70</td>
<td align="left">0.00020</td>
</tr>
<tr class="even">
<td>152</td>
<td align="left"><a href="GO:0072078" class="uri">GO:0072078</a></td>
<td align="left">nephron tubule morphogenesis</td>
<td align="right">30</td>
<td align="right">16</td>
<td align="right">6.70</td>
<td align="left">0.00020</td>
</tr>
<tr class="odd">
<td>153</td>
<td align="left"><a href="GO:0072088" class="uri">GO:0072088</a></td>
<td align="left">nephron epithelium morphogenesis</td>
<td align="right">30</td>
<td align="right">16</td>
<td align="right">6.70</td>
<td align="left">0.00020</td>
</tr>
<tr class="even">
<td>154</td>
<td align="left"><a href="GO:0007129" class="uri">GO:0007129</a></td>
<td align="left">synapsis</td>
<td align="right">17</td>
<td align="right">11</td>
<td align="right">3.79</td>
<td align="left">0.00021</td>
</tr>
<tr class="odd">
<td>155</td>
<td align="left"><a href="GO:0007281" class="uri">GO:0007281</a></td>
<td align="left">germ cell development</td>
<td align="right">113</td>
<td align="right">42</td>
<td align="right">25.22</td>
<td align="left">0.00022</td>
</tr>
<tr class="even">
<td>156</td>
<td align="left"><a href="GO:0072080" class="uri">GO:0072080</a></td>
<td align="left">nephron tubule development</td>
<td align="right">36</td>
<td align="right">18</td>
<td align="right">8.04</td>
<td align="left">0.00024</td>
</tr>
<tr class="odd">
<td>157</td>
<td align="left"><a href="GO:0044706" class="uri">GO:0044706</a></td>
<td align="left">multi-multicellular organism process</td>
<td align="right">54</td>
<td align="right">24</td>
<td align="right">12.05</td>
<td align="left">0.00024</td>
</tr>
<tr class="even">
<td>158</td>
<td align="left"><a href="GO:0042455" class="uri">GO:0042455</a></td>
<td align="left">ribonucleoside biosynthetic process</td>
<td align="right">45</td>
<td align="right">21</td>
<td align="right">10.05</td>
<td align="left">0.00025</td>
</tr>
<tr class="odd">
<td>159</td>
<td align="left"><a href="GO:0007067" class="uri">GO:0007067</a></td>
<td align="left">mitotic nuclear division</td>
<td align="right">222</td>
<td align="right">72</td>
<td align="right">49.56</td>
<td align="left">0.00027</td>
</tr>
<tr class="even">
<td>160</td>
<td align="left"><a href="GO:0030212" class="uri">GO:0030212</a></td>
<td align="left">hyaluronan metabolic process</td>
<td align="right">15</td>
<td align="right">10</td>
<td align="right">3.35</td>
<td align="left">0.00029</td>
</tr>
<tr class="odd">
<td>161</td>
<td align="left"><a href="GO:0060675" class="uri">GO:0060675</a></td>
<td align="left">ureteric bud morphogenesis</td>
<td align="right">28</td>
<td align="right">15</td>
<td align="right">6.25</td>
<td align="left">0.00030</td>
</tr>
<tr class="even">
<td>162</td>
<td align="left"><a href="GO:0072171" class="uri">GO:0072171</a></td>
<td align="left">mesonephric tubule morphogenesis</td>
<td align="right">28</td>
<td align="right">15</td>
<td align="right">6.25</td>
<td align="left">0.00030</td>
</tr>
<tr class="odd">
<td>163</td>
<td align="left"><a href="GO:0010466" class="uri">GO:0010466</a></td>
<td align="left">negative regulation of peptidase activit...</td>
<td align="right">61</td>
<td align="right">26</td>
<td align="right">13.62</td>
<td align="left">0.00030</td>
</tr>
<tr class="even">
<td>164</td>
<td align="left"><a href="GO:0072330" class="uri">GO:0072330</a></td>
<td align="left">monocarboxylic acid biosynthetic process</td>
<td align="right">74</td>
<td align="right">30</td>
<td align="right">16.52</td>
<td align="left">0.00031</td>
</tr>
<tr class="odd">
<td>165</td>
<td align="left"><a href="GO:0007423" class="uri">GO:0007423</a></td>
<td align="left">sensory organ development</td>
<td align="right">238</td>
<td align="right">76</td>
<td align="right">53.13</td>
<td align="left">0.00031</td>
</tr>
<tr class="even">
<td>166</td>
<td align="left"><a href="GO:0060173" class="uri">GO:0060173</a></td>
<td align="left">limb development</td>
<td align="right">58</td>
<td align="right">25</td>
<td align="right">12.95</td>
<td align="left">0.00032</td>
</tr>
<tr class="odd">
<td>167</td>
<td align="left"><a href="GO:0001708" class="uri">GO:0001708</a></td>
<td align="left">cell fate specification</td>
<td align="right">31</td>
<td align="right">16</td>
<td align="right">6.92</td>
<td align="left">0.00033</td>
</tr>
<tr class="even">
<td>168</td>
<td align="left"><a href="GO:0060119" class="uri">GO:0060119</a></td>
<td align="left">inner ear receptor cell development</td>
<td align="right">31</td>
<td align="right">16</td>
<td align="right">6.92</td>
<td align="left">0.00033</td>
</tr>
<tr class="odd">
<td>169</td>
<td align="left"><a href="GO:0046165" class="uri">GO:0046165</a></td>
<td align="left">alcohol biosynthetic process</td>
<td align="right">49</td>
<td align="right">22</td>
<td align="right">10.94</td>
<td align="left">0.00036</td>
</tr>
<tr class="even">
<td>170</td>
<td align="left"><a href="GO:0051297" class="uri">GO:0051297</a></td>
<td align="left">centrosome organization</td>
<td align="right">46</td>
<td align="right">21</td>
<td align="right">10.27</td>
<td align="left">0.00037</td>
</tr>
<tr class="odd">
<td>171</td>
<td align="left"><a href="GO:0042451" class="uri">GO:0042451</a></td>
<td align="left">purine nucleoside biosynthetic process</td>
<td align="right">40</td>
<td align="right">19</td>
<td align="right">8.93</td>
<td align="left">0.00037</td>
</tr>
<tr class="even">
<td>172</td>
<td align="left"><a href="GO:0061572" class="uri">GO:0061572</a></td>
<td align="left">actin filament bundle organization</td>
<td align="right">43</td>
<td align="right">20</td>
<td align="right">9.60</td>
<td align="left">0.00037</td>
</tr>
<tr class="odd">
<td>173</td>
<td align="left"><a href="GO:1901617" class="uri">GO:1901617</a></td>
<td align="left">organic hydroxy compound biosynthetic pr...</td>
<td align="right">78</td>
<td align="right">31</td>
<td align="right">17.41</td>
<td align="left">0.00037</td>
</tr>
<tr class="even">
<td>174</td>
<td align="left"><a href="GO:0051963" class="uri">GO:0051963</a></td>
<td align="left">regulation of synapse assembly</td>
<td align="right">65</td>
<td align="right">27</td>
<td align="right">14.51</td>
<td align="left">0.00038</td>
</tr>
<tr class="odd">
<td>175</td>
<td align="left"><a href="GO:0021517" class="uri">GO:0021517</a></td>
<td align="left">ventral spinal cord development</td>
<td align="right">13</td>
<td align="right">9</td>
<td align="right">2.90</td>
<td align="left">0.00040</td>
</tr>
<tr class="even">
<td>176</td>
<td align="left"><a href="GO:0051301" class="uri">GO:0051301</a></td>
<td align="left">cell division</td>
<td align="right">320</td>
<td align="right">97</td>
<td align="right">71.43</td>
<td align="left">0.00041</td>
</tr>
<tr class="odd">
<td>177</td>
<td align="left"><a href="GO:0048565" class="uri">GO:0048565</a></td>
<td align="left">digestive tract development</td>
<td align="right">59</td>
<td align="right">25</td>
<td align="right">13.17</td>
<td align="left">0.00044</td>
</tr>
<tr class="even">
<td>178</td>
<td align="left"><a href="GO:0009612" class="uri">GO:0009612</a></td>
<td align="left">response to mechanical stimulus</td>
<td align="right">82</td>
<td align="right">32</td>
<td align="right">18.30</td>
<td align="left">0.00045</td>
</tr>
<tr class="odd">
<td>179</td>
<td align="left"><a href="GO:0006631" class="uri">GO:0006631</a></td>
<td align="left">fatty acid metabolic process</td>
<td align="right">134</td>
<td align="right">47</td>
<td align="right">29.91</td>
<td align="left">0.00045</td>
</tr>
<tr class="even">
<td>180</td>
<td align="left"><a href="GO:0022412" class="uri">GO:0022412</a></td>
<td align="left">cellular process involved in reproductio...</td>
<td align="right">134</td>
<td align="right">47</td>
<td align="right">29.91</td>
<td align="left">0.00045</td>
</tr>
<tr class="odd">
<td>181</td>
<td align="left"><a href="GO:0072001" class="uri">GO:0072001</a></td>
<td align="left">renal system development</td>
<td align="right">145</td>
<td align="right">50</td>
<td align="right">32.37</td>
<td align="left">0.00048</td>
</tr>
<tr class="even">
<td>182</td>
<td align="left"><a href="GO:0015893" class="uri">GO:0015893</a></td>
<td align="left">drug transport</td>
<td align="right">29</td>
<td align="right">15</td>
<td align="right">6.47</td>
<td align="left">0.00049</td>
</tr>
<tr class="odd">
<td>183</td>
<td align="left"><a href="GO:0009698" class="uri">GO:0009698</a></td>
<td align="left">phenylpropanoid metabolic process</td>
<td align="right">11</td>
<td align="right">8</td>
<td align="right">2.46</td>
<td align="left">0.00052</td>
</tr>
<tr class="even">
<td>184</td>
<td align="left"><a href="GO:0030705" class="uri">GO:0030705</a></td>
<td align="left">cytoskeleton-dependent intracellular tra...</td>
<td align="right">47</td>
<td align="right">21</td>
<td align="right">10.49</td>
<td align="left">0.00053</td>
</tr>
<tr class="odd">
<td>185</td>
<td align="left"><a href="GO:0061005" class="uri">GO:0061005</a></td>
<td align="left">cell differentiation involved in kidney ...</td>
<td align="right">35</td>
<td align="right">17</td>
<td align="right">7.81</td>
<td align="left">0.00054</td>
</tr>
<tr class="even">
<td>186</td>
<td align="left"><a href="GO:0071695" class="uri">GO:0071695</a></td>
<td align="left">anatomical structure maturation</td>
<td align="right">35</td>
<td align="right">17</td>
<td align="right">7.81</td>
<td align="left">0.00054</td>
</tr>
<tr class="odd">
<td>187</td>
<td align="left"><a href="GO:0030326" class="uri">GO:0030326</a></td>
<td align="left">embryonic limb morphogenesis</td>
<td align="right">41</td>
<td align="right">19</td>
<td align="right">9.15</td>
<td align="left">0.00055</td>
</tr>
<tr class="even">
<td>188</td>
<td align="left"><a href="GO:0060993" class="uri">GO:0060993</a></td>
<td align="left">kidney morphogenesis</td>
<td align="right">63</td>
<td align="right">26</td>
<td align="right">14.06</td>
<td align="left">0.00055</td>
</tr>
<tr class="odd">
<td>189</td>
<td align="left"><a href="GO:0051026" class="uri">GO:0051026</a></td>
<td align="left">chiasma assembly</td>
<td align="right">5</td>
<td align="right">5</td>
<td align="right">1.12</td>
<td align="left">0.00055</td>
</tr>
<tr class="even">
<td>190</td>
<td align="left"><a href="GO:1904158" class="uri">GO:1904158</a></td>
<td align="left">axonemal central apparatus assembly</td>
<td align="right">5</td>
<td align="right">5</td>
<td align="right">1.12</td>
<td align="left">0.00055</td>
</tr>
<tr class="odd">
<td>191</td>
<td align="left"><a href="GO:0048839" class="uri">GO:0048839</a></td>
<td align="left">inner ear development</td>
<td align="right">107</td>
<td align="right">39</td>
<td align="right">23.89</td>
<td align="left">0.00057</td>
</tr>
<tr class="even">
<td>192</td>
<td align="left"><a href="GO:0032504" class="uri">GO:0032504</a></td>
<td align="left">multicellular organism reproduction</td>
<td align="right">273</td>
<td align="right">84</td>
<td align="right">60.94</td>
<td align="left">0.00059</td>
</tr>
<tr class="odd">
<td>193</td>
<td align="left"><a href="GO:0000038" class="uri">GO:0000038</a></td>
<td align="left">very long-chain fatty acid metabolic pro...</td>
<td align="right">16</td>
<td align="right">10</td>
<td align="right">3.57</td>
<td align="left">0.00063</td>
</tr>
<tr class="even">
<td>194</td>
<td align="left"><a href="GO:0043583" class="uri">GO:0043583</a></td>
<td align="left">ear development</td>
<td align="right">111</td>
<td align="right">40</td>
<td align="right">24.78</td>
<td align="left">0.00064</td>
</tr>
<tr class="odd">
<td>195</td>
<td align="left"><a href="GO:0021511" class="uri">GO:0021511</a></td>
<td align="left">spinal cord patterning</td>
<td align="right">9</td>
<td align="right">7</td>
<td align="right">2.01</td>
<td align="left">0.00064</td>
</tr>
<tr class="even">
<td>196</td>
<td align="left"><a href="GO:0035222" class="uri">GO:0035222</a></td>
<td align="left">wing disc pattern formation</td>
<td align="right">9</td>
<td align="right">7</td>
<td align="right">2.01</td>
<td align="left">0.00064</td>
</tr>
<tr class="odd">
<td>197</td>
<td align="left"><a href="GO:0044458" class="uri">GO:0044458</a></td>
<td align="left">motile cilium assembly</td>
<td align="right">9</td>
<td align="right">7</td>
<td align="right">2.01</td>
<td align="left">0.00064</td>
</tr>
<tr class="even">
<td>198</td>
<td align="left"><a href="GO:0048609" class="uri">GO:0048609</a></td>
<td align="left">multicellular organismal reproductive pr...</td>
<td align="right">270</td>
<td align="right">83</td>
<td align="right">60.27</td>
<td align="left">0.00066</td>
</tr>
<tr class="odd">
<td>199</td>
<td align="left"><a href="GO:0007052" class="uri">GO:0007052</a></td>
<td align="left">mitotic spindle organization</td>
<td align="right">54</td>
<td align="right">23</td>
<td align="right">12.05</td>
<td align="left">0.00067</td>
</tr>
<tr class="even">
<td>200</td>
<td align="left"><a href="GO:0007051" class="uri">GO:0007051</a></td>
<td align="left">spindle organization</td>
<td align="right">77</td>
<td align="right">30</td>
<td align="right">17.19</td>
<td align="left">0.00069</td>
</tr>
<tr class="odd">
<td>201</td>
<td align="left"><a href="GO:0006561" class="uri">GO:0006561</a></td>
<td align="left">proline biosynthetic process</td>
<td align="right">7</td>
<td align="right">6</td>
<td align="right">1.56</td>
<td align="left">0.00070</td>
</tr>
<tr class="even">
<td>202</td>
<td align="left"><a href="GO:0007450" class="uri">GO:0007450</a></td>
<td align="left">dorsal/ventral pattern formation, imagin...</td>
<td align="right">7</td>
<td align="right">6</td>
<td align="right">1.56</td>
<td align="left">0.00070</td>
</tr>
<tr class="odd">
<td>203</td>
<td align="left"><a href="GO:0007451" class="uri">GO:0007451</a></td>
<td align="left">dorsal/ventral lineage restriction, imag...</td>
<td align="right">7</td>
<td align="right">6</td>
<td align="right">1.56</td>
<td align="left">0.00070</td>
</tr>
<tr class="even">
<td>204</td>
<td align="left"><a href="GO:0010457" class="uri">GO:0010457</a></td>
<td align="left">centriole-centriole cohesion</td>
<td align="right">7</td>
<td align="right">6</td>
<td align="right">1.56</td>
<td align="left">0.00070</td>
</tr>
<tr class="odd">
<td>205</td>
<td align="left"><a href="GO:0035285" class="uri">GO:0035285</a></td>
<td align="left">appendage segmentation</td>
<td align="right">7</td>
<td align="right">6</td>
<td align="right">1.56</td>
<td align="left">0.00070</td>
</tr>
<tr class="even">
<td>206</td>
<td align="left"><a href="GO:0036011" class="uri">GO:0036011</a></td>
<td align="left">imaginal disc-derived leg segmentation</td>
<td align="right">7</td>
<td align="right">6</td>
<td align="right">1.56</td>
<td align="left">0.00070</td>
</tr>
<tr class="odd">
<td>207</td>
<td align="left"><a href="GO:0048190" class="uri">GO:0048190</a></td>
<td align="left">wing disc dorsal/ventral pattern formati...</td>
<td align="right">7</td>
<td align="right">6</td>
<td align="right">1.56</td>
<td align="left">0.00070</td>
</tr>
<tr class="even">
<td>208</td>
<td align="left"><a href="GO:0061525" class="uri">GO:0061525</a></td>
<td align="left">hindgut development</td>
<td align="right">7</td>
<td align="right">6</td>
<td align="right">1.56</td>
<td align="left">0.00070</td>
</tr>
<tr class="odd">
<td>209</td>
<td align="left"><a href="GO:0003002" class="uri">GO:0003002</a></td>
<td align="left">regionalization</td>
<td align="right">133</td>
<td align="right">46</td>
<td align="right">29.69</td>
<td align="left">0.00074</td>
</tr>
<tr class="even">
<td>210</td>
<td align="left"><a href="GO:0007157" class="uri">GO:0007157</a></td>
<td align="left">heterophilic cell-cell adhesion via plas...</td>
<td align="right">42</td>
<td align="right">19</td>
<td align="right">9.38</td>
<td align="left">0.00080</td>
</tr>
<tr class="odd">
<td>211</td>
<td align="left"><a href="GO:0051225" class="uri">GO:0051225</a></td>
<td align="left">spindle assembly</td>
<td align="right">42</td>
<td align="right">19</td>
<td align="right">9.38</td>
<td align="left">0.00080</td>
</tr>
<tr class="even">
<td>212</td>
<td align="left"><a href="GO:0045132" class="uri">GO:0045132</a></td>
<td align="left">meiotic chromosome segregation</td>
<td align="right">33</td>
<td align="right">16</td>
<td align="right">7.37</td>
<td align="left">0.00081</td>
</tr>
<tr class="odd">
<td>213</td>
<td align="left"><a href="GO:0001657" class="uri">GO:0001657</a></td>
<td align="left">ureteric bud development</td>
<td align="right">39</td>
<td align="right">18</td>
<td align="right">8.71</td>
<td align="left">0.00081</td>
</tr>
<tr class="even">
<td>214</td>
<td align="left"><a href="GO:0046129" class="uri">GO:0046129</a></td>
<td align="left">purine ribonucleoside biosynthetic proce...</td>
<td align="right">39</td>
<td align="right">18</td>
<td align="right">8.71</td>
<td align="left">0.00081</td>
</tr>
<tr class="odd">
<td>215</td>
<td align="left"><a href="GO:0016525" class="uri">GO:0016525</a></td>
<td align="left">negative regulation of angiogenesis</td>
<td align="right">14</td>
<td align="right">9</td>
<td align="right">3.13</td>
<td align="left">0.00089</td>
</tr>
<tr class="even">
<td>216</td>
<td align="left"><a href="GO:0050793" class="uri">GO:0050793</a></td>
<td align="left">regulation of developmental process</td>
<td align="right">726</td>
<td align="right">196</td>
<td align="right">162.06</td>
<td align="left">0.00094</td>
</tr>
<tr class="odd">
<td>217</td>
<td align="left"><a href="GO:0050953" class="uri">GO:0050953</a></td>
<td align="left">sensory perception of light stimulus</td>
<td align="right">75</td>
<td align="right">29</td>
<td align="right">16.74</td>
<td align="left">0.00097</td>
</tr>
<tr class="even">
<td>218</td>
<td align="left"><a href="GO:0010842" class="uri">GO:0010842</a></td>
<td align="left">retina layer formation</td>
<td align="right">22</td>
<td align="right">12</td>
<td align="right">4.91</td>
<td align="left">0.00098</td>
</tr>
<tr class="odd">
<td>219</td>
<td align="left"><a href="GO:0060291" class="uri">GO:0060291</a></td>
<td align="left">long-term synaptic potentiation</td>
<td align="right">22</td>
<td align="right">12</td>
<td align="right">4.91</td>
<td align="left">0.00098</td>
</tr>
<tr class="even">
<td>220</td>
<td align="left"><a href="GO:0031023" class="uri">GO:0031023</a></td>
<td align="left">microtubule organizing center organizati...</td>
<td align="right">49</td>
<td align="right">21</td>
<td align="right">10.94</td>
<td align="left">0.00103</td>
</tr>
<tr class="odd">
<td>221</td>
<td align="left"><a href="GO:0030856" class="uri">GO:0030856</a></td>
<td align="left">regulation of epithelial cell differenti...</td>
<td align="right">62</td>
<td align="right">25</td>
<td align="right">13.84</td>
<td align="left">0.00105</td>
</tr>
<tr class="even">
<td>222</td>
<td align="left"><a href="GO:0035113" class="uri">GO:0035113</a></td>
<td align="left">embryonic appendage morphogenesis</td>
<td align="right">46</td>
<td align="right">20</td>
<td align="right">10.27</td>
<td align="left">0.00108</td>
</tr>
<tr class="odd">
<td>223</td>
<td align="left"><a href="GO:0007613" class="uri">GO:0007613</a></td>
<td align="left">memory</td>
<td align="right">43</td>
<td align="right">19</td>
<td align="right">9.60</td>
<td align="left">0.00113</td>
</tr>
<tr class="even">
<td>224</td>
<td align="left"><a href="GO:0001823" class="uri">GO:0001823</a></td>
<td align="left">mesonephros development</td>
<td align="right">40</td>
<td align="right">18</td>
<td align="right">8.93</td>
<td align="left">0.00118</td>
</tr>
<tr class="odd">
<td>225</td>
<td align="left"><a href="GO:0072163" class="uri">GO:0072163</a></td>
<td align="left">mesonephric epithelium development</td>
<td align="right">40</td>
<td align="right">18</td>
<td align="right">8.93</td>
<td align="left">0.00118</td>
</tr>
<tr class="even">
<td>226</td>
<td align="left"><a href="GO:0072164" class="uri">GO:0072164</a></td>
<td align="left">mesonephric tubule development</td>
<td align="right">40</td>
<td align="right">18</td>
<td align="right">8.93</td>
<td align="left">0.00118</td>
</tr>
<tr class="odd">
<td>227</td>
<td align="left"><a href="GO:0007098" class="uri">GO:0007098</a></td>
<td align="left">centrosome cycle</td>
<td align="right">34</td>
<td align="right">16</td>
<td align="right">7.59</td>
<td align="left">0.00122</td>
</tr>
<tr class="even">
<td>228</td>
<td align="left"><a href="GO:0009142" class="uri">GO:0009142</a></td>
<td align="left">nucleoside triphosphate biosynthetic pro...</td>
<td align="right">34</td>
<td align="right">16</td>
<td align="right">7.59</td>
<td align="left">0.00122</td>
</tr>
<tr class="odd">
<td>229</td>
<td align="left"><a href="GO:0060041" class="uri">GO:0060041</a></td>
<td align="left">retina development in camera-type eye</td>
<td align="right">56</td>
<td align="right">23</td>
<td align="right">12.50</td>
<td align="left">0.00122</td>
</tr>
<tr class="even">
<td>230</td>
<td align="left"><a href="GO:0021532" class="uri">GO:0021532</a></td>
<td align="left">neural tube patterning</td>
<td align="right">12</td>
<td align="right">8</td>
<td align="right">2.68</td>
<td align="left">0.00125</td>
</tr>
<tr class="odd">
<td>231</td>
<td align="left"><a href="GO:0046605" class="uri">GO:0046605</a></td>
<td align="left">regulation of centrosome cycle</td>
<td align="right">12</td>
<td align="right">8</td>
<td align="right">2.68</td>
<td align="left">0.00125</td>
</tr>
<tr class="even">
<td>232</td>
<td align="left"><a href="GO:0001822" class="uri">GO:0001822</a></td>
<td align="left">kidney development</td>
<td align="right">136</td>
<td align="right">46</td>
<td align="right">30.36</td>
<td align="left">0.00126</td>
</tr>
<tr class="odd">
<td>233</td>
<td align="left"><a href="GO:0055123" class="uri">GO:0055123</a></td>
<td align="left">digestive system development</td>
<td align="right">63</td>
<td align="right">25</td>
<td align="right">14.06</td>
<td align="left">0.00138</td>
</tr>
<tr class="even">
<td>234</td>
<td align="left"><a href="GO:0001655" class="uri">GO:0001655</a></td>
<td align="left">urogenital system development</td>
<td align="right">155</td>
<td align="right">51</td>
<td align="right">34.60</td>
<td align="left">0.00142</td>
</tr>
<tr class="odd">
<td>235</td>
<td align="left"><a href="GO:0031987" class="uri">GO:0031987</a></td>
<td align="left">locomotion involved in locomotory behavi...</td>
<td align="right">20</td>
<td align="right">11</td>
<td align="right">4.46</td>
<td align="left">0.00147</td>
</tr>
<tr class="even">
<td>236</td>
<td align="left"><a href="GO:0072073" class="uri">GO:0072073</a></td>
<td align="left">kidney epithelium development</td>
<td align="right">60</td>
<td align="right">24</td>
<td align="right">13.39</td>
<td align="left">0.00150</td>
</tr>
<tr class="odd">
<td>237</td>
<td align="left"><a href="GO:1903047" class="uri">GO:1903047</a></td>
<td align="left">mitotic cell cycle process</td>
<td align="right">395</td>
<td align="right">113</td>
<td align="right">88.17</td>
<td align="left">0.00153</td>
</tr>
<tr class="even">
<td>238</td>
<td align="left"><a href="GO:0007411" class="uri">GO:0007411</a></td>
<td align="left">axon guidance</td>
<td align="right">178</td>
<td align="right">57</td>
<td align="right">39.73</td>
<td align="left">0.00158</td>
</tr>
<tr class="odd">
<td>239</td>
<td align="left"><a href="GO:0097485" class="uri">GO:0097485</a></td>
<td align="left">neuron projection guidance</td>
<td align="right">178</td>
<td align="right">57</td>
<td align="right">39.73</td>
<td align="left">0.00158</td>
</tr>
<tr class="even">
<td>240</td>
<td align="left"><a href="GO:0009953" class="uri">GO:0009953</a></td>
<td align="left">dorsal/ventral pattern formation</td>
<td align="right">44</td>
<td align="right">19</td>
<td align="right">9.82</td>
<td align="left">0.00158</td>
</tr>
<tr class="odd">
<td>241</td>
<td align="left"><a href="GO:0016053" class="uri">GO:0016053</a></td>
<td align="left">organic acid biosynthetic process</td>
<td align="right">141</td>
<td align="right">47</td>
<td align="right">31.47</td>
<td align="left">0.00159</td>
</tr>
<tr class="even">
<td>242</td>
<td align="left"><a href="GO:0046394" class="uri">GO:0046394</a></td>
<td align="left">carboxylic acid biosynthetic process</td>
<td align="right">141</td>
<td align="right">47</td>
<td align="right">31.47</td>
<td align="left">0.00159</td>
</tr>
<tr class="odd">
<td>243</td>
<td align="left"><a href="GO:0060541" class="uri">GO:0060541</a></td>
<td align="left">respiratory system development</td>
<td align="right">84</td>
<td align="right">31</td>
<td align="right">18.75</td>
<td align="left">0.00161</td>
</tr>
<tr class="even">
<td>244</td>
<td align="left"><a href="GO:0042220" class="uri">GO:0042220</a></td>
<td align="left">response to cocaine</td>
<td align="right">23</td>
<td align="right">12</td>
<td align="right">5.13</td>
<td align="left">0.00164</td>
</tr>
<tr class="odd">
<td>245</td>
<td align="left"><a href="GO:0010970" class="uri">GO:0010970</a></td>
<td align="left">microtubule-based transport</td>
<td align="right">41</td>
<td align="right">18</td>
<td align="right">9.15</td>
<td align="left">0.00166</td>
</tr>
<tr class="even">
<td>246</td>
<td align="left"><a href="GO:0042490" class="uri">GO:0042490</a></td>
<td align="left">mechanoreceptor differentiation</td>
<td align="right">41</td>
<td align="right">18</td>
<td align="right">9.15</td>
<td align="left">0.00166</td>
</tr>
<tr class="odd">
<td>247</td>
<td align="left"><a href="GO:0060113" class="uri">GO:0060113</a></td>
<td align="left">inner ear receptor cell differentiation</td>
<td align="right">41</td>
<td align="right">18</td>
<td align="right">9.15</td>
<td align="left">0.00166</td>
</tr>
<tr class="even">
<td>248</td>
<td align="left"><a href="GO:0006883" class="uri">GO:0006883</a></td>
<td align="left">cellular sodium ion homeostasis</td>
<td align="right">10</td>
<td align="right">7</td>
<td align="right">2.23</td>
<td align="left">0.00172</td>
</tr>
<tr class="odd">
<td>249</td>
<td align="left"><a href="GO:0009220" class="uri">GO:0009220</a></td>
<td align="left">pyrimidine ribonucleotide biosynthetic p...</td>
<td align="right">10</td>
<td align="right">7</td>
<td align="right">2.23</td>
<td align="left">0.00172</td>
</tr>
<tr class="even">
<td>250</td>
<td align="left"><a href="GO:0035115" class="uri">GO:0035115</a></td>
<td align="left">embryonic forelimb morphogenesis</td>
<td align="right">10</td>
<td align="right">7</td>
<td align="right">2.23</td>
<td align="left">0.00172</td>
</tr>
<tr class="odd">
<td>251</td>
<td align="left"><a href="GO:0060438" class="uri">GO:0060438</a></td>
<td align="left">trachea development</td>
<td align="right">10</td>
<td align="right">7</td>
<td align="right">2.23</td>
<td align="left">0.00172</td>
</tr>
<tr class="even">
<td>252</td>
<td align="left"><a href="GO:0048016" class="uri">GO:0048016</a></td>
<td align="left">inositol phosphate-mediated signaling</td>
<td align="right">15</td>
<td align="right">9</td>
<td align="right">3.35</td>
<td align="left">0.00179</td>
</tr>
<tr class="odd">
<td>253</td>
<td align="left"><a href="GO:1901343" class="uri">GO:1901343</a></td>
<td align="left">negative regulation of vasculature devel...</td>
<td align="right">15</td>
<td align="right">9</td>
<td align="right">3.35</td>
<td align="left">0.00179</td>
</tr>
<tr class="even">
<td>254</td>
<td align="left"><a href="GO:2000181" class="uri">GO:2000181</a></td>
<td align="left">negative regulation of blood vessel morp...</td>
<td align="right">15</td>
<td align="right">9</td>
<td align="right">3.35</td>
<td align="left">0.00179</td>
</tr>
<tr class="odd">
<td>255</td>
<td align="left"><a href="GO:0006637" class="uri">GO:0006637</a></td>
<td align="left">acyl-CoA metabolic process</td>
<td align="right">29</td>
<td align="right">14</td>
<td align="right">6.47</td>
<td align="left">0.00180</td>
</tr>
<tr class="even">
<td>256</td>
<td align="left"><a href="GO:0035383" class="uri">GO:0035383</a></td>
<td align="left">thioester metabolic process</td>
<td align="right">29</td>
<td align="right">14</td>
<td align="right">6.47</td>
<td align="left">0.00180</td>
</tr>
<tr class="odd">
<td>257</td>
<td align="left"><a href="GO:0035850" class="uri">GO:0035850</a></td>
<td align="left">epithelial cell differentiation involved...</td>
<td align="right">32</td>
<td align="right">15</td>
<td align="right">7.14</td>
<td align="left">0.00181</td>
</tr>
<tr class="even">
<td>258</td>
<td align="left"><a href="GO:0042493" class="uri">GO:0042493</a></td>
<td align="left">response to drug</td>
<td align="right">131</td>
<td align="right">44</td>
<td align="right">29.24</td>
<td align="left">0.00186</td>
</tr>
<tr class="odd">
<td>259</td>
<td align="left"><a href="GO:1902652" class="uri">GO:1902652</a></td>
<td align="left">secondary alcohol metabolic process</td>
<td align="right">61</td>
<td align="right">24</td>
<td align="right">13.62</td>
<td align="left">0.00196</td>
</tr>
<tr class="even">
<td>260</td>
<td align="left"><a href="GO:0044283" class="uri">GO:0044283</a></td>
<td align="left">small molecule biosynthetic process</td>
<td align="right">233</td>
<td align="right">71</td>
<td align="right">52.01</td>
<td align="left">0.00203</td>
</tr>
<tr class="odd">
<td>261</td>
<td align="left"><a href="GO:0000083" class="uri">GO:0000083</a></td>
<td align="left">regulation of transcription involved in ...</td>
<td align="right">8</td>
<td align="right">6</td>
<td align="right">1.79</td>
<td align="left">0.00225</td>
</tr>
<tr class="even">
<td>262</td>
<td align="left"><a href="GO:0007478" class="uri">GO:0007478</a></td>
<td align="left">leg disc morphogenesis</td>
<td align="right">8</td>
<td align="right">6</td>
<td align="right">1.79</td>
<td align="left">0.00225</td>
</tr>
<tr class="odd">
<td>263</td>
<td align="left"><a href="GO:0007480" class="uri">GO:0007480</a></td>
<td align="left">imaginal disc-derived leg morphogenesis</td>
<td align="right">8</td>
<td align="right">6</td>
<td align="right">1.79</td>
<td align="left">0.00225</td>
</tr>
<tr class="even">
<td>264</td>
<td align="left"><a href="GO:0008587" class="uri">GO:0008587</a></td>
<td align="left">imaginal disc-derived wing margin morpho...</td>
<td align="right">8</td>
<td align="right">6</td>
<td align="right">1.79</td>
<td align="left">0.00225</td>
</tr>
<tr class="odd">
<td>265</td>
<td align="left"><a href="GO:0009084" class="uri">GO:0009084</a></td>
<td align="left">glutamine family amino acid biosynthetic...</td>
<td align="right">8</td>
<td align="right">6</td>
<td align="right">1.79</td>
<td align="left">0.00225</td>
</tr>
<tr class="even">
<td>266</td>
<td align="left"><a href="GO:0009808" class="uri">GO:0009808</a></td>
<td align="left">lignin metabolic process</td>
<td align="right">8</td>
<td align="right">6</td>
<td align="right">1.79</td>
<td align="left">0.00225</td>
</tr>
<tr class="odd">
<td>267</td>
<td align="left"><a href="GO:0035218" class="uri">GO:0035218</a></td>
<td align="left">leg disc development</td>
<td align="right">8</td>
<td align="right">6</td>
<td align="right">1.79</td>
<td align="left">0.00225</td>
</tr>
<tr class="even">
<td>268</td>
<td align="left"><a href="GO:0035721" class="uri">GO:0035721</a></td>
<td align="left">intraciliary retrograde transport</td>
<td align="right">8</td>
<td align="right">6</td>
<td align="right">1.79</td>
<td align="left">0.00225</td>
</tr>
<tr class="odd">
<td>269</td>
<td align="left"><a href="GO:0046132" class="uri">GO:0046132</a></td>
<td align="left">pyrimidine ribonucleoside biosynthetic p...</td>
<td align="right">8</td>
<td align="right">6</td>
<td align="right">1.79</td>
<td align="left">0.00225</td>
</tr>
<tr class="even">
<td>270</td>
<td align="left"><a href="GO:0051593" class="uri">GO:0051593</a></td>
<td align="left">response to folic acid</td>
<td align="right">8</td>
<td align="right">6</td>
<td align="right">1.79</td>
<td align="left">0.00225</td>
</tr>
<tr class="odd">
<td>271</td>
<td align="left"><a href="GO:0071436" class="uri">GO:0071436</a></td>
<td align="left">sodium ion export</td>
<td align="right">8</td>
<td align="right">6</td>
<td align="right">1.79</td>
<td align="left">0.00225</td>
</tr>
<tr class="even">
<td>272</td>
<td align="left"><a href="GO:1901660" class="uri">GO:1901660</a></td>
<td align="left">calcium ion export</td>
<td align="right">8</td>
<td align="right">6</td>
<td align="right">1.79</td>
<td align="left">0.00225</td>
</tr>
<tr class="odd">
<td>273</td>
<td align="left"><a href="GO:0061371" class="uri">GO:0061371</a></td>
<td align="left">determination of heart left/right asymme...</td>
<td align="right">21</td>
<td align="right">11</td>
<td align="right">4.69</td>
<td align="left">0.00247</td>
</tr>
<tr class="even">
<td>274</td>
<td align="left"><a href="GO:0072207" class="uri">GO:0072207</a></td>
<td align="left">metanephric epithelium development</td>
<td align="right">21</td>
<td align="right">11</td>
<td align="right">4.69</td>
<td align="left">0.00247</td>
</tr>
<tr class="odd">
<td>275</td>
<td align="left"><a href="GO:0072243" class="uri">GO:0072243</a></td>
<td align="left">metanephric nephron epithelium developme...</td>
<td align="right">21</td>
<td align="right">11</td>
<td align="right">4.69</td>
<td align="left">0.00247</td>
</tr>
<tr class="even">
<td>276</td>
<td align="left"><a href="GO:0005979" class="uri">GO:0005979</a></td>
<td align="left">regulation of glycogen biosynthetic proc...</td>
<td align="right">4</td>
<td align="right">4</td>
<td align="right">0.89</td>
<td align="left">0.00247</td>
</tr>
<tr class="odd">
<td>277</td>
<td align="left"><a href="GO:0006183" class="uri">GO:0006183</a></td>
<td align="left">GTP biosynthetic process</td>
<td align="right">4</td>
<td align="right">4</td>
<td align="right">0.89</td>
<td align="left">0.00247</td>
</tr>
<tr class="even">
<td>278</td>
<td align="left"><a href="GO:0007144" class="uri">GO:0007144</a></td>
<td align="left">female meiosis I</td>
<td align="right">4</td>
<td align="right">4</td>
<td align="right">0.89</td>
<td align="left">0.00247</td>
</tr>
<tr class="odd">
<td>279</td>
<td align="left"><a href="GO:0009157" class="uri">GO:0009157</a></td>
<td align="left">deoxyribonucleoside monophosphate biosyn...</td>
<td align="right">4</td>
<td align="right">4</td>
<td align="right">0.89</td>
<td align="left">0.00247</td>
</tr>
<tr class="even">
<td>280</td>
<td align="left"><a href="GO:0009162" class="uri">GO:0009162</a></td>
<td align="left">deoxyribonucleoside monophosphate metabo...</td>
<td align="right">4</td>
<td align="right">4</td>
<td align="right">0.89</td>
<td align="left">0.00247</td>
</tr>
<tr class="odd">
<td>281</td>
<td align="left"><a href="GO:0010962" class="uri">GO:0010962</a></td>
<td align="left">regulation of glucan biosynthetic proces...</td>
<td align="right">4</td>
<td align="right">4</td>
<td align="right">0.89</td>
<td align="left">0.00247</td>
</tr>
<tr class="even">
<td>282</td>
<td align="left"><a href="GO:0035215" class="uri">GO:0035215</a></td>
<td align="left">genital disc development</td>
<td align="right">4</td>
<td align="right">4</td>
<td align="right">0.89</td>
<td align="left">0.00247</td>
</tr>
<tr class="odd">
<td>283</td>
<td align="left"><a href="GO:0035845" class="uri">GO:0035845</a></td>
<td align="left">photoreceptor cell outer segment organiz...</td>
<td align="right">4</td>
<td align="right">4</td>
<td align="right">0.89</td>
<td align="left">0.00247</td>
</tr>
<tr class="even">
<td>284</td>
<td align="left"><a href="GO:0045463" class="uri">GO:0045463</a></td>
<td align="left">R8 cell development</td>
<td align="right">4</td>
<td align="right">4</td>
<td align="right">0.89</td>
<td align="left">0.00247</td>
</tr>
<tr class="odd">
<td>285</td>
<td align="left"><a href="GO:0045465" class="uri">GO:0045465</a></td>
<td align="left">R8 cell differentiation</td>
<td align="right">4</td>
<td align="right">4</td>
<td align="right">0.89</td>
<td align="left">0.00247</td>
</tr>
<tr class="even">
<td>286</td>
<td align="left"><a href="GO:0051988" class="uri">GO:0051988</a></td>
<td align="left">regulation of attachment of spindle micr...</td>
<td align="right">4</td>
<td align="right">4</td>
<td align="right">0.89</td>
<td align="left">0.00247</td>
</tr>
<tr class="odd">
<td>287</td>
<td align="left"><a href="GO:0070873" class="uri">GO:0070873</a></td>
<td align="left">regulation of glycogen metabolic process</td>
<td align="right">4</td>
<td align="right">4</td>
<td align="right">0.89</td>
<td align="left">0.00247</td>
</tr>
<tr class="even">
<td>288</td>
<td align="left"><a href="GO:1901070" class="uri">GO:1901070</a></td>
<td align="left">guanosine-containing compound biosynthet...</td>
<td align="right">4</td>
<td align="right">4</td>
<td align="right">0.89</td>
<td align="left">0.00247</td>
</tr>
<tr class="odd">
<td>289</td>
<td align="left"><a href="GO:2001013" class="uri">GO:2001013</a></td>
<td align="left">epithelial cell proliferation involved i...</td>
<td align="right">4</td>
<td align="right">4</td>
<td align="right">0.89</td>
<td align="left">0.00247</td>
</tr>
<tr class="even">
<td>290</td>
<td align="left"><a href="GO:2001293" class="uri">GO:2001293</a></td>
<td align="left">malonyl-CoA metabolic process</td>
<td align="right">4</td>
<td align="right">4</td>
<td align="right">0.89</td>
<td align="left">0.00247</td>
</tr>
<tr class="odd">
<td>291</td>
<td align="left"><a href="GO:2001295" class="uri">GO:2001295</a></td>
<td align="left">malonyl-CoA biosynthetic process</td>
<td align="right">4</td>
<td align="right">4</td>
<td align="right">0.89</td>
<td align="left">0.00247</td>
</tr>
<tr class="even">
<td>292</td>
<td align="left"><a href="GO:0001947" class="uri">GO:0001947</a></td>
<td align="left">heart looping</td>
<td align="right">13</td>
<td align="right">8</td>
<td align="right">2.90</td>
<td align="left">0.00262</td>
</tr>
<tr class="odd">
<td>293</td>
<td align="left"><a href="GO:0006221" class="uri">GO:0006221</a></td>
<td align="left">pyrimidine nucleotide biosynthetic proce...</td>
<td align="right">13</td>
<td align="right">8</td>
<td align="right">2.90</td>
<td align="left">0.00262</td>
</tr>
<tr class="even">
<td>294</td>
<td align="left"><a href="GO:0055078" class="uri">GO:0055078</a></td>
<td align="left">sodium ion homeostasis</td>
<td align="right">13</td>
<td align="right">8</td>
<td align="right">2.90</td>
<td align="left">0.00262</td>
</tr>
<tr class="odd">
<td>295</td>
<td align="left"><a href="GO:0030199" class="uri">GO:0030199</a></td>
<td align="left">collagen fibril organization</td>
<td align="right">24</td>
<td align="right">12</td>
<td align="right">5.36</td>
<td align="left">0.00263</td>
</tr>
<tr class="even">
<td>296</td>
<td align="left"><a href="GO:0030858" class="uri">GO:0030858</a></td>
<td align="left">positive regulation of epithelial cell d...</td>
<td align="right">30</td>
<td align="right">14</td>
<td align="right">6.70</td>
<td align="left">0.00269</td>
</tr>
<tr class="odd">
<td>297</td>
<td align="left"><a href="GO:0033273" class="uri">GO:0033273</a></td>
<td align="left">response to vitamin</td>
<td align="right">30</td>
<td align="right">14</td>
<td align="right">6.70</td>
<td align="left">0.00269</td>
</tr>
<tr class="even">
<td>298</td>
<td align="left"><a href="GO:0008589" class="uri">GO:0008589</a></td>
<td align="left">regulation of smoothened signaling pathw...</td>
<td align="right">27</td>
<td align="right">13</td>
<td align="right">6.03</td>
<td align="left">0.00269</td>
</tr>
<tr class="odd">
<td>299</td>
<td align="left"><a href="GO:0018298" class="uri">GO:0018298</a></td>
<td align="left">protein-chromophore linkage</td>
<td align="right">27</td>
<td align="right">13</td>
<td align="right">6.03</td>
<td align="left">0.00269</td>
</tr>
<tr class="even">
<td>300</td>
<td align="left"><a href="GO:0006241" class="uri">GO:0006241</a></td>
<td align="left">CTP biosynthetic process</td>
<td align="right">6</td>
<td align="right">5</td>
<td align="right">1.34</td>
<td align="left">0.00269</td>
</tr>
<tr class="odd">
<td>301</td>
<td align="left"><a href="GO:0007442" class="uri">GO:0007442</a></td>
<td align="left">hindgut morphogenesis</td>
<td align="right">6</td>
<td align="right">5</td>
<td align="right">1.34</td>
<td align="left">0.00269</td>
</tr>
<tr class="even">
<td>302</td>
<td align="left"><a href="GO:0033383" class="uri">GO:0033383</a></td>
<td align="left">geranyl diphosphate metabolic process</td>
<td align="right">6</td>
<td align="right">5</td>
<td align="right">1.34</td>
<td align="left">0.00269</td>
</tr>
<tr class="odd">
<td>303</td>
<td align="left"><a href="GO:0033384" class="uri">GO:0033384</a></td>
<td align="left">geranyl diphosphate biosynthetic process</td>
<td align="right">6</td>
<td align="right">5</td>
<td align="right">1.34</td>
<td align="left">0.00269</td>
</tr>
<tr class="even">
<td>304</td>
<td align="left"><a href="GO:0042761" class="uri">GO:0042761</a></td>
<td align="left">very long-chain fatty acid biosynthetic ...</td>
<td align="right">6</td>
<td align="right">5</td>
<td align="right">1.34</td>
<td align="left">0.00269</td>
</tr>
<tr class="odd">
<td>305</td>
<td align="left"><a href="GO:0045337" class="uri">GO:0045337</a></td>
<td align="left">farnesyl diphosphate biosynthetic proces...</td>
<td align="right">6</td>
<td align="right">5</td>
<td align="right">1.34</td>
<td align="left">0.00269</td>
</tr>
<tr class="even">
<td>306</td>
<td align="left"><a href="GO:0045338" class="uri">GO:0045338</a></td>
<td align="left">farnesyl diphosphate metabolic process</td>
<td align="right">6</td>
<td align="right">5</td>
<td align="right">1.34</td>
<td align="left">0.00269</td>
</tr>
<tr class="odd">
<td>307</td>
<td align="left"><a href="GO:0046036" class="uri">GO:0046036</a></td>
<td align="left">CTP metabolic process</td>
<td align="right">6</td>
<td align="right">5</td>
<td align="right">1.34</td>
<td align="left">0.00269</td>
</tr>
<tr class="even">
<td>308</td>
<td align="left"><a href="GO:0035051" class="uri">GO:0035051</a></td>
<td align="left">cardiocyte differentiation</td>
<td align="right">49</td>
<td align="right">20</td>
<td align="right">10.94</td>
<td align="left">0.00273</td>
</tr>
<tr class="odd">
<td>309</td>
<td align="left"><a href="GO:0008203" class="uri">GO:0008203</a></td>
<td align="left">cholesterol metabolic process</td>
<td align="right">59</td>
<td align="right">23</td>
<td align="right">13.17</td>
<td align="left">0.00278</td>
</tr>
<tr class="even">
<td>310</td>
<td align="left"><a href="GO:0016125" class="uri">GO:0016125</a></td>
<td align="left">sterol metabolic process</td>
<td align="right">66</td>
<td align="right">25</td>
<td align="right">14.73</td>
<td align="left">0.00294</td>
</tr>
<tr class="odd">
<td>311</td>
<td align="left"><a href="GO:0006694" class="uri">GO:0006694</a></td>
<td align="left">steroid biosynthetic process</td>
<td align="right">46</td>
<td align="right">19</td>
<td align="right">10.27</td>
<td align="left">0.00295</td>
</tr>
<tr class="even">
<td>312</td>
<td align="left"><a href="GO:0071840" class="uri">GO:0071840</a></td>
<td align="left">cellular component organization or bioge...</td>
<td align="right">2238</td>
<td align="right">544</td>
<td align="right">499.58</td>
<td align="left">0.00296</td>
</tr>
<tr class="odd">
<td>313</td>
<td align="left"><a href="GO:0050807" class="uri">GO:0050807</a></td>
<td align="left">regulation of synapse organization</td>
<td align="right">87</td>
<td align="right">31</td>
<td align="right">19.42</td>
<td align="left">0.00305</td>
</tr>
<tr class="even">
<td>314</td>
<td align="left"><a href="GO:0050890" class="uri">GO:0050890</a></td>
<td align="left">cognition</td>
<td align="right">80</td>
<td align="right">29</td>
<td align="right">17.86</td>
<td align="left">0.00306</td>
</tr>
<tr class="odd">
<td>315</td>
<td align="left"><a href="GO:0090596" class="uri">GO:0090596</a></td>
<td align="left">sensory organ morphogenesis</td>
<td align="right">127</td>
<td align="right">42</td>
<td align="right">28.35</td>
<td align="left">0.00324</td>
</tr>
<tr class="even">
<td>316</td>
<td align="left"><a href="GO:0045165" class="uri">GO:0045165</a></td>
<td align="left">cell fate commitment</td>
<td align="right">77</td>
<td align="right">28</td>
<td align="right">17.19</td>
<td align="left">0.00340</td>
</tr>
<tr class="odd">
<td>317</td>
<td align="left"><a href="GO:0031589" class="uri">GO:0031589</a></td>
<td align="left">cell-substrate adhesion</td>
<td align="right">120</td>
<td align="right">40</td>
<td align="right">26.79</td>
<td align="left">0.00342</td>
</tr>
<tr class="even">
<td>318</td>
<td align="left"><a href="GO:0010564" class="uri">GO:0010564</a></td>
<td align="left">regulation of cell cycle process</td>
<td align="right">165</td>
<td align="right">52</td>
<td align="right">36.83</td>
<td align="left">0.00361</td>
</tr>
<tr class="odd">
<td>319</td>
<td align="left"><a href="GO:0044272" class="uri">GO:0044272</a></td>
<td align="left">sulfur compound biosynthetic process</td>
<td align="right">50</td>
<td align="right">20</td>
<td align="right">11.16</td>
<td align="left">0.00362</td>
</tr>
<tr class="even">
<td>320</td>
<td align="left"><a href="GO:0072202" class="uri">GO:0072202</a></td>
<td align="left">cell differentiation involved in metanep...</td>
<td align="right">19</td>
<td align="right">10</td>
<td align="right">4.24</td>
<td align="left">0.00371</td>
</tr>
<tr class="odd">
<td>321</td>
<td align="left"><a href="GO:0048568" class="uri">GO:0048568</a></td>
<td align="left">embryonic organ development</td>
<td align="right">154</td>
<td align="right">49</td>
<td align="right">34.38</td>
<td align="left">0.00377</td>
</tr>
<tr class="even">
<td>322</td>
<td align="left"><a href="GO:0044259" class="uri">GO:0044259</a></td>
<td align="left">multicellular organismal macromolecule m...</td>
<td align="right">34</td>
<td align="right">15</td>
<td align="right">7.59</td>
<td align="left">0.00377</td>
</tr>
<tr class="odd">
<td>323</td>
<td align="left"><a href="GO:0048545" class="uri">GO:0048545</a></td>
<td align="left">response to steroid hormone</td>
<td align="right">81</td>
<td align="right">29</td>
<td align="right">18.08</td>
<td align="left">0.00378</td>
</tr>
<tr class="even">
<td>324</td>
<td align="left"><a href="GO:0002026" class="uri">GO:0002026</a></td>
<td align="left">regulation of the force of heart contrac...</td>
<td align="right">11</td>
<td align="right">7</td>
<td align="right">2.46</td>
<td align="left">0.00381</td>
</tr>
<tr class="odd">
<td>325</td>
<td align="left"><a href="GO:0007447" class="uri">GO:0007447</a></td>
<td align="left">imaginal disc pattern formation</td>
<td align="right">11</td>
<td align="right">7</td>
<td align="right">2.46</td>
<td align="left">0.00381</td>
</tr>
<tr class="even">
<td>326</td>
<td align="left"><a href="GO:0035136" class="uri">GO:0035136</a></td>
<td align="left">forelimb morphogenesis</td>
<td align="right">11</td>
<td align="right">7</td>
<td align="right">2.46</td>
<td align="left">0.00381</td>
</tr>
<tr class="odd">
<td>327</td>
<td align="left"><a href="GO:0046134" class="uri">GO:0046134</a></td>
<td align="left">pyrimidine nucleoside biosynthetic proce...</td>
<td align="right">11</td>
<td align="right">7</td>
<td align="right">2.46</td>
<td align="left">0.00381</td>
</tr>
<tr class="even">
<td>328</td>
<td align="left"><a href="GO:0003407" class="uri">GO:0003407</a></td>
<td align="left">neural retina development</td>
<td align="right">28</td>
<td align="right">13</td>
<td align="right">6.25</td>
<td align="left">0.00402</td>
</tr>
<tr class="odd">
<td>329</td>
<td align="left"><a href="GO:0019748" class="uri">GO:0019748</a></td>
<td align="left">secondary metabolic process</td>
<td align="right">28</td>
<td align="right">13</td>
<td align="right">6.25</td>
<td align="left">0.00402</td>
</tr>
<tr class="even">
<td>330</td>
<td align="left"><a href="GO:0007611" class="uri">GO:0007611</a></td>
<td align="left">learning or memory</td>
<td align="right">78</td>
<td align="right">28</td>
<td align="right">17.41</td>
<td align="left">0.00420</td>
</tr>
<tr class="odd">
<td>331</td>
<td align="left"><a href="GO:0007420" class="uri">GO:0007420</a></td>
<td align="left">brain development</td>
<td align="right">295</td>
<td align="right">85</td>
<td align="right">65.85</td>
<td align="left">0.00459</td>
</tr>
<tr class="even">
<td>332</td>
<td align="left"><a href="GO:0051726" class="uri">GO:0051726</a></td>
<td align="left">regulation of cell cycle</td>
<td align="right">307</td>
<td align="right">88</td>
<td align="right">68.53</td>
<td align="left">0.00463</td>
</tr>
<tr class="odd">
<td>333</td>
<td align="left"><a href="GO:0003013" class="uri">GO:0003013</a></td>
<td align="left">circulatory system process</td>
<td align="right">100</td>
<td align="right">34</td>
<td align="right">22.32</td>
<td align="left">0.00470</td>
</tr>
<tr class="even">
<td>334</td>
<td align="left"><a href="GO:0044089" class="uri">GO:0044089</a></td>
<td align="left">positive regulation of cellular componen...</td>
<td align="right">122</td>
<td align="right">40</td>
<td align="right">27.23</td>
<td align="left">0.00474</td>
</tr>
<tr class="odd">
<td>335</td>
<td align="left"><a href="GO:0050803" class="uri">GO:0050803</a></td>
<td align="left">regulation of synapse structure or activ...</td>
<td align="right">122</td>
<td align="right">40</td>
<td align="right">27.23</td>
<td align="left">0.00474</td>
</tr>
<tr class="even">
<td>336</td>
<td align="left"><a href="GO:0007049" class="uri">GO:0007049</a></td>
<td align="left">cell cycle</td>
<td align="right">692</td>
<td align="right">182</td>
<td align="right">154.47</td>
<td align="left">0.00501</td>
</tr>
<tr class="odd">
<td>337</td>
<td align="left"><a href="GO:0048754" class="uri">GO:0048754</a></td>
<td align="left">branching morphogenesis of an epithelial...</td>
<td align="right">58</td>
<td align="right">22</td>
<td align="right">12.95</td>
<td align="left">0.00502</td>
</tr>
<tr class="even">
<td>338</td>
<td align="left"><a href="GO:0006996" class="uri">GO:0006996</a></td>
<td align="left">organelle organization</td>
<td align="right">1281</td>
<td align="right">321</td>
<td align="right">285.95</td>
<td align="left">0.00510</td>
</tr>
<tr class="odd">
<td>339</td>
<td align="left"><a href="GO:0007601" class="uri">GO:0007601</a></td>
<td align="left">visual perception</td>
<td align="right">72</td>
<td align="right">26</td>
<td align="right">16.07</td>
<td align="left">0.00520</td>
</tr>
<tr class="even">
<td>340</td>
<td align="left"><a href="GO:0031960" class="uri">GO:0031960</a></td>
<td align="left">response to corticosteroid</td>
<td align="right">35</td>
<td align="right">15</td>
<td align="right">7.81</td>
<td align="left">0.00525</td>
</tr>
<tr class="odd">
<td>341</td>
<td align="left"><a href="GO:0007141" class="uri">GO:0007141</a></td>
<td align="left">male meiosis I</td>
<td align="right">9</td>
<td align="right">6</td>
<td align="right">2.01</td>
<td align="left">0.00548</td>
</tr>
<tr class="even">
<td>342</td>
<td align="left"><a href="GO:0007391" class="uri">GO:0007391</a></td>
<td align="left">dorsal closure</td>
<td align="right">9</td>
<td align="right">6</td>
<td align="right">2.01</td>
<td align="left">0.00548</td>
</tr>
<tr class="odd">
<td>343</td>
<td align="left"><a href="GO:0008608" class="uri">GO:0008608</a></td>
<td align="left">attachment of spindle microtubules to ki...</td>
<td align="right">9</td>
<td align="right">6</td>
<td align="right">2.01</td>
<td align="left">0.00548</td>
</tr>
<tr class="even">
<td>344</td>
<td align="left"><a href="GO:0009148" class="uri">GO:0009148</a></td>
<td align="left">pyrimidine nucleoside triphosphate biosy...</td>
<td align="right">9</td>
<td align="right">6</td>
<td align="right">2.01</td>
<td align="left">0.00548</td>
</tr>
<tr class="odd">
<td>345</td>
<td align="left"><a href="GO:0021522" class="uri">GO:0021522</a></td>
<td align="left">spinal cord motor neuron differentiation</td>
<td align="right">9</td>
<td align="right">6</td>
<td align="right">2.01</td>
<td align="left">0.00548</td>
</tr>
<tr class="even">
<td>346</td>
<td align="left"><a href="GO:0036099" class="uri">GO:0036099</a></td>
<td align="left">female germ-line stem cell population ma...</td>
<td align="right">9</td>
<td align="right">6</td>
<td align="right">2.01</td>
<td align="left">0.00548</td>
</tr>
<tr class="odd">
<td>347</td>
<td align="left"><a href="GO:0090128" class="uri">GO:0090128</a></td>
<td align="left">regulation of synapse maturation</td>
<td align="right">9</td>
<td align="right">6</td>
<td align="right">2.01</td>
<td align="left">0.00548</td>
</tr>
<tr class="even">
<td>348</td>
<td align="left"><a href="GO:0090129" class="uri">GO:0090129</a></td>
<td align="left">positive regulation of synapse maturatio...</td>
<td align="right">9</td>
<td align="right">6</td>
<td align="right">2.01</td>
<td align="left">0.00548</td>
</tr>
<tr class="odd">
<td>349</td>
<td align="left"><a href="GO:0032963" class="uri">GO:0032963</a></td>
<td align="left">collagen metabolic process</td>
<td align="right">32</td>
<td align="right">14</td>
<td align="right">7.14</td>
<td align="left">0.00556</td>
</tr>
<tr class="even">
<td>350</td>
<td align="left"><a href="GO:0042461" class="uri">GO:0042461</a></td>
<td align="left">photoreceptor cell development</td>
<td align="right">32</td>
<td align="right">14</td>
<td align="right">7.14</td>
<td align="left">0.00556</td>
</tr>
<tr class="odd">
<td>351</td>
<td align="left"><a href="GO:0051384" class="uri">GO:0051384</a></td>
<td align="left">response to glucocorticoid</td>
<td align="right">32</td>
<td align="right">14</td>
<td align="right">7.14</td>
<td align="left">0.00556</td>
</tr>
<tr class="even">
<td>352</td>
<td align="left"><a href="GO:0061138" class="uri">GO:0061138</a></td>
<td align="left">morphogenesis of a branching epithelium</td>
<td align="right">62</td>
<td align="right">23</td>
<td align="right">13.84</td>
<td align="left">0.00573</td>
</tr>
<tr class="odd">
<td>353</td>
<td align="left"><a href="GO:0007160" class="uri">GO:0007160</a></td>
<td align="left">cell-matrix adhesion</td>
<td align="right">69</td>
<td align="right">25</td>
<td align="right">15.40</td>
<td align="left">0.00579</td>
</tr>
<tr class="even">
<td>354</td>
<td align="left"><a href="GO:2000026" class="uri">GO:2000026</a></td>
<td align="left">regulation of multicellular organismal d...</td>
<td align="right">578</td>
<td align="right">154</td>
<td align="right">129.03</td>
<td align="left">0.00582</td>
</tr>
<tr class="odd">
<td>355</td>
<td align="left"><a href="GO:0003143" class="uri">GO:0003143</a></td>
<td align="left">embryonic heart tube morphogenesis</td>
<td align="right">20</td>
<td align="right">10</td>
<td align="right">4.46</td>
<td align="left">0.00596</td>
</tr>
<tr class="even">
<td>356</td>
<td align="left"><a href="GO:0006730" class="uri">GO:0006730</a></td>
<td align="left">one-carbon metabolic process</td>
<td align="right">20</td>
<td align="right">10</td>
<td align="right">4.46</td>
<td align="left">0.00596</td>
</tr>
<tr class="odd">
<td>357</td>
<td align="left"><a href="GO:0072170" class="uri">GO:0072170</a></td>
<td align="left">metanephric tubule development</td>
<td align="right">20</td>
<td align="right">10</td>
<td align="right">4.46</td>
<td align="left">0.00596</td>
</tr>
<tr class="even">
<td>358</td>
<td align="left"><a href="GO:0072215" class="uri">GO:0072215</a></td>
<td align="left">regulation of metanephros development</td>
<td align="right">20</td>
<td align="right">10</td>
<td align="right">4.46</td>
<td align="left">0.00596</td>
</tr>
<tr class="odd">
<td>359</td>
<td align="left"><a href="GO:0072234" class="uri">GO:0072234</a></td>
<td align="left">metanephric nephron tubule development</td>
<td align="right">20</td>
<td align="right">10</td>
<td align="right">4.46</td>
<td align="left">0.00596</td>
</tr>
<tr class="even">
<td>360</td>
<td align="left"><a href="GO:0030574" class="uri">GO:0030574</a></td>
<td align="left">collagen catabolic process</td>
<td align="right">23</td>
<td align="right">11</td>
<td align="right">5.13</td>
<td align="left">0.00606</td>
</tr>
<tr class="odd">
<td>361</td>
<td align="left"><a href="GO:0042398" class="uri">GO:0042398</a></td>
<td align="left">cellular modified amino acid biosyntheti...</td>
<td align="right">23</td>
<td align="right">11</td>
<td align="right">5.13</td>
<td align="left">0.00606</td>
</tr>
<tr class="even">
<td>362</td>
<td align="left"><a href="GO:0008015" class="uri">GO:0008015</a></td>
<td align="left">blood circulation</td>
<td align="right">98</td>
<td align="right">33</td>
<td align="right">21.88</td>
<td align="left">0.00626</td>
</tr>
<tr class="odd">
<td>363</td>
<td align="left"><a href="GO:0032886" class="uri">GO:0032886</a></td>
<td align="left">regulation of microtubule-based process</td>
<td align="right">66</td>
<td align="right">24</td>
<td align="right">14.73</td>
<td align="left">0.00644</td>
</tr>
<tr class="even">
<td>364</td>
<td align="left"><a href="GO:0035335" class="uri">GO:0035335</a></td>
<td align="left">peptidyl-tyrosine dephosphorylation</td>
<td align="right">70</td>
<td align="right">25</td>
<td align="right">15.63</td>
<td align="left">0.00714</td>
</tr>
<tr class="odd">
<td>365</td>
<td align="left"><a href="GO:0051271" class="uri">GO:0051271</a></td>
<td align="left">negative regulation of cellular componen...</td>
<td align="right">36</td>
<td align="right">15</td>
<td align="right">8.04</td>
<td align="left">0.00718</td>
</tr>
<tr class="even">
<td>366</td>
<td align="left"><a href="GO:0001700" class="uri">GO:0001700</a></td>
<td align="left">embryonic development via the syncytial ...</td>
<td align="right">12</td>
<td align="right">7</td>
<td align="right">2.68</td>
<td align="left">0.00740</td>
</tr>
<tr class="odd">
<td>367</td>
<td align="left"><a href="GO:0007143" class="uri">GO:0007143</a></td>
<td align="left">female meiotic division</td>
<td align="right">12</td>
<td align="right">7</td>
<td align="right">2.68</td>
<td align="left">0.00740</td>
</tr>
<tr class="even">
<td>368</td>
<td align="left"><a href="GO:0009218" class="uri">GO:0009218</a></td>
<td align="left">pyrimidine ribonucleotide metabolic proc...</td>
<td align="right">12</td>
<td align="right">7</td>
<td align="right">2.68</td>
<td align="left">0.00740</td>
</tr>
<tr class="odd">
<td>369</td>
<td align="left"><a href="GO:0021515" class="uri">GO:0021515</a></td>
<td align="left">cell differentiation in spinal cord</td>
<td align="right">12</td>
<td align="right">7</td>
<td align="right">2.68</td>
<td align="left">0.00740</td>
</tr>
<tr class="even">
<td>370</td>
<td align="left"><a href="GO:0030718" class="uri">GO:0030718</a></td>
<td align="left">germ-line stem cell population maintenan...</td>
<td align="right">12</td>
<td align="right">7</td>
<td align="right">2.68</td>
<td align="left">0.00740</td>
</tr>
<tr class="odd">
<td>371</td>
<td align="left"><a href="GO:0034453" class="uri">GO:0034453</a></td>
<td align="left">microtubule anchoring</td>
<td align="right">12</td>
<td align="right">7</td>
<td align="right">2.68</td>
<td align="left">0.00740</td>
</tr>
<tr class="even">
<td>372</td>
<td align="left"><a href="GO:0060135" class="uri">GO:0060135</a></td>
<td align="left">maternal process involved in female preg...</td>
<td align="right">12</td>
<td align="right">7</td>
<td align="right">2.68</td>
<td align="left">0.00740</td>
</tr>
<tr class="odd">
<td>373</td>
<td align="left"><a href="GO:0007288" class="uri">GO:0007288</a></td>
<td align="left">sperm axoneme assembly</td>
<td align="right">7</td>
<td align="right">5</td>
<td align="right">1.56</td>
<td align="left">0.00769</td>
</tr>
<tr class="even">
<td>374</td>
<td align="left"><a href="GO:0009209" class="uri">GO:0009209</a></td>
<td align="left">pyrimidine ribonucleoside triphosphate b...</td>
<td align="right">7</td>
<td align="right">5</td>
<td align="right">1.56</td>
<td align="left">0.00769</td>
</tr>
<tr class="odd">
<td>375</td>
<td align="left"><a href="GO:0009263" class="uri">GO:0009263</a></td>
<td align="left">deoxyribonucleotide biosynthetic process</td>
<td align="right">7</td>
<td align="right">5</td>
<td align="right">1.56</td>
<td align="left">0.00769</td>
</tr>
<tr class="even">
<td>376</td>
<td align="left"><a href="GO:0015949" class="uri">GO:0015949</a></td>
<td align="left">nucleobase-containing small molecule int...</td>
<td align="right">7</td>
<td align="right">5</td>
<td align="right">1.56</td>
<td align="left">0.00769</td>
</tr>
<tr class="odd">
<td>377</td>
<td align="left"><a href="GO:0021513" class="uri">GO:0021513</a></td>
<td align="left">spinal cord dorsal/ventral patterning</td>
<td align="right">7</td>
<td align="right">5</td>
<td align="right">1.56</td>
<td align="left">0.00769</td>
</tr>
<tr class="even">
<td>378</td>
<td align="left"><a href="GO:0045792" class="uri">GO:0045792</a></td>
<td align="left">negative regulation of cell size</td>
<td align="right">7</td>
<td align="right">5</td>
<td align="right">1.56</td>
<td align="left">0.00769</td>
</tr>
<tr class="odd">
<td>379</td>
<td align="left"><a href="GO:0046271" class="uri">GO:0046271</a></td>
<td align="left">phenylpropanoid catabolic process</td>
<td align="right">7</td>
<td align="right">5</td>
<td align="right">1.56</td>
<td align="left">0.00769</td>
</tr>
<tr class="even">
<td>380</td>
<td align="left"><a href="GO:0046274" class="uri">GO:0046274</a></td>
<td align="left">lignin catabolic process</td>
<td align="right">7</td>
<td align="right">5</td>
<td align="right">1.56</td>
<td align="left">0.00769</td>
</tr>
<tr class="odd">
<td>381</td>
<td align="left"><a href="GO:0046331" class="uri">GO:0046331</a></td>
<td align="left">lateral inhibition</td>
<td align="right">7</td>
<td align="right">5</td>
<td align="right">1.56</td>
<td align="left">0.00769</td>
</tr>
<tr class="even">
<td>382</td>
<td align="left"><a href="GO:0097369" class="uri">GO:0097369</a></td>
<td align="left">sodium ion import</td>
<td align="right">7</td>
<td align="right">5</td>
<td align="right">1.56</td>
<td align="left">0.00769</td>
</tr>
<tr class="odd">
<td>383</td>
<td align="left"><a href="GO:0032465" class="uri">GO:0032465</a></td>
<td align="left">regulation of cytokinesis</td>
<td align="right">33</td>
<td align="right">14</td>
<td align="right">7.37</td>
<td align="left">0.00770</td>
</tr>
<tr class="even">
<td>384</td>
<td align="left"><a href="GO:0042330" class="uri">GO:0042330</a></td>
<td align="left">taxis</td>
<td align="right">233</td>
<td align="right">68</td>
<td align="right">52.01</td>
<td align="left">0.00774</td>
</tr>
<tr class="odd">
<td>385</td>
<td align="left"><a href="GO:0044839" class="uri">GO:0044839</a></td>
<td align="left">cell cycle G2/M phase transition</td>
<td align="right">53</td>
<td align="right">20</td>
<td align="right">11.83</td>
<td align="left">0.00779</td>
</tr>
<tr class="even">
<td>386</td>
<td align="left"><a href="GO:0007417" class="uri">GO:0007417</a></td>
<td align="left">central nervous system development</td>
<td align="right">349</td>
<td align="right">97</td>
<td align="right">77.91</td>
<td align="left">0.00799</td>
</tr>
<tr class="odd">
<td>387</td>
<td align="left"><a href="GO:0055007" class="uri">GO:0055007</a></td>
<td align="left">cardiac muscle cell differentiation</td>
<td align="right">43</td>
<td align="right">17</td>
<td align="right">9.60</td>
<td align="left">0.00809</td>
</tr>
<tr class="even">
<td>388</td>
<td align="left"><a href="GO:0060322" class="uri">GO:0060322</a></td>
<td align="left">head development</td>
<td align="right">309</td>
<td align="right">87</td>
<td align="right">68.98</td>
<td align="left">0.00815</td>
</tr>
<tr class="odd">
<td>389</td>
<td align="left"><a href="GO:0050806" class="uri">GO:0050806</a></td>
<td align="left">positive regulation of synaptic transmis...</td>
<td align="right">30</td>
<td align="right">13</td>
<td align="right">6.70</td>
<td align="left">0.00820</td>
</tr>
<tr class="even">
<td>390</td>
<td align="left"><a href="GO:0018212" class="uri">GO:0018212</a></td>
<td align="left">peptidyl-tyrosine modification</td>
<td align="right">179</td>
<td align="right">54</td>
<td align="right">39.96</td>
<td align="left">0.00826</td>
</tr>
<tr class="odd">
<td>391</td>
<td align="left"><a href="GO:0007616" class="uri">GO:0007616</a></td>
<td align="left">long-term memory</td>
<td align="right">15</td>
<td align="right">8</td>
<td align="right">3.35</td>
<td align="left">0.00851</td>
</tr>
<tr class="even">
<td>392</td>
<td align="left"><a href="GO:0033865" class="uri">GO:0033865</a></td>
<td align="left">nucleoside bisphosphate metabolic proces...</td>
<td align="right">15</td>
<td align="right">8</td>
<td align="right">3.35</td>
<td align="left">0.00851</td>
</tr>
<tr class="odd">
<td>393</td>
<td align="left"><a href="GO:0033875" class="uri">GO:0033875</a></td>
<td align="left">ribonucleoside bisphosphate metabolic pr...</td>
<td align="right">15</td>
<td align="right">8</td>
<td align="right">3.35</td>
<td align="left">0.00851</td>
</tr>
<tr class="even">
<td>394</td>
<td align="left"><a href="GO:0034032" class="uri">GO:0034032</a></td>
<td align="left">purine nucleoside bisphosphate metabolic...</td>
<td align="right">15</td>
<td align="right">8</td>
<td align="right">3.35</td>
<td align="left">0.00851</td>
</tr>
<tr class="odd">
<td>395</td>
<td align="left"><a href="GO:0072528" class="uri">GO:0072528</a></td>
<td align="left">pyrimidine-containing compound biosynthe...</td>
<td align="right">15</td>
<td align="right">8</td>
<td align="right">3.35</td>
<td align="left">0.00851</td>
</tr>
<tr class="even">
<td>396</td>
<td align="left"><a href="GO:0000086" class="uri">GO:0000086</a></td>
<td align="left">G2/M transition of mitotic cell cycle</td>
<td align="right">50</td>
<td align="right">19</td>
<td align="right">11.16</td>
<td align="left">0.00861</td>
</tr>
<tr class="odd">
<td>397</td>
<td align="left"><a href="GO:0051346" class="uri">GO:0051346</a></td>
<td align="left">negative regulation of hydrolase activit...</td>
<td align="right">111</td>
<td align="right">36</td>
<td align="right">24.78</td>
<td align="left">0.00863</td>
</tr>
<tr class="even">
<td>398</td>
<td align="left"><a href="GO:0051298" class="uri">GO:0051298</a></td>
<td align="left">centrosome duplication</td>
<td align="right">27</td>
<td align="right">12</td>
<td align="right">6.03</td>
<td align="left">0.00864</td>
</tr>
<tr class="odd">
<td>399</td>
<td align="left"><a href="GO:0050680" class="uri">GO:0050680</a></td>
<td align="left">negative regulation of epithelial cell p...</td>
<td align="right">40</td>
<td align="right">16</td>
<td align="right">8.93</td>
<td align="left">0.00884</td>
</tr>
<tr class="even">
<td>400</td>
<td align="left"><a href="GO:0007424" class="uri">GO:0007424</a></td>
<td align="left">open tracheal system development</td>
<td align="right">24</td>
<td align="right">11</td>
<td align="right">5.36</td>
<td align="left">0.00897</td>
</tr>
<tr class="odd">
<td>401</td>
<td align="left"><a href="GO:0043931" class="uri">GO:0043931</a></td>
<td align="left">ossification involved in bone maturation</td>
<td align="right">24</td>
<td align="right">11</td>
<td align="right">5.36</td>
<td align="left">0.00897</td>
</tr>
<tr class="even">
<td>402</td>
<td align="left"><a href="GO:0006898" class="uri">GO:0006898</a></td>
<td align="left">receptor-mediated endocytosis</td>
<td align="right">130</td>
<td align="right">41</td>
<td align="right">29.02</td>
<td align="left">0.00900</td>
</tr>
<tr class="odd">
<td>403</td>
<td align="left"><a href="GO:0001709" class="uri">GO:0001709</a></td>
<td align="left">cell fate determination</td>
<td align="right">18</td>
<td align="right">9</td>
<td align="right">4.02</td>
<td align="left">0.00902</td>
</tr>
<tr class="even">
<td>404</td>
<td align="left"><a href="GO:0006636" class="uri">GO:0006636</a></td>
<td align="left">unsaturated fatty acid biosynthetic proc...</td>
<td align="right">18</td>
<td align="right">9</td>
<td align="right">4.02</td>
<td align="left">0.00902</td>
</tr>
<tr class="odd">
<td>405</td>
<td align="left"><a href="GO:0072257" class="uri">GO:0072257</a></td>
<td align="left">metanephric nephron tubule epithelial ce...</td>
<td align="right">18</td>
<td align="right">9</td>
<td align="right">4.02</td>
<td align="left">0.00902</td>
</tr>
<tr class="even">
<td>406</td>
<td align="left"><a href="GO:0072307" class="uri">GO:0072307</a></td>
<td align="left">regulation of metanephric nephron tubule...</td>
<td align="right">18</td>
<td align="right">9</td>
<td align="right">4.02</td>
<td align="left">0.00902</td>
</tr>
<tr class="odd">
<td>407</td>
<td align="left"><a href="GO:0000904" class="uri">GO:0000904</a></td>
<td align="left">cell morphogenesis involved in different...</td>
<td align="right">355</td>
<td align="right">98</td>
<td align="right">79.25</td>
<td align="left">0.00942</td>
</tr>
<tr class="even">
<td>408</td>
<td align="left"><a href="GO:0021510" class="uri">GO:0021510</a></td>
<td align="left">spinal cord development</td>
<td align="right">37</td>
<td align="right">15</td>
<td align="right">8.26</td>
<td align="left">0.00964</td>
</tr>
<tr class="odd">
<td>409</td>
<td align="left"><a href="GO:0044236" class="uri">GO:0044236</a></td>
<td align="left">multicellular organismal metabolic proce...</td>
<td align="right">37</td>
<td align="right">15</td>
<td align="right">8.26</td>
<td align="left">0.00964</td>
</tr>
<tr class="even">
<td>410</td>
<td align="left"><a href="GO:0051270" class="uri">GO:0051270</a></td>
<td align="left">regulation of cellular component movemen...</td>
<td align="right">192</td>
<td align="right">57</td>
<td align="right">42.86</td>
<td align="left">0.00969</td>
</tr>
<tr class="odd">
<td>411</td>
<td align="left"><a href="GO:0000278" class="uri">GO:0000278</a></td>
<td align="left">mitotic cell cycle</td>
<td align="right">441</td>
<td align="right">119</td>
<td align="right">98.44</td>
<td align="left">0.00971</td>
</tr>
<tr class="even">
<td>412</td>
<td align="left"><a href="GO:0003006" class="uri">GO:0003006</a></td>
<td align="left">developmental process involved in reprod...</td>
<td align="right">263</td>
<td align="right">75</td>
<td align="right">58.71</td>
<td align="left">0.00976</td>
</tr>
<tr class="odd">
<td>413</td>
<td align="left"><a href="GO:0006228" class="uri">GO:0006228</a></td>
<td align="left">UTP biosynthetic process</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.12</td>
<td align="left">0.01017</td>
</tr>
<tr class="even">
<td>414</td>
<td align="left"><a href="GO:0009453" class="uri">GO:0009453</a></td>
<td align="left">energy taxis</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.12</td>
<td align="left">0.01017</td>
</tr>
<tr class="odd">
<td>415</td>
<td align="left"><a href="GO:0010453" class="uri">GO:0010453</a></td>
<td align="left">regulation of cell fate commitment</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.12</td>
<td align="left">0.01017</td>
</tr>
<tr class="even">
<td>416</td>
<td align="left"><a href="GO:0010596" class="uri">GO:0010596</a></td>
<td align="left">negative regulation of endothelial cell ...</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.12</td>
<td align="left">0.01017</td>
</tr>
<tr class="odd">
<td>417</td>
<td align="left"><a href="GO:0015937" class="uri">GO:0015937</a></td>
<td align="left">coenzyme A biosynthetic process</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.12</td>
<td align="left">0.01017</td>
</tr>
<tr class="even">
<td>418</td>
<td align="left"><a href="GO:0032466" class="uri">GO:0032466</a></td>
<td align="left">negative regulation of cytokinesis</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.12</td>
<td align="left">0.01017</td>
</tr>
<tr class="odd">
<td>419</td>
<td align="left"><a href="GO:0032881" class="uri">GO:0032881</a></td>
<td align="left">regulation of polysaccharide metabolic p...</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.12</td>
<td align="left">0.01017</td>
</tr>
<tr class="even">
<td>420</td>
<td align="left"><a href="GO:0032885" class="uri">GO:0032885</a></td>
<td align="left">regulation of polysaccharide biosyntheti...</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.12</td>
<td align="left">0.01017</td>
</tr>
<tr class="odd">
<td>421</td>
<td align="left"><a href="GO:0033866" class="uri">GO:0033866</a></td>
<td align="left">nucleoside bisphosphate biosynthetic pro...</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.12</td>
<td align="left">0.01017</td>
</tr>
<tr class="even">
<td>422</td>
<td align="left"><a href="GO:0034030" class="uri">GO:0034030</a></td>
<td align="left">ribonucleoside bisphosphate biosynthetic...</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.12</td>
<td align="left">0.01017</td>
</tr>
<tr class="odd">
<td>423</td>
<td align="left"><a href="GO:0034033" class="uri">GO:0034033</a></td>
<td align="left">purine nucleoside bisphosphate biosynthe...</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.12</td>
<td align="left">0.01017</td>
</tr>
<tr class="even">
<td>424</td>
<td align="left"><a href="GO:0042331" class="uri">GO:0042331</a></td>
<td align="left">phototaxis</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.12</td>
<td align="left">0.01017</td>
</tr>
<tr class="odd">
<td>425</td>
<td align="left"><a href="GO:0042659" class="uri">GO:0042659</a></td>
<td align="left">regulation of cell fate specification</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.12</td>
<td align="left">0.01017</td>
</tr>
<tr class="even">
<td>426</td>
<td align="left"><a href="GO:0060287" class="uri">GO:0060287</a></td>
<td align="left">epithelial cilium movement involved in d...</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.12</td>
<td align="left">0.01017</td>
</tr>
<tr class="odd">
<td>427</td>
<td align="left"><a href="GO:1900273" class="uri">GO:1900273</a></td>
<td align="left">positive regulation of long-term synapti...</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.12</td>
<td align="left">0.01017</td>
</tr>
<tr class="even">
<td>428</td>
<td align="left"><a href="GO:2000463" class="uri">GO:2000463</a></td>
<td align="left">positive regulation of excitatory postsy...</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.12</td>
<td align="left">0.01017</td>
</tr>
<tr class="odd">
<td>429</td>
<td align="left"><a href="GO:0010951" class="uri">GO:0010951</a></td>
<td align="left">negative regulation of endopeptidase act...</td>
<td align="right">51</td>
<td align="right">19</td>
<td align="right">11.38</td>
<td align="left">0.01092</td>
</tr>
<tr class="even">
<td>430</td>
<td align="left"><a href="GO:0000727" class="uri">GO:0000727</a></td>
<td align="left">double-strand break repair via break-ind...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="odd">
<td>431</td>
<td align="left"><a href="GO:0002315" class="uri">GO:0002315</a></td>
<td align="left">marginal zone B cell differentiation</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="even">
<td>432</td>
<td align="left"><a href="GO:0003356" class="uri">GO:0003356</a></td>
<td align="left">regulation of cilium beat frequency</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="odd">
<td>433</td>
<td align="left"><a href="GO:0006348" class="uri">GO:0006348</a></td>
<td align="left">chromatin silencing at telomere</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="even">
<td>434</td>
<td align="left"><a href="GO:0007443" class="uri">GO:0007443</a></td>
<td align="left">Malpighian tubule morphogenesis</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="odd">
<td>435</td>
<td align="left"><a href="GO:0007449" class="uri">GO:0007449</a></td>
<td align="left">proximal/distal pattern formation, imagi...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="even">
<td>436</td>
<td align="left"><a href="GO:0007460" class="uri">GO:0007460</a></td>
<td align="left">R8 cell fate commitment</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="odd">
<td>437</td>
<td align="left"><a href="GO:0007473" class="uri">GO:0007473</a></td>
<td align="left">wing disc proximal/distal pattern format...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="even">
<td>438</td>
<td align="left"><a href="GO:0007483" class="uri">GO:0007483</a></td>
<td align="left">genital disc morphogenesis</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="odd">
<td>439</td>
<td align="left"><a href="GO:0007484" class="uri">GO:0007484</a></td>
<td align="left">imaginal disc-derived genitalia developm...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="even">
<td>440</td>
<td align="left"><a href="GO:0007485" class="uri">GO:0007485</a></td>
<td align="left">imaginal disc-derived male genitalia dev...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="odd">
<td>441</td>
<td align="left"><a href="GO:0007521" class="uri">GO:0007521</a></td>
<td align="left">muscle cell fate determination</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="even">
<td>442</td>
<td align="left"><a href="GO:0007615" class="uri">GO:0007615</a></td>
<td align="left">anesthesia-resistant memory</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="odd">
<td>443</td>
<td align="left"><a href="GO:0009996" class="uri">GO:0009996</a></td>
<td align="left">negative regulation of cell fate specifi...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="even">
<td>444</td>
<td align="left"><a href="GO:0010002" class="uri">GO:0010002</a></td>
<td align="left">cardioblast differentiation</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="odd">
<td>445</td>
<td align="left"><a href="GO:0010454" class="uri">GO:0010454</a></td>
<td align="left">negative regulation of cell fate commitm...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="even">
<td>446</td>
<td align="left"><a href="GO:0016330" class="uri">GO:0016330</a></td>
<td align="left">second mitotic wave involved in compound...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="odd">
<td>447</td>
<td align="left"><a href="GO:0016348" class="uri">GO:0016348</a></td>
<td align="left">imaginal disc-derived leg joint morphoge...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="even">
<td>448</td>
<td align="left"><a href="GO:0018917" class="uri">GO:0018917</a></td>
<td align="left">fluorene metabolic process</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="odd">
<td>449</td>
<td align="left"><a href="GO:0021910" class="uri">GO:0021910</a></td>
<td align="left">smoothened signaling pathway involved in...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="even">
<td>450</td>
<td align="left"><a href="GO:0021943" class="uri">GO:0021943</a></td>
<td align="left">formation of radial glial scaffolds</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="odd">
<td>451</td>
<td align="left"><a href="GO:0030046" class="uri">GO:0030046</a></td>
<td align="left">parallel actin filament bundle assembly</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="even">
<td>452</td>
<td align="left"><a href="GO:0030713" class="uri">GO:0030713</a></td>
<td align="left">ovarian follicle cell stalk formation</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="odd">
<td>453</td>
<td align="left"><a href="GO:0032211" class="uri">GO:0032211</a></td>
<td align="left">negative regulation of telomere maintena...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="even">
<td>454</td>
<td align="left"><a href="GO:0035126" class="uri">GO:0035126</a></td>
<td align="left">post-embryonic genitalia morphogenesis</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="odd">
<td>455</td>
<td align="left"><a href="GO:0035155" class="uri">GO:0035155</a></td>
<td align="left">negative regulation of terminal cell fat...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="even">
<td>456</td>
<td align="left"><a href="GO:0035156" class="uri">GO:0035156</a></td>
<td align="left">fusion cell fate specification</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="odd">
<td>457</td>
<td align="left"><a href="GO:0035157" class="uri">GO:0035157</a></td>
<td align="left">negative regulation of fusion cell fate ...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="even">
<td>458</td>
<td align="left"><a href="GO:0035163" class="uri">GO:0035163</a></td>
<td align="left">embryonic hemocyte differentiation</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="odd">
<td>459</td>
<td align="left"><a href="GO:0035165" class="uri">GO:0035165</a></td>
<td align="left">embryonic crystal cell differentiation</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="even">
<td>460</td>
<td align="left"><a href="GO:0035170" class="uri">GO:0035170</a></td>
<td align="left">lymph gland crystal cell differentiation</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="odd">
<td>461</td>
<td align="left"><a href="GO:0035336" class="uri">GO:0035336</a></td>
<td align="left">long-chain fatty-acyl-CoA metabolic proc...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="even">
<td>462</td>
<td align="left"><a href="GO:0035338" class="uri">GO:0035338</a></td>
<td align="left">long-chain fatty-acyl-CoA biosynthetic p...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="odd">
<td>463</td>
<td align="left"><a href="GO:0042480" class="uri">GO:0042480</a></td>
<td align="left">negative regulation of eye photoreceptor...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="even">
<td>464</td>
<td align="left"><a href="GO:0042676" class="uri">GO:0042676</a></td>
<td align="left">compound eye cone cell fate commitment</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="odd">
<td>465</td>
<td align="left"><a href="GO:0042684" class="uri">GO:0042684</a></td>
<td align="left">cardioblast cell fate commitment</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="even">
<td>466</td>
<td align="left"><a href="GO:0042685" class="uri">GO:0042685</a></td>
<td align="left">cardioblast cell fate specification</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="odd">
<td>467</td>
<td align="left"><a href="GO:0042686" class="uri">GO:0042686</a></td>
<td align="left">regulation of cardioblast cell fate spec...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="even">
<td>468</td>
<td align="left"><a href="GO:0042688" class="uri">GO:0042688</a></td>
<td align="left">crystal cell differentiation</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="odd">
<td>469</td>
<td align="left"><a href="GO:0042689" class="uri">GO:0042689</a></td>
<td align="left">regulation of crystal cell differentiati...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="even">
<td>470</td>
<td align="left"><a href="GO:0042738" class="uri">GO:0042738</a></td>
<td align="left">exogenous drug catabolic process</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="odd">
<td>471</td>
<td align="left"><a href="GO:0043537" class="uri">GO:0043537</a></td>
<td align="left">negative regulation of blood vessel endo...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="even">
<td>472</td>
<td align="left"><a href="GO:0043576" class="uri">GO:0043576</a></td>
<td align="left">regulation of respiratory gaseous exchan...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="odd">
<td>473</td>
<td align="left"><a href="GO:0044806" class="uri">GO:0044806</a></td>
<td align="left">G-quadruplex DNA unwinding</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="even">
<td>474</td>
<td align="left"><a href="GO:0045316" class="uri">GO:0045316</a></td>
<td align="left">negative regulation of compound eye phot...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="odd">
<td>475</td>
<td align="left"><a href="GO:0045468" class="uri">GO:0045468</a></td>
<td align="left">regulation of R8 cell spacing in compoun...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="even">
<td>476</td>
<td align="left"><a href="GO:0046533" class="uri">GO:0046533</a></td>
<td align="left">negative regulation of photoreceptor cel...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="odd">
<td>477</td>
<td align="left"><a href="GO:0046602" class="uri">GO:0046602</a></td>
<td align="left">regulation of mitotic centrosome separat...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="even">
<td>478</td>
<td align="left"><a href="GO:0048052" class="uri">GO:0048052</a></td>
<td align="left">R1/R6 cell differentiation</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="odd">
<td>479</td>
<td align="left"><a href="GO:0048619" class="uri">GO:0048619</a></td>
<td align="left">embryonic hindgut morphogenesis</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="even">
<td>480</td>
<td align="left"><a href="GO:0048803" class="uri">GO:0048803</a></td>
<td align="left">imaginal disc-derived male genitalia mor...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="odd">
<td>481</td>
<td align="left"><a href="GO:0048805" class="uri">GO:0048805</a></td>
<td align="left">imaginal disc-derived genitalia morphoge...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="even">
<td>482</td>
<td align="left"><a href="GO:0048808" class="uri">GO:0048808</a></td>
<td align="left">male genitalia morphogenesis</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="odd">
<td>483</td>
<td align="left"><a href="GO:0051256" class="uri">GO:0051256</a></td>
<td align="left">mitotic spindle midzone assembly</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="even">
<td>484</td>
<td align="left"><a href="GO:0051890" class="uri">GO:0051890</a></td>
<td align="left">regulation of cardioblast differentiatio...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="odd">
<td>485</td>
<td align="left"><a href="GO:0051974" class="uri">GO:0051974</a></td>
<td align="left">negative regulation of telomerase activi...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="even">
<td>486</td>
<td align="left"><a href="GO:0060250" class="uri">GO:0060250</a></td>
<td align="left">germ-line stem-cell niche homeostasis</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="odd">
<td>487</td>
<td align="left"><a href="GO:0060288" class="uri">GO:0060288</a></td>
<td align="left">formation of a compartment boundary</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="even">
<td>488</td>
<td align="left"><a href="GO:0060911" class="uri">GO:0060911</a></td>
<td align="left">cardiac cell fate commitment</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="odd">
<td>489</td>
<td align="left"><a href="GO:0060912" class="uri">GO:0060912</a></td>
<td align="left">cardiac cell fate specification</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="even">
<td>490</td>
<td align="left"><a href="GO:0061042" class="uri">GO:0061042</a></td>
<td align="left">vascular wound healing</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="odd">
<td>491</td>
<td align="left"><a href="GO:0061331" class="uri">GO:0061331</a></td>
<td align="left">epithelial cell proliferation involved i...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="even">
<td>492</td>
<td align="left"><a href="GO:0070986" class="uri">GO:0070986</a></td>
<td align="left">left/right axis specification</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="odd">
<td>493</td>
<td align="left"><a href="GO:0072499" class="uri">GO:0072499</a></td>
<td align="left">photoreceptor cell axon guidance</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="even">
<td>494</td>
<td align="left"><a href="GO:0072578" class="uri">GO:0072578</a></td>
<td align="left">neurotransmitter-gated ion channel clust...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="odd">
<td>495</td>
<td align="left"><a href="GO:0097114" class="uri">GO:0097114</a></td>
<td align="left">NMDA glutamate receptor clustering</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="even">
<td>496</td>
<td align="left"><a href="GO:0097475" class="uri">GO:0097475</a></td>
<td align="left">motor neuron migration</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="odd">
<td>497</td>
<td align="left"><a href="GO:1900451" class="uri">GO:1900451</a></td>
<td align="left">positive regulation of glutamate recepto...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="even">
<td>498</td>
<td align="left"><a href="GO:2000043" class="uri">GO:2000043</a></td>
<td align="left">regulation of cardiac cell fate specific...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="odd">
<td>499</td>
<td align="left"><a href="GO:2000047" class="uri">GO:2000047</a></td>
<td align="left">regulation of cell-cell adhesion mediate...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="even">
<td>500</td>
<td align="left"><a href="GO:2000048" class="uri">GO:2000048</a></td>
<td align="left">negative regulation of cell-cell adhesio...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="odd">
<td>501</td>
<td align="left"><a href="GO:2000969" class="uri">GO:2000969</a></td>
<td align="left">positive regulation of alpha-amino-3-hyd...</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.67</td>
<td align="left">0.01111</td>
</tr>
<tr class="even">
<td>502</td>
<td align="left"><a href="GO:0032467" class="uri">GO:0032467</a></td>
<td align="left">positive regulation of cytokinesis</td>
<td align="right">10</td>
<td align="right">6</td>
<td align="right">2.23</td>
<td align="left">0.01113</td>
</tr>
<tr class="odd">
<td>503</td>
<td align="left"><a href="GO:0045466" class="uri">GO:0045466</a></td>
<td align="left">R7 cell differentiation</td>
<td align="right">10</td>
<td align="right">6</td>
<td align="right">2.23</td>
<td align="left">0.01113</td>
</tr>
<tr class="even">
<td>504</td>
<td align="left"><a href="GO:0055094" class="uri">GO:0055094</a></td>
<td align="left">response to lipoprotein particle</td>
<td align="right">10</td>
<td align="right">6</td>
<td align="right">2.23</td>
<td align="left">0.01113</td>
</tr>
<tr class="odd">
<td>505</td>
<td align="left"><a href="GO:0055098" class="uri">GO:0055098</a></td>
<td align="left">response to low-density lipoprotein part...</td>
<td align="right">10</td>
<td align="right">6</td>
<td align="right">2.23</td>
<td align="left">0.01113</td>
</tr>
<tr class="even">
<td>506</td>
<td align="left"><a href="GO:0060074" class="uri">GO:0060074</a></td>
<td align="left">synapse maturation</td>
<td align="right">10</td>
<td align="right">6</td>
<td align="right">2.23</td>
<td align="left">0.01113</td>
</tr>
<tr class="odd">
<td>507</td>
<td align="left"><a href="GO:0007612" class="uri">GO:0007612</a></td>
<td align="left">learning</td>
<td align="right">41</td>
<td align="right">16</td>
<td align="right">9.15</td>
<td align="left">0.01156</td>
</tr>
<tr class="even">
<td>508</td>
<td align="left"><a href="GO:0033559" class="uri">GO:0033559</a></td>
<td align="left">unsaturated fatty acid metabolic process</td>
<td align="right">41</td>
<td align="right">16</td>
<td align="right">9.15</td>
<td align="left">0.01156</td>
</tr>
<tr class="odd">
<td>509</td>
<td align="left"><a href="GO:0045765" class="uri">GO:0045765</a></td>
<td align="left">regulation of angiogenesis</td>
<td align="right">41</td>
<td align="right">16</td>
<td align="right">9.15</td>
<td align="left">0.01156</td>
</tr>
<tr class="even">
<td>510</td>
<td align="left"><a href="GO:0018108" class="uri">GO:0018108</a></td>
<td align="left">peptidyl-tyrosine phosphorylation</td>
<td align="right">178</td>
<td align="right">53</td>
<td align="right">39.73</td>
<td align="left">0.01159</td>
</tr>
<tr class="odd">
<td>511</td>
<td align="left"><a href="GO:0046661" class="uri">GO:0046661</a></td>
<td align="left">male sex differentiation</td>
<td align="right">28</td>
<td align="right">12</td>
<td align="right">6.25</td>
<td align="left">0.01210</td>
</tr>
<tr class="even">
<td>512</td>
<td align="left"><a href="GO:0060042" class="uri">GO:0060042</a></td>
<td align="left">retina morphogenesis in camera-type eye</td>
<td align="right">28</td>
<td align="right">12</td>
<td align="right">6.25</td>
<td align="left">0.01210</td>
</tr>
<tr class="odd">
<td>513</td>
<td align="left"><a href="GO:0051493" class="uri">GO:0051493</a></td>
<td align="left">regulation of cytoskeleton organization</td>
<td align="right">121</td>
<td align="right">38</td>
<td align="right">27.01</td>
<td align="left">0.01249</td>
</tr>
<tr class="even">
<td>514</td>
<td align="left"><a href="GO:0007552" class="uri">GO:0007552</a></td>
<td align="left">metamorphosis</td>
<td align="right">25</td>
<td align="right">11</td>
<td align="right">5.58</td>
<td align="left">0.01285</td>
</tr>
<tr class="odd">
<td>515</td>
<td align="left"><a href="GO:0007560" class="uri">GO:0007560</a></td>
<td align="left">imaginal disc morphogenesis</td>
<td align="right">25</td>
<td align="right">11</td>
<td align="right">5.58</td>
<td align="left">0.01285</td>
</tr>
<tr class="even">
<td>516</td>
<td align="left"><a href="GO:0044243" class="uri">GO:0044243</a></td>
<td align="left">multicellular organismal catabolic proce...</td>
<td align="right">25</td>
<td align="right">11</td>
<td align="right">5.58</td>
<td align="left">0.01285</td>
</tr>
<tr class="odd">
<td>517</td>
<td align="left"><a href="GO:0048707" class="uri">GO:0048707</a></td>
<td align="left">instar larval or pupal morphogenesis</td>
<td align="right">25</td>
<td align="right">11</td>
<td align="right">5.58</td>
<td align="left">0.01285</td>
</tr>
<tr class="even">
<td>518</td>
<td align="left"><a href="GO:0048799" class="uri">GO:0048799</a></td>
<td align="left">organ maturation</td>
<td align="right">25</td>
<td align="right">11</td>
<td align="right">5.58</td>
<td align="left">0.01285</td>
</tr>
<tr class="odd">
<td>519</td>
<td align="left"><a href="GO:0070977" class="uri">GO:0070977</a></td>
<td align="left">bone maturation</td>
<td align="right">25</td>
<td align="right">11</td>
<td align="right">5.58</td>
<td align="left">0.01285</td>
</tr>
<tr class="even">
<td>520</td>
<td align="left"><a href="GO:0090183" class="uri">GO:0090183</a></td>
<td align="left">regulation of kidney development</td>
<td align="right">25</td>
<td align="right">11</td>
<td align="right">5.58</td>
<td align="left">0.01285</td>
</tr>
<tr class="odd">
<td>521</td>
<td align="left"><a href="GO:0007566" class="uri">GO:0007566</a></td>
<td align="left">embryo implantation</td>
<td align="right">13</td>
<td align="right">7</td>
<td align="right">2.90</td>
<td align="left">0.01297</td>
</tr>
<tr class="even">
<td>522</td>
<td align="left"><a href="GO:0009954" class="uri">GO:0009954</a></td>
<td align="left">proximal/distal pattern formation</td>
<td align="right">13</td>
<td align="right">7</td>
<td align="right">2.90</td>
<td align="left">0.01297</td>
</tr>
<tr class="odd">
<td>523</td>
<td align="left"><a href="GO:0061512" class="uri">GO:0061512</a></td>
<td align="left">protein localization to cilium</td>
<td align="right">13</td>
<td align="right">7</td>
<td align="right">2.90</td>
<td align="left">0.01297</td>
</tr>
<tr class="even">
<td>524</td>
<td align="left"><a href="GO:0071260" class="uri">GO:0071260</a></td>
<td align="left">cellular response to mechanical stimulus</td>
<td align="right">13</td>
<td align="right">7</td>
<td align="right">2.90</td>
<td align="left">0.01297</td>
</tr>
<tr class="odd">
<td>525</td>
<td align="left"><a href="GO:0044770" class="uri">GO:0044770</a></td>
<td align="left">cell cycle phase transition</td>
<td align="right">179</td>
<td align="right">53</td>
<td align="right">39.96</td>
<td align="left">0.01300</td>
</tr>
<tr class="even">
<td>526</td>
<td align="left"><a href="GO:0050673" class="uri">GO:0050673</a></td>
<td align="left">epithelial cell proliferation</td>
<td align="right">84</td>
<td align="right">28</td>
<td align="right">18.75</td>
<td align="left">0.01301</td>
</tr>
<tr class="odd">
<td>527</td>
<td align="left"><a href="GO:0031214" class="uri">GO:0031214</a></td>
<td align="left">biomineral tissue development</td>
<td align="right">22</td>
<td align="right">10</td>
<td align="right">4.91</td>
<td align="left">0.01345</td>
</tr>
<tr class="even">
<td>528</td>
<td align="left"><a href="GO:0072160" class="uri">GO:0072160</a></td>
<td align="left">nephron tubule epithelial cell different...</td>
<td align="right">22</td>
<td align="right">10</td>
<td align="right">4.91</td>
<td align="left">0.01345</td>
</tr>
<tr class="odd">
<td>529</td>
<td align="left"><a href="GO:0051962" class="uri">GO:0051962</a></td>
<td align="left">positive regulation of nervous system de...</td>
<td align="right">156</td>
<td align="right">47</td>
<td align="right">34.82</td>
<td align="left">0.01346</td>
</tr>
<tr class="even">
<td>530</td>
<td align="left"><a href="GO:0006213" class="uri">GO:0006213</a></td>
<td align="left">pyrimidine nucleoside metabolic process</td>
<td align="right">16</td>
<td align="right">8</td>
<td align="right">3.57</td>
<td align="left">0.01373</td>
</tr>
<tr class="odd">
<td>531</td>
<td align="left"><a href="GO:0006692" class="uri">GO:0006692</a></td>
<td align="left">prostanoid metabolic process</td>
<td align="right">16</td>
<td align="right">8</td>
<td align="right">3.57</td>
<td align="left">0.01373</td>
</tr>
<tr class="even">
<td>532</td>
<td align="left"><a href="GO:0006693" class="uri">GO:0006693</a></td>
<td align="left">prostaglandin metabolic process</td>
<td align="right">16</td>
<td align="right">8</td>
<td align="right">3.57</td>
<td align="left">0.01373</td>
</tr>
<tr class="odd">
<td>533</td>
<td align="left"><a href="GO:0003382" class="uri">GO:0003382</a></td>
<td align="left">epithelial cell morphogenesis</td>
<td align="right">19</td>
<td align="right">9</td>
<td align="right">4.24</td>
<td align="left">0.01380</td>
</tr>
<tr class="even">
<td>534</td>
<td align="left"><a href="GO:0008584" class="uri">GO:0008584</a></td>
<td align="left">male gonad development</td>
<td align="right">19</td>
<td align="right">9</td>
<td align="right">4.24</td>
<td align="left">0.01380</td>
</tr>
<tr class="odd">
<td>535</td>
<td align="left"><a href="GO:0034394" class="uri">GO:0034394</a></td>
<td align="left">protein localization to cell surface</td>
<td align="right">19</td>
<td align="right">9</td>
<td align="right">4.24</td>
<td align="left">0.01380</td>
</tr>
<tr class="even">
<td>536</td>
<td align="left"><a href="GO:0046546" class="uri">GO:0046546</a></td>
<td align="left">development of primary male sexual chara...</td>
<td align="right">19</td>
<td align="right">9</td>
<td align="right">4.24</td>
<td align="left">0.01380</td>
</tr>
<tr class="odd">
<td>537</td>
<td align="left"><a href="GO:0009798" class="uri">GO:0009798</a></td>
<td align="left">axis specification</td>
<td align="right">32</td>
<td align="right">13</td>
<td align="right">7.14</td>
<td align="left">0.01522</td>
</tr>
<tr class="even">
<td>538</td>
<td align="left"><a href="GO:0030336" class="uri">GO:0030336</a></td>
<td align="left">negative regulation of cell migration</td>
<td align="right">32</td>
<td align="right">13</td>
<td align="right">7.14</td>
<td align="left">0.01522</td>
</tr>
<tr class="odd">
<td>539</td>
<td align="left"><a href="GO:0008202" class="uri">GO:0008202</a></td>
<td align="left">steroid metabolic process</td>
<td align="right">89</td>
<td align="right">29</td>
<td align="right">19.87</td>
<td align="left">0.01615</td>
</tr>
<tr class="even">
<td>540</td>
<td align="left"><a href="GO:0006690" class="uri">GO:0006690</a></td>
<td align="left">icosanoid metabolic process</td>
<td align="right">39</td>
<td align="right">15</td>
<td align="right">8.71</td>
<td align="left">0.01649</td>
</tr>
<tr class="odd">
<td>541</td>
<td align="left"><a href="GO:1901568" class="uri">GO:1901568</a></td>
<td align="left">fatty acid derivative metabolic process</td>
<td align="right">39</td>
<td align="right">15</td>
<td align="right">8.71</td>
<td align="left">0.01649</td>
</tr>
<tr class="even">
<td>542</td>
<td align="left"><a href="GO:0048563" class="uri">GO:0048563</a></td>
<td align="left">post-embryonic organ morphogenesis</td>
<td align="right">29</td>
<td align="right">12</td>
<td align="right">6.47</td>
<td align="left">0.01655</td>
</tr>
<tr class="odd">
<td>543</td>
<td align="left"><a href="GO:0003015" class="uri">GO:0003015</a></td>
<td align="left">heart process</td>
<td align="right">60</td>
<td align="right">21</td>
<td align="right">13.39</td>
<td align="left">0.01665</td>
</tr>
<tr class="even">
<td>544</td>
<td align="left"><a href="GO:0000712" class="uri">GO:0000712</a></td>
<td align="left">resolution of meiotic recombination inte...</td>
<td align="right">8</td>
<td align="right">5</td>
<td align="right">1.79</td>
<td align="left">0.01675</td>
</tr>
<tr class="odd">
<td>545</td>
<td align="left"><a href="GO:0001955" class="uri">GO:0001955</a></td>
<td align="left">blood vessel maturation</td>
<td align="right">8</td>
<td align="right">5</td>
<td align="right">1.79</td>
<td align="left">0.01675</td>
</tr>
<tr class="even">
<td>546</td>
<td align="left"><a href="GO:0030913" class="uri">GO:0030913</a></td>
<td align="left">paranodal junction assembly</td>
<td align="right">8</td>
<td align="right">5</td>
<td align="right">1.79</td>
<td align="left">0.01675</td>
</tr>
<tr class="odd">
<td>547</td>
<td align="left"><a href="GO:0035418" class="uri">GO:0035418</a></td>
<td align="left">protein localization to synapse</td>
<td align="right">8</td>
<td align="right">5</td>
<td align="right">1.79</td>
<td align="left">0.01675</td>
</tr>
<tr class="even">
<td>548</td>
<td align="left"><a href="GO:0051764" class="uri">GO:0051764</a></td>
<td align="left">actin crosslink formation</td>
<td align="right">8</td>
<td align="right">5</td>
<td align="right">1.79</td>
<td align="left">0.01675</td>
</tr>
<tr class="odd">
<td>549</td>
<td align="left"><a href="GO:0098735" class="uri">GO:0098735</a></td>
<td align="left">positive regulation of the force of hear...</td>
<td align="right">8</td>
<td align="right">5</td>
<td align="right">1.79</td>
<td align="left">0.01675</td>
</tr>
<tr class="even">
<td>550</td>
<td align="left"><a href="GO:0070507" class="uri">GO:0070507</a></td>
<td align="left">regulation of microtubule cytoskeleton o...</td>
<td align="right">53</td>
<td align="right">19</td>
<td align="right">11.83</td>
<td align="left">0.01700</td>
</tr>
<tr class="odd">
<td>551</td>
<td align="left"><a href="GO:0006935" class="uri">GO:0006935</a></td>
<td align="left">chemotaxis</td>
<td align="right">229</td>
<td align="right">65</td>
<td align="right">51.12</td>
<td align="left">0.01709</td>
</tr>
<tr class="even">
<td>552</td>
<td align="left"><a href="GO:0072661" class="uri">GO:0072661</a></td>
<td align="left">protein targeting to plasma membrane</td>
<td align="right">26</td>
<td align="right">11</td>
<td align="right">5.80</td>
<td align="left">0.01787</td>
</tr>
<tr class="odd">
<td>553</td>
<td align="left"><a href="GO:0097306" class="uri">GO:0097306</a></td>
<td align="left">cellular response to alcohol</td>
<td align="right">26</td>
<td align="right">11</td>
<td align="right">5.80</td>
<td align="left">0.01787</td>
</tr>
<tr class="even">
<td>554</td>
<td align="left"><a href="GO:0098609" class="uri">GO:0098609</a></td>
<td align="left">cell-cell adhesion</td>
<td align="right">351</td>
<td align="right">95</td>
<td align="right">78.35</td>
<td align="left">0.01811</td>
</tr>
<tr class="odd">
<td>555</td>
<td align="left"><a href="GO:0022617" class="uri">GO:0022617</a></td>
<td align="left">extracellular matrix disassembly</td>
<td align="right">36</td>
<td align="right">14</td>
<td align="right">8.04</td>
<td align="left">0.01823</td>
</tr>
<tr class="even">
<td>556</td>
<td align="left"><a href="GO:0008610" class="uri">GO:0008610</a></td>
<td align="left">lipid biosynthetic process</td>
<td align="right">202</td>
<td align="right">58</td>
<td align="right">45.09</td>
<td align="left">0.01854</td>
</tr>
<tr class="odd">
<td>557</td>
<td align="left"><a href="GO:0001738" class="uri">GO:0001738</a></td>
<td align="left">morphogenesis of a polarized epithelium</td>
<td align="right">43</td>
<td align="right">16</td>
<td align="right">9.60</td>
<td align="left">0.01894</td>
</tr>
<tr class="even">
<td>558</td>
<td align="left"><a href="GO:0007584" class="uri">GO:0007584</a></td>
<td align="left">response to nutrient</td>
<td align="right">50</td>
<td align="right">18</td>
<td align="right">11.16</td>
<td align="left">0.01900</td>
</tr>
<tr class="odd">
<td>559</td>
<td align="left"><a href="GO:0007009" class="uri">GO:0007009</a></td>
<td align="left">plasma membrane organization</td>
<td align="right">128</td>
<td align="right">39</td>
<td align="right">28.57</td>
<td align="left">0.01912</td>
</tr>
<tr class="even">
<td>560</td>
<td align="left"><a href="GO:0048699" class="uri">GO:0048699</a></td>
<td align="left">generation of neurons</td>
<td align="right">576</td>
<td align="right">149</td>
<td align="right">128.58</td>
<td align="left">0.01947</td>
</tr>
<tr class="odd">
<td>561</td>
<td align="left"><a href="GO:0021700" class="uri">GO:0021700</a></td>
<td align="left">developmental maturation</td>
<td align="right">113</td>
<td align="right">35</td>
<td align="right">25.22</td>
<td align="left">0.01991</td>
</tr>
<tr class="even">
<td>562</td>
<td align="left"><a href="GO:0009147" class="uri">GO:0009147</a></td>
<td align="left">pyrimidine nucleoside triphosphate metab...</td>
<td align="right">11</td>
<td align="right">6</td>
<td align="right">2.46</td>
<td align="left">0.01992</td>
</tr>
<tr class="odd">
<td>563</td>
<td align="left"><a href="GO:0042908" class="uri">GO:0042908</a></td>
<td align="left">xenobiotic transport</td>
<td align="right">11</td>
<td align="right">6</td>
<td align="right">2.46</td>
<td align="left">0.01992</td>
</tr>
<tr class="even">
<td>564</td>
<td align="left"><a href="GO:0045747" class="uri">GO:0045747</a></td>
<td align="left">positive regulation of Notch signaling p...</td>
<td align="right">11</td>
<td align="right">6</td>
<td align="right">2.46</td>
<td align="left">0.01992</td>
</tr>
<tr class="odd">
<td>565</td>
<td align="left"><a href="GO:0055013" class="uri">GO:0055013</a></td>
<td align="left">cardiac muscle cell development</td>
<td align="right">33</td>
<td align="right">13</td>
<td align="right">7.37</td>
<td align="left">0.02010</td>
</tr>
<tr class="even">
<td>566</td>
<td align="left"><a href="GO:0022008" class="uri">GO:0022008</a></td>
<td align="left">neurogenesis</td>
<td align="right">606</td>
<td align="right">156</td>
<td align="right">135.28</td>
<td align="left">0.02018</td>
</tr>
<tr class="odd">
<td>567</td>
<td align="left"><a href="GO:0042060" class="uri">GO:0042060</a></td>
<td align="left">wound healing</td>
<td align="right">148</td>
<td align="right">44</td>
<td align="right">33.04</td>
<td align="left">0.02071</td>
</tr>
<tr class="even">
<td>568</td>
<td align="left"><a href="GO:0003338" class="uri">GO:0003338</a></td>
<td align="left">metanephros morphogenesis</td>
<td align="right">17</td>
<td align="right">8</td>
<td align="right">3.79</td>
<td align="left">0.02095</td>
</tr>
<tr class="odd">
<td>569</td>
<td align="left"><a href="GO:0048546" class="uri">GO:0048546</a></td>
<td align="left">digestive tract morphogenesis</td>
<td align="right">17</td>
<td align="right">8</td>
<td align="right">3.79</td>
<td align="left">0.02095</td>
</tr>
<tr class="even">
<td>570</td>
<td align="left"><a href="GO:0001516" class="uri">GO:0001516</a></td>
<td align="left">prostaglandin biosynthetic process</td>
<td align="right">14</td>
<td align="right">7</td>
<td align="right">3.13</td>
<td align="left">0.02102</td>
</tr>
<tr class="odd">
<td>571</td>
<td align="left"><a href="GO:0001941" class="uri">GO:0001941</a></td>
<td align="left">postsynaptic membrane organization</td>
<td align="right">14</td>
<td align="right">7</td>
<td align="right">3.13</td>
<td align="left">0.02102</td>
</tr>
<tr class="even">
<td>572</td>
<td align="left"><a href="GO:0015936" class="uri">GO:0015936</a></td>
<td align="left">coenzyme A metabolic process</td>
<td align="right">14</td>
<td align="right">7</td>
<td align="right">3.13</td>
<td align="left">0.02102</td>
</tr>
<tr class="odd">
<td>573</td>
<td align="left"><a href="GO:0046456" class="uri">GO:0046456</a></td>
<td align="left">icosanoid biosynthetic process</td>
<td align="right">14</td>
<td align="right">7</td>
<td align="right">3.13</td>
<td align="left">0.02102</td>
</tr>
<tr class="even">
<td>574</td>
<td align="left"><a href="GO:0046457" class="uri">GO:0046457</a></td>
<td align="left">prostanoid biosynthetic process</td>
<td align="right">14</td>
<td align="right">7</td>
<td align="right">3.13</td>
<td align="left">0.02102</td>
</tr>
<tr class="odd">
<td>575</td>
<td align="left"><a href="GO:1901570" class="uri">GO:1901570</a></td>
<td align="left">fatty acid derivative biosynthetic proce...</td>
<td align="right">14</td>
<td align="right">7</td>
<td align="right">3.13</td>
<td align="left">0.02102</td>
</tr>
<tr class="even">
<td>576</td>
<td align="left"><a href="GO:0072522" class="uri">GO:0072522</a></td>
<td align="left">purine-containing compound biosynthetic ...</td>
<td align="right">87</td>
<td align="right">28</td>
<td align="right">19.42</td>
<td align="left">0.02116</td>
</tr>
<tr class="odd">
<td>577</td>
<td align="left"><a href="GO:2000145" class="uri">GO:2000145</a></td>
<td align="left">regulation of cell motility</td>
<td align="right">160</td>
<td align="right">47</td>
<td align="right">35.72</td>
<td align="left">0.02139</td>
</tr>
<tr class="even">
<td>578</td>
<td align="left"><a href="GO:0044772" class="uri">GO:0044772</a></td>
<td align="left">mitotic cell cycle phase transition</td>
<td align="right">168</td>
<td align="right">49</td>
<td align="right">37.50</td>
<td align="left">0.02175</td>
</tr>
<tr class="odd">
<td>579</td>
<td align="left"><a href="GO:0031623" class="uri">GO:0031623</a></td>
<td align="left">receptor internalization</td>
<td align="right">30</td>
<td align="right">12</td>
<td align="right">6.70</td>
<td align="left">0.02211</td>
</tr>
<tr class="even">
<td>580</td>
<td align="left"><a href="GO:0060047" class="uri">GO:0060047</a></td>
<td align="left">heart contraction</td>
<td align="right">58</td>
<td align="right">20</td>
<td align="right">12.95</td>
<td align="left">0.02266</td>
</tr>
<tr class="odd">
<td>581</td>
<td align="left"><a href="GO:0007219" class="uri">GO:0007219</a></td>
<td align="left">Notch signaling pathway</td>
<td align="right">69</td>
<td align="right">23</td>
<td align="right">15.40</td>
<td align="left">0.02302</td>
</tr>
<tr class="even">
<td>582</td>
<td align="left"><a href="GO:0050795" class="uri">GO:0050795</a></td>
<td align="left">regulation of behavior</td>
<td align="right">69</td>
<td align="right">23</td>
<td align="right">15.40</td>
<td align="left">0.02302</td>
</tr>
<tr class="odd">
<td>583</td>
<td align="left"><a href="GO:0007585" class="uri">GO:0007585</a></td>
<td align="left">respiratory gaseous exchange</td>
<td align="right">37</td>
<td align="right">14</td>
<td align="right">8.26</td>
<td align="left">0.02346</td>
</tr>
<tr class="even">
<td>584</td>
<td align="left"><a href="GO:0050808" class="uri">GO:0050808</a></td>
<td align="left">synapse organization</td>
<td align="right">153</td>
<td align="right">45</td>
<td align="right">34.15</td>
<td align="left">0.02347</td>
</tr>
<tr class="odd">
<td>585</td>
<td align="left"><a href="GO:0001503" class="uri">GO:0001503</a></td>
<td align="left">ossification</td>
<td align="right">122</td>
<td align="right">37</td>
<td align="right">27.23</td>
<td align="left">0.02370</td>
</tr>
<tr class="even">
<td>586</td>
<td align="left"><a href="GO:0097305" class="uri">GO:0097305</a></td>
<td align="left">response to alcohol</td>
<td align="right">84</td>
<td align="right">27</td>
<td align="right">18.75</td>
<td align="left">0.02370</td>
</tr>
<tr class="odd">
<td>587</td>
<td align="left"><a href="GO:0055001" class="uri">GO:0055001</a></td>
<td align="left">muscle cell development</td>
<td align="right">73</td>
<td align="right">24</td>
<td align="right">16.30</td>
<td align="left">0.02428</td>
</tr>
<tr class="even">
<td>588</td>
<td align="left"><a href="GO:0048562" class="uri">GO:0048562</a></td>
<td align="left">embryonic organ morphogenesis</td>
<td align="right">107</td>
<td align="right">33</td>
<td align="right">23.89</td>
<td align="left">0.02479</td>
</tr>
<tr class="odd">
<td>589</td>
<td align="left"><a href="GO:0000912" class="uri">GO:0000912</a></td>
<td align="left">assembly of actomyosin apparatus involve...</td>
<td align="right">6</td>
<td align="right">4</td>
<td align="right">1.34</td>
<td align="left">0.02512</td>
</tr>
<tr class="even">
<td>590</td>
<td align="left"><a href="GO:0000915" class="uri">GO:0000915</a></td>
<td align="left">actomyosin contractile ring assembly</td>
<td align="right">6</td>
<td align="right">4</td>
<td align="right">1.34</td>
<td align="left">0.02512</td>
</tr>
<tr class="odd">
<td>591</td>
<td align="left"><a href="GO:0001556" class="uri">GO:0001556</a></td>
<td align="left">oocyte maturation</td>
<td align="right">6</td>
<td align="right">4</td>
<td align="right">1.34</td>
<td align="left">0.02512</td>
</tr>
<tr class="even">
<td>592</td>
<td align="left"><a href="GO:0003352" class="uri">GO:0003352</a></td>
<td align="left">regulation of cilium movement</td>
<td align="right">6</td>
<td align="right">4</td>
<td align="right">1.34</td>
<td align="left">0.02512</td>
</tr>
<tr class="odd">
<td>593</td>
<td align="left"><a href="GO:0006085" class="uri">GO:0006085</a></td>
<td align="left">acetyl-CoA biosynthetic process</td>
<td align="right">6</td>
<td align="right">4</td>
<td align="right">1.34</td>
<td align="left">0.02512</td>
</tr>
<tr class="even">
<td>594</td>
<td align="left"><a href="GO:0007440" class="uri">GO:0007440</a></td>
<td align="left">foregut morphogenesis</td>
<td align="right">6</td>
<td align="right">4</td>
<td align="right">1.34</td>
<td align="left">0.02512</td>
</tr>
<tr class="odd">
<td>595</td>
<td align="left"><a href="GO:0021904" class="uri">GO:0021904</a></td>
<td align="left">dorsal/ventral neural tube patterning</td>
<td align="right">6</td>
<td align="right">4</td>
<td align="right">1.34</td>
<td align="left">0.02512</td>
</tr>
<tr class="even">
<td>596</td>
<td align="left"><a href="GO:0035214" class="uri">GO:0035214</a></td>
<td align="left">eye-antennal disc development</td>
<td align="right">6</td>
<td align="right">4</td>
<td align="right">1.34</td>
<td align="left">0.02512</td>
</tr>
<tr class="odd">
<td>597</td>
<td align="left"><a href="GO:0035735" class="uri">GO:0035735</a></td>
<td align="left">intraciliary transport involved in ciliu...</td>
<td align="right">6</td>
<td align="right">4</td>
<td align="right">1.34</td>
<td align="left">0.02512</td>
</tr>
<tr class="even">
<td>598</td>
<td align="left"><a href="GO:0035999" class="uri">GO:0035999</a></td>
<td align="left">tetrahydrofolate interconversion</td>
<td align="right">6</td>
<td align="right">4</td>
<td align="right">1.34</td>
<td align="left">0.02512</td>
</tr>
<tr class="odd">
<td>599</td>
<td align="left"><a href="GO:0042693" class="uri">GO:0042693</a></td>
<td align="left">muscle cell fate commitment</td>
<td align="right">6</td>
<td align="right">4</td>
<td align="right">1.34</td>
<td align="left">0.02512</td>
</tr>
<tr class="even">
<td>600</td>
<td align="left"><a href="GO:0042771" class="uri">GO:0042771</a></td>
<td align="left">intrinsic apoptotic signaling pathway in...</td>
<td align="right">6</td>
<td align="right">4</td>
<td align="right">1.34</td>
<td align="left">0.02512</td>
</tr>
<tr class="odd">
<td>601</td>
<td align="left"><a href="GO:0044837" class="uri">GO:0044837</a></td>
<td align="left">actomyosin contractile ring organization</td>
<td align="right">6</td>
<td align="right">4</td>
<td align="right">1.34</td>
<td align="left">0.02512</td>
</tr>
<tr class="even">
<td>602</td>
<td align="left"><a href="GO:0045940" class="uri">GO:0045940</a></td>
<td align="left">positive regulation of steroid metabolic...</td>
<td align="right">6</td>
<td align="right">4</td>
<td align="right">1.34</td>
<td align="left">0.02512</td>
</tr>
<tr class="odd">
<td>603</td>
<td align="left"><a href="GO:0060831" class="uri">GO:0060831</a></td>
<td align="left">smoothened signaling pathway involved in...</td>
<td align="right">6</td>
<td align="right">4</td>
<td align="right">1.34</td>
<td align="left">0.02512</td>
</tr>
<tr class="even">
<td>604</td>
<td align="left"><a href="GO:0090205" class="uri">GO:0090205</a></td>
<td align="left">positive regulation of cholesterol metab...</td>
<td align="right">6</td>
<td align="right">4</td>
<td align="right">1.34</td>
<td align="left">0.02512</td>
</tr>
<tr class="odd">
<td>605</td>
<td align="left"><a href="GO:1900271" class="uri">GO:1900271</a></td>
<td align="left">regulation of long-term synaptic potenti...</td>
<td align="right">6</td>
<td align="right">4</td>
<td align="right">1.34</td>
<td align="left">0.02512</td>
</tr>
<tr class="even">
<td>606</td>
<td align="left"><a href="GO:1902001" class="uri">GO:1902001</a></td>
<td align="left">fatty acid transmembrane transport</td>
<td align="right">6</td>
<td align="right">4</td>
<td align="right">1.34</td>
<td align="left">0.02512</td>
</tr>
<tr class="odd">
<td>607</td>
<td align="left"><a href="GO:0016331" class="uri">GO:0016331</a></td>
<td align="left">morphogenesis of embryonic epithelium</td>
<td align="right">55</td>
<td align="right">19</td>
<td align="right">12.28</td>
<td align="left">0.02543</td>
</tr>
<tr class="even">
<td>608</td>
<td align="left"><a href="GO:0044843" class="uri">GO:0044843</a></td>
<td align="left">cell cycle G1/S phase transition</td>
<td align="right">92</td>
<td align="right">29</td>
<td align="right">20.54</td>
<td align="left">0.02547</td>
</tr>
<tr class="odd">
<td>609</td>
<td align="left"><a href="GO:0006869" class="uri">GO:0006869</a></td>
<td align="left">lipid transport</td>
<td align="right">142</td>
<td align="right">42</td>
<td align="right">31.70</td>
<td align="left">0.02548</td>
</tr>
<tr class="even">
<td>610</td>
<td align="left"><a href="GO:0007166" class="uri">GO:0007166</a></td>
<td align="left">cell surface receptor signaling pathway</td>
<td align="right">801</td>
<td align="right">201</td>
<td align="right">178.80</td>
<td align="left">0.02553</td>
</tr>
<tr class="odd">
<td>611</td>
<td align="left"><a href="GO:0055006" class="uri">GO:0055006</a></td>
<td align="left">cardiac cell development</td>
<td align="right">34</td>
<td align="right">13</td>
<td align="right">7.59</td>
<td align="left">0.02607</td>
</tr>
<tr class="even">
<td>612</td>
<td align="left"><a href="GO:2000146" class="uri">GO:2000146</a></td>
<td align="left">negative regulation of cell motility</td>
<td align="right">34</td>
<td align="right">13</td>
<td align="right">7.59</td>
<td align="left">0.02607</td>
</tr>
<tr class="odd">
<td>613</td>
<td align="left"><a href="GO:0030334" class="uri">GO:0030334</a></td>
<td align="left">regulation of cell migration</td>
<td align="right">154</td>
<td align="right">45</td>
<td align="right">34.38</td>
<td align="left">0.02621</td>
</tr>
<tr class="even">
<td>614</td>
<td align="left"><a href="GO:0060348" class="uri">GO:0060348</a></td>
<td align="left">bone development</td>
<td align="right">81</td>
<td align="right">26</td>
<td align="right">18.08</td>
<td align="left">0.02656</td>
</tr>
<tr class="odd">
<td>615</td>
<td align="left"><a href="GO:1901607" class="uri">GO:1901607</a></td>
<td align="left">alpha-amino acid biosynthetic process</td>
<td align="right">41</td>
<td align="right">15</td>
<td align="right">9.15</td>
<td align="left">0.02658</td>
</tr>
<tr class="even">
<td>616</td>
<td align="left"><a href="GO:0016477" class="uri">GO:0016477</a></td>
<td align="left">cell migration</td>
<td align="right">398</td>
<td align="right">105</td>
<td align="right">88.84</td>
<td align="left">0.02740</td>
</tr>
<tr class="odd">
<td>617</td>
<td align="left"><a href="GO:0001654" class="uri">GO:0001654</a></td>
<td align="left">eye development</td>
<td align="right">143</td>
<td align="right">42</td>
<td align="right">31.92</td>
<td align="left">0.02853</td>
</tr>
<tr class="even">
<td>618</td>
<td align="left"><a href="GO:0006084" class="uri">GO:0006084</a></td>
<td align="left">acetyl-CoA metabolic process</td>
<td align="right">21</td>
<td align="right">9</td>
<td align="right">4.69</td>
<td align="left">0.02854</td>
</tr>
<tr class="odd">
<td>619</td>
<td align="left"><a href="GO:0009064" class="uri">GO:0009064</a></td>
<td align="left">glutamine family amino acid metabolic pr...</td>
<td align="right">21</td>
<td align="right">9</td>
<td align="right">4.69</td>
<td align="left">0.02854</td>
</tr>
<tr class="even">
<td>620</td>
<td align="left"><a href="GO:0030282" class="uri">GO:0030282</a></td>
<td align="left">bone mineralization</td>
<td align="right">21</td>
<td align="right">9</td>
<td align="right">4.69</td>
<td align="left">0.02854</td>
</tr>
<tr class="odd">
<td>621</td>
<td align="left"><a href="GO:0043534" class="uri">GO:0043534</a></td>
<td align="left">blood vessel endothelial cell migration</td>
<td align="right">21</td>
<td align="right">9</td>
<td align="right">4.69</td>
<td align="left">0.02854</td>
</tr>
<tr class="even">
<td>622</td>
<td align="left"><a href="GO:0072182" class="uri">GO:0072182</a></td>
<td align="left">regulation of nephron tubule epithelial ...</td>
<td align="right">21</td>
<td align="right">9</td>
<td align="right">4.69</td>
<td align="left">0.02854</td>
</tr>
<tr class="odd">
<td>623</td>
<td align="left"><a href="GO:2000696" class="uri">GO:2000696</a></td>
<td align="left">regulation of epithelial cell differenti...</td>
<td align="right">21</td>
<td align="right">9</td>
<td align="right">4.69</td>
<td align="left">0.02854</td>
</tr>
<tr class="even">
<td>624</td>
<td align="left"><a href="GO:0046530" class="uri">GO:0046530</a></td>
<td align="left">photoreceptor cell differentiation</td>
<td align="right">52</td>
<td align="right">18</td>
<td align="right">11.61</td>
<td align="left">0.02856</td>
</tr>
<tr class="odd">
<td>625</td>
<td align="left"><a href="GO:0008544" class="uri">GO:0008544</a></td>
<td align="left">epidermis development</td>
<td align="right">74</td>
<td align="right">24</td>
<td align="right">16.52</td>
<td align="left">0.02857</td>
</tr>
<tr class="even">
<td>626</td>
<td align="left"><a href="GO:0002165" class="uri">GO:0002165</a></td>
<td align="left">instar larval or pupal development</td>
<td align="right">31</td>
<td align="right">12</td>
<td align="right">6.92</td>
<td align="left">0.02894</td>
</tr>
<tr class="odd">
<td>627</td>
<td align="left"><a href="GO:0008299" class="uri">GO:0008299</a></td>
<td align="left">isoprenoid biosynthetic process</td>
<td align="right">31</td>
<td align="right">12</td>
<td align="right">6.92</td>
<td align="left">0.02894</td>
</tr>
<tr class="even">
<td>628</td>
<td align="left"><a href="GO:0010810" class="uri">GO:0010810</a></td>
<td align="left">regulation of cell-substrate adhesion</td>
<td align="right">45</td>
<td align="right">16</td>
<td align="right">10.05</td>
<td align="left">0.02945</td>
</tr>
<tr class="odd">
<td>629</td>
<td align="left"><a href="GO:1901342" class="uri">GO:1901342</a></td>
<td align="left">regulation of vasculature development</td>
<td align="right">45</td>
<td align="right">16</td>
<td align="right">10.05</td>
<td align="left">0.02945</td>
</tr>
<tr class="even">
<td>630</td>
<td align="left"><a href="GO:0001501" class="uri">GO:0001501</a></td>
<td align="left">skeletal system development</td>
<td align="right">163</td>
<td align="right">47</td>
<td align="right">36.39</td>
<td align="left">0.02951</td>
</tr>
<tr class="odd">
<td>631</td>
<td align="left"><a href="GO:0035725" class="uri">GO:0035725</a></td>
<td align="left">sodium ion transmembrane transport</td>
<td align="right">38</td>
<td align="right">14</td>
<td align="right">8.48</td>
<td align="left">0.02973</td>
</tr>
<tr class="even">
<td>632</td>
<td align="left"><a href="GO:0042472" class="uri">GO:0042472</a></td>
<td align="left">inner ear morphogenesis</td>
<td align="right">38</td>
<td align="right">14</td>
<td align="right">8.48</td>
<td align="left">0.02973</td>
</tr>
<tr class="odd">
<td>633</td>
<td align="left"><a href="GO:0030500" class="uri">GO:0030500</a></td>
<td align="left">regulation of bone mineralization</td>
<td align="right">18</td>
<td align="right">8</td>
<td align="right">4.02</td>
<td align="left">0.03050</td>
</tr>
<tr class="even">
<td>634</td>
<td align="left"><a href="GO:0035994" class="uri">GO:0035994</a></td>
<td align="left">response to muscle stretch</td>
<td align="right">18</td>
<td align="right">8</td>
<td align="right">4.02</td>
<td align="left">0.03050</td>
</tr>
<tr class="odd">
<td>635</td>
<td align="left"><a href="GO:0070167" class="uri">GO:0070167</a></td>
<td align="left">regulation of biomineral tissue developm...</td>
<td align="right">18</td>
<td align="right">8</td>
<td align="right">4.02</td>
<td align="left">0.03050</td>
</tr>
<tr class="even">
<td>636</td>
<td align="left"><a href="GO:0040013" class="uri">GO:0040013</a></td>
<td align="left">negative regulation of locomotion</td>
<td align="right">56</td>
<td align="right">19</td>
<td align="right">12.50</td>
<td align="left">0.03067</td>
</tr>
<tr class="odd">
<td>637</td>
<td align="left"><a href="GO:0003334" class="uri">GO:0003334</a></td>
<td align="left">keratinocyte development</td>
<td align="right">9</td>
<td align="right">5</td>
<td align="right">2.01</td>
<td align="left">0.03084</td>
</tr>
<tr class="even">
<td>638</td>
<td align="left"><a href="GO:0009208" class="uri">GO:0009208</a></td>
<td align="left">pyrimidine ribonucleoside triphosphate m...</td>
<td align="right">9</td>
<td align="right">5</td>
<td align="right">2.01</td>
<td align="left">0.03084</td>
</tr>
<tr class="odd">
<td>639</td>
<td align="left"><a href="GO:0042335" class="uri">GO:0042335</a></td>
<td align="left">cuticle development</td>
<td align="right">9</td>
<td align="right">5</td>
<td align="right">2.01</td>
<td align="left">0.03084</td>
</tr>
<tr class="even">
<td>640</td>
<td align="left"><a href="GO:0044331" class="uri">GO:0044331</a></td>
<td align="left">cell-cell adhesion mediated by cadherin</td>
<td align="right">9</td>
<td align="right">5</td>
<td align="right">2.01</td>
<td align="left">0.03084</td>
</tr>
<tr class="odd">
<td>641</td>
<td align="left"><a href="GO:0071402" class="uri">GO:0071402</a></td>
<td align="left">cellular response to lipoprotein particl...</td>
<td align="right">9</td>
<td align="right">5</td>
<td align="right">2.01</td>
<td align="left">0.03084</td>
</tr>
<tr class="even">
<td>642</td>
<td align="left"><a href="GO:0071404" class="uri">GO:0071404</a></td>
<td align="left">cellular response to low-density lipopro...</td>
<td align="right">9</td>
<td align="right">5</td>
<td align="right">2.01</td>
<td align="left">0.03084</td>
</tr>
<tr class="odd">
<td>643</td>
<td align="left"><a href="GO:0001667" class="uri">GO:0001667</a></td>
<td align="left">ameboidal-type cell migration</td>
<td align="right">82</td>
<td align="right">26</td>
<td align="right">18.30</td>
<td align="left">0.03089</td>
</tr>
<tr class="even">
<td>644</td>
<td align="left"><a href="GO:0032365" class="uri">GO:0032365</a></td>
<td align="left">intracellular lipid transport</td>
<td align="right">15</td>
<td align="right">7</td>
<td align="right">3.35</td>
<td align="left">0.03197</td>
</tr>
<tr class="odd">
<td>645</td>
<td align="left"><a href="GO:0050922" class="uri">GO:0050922</a></td>
<td align="left">negative regulation of chemotaxis</td>
<td align="right">15</td>
<td align="right">7</td>
<td align="right">3.35</td>
<td align="left">0.03197</td>
</tr>
<tr class="even">
<td>646</td>
<td align="left"><a href="GO:0010811" class="uri">GO:0010811</a></td>
<td align="left">positive regulation of cell-substrate ad...</td>
<td align="right">28</td>
<td align="right">11</td>
<td align="right">6.25</td>
<td align="left">0.03207</td>
</tr>
<tr class="odd">
<td>647</td>
<td align="left"><a href="GO:0008652" class="uri">GO:0008652</a></td>
<td align="left">cellular amino acid biosynthetic process</td>
<td align="right">49</td>
<td align="right">17</td>
<td align="right">10.94</td>
<td align="left">0.03209</td>
</tr>
<tr class="even">
<td>648</td>
<td align="left"><a href="GO:0046131" class="uri">GO:0046131</a></td>
<td align="left">pyrimidine ribonucleoside metabolic proc...</td>
<td align="right">12</td>
<td align="right">6</td>
<td align="right">2.68</td>
<td align="left">0.03243</td>
</tr>
<tr class="odd">
<td>649</td>
<td align="left"><a href="GO:0007186" class="uri">GO:0007186</a></td>
<td align="left">G-protein coupled receptor signaling pat...</td>
<td align="right">208</td>
<td align="right">58</td>
<td align="right">46.43</td>
<td align="left">0.03278</td>
</tr>
<tr class="even">
<td>650</td>
<td align="left"><a href="GO:0006959" class="uri">GO:0006959</a></td>
<td align="left">humoral immune response</td>
<td align="right">42</td>
<td align="right">15</td>
<td align="right">9.38</td>
<td align="left">0.03306</td>
</tr>
<tr class="odd">
<td>651</td>
<td align="left"><a href="GO:0002040" class="uri">GO:0002040</a></td>
<td align="left">sprouting angiogenesis</td>
<td align="right">35</td>
<td align="right">13</td>
<td align="right">7.81</td>
<td align="left">0.03324</td>
</tr>
<tr class="even">
<td>652</td>
<td align="left"><a href="GO:0007444" class="uri">GO:0007444</a></td>
<td align="left">imaginal disc development</td>
<td align="right">35</td>
<td align="right">13</td>
<td align="right">7.81</td>
<td align="left">0.03324</td>
</tr>
<tr class="odd">
<td>653</td>
<td align="left"><a href="GO:0050678" class="uri">GO:0050678</a></td>
<td align="left">regulation of epithelial cell proliferat...</td>
<td align="right">64</td>
<td align="right">21</td>
<td align="right">14.29</td>
<td align="left">0.03439</td>
</tr>
<tr class="even">
<td>654</td>
<td align="left"><a href="GO:0007610" class="uri">GO:0007610</a></td>
<td align="left">behavior</td>
<td align="right">245</td>
<td align="right">67</td>
<td align="right">54.69</td>
<td align="left">0.03446</td>
</tr>
<tr class="odd">
<td>655</td>
<td align="left"><a href="GO:0014070" class="uri">GO:0014070</a></td>
<td align="left">response to organic cyclic compound</td>
<td align="right">245</td>
<td align="right">67</td>
<td align="right">54.69</td>
<td align="left">0.03446</td>
</tr>
<tr class="even">
<td>656</td>
<td align="left"><a href="GO:0021987" class="uri">GO:0021987</a></td>
<td align="left">cerebral cortex development</td>
<td align="right">53</td>
<td align="right">18</td>
<td align="right">11.83</td>
<td align="left">0.03450</td>
</tr>
<tr class="odd">
<td>657</td>
<td align="left"><a href="GO:0006801" class="uri">GO:0006801</a></td>
<td align="left">superoxide metabolic process</td>
<td align="right">25</td>
<td align="right">10</td>
<td align="right">5.58</td>
<td align="left">0.03543</td>
</tr>
<tr class="even">
<td>658</td>
<td align="left"><a href="GO:0008045" class="uri">GO:0008045</a></td>
<td align="left">motor neuron axon guidance</td>
<td align="right">25</td>
<td align="right">10</td>
<td align="right">5.58</td>
<td align="left">0.03543</td>
</tr>
<tr class="odd">
<td>659</td>
<td align="left"><a href="GO:0055002" class="uri">GO:0055002</a></td>
<td align="left">striated muscle cell development</td>
<td align="right">68</td>
<td align="right">22</td>
<td align="right">15.18</td>
<td align="left">0.03601</td>
</tr>
<tr class="even">
<td>660</td>
<td align="left"><a href="GO:0060401" class="uri">GO:0060401</a></td>
<td align="left">cytosolic calcium ion transport</td>
<td align="right">46</td>
<td align="right">16</td>
<td align="right">10.27</td>
<td align="left">0.03609</td>
</tr>
<tr class="odd">
<td>661</td>
<td align="left"><a href="GO:0060402" class="uri">GO:0060402</a></td>
<td align="left">calcium ion transport into cytosol</td>
<td align="right">46</td>
<td align="right">16</td>
<td align="right">10.27</td>
<td align="left">0.03609</td>
</tr>
<tr class="even">
<td>662</td>
<td align="left"><a href="GO:1901615" class="uri">GO:1901615</a></td>
<td align="left">organic hydroxy compound metabolic proce...</td>
<td align="right">181</td>
<td align="right">51</td>
<td align="right">40.40</td>
<td align="left">0.03626</td>
</tr>
<tr class="odd">
<td>663</td>
<td align="left"><a href="GO:0000082" class="uri">GO:0000082</a></td>
<td align="left">G1/S transition of mitotic cell cycle</td>
<td align="right">87</td>
<td align="right">27</td>
<td align="right">19.42</td>
<td align="left">0.03675</td>
</tr>
<tr class="even">
<td>664</td>
<td align="left"><a href="GO:0000730" class="uri">GO:0000730</a></td>
<td align="left">DNA recombinase assembly</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="odd">
<td>665</td>
<td align="left"><a href="GO:0001755" class="uri">GO:0001755</a></td>
<td align="left">neural crest cell migration</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="even">
<td>666</td>
<td align="left"><a href="GO:0001893" class="uri">GO:0001893</a></td>
<td align="left">maternal placenta development</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="odd">
<td>667</td>
<td align="left"><a href="GO:0002052" class="uri">GO:0002052</a></td>
<td align="left">positive regulation of neuroblast prolif...</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="even">
<td>668</td>
<td align="left"><a href="GO:0002313" class="uri">GO:0002313</a></td>
<td align="left">mature B cell differentiation involved i...</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="odd">
<td>669</td>
<td align="left"><a href="GO:0002921" class="uri">GO:0002921</a></td>
<td align="left">negative regulation of humoral immune re...</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="even">
<td>670</td>
<td align="left"><a href="GO:0006968" class="uri">GO:0006968</a></td>
<td align="left">cellular defense response</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="odd">
<td>671</td>
<td align="left"><a href="GO:0007304" class="uri">GO:0007304</a></td>
<td align="left">chorion-containing eggshell formation</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="even">
<td>672</td>
<td align="left"><a href="GO:0007306" class="uri">GO:0007306</a></td>
<td align="left">eggshell chorion assembly</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="odd">
<td>673</td>
<td align="left"><a href="GO:0007386" class="uri">GO:0007386</a></td>
<td align="left">compartment pattern specification</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="even">
<td>674</td>
<td align="left"><a href="GO:0007400" class="uri">GO:0007400</a></td>
<td align="left">neuroblast fate determination</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="odd">
<td>675</td>
<td align="left"><a href="GO:0007403" class="uri">GO:0007403</a></td>
<td align="left">glial cell fate determination</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="even">
<td>676</td>
<td align="left"><a href="GO:0007474" class="uri">GO:0007474</a></td>
<td align="left">imaginal disc-derived wing vein specific...</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="odd">
<td>677</td>
<td align="left"><a href="GO:0009186" class="uri">GO:0009186</a></td>
<td align="left">deoxyribonucleoside diphosphate metaboli...</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="even">
<td>678</td>
<td align="left"><a href="GO:0009265" class="uri">GO:0009265</a></td>
<td align="left">2'-deoxyribonucleotide biosynthetic proc...</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="odd">
<td>679</td>
<td align="left"><a href="GO:0009608" class="uri">GO:0009608</a></td>
<td align="left">response to symbiont</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="even">
<td>680</td>
<td align="left"><a href="GO:0010710" class="uri">GO:0010710</a></td>
<td align="left">regulation of collagen catabolic process</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="odd">
<td>681</td>
<td align="left"><a href="GO:0016128" class="uri">GO:0016128</a></td>
<td align="left">phytosteroid metabolic process</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="even">
<td>682</td>
<td align="left"><a href="GO:0016129" class="uri">GO:0016129</a></td>
<td align="left">phytosteroid biosynthetic process</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="odd">
<td>683</td>
<td align="left"><a href="GO:0018879" class="uri">GO:0018879</a></td>
<td align="left">biphenyl metabolic process</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="even">
<td>684</td>
<td align="left"><a href="GO:0021781" class="uri">GO:0021781</a></td>
<td align="left">glial cell fate commitment</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="odd">
<td>685</td>
<td align="left"><a href="GO:0021932" class="uri">GO:0021932</a></td>
<td align="left">hindbrain radial glia guided cell migrat...</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="even">
<td>686</td>
<td align="left"><a href="GO:0030703" class="uri">GO:0030703</a></td>
<td align="left">eggshell formation</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="odd">
<td>687</td>
<td align="left"><a href="GO:0030720" class="uri">GO:0030720</a></td>
<td align="left">oocyte localization involved in germariu...</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="even">
<td>688</td>
<td align="left"><a href="GO:0035017" class="uri">GO:0035017</a></td>
<td align="left">cuticle pattern formation</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="odd">
<td>689</td>
<td align="left"><a href="GO:0035153" class="uri">GO:0035153</a></td>
<td align="left">epithelial cell type specification, open...</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="even">
<td>690</td>
<td align="left"><a href="GO:0035154" class="uri">GO:0035154</a></td>
<td align="left">terminal cell fate specification, open t...</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="odd">
<td>691</td>
<td align="left"><a href="GO:0035608" class="uri">GO:0035608</a></td>
<td align="left">protein deglutamylation</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="even">
<td>692</td>
<td align="left"><a href="GO:0036335" class="uri">GO:0036335</a></td>
<td align="left">intestinal stem cell homeostasis</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="odd">
<td>693</td>
<td align="left"><a href="GO:0039019" class="uri">GO:0039019</a></td>
<td align="left">pronephric nephron development</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="even">
<td>694</td>
<td align="left"><a href="GO:0039020" class="uri">GO:0039020</a></td>
<td align="left">pronephric nephron tubule development</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="odd">
<td>695</td>
<td align="left"><a href="GO:0040036" class="uri">GO:0040036</a></td>
<td align="left">regulation of fibroblast growth factor r...</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="even">
<td>696</td>
<td align="left"><a href="GO:0042148" class="uri">GO:0042148</a></td>
<td align="left">strand invasion</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="odd">
<td>697</td>
<td align="left"><a href="GO:0042478" class="uri">GO:0042478</a></td>
<td align="left">regulation of eye photoreceptor cell dev...</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="even">
<td>698</td>
<td align="left"><a href="GO:0042675" class="uri">GO:0042675</a></td>
<td align="left">compound eye cone cell differentiation</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="odd">
<td>699</td>
<td align="left"><a href="GO:0042737" class="uri">GO:0042737</a></td>
<td align="left">drug catabolic process</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="even">
<td>700</td>
<td align="left"><a href="GO:0043052" class="uri">GO:0043052</a></td>
<td align="left">thermotaxis</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="odd">
<td>701</td>
<td align="left"><a href="GO:0045314" class="uri">GO:0045314</a></td>
<td align="left">regulation of compound eye photoreceptor...</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="even">
<td>702</td>
<td align="left"><a href="GO:0046385" class="uri">GO:0046385</a></td>
<td align="left">deoxyribose phosphate biosynthetic proce...</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="odd">
<td>703</td>
<td align="left"><a href="GO:0046618" class="uri">GO:0046618</a></td>
<td align="left">drug export</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="even">
<td>704</td>
<td align="left"><a href="GO:0046843" class="uri">GO:0046843</a></td>
<td align="left">dorsal appendage formation</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="odd">
<td>705</td>
<td align="left"><a href="GO:0048867" class="uri">GO:0048867</a></td>
<td align="left">stem cell fate determination</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="even">
<td>706</td>
<td align="left"><a href="GO:0051255" class="uri">GO:0051255</a></td>
<td align="left">spindle midzone assembly</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="odd">
<td>707</td>
<td align="left"><a href="GO:0051481" class="uri">GO:0051481</a></td>
<td align="left">negative regulation of cytosolic calcium...</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="even">
<td>708</td>
<td align="left"><a href="GO:0055129" class="uri">GO:0055129</a></td>
<td align="left">L-proline biosynthetic process</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="odd">
<td>709</td>
<td align="left"><a href="GO:0060055" class="uri">GO:0060055</a></td>
<td align="left">angiogenesis involved in wound healing</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="even">
<td>710</td>
<td align="left"><a href="GO:0060136" class="uri">GO:0060136</a></td>
<td align="left">embryonic process involved in female pre...</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="odd">
<td>711</td>
<td align="left"><a href="GO:0061003" class="uri">GO:0061003</a></td>
<td align="left">positive regulation of dendritic spine m...</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="even">
<td>712</td>
<td align="left"><a href="GO:0061382" class="uri">GO:0061382</a></td>
<td align="left">Malpighian tubule tip cell differentiati...</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="odd">
<td>713</td>
<td align="left"><a href="GO:0071205" class="uri">GO:0071205</a></td>
<td align="left">protein localization to juxtaparanode re...</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="even">
<td>714</td>
<td align="left"><a href="GO:0071361" class="uri">GO:0071361</a></td>
<td align="left">cellular response to ethanol</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="odd">
<td>715</td>
<td align="left"><a href="GO:0071398" class="uri">GO:0071398</a></td>
<td align="left">cellular response to fatty acid</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="even">
<td>716</td>
<td align="left"><a href="GO:0071625" class="uri">GO:0071625</a></td>
<td align="left">vocalization behavior</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="odd">
<td>717</td>
<td align="left"><a href="GO:2000311" class="uri">GO:2000311</a></td>
<td align="left">regulation of alpha-amino-3-hydroxy-5-me...</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.89</td>
<td align="left">0.03700</td>
</tr>
<tr class="even">
<td>718</td>
<td align="left"><a href="GO:0042471" class="uri">GO:0042471</a></td>
<td align="left">ear morphogenesis</td>
<td align="right">39</td>
<td align="right">14</td>
<td align="right">8.71</td>
<td align="left">0.03713</td>
</tr>
<tr class="odd">
<td>719</td>
<td align="left"><a href="GO:0010923" class="uri">GO:0010923</a></td>
<td align="left">negative regulation of phosphatase activ...</td>
<td align="right">22</td>
<td align="right">9</td>
<td align="right">4.91</td>
<td align="left">0.03899</td>
</tr>
<tr class="even">
<td>720</td>
<td align="left"><a href="GO:0009792" class="uri">GO:0009792</a></td>
<td align="left">embryo development ending in birth or eg...</td>
<td align="right">255</td>
<td align="right">69</td>
<td align="right">56.92</td>
<td align="left">0.03992</td>
</tr>
<tr class="odd">
<td>721</td>
<td align="left"><a href="GO:0006164" class="uri">GO:0006164</a></td>
<td align="left">purine nucleotide biosynthetic process</td>
<td align="right">80</td>
<td align="right">25</td>
<td align="right">17.86</td>
<td align="left">0.04004</td>
</tr>
<tr class="even">
<td>722</td>
<td align="left"><a href="GO:0090068" class="uri">GO:0090068</a></td>
<td align="left">positive regulation of cell cycle proces...</td>
<td align="right">65</td>
<td align="right">21</td>
<td align="right">14.51</td>
<td align="left">0.04047</td>
</tr>
<tr class="odd">
<td>723</td>
<td align="left"><a href="GO:1903510" class="uri">GO:1903510</a></td>
<td align="left">mucopolysaccharide metabolic process</td>
<td align="right">43</td>
<td align="right">15</td>
<td align="right">9.60</td>
<td align="left">0.04063</td>
</tr>
<tr class="even">
<td>724</td>
<td align="left"><a href="GO:0030048" class="uri">GO:0030048</a></td>
<td align="left">actin filament-based movement</td>
<td align="right">29</td>
<td align="right">11</td>
<td align="right">6.47</td>
<td align="left">0.04155</td>
</tr>
<tr class="odd">
<td>725</td>
<td align="left"><a href="GO:0009124" class="uri">GO:0009124</a></td>
<td align="left">nucleoside monophosphate biosynthetic pr...</td>
<td align="right">36</td>
<td align="right">13</td>
<td align="right">8.04</td>
<td align="left">0.04173</td>
</tr>
<tr class="even">
<td>726</td>
<td align="left"><a href="GO:0006066" class="uri">GO:0006066</a></td>
<td align="left">alcohol metabolic process</td>
<td align="right">123</td>
<td align="right">36</td>
<td align="right">27.46</td>
<td align="left">0.04236</td>
</tr>
<tr class="odd">
<td>727</td>
<td align="left"><a href="GO:0006220" class="uri">GO:0006220</a></td>
<td align="left">pyrimidine nucleotide metabolic process</td>
<td align="right">19</td>
<td align="right">8</td>
<td align="right">4.24</td>
<td align="left">0.04264</td>
</tr>
<tr class="even">
<td>728</td>
<td align="left"><a href="GO:0048592" class="uri">GO:0048592</a></td>
<td align="left">eye morphogenesis</td>
<td align="right">92</td>
<td align="right">28</td>
<td align="right">20.54</td>
<td align="left">0.04300</td>
</tr>
<tr class="odd">
<td>729</td>
<td align="left"><a href="GO:0009913" class="uri">GO:0009913</a></td>
<td align="left">epidermal cell differentiation</td>
<td align="right">58</td>
<td align="right">19</td>
<td align="right">12.95</td>
<td align="left">0.04353</td>
</tr>
<tr class="even">
<td>730</td>
<td align="left"><a href="GO:0035282" class="uri">GO:0035282</a></td>
<td align="left">segmentation</td>
<td align="right">47</td>
<td align="right">16</td>
<td align="right">10.49</td>
<td align="left">0.04376</td>
</tr>
<tr class="odd">
<td>731</td>
<td align="left"><a href="GO:0030029" class="uri">GO:0030029</a></td>
<td align="left">actin filament-based process</td>
<td align="right">232</td>
<td align="right">63</td>
<td align="right">51.79</td>
<td align="left">0.04492</td>
</tr>
<tr class="even">
<td>732</td>
<td align="left"><a href="GO:0001101" class="uri">GO:0001101</a></td>
<td align="left">response to acid chemical</td>
<td align="right">62</td>
<td align="right">20</td>
<td align="right">13.84</td>
<td align="left">0.04552</td>
</tr>
<tr class="odd">
<td>733</td>
<td align="left"><a href="GO:0002011" class="uri">GO:0002011</a></td>
<td align="left">morphogenesis of an epithelial sheet</td>
<td align="right">16</td>
<td align="right">7</td>
<td align="right">3.57</td>
<td align="left">0.04616</td>
</tr>
<tr class="even">
<td>734</td>
<td align="left"><a href="GO:0010623" class="uri">GO:0010623</a></td>
<td align="left">developmental programmed cell death</td>
<td align="right">16</td>
<td align="right">7</td>
<td align="right">3.57</td>
<td align="left">0.04616</td>
</tr>
<tr class="odd">
<td>735</td>
<td align="left"><a href="GO:0032835" class="uri">GO:0032835</a></td>
<td align="left">glomerulus development</td>
<td align="right">16</td>
<td align="right">7</td>
<td align="right">3.57</td>
<td align="left">0.04616</td>
</tr>
<tr class="even">
<td>736</td>
<td align="left"><a href="GO:0071320" class="uri">GO:0071320</a></td>
<td align="left">cellular response to cAMP</td>
<td align="right">16</td>
<td align="right">7</td>
<td align="right">3.57</td>
<td align="left">0.04616</td>
</tr>
<tr class="odd">
<td>737</td>
<td align="left"><a href="GO:0009201" class="uri">GO:0009201</a></td>
<td align="left">ribonucleoside triphosphate biosynthetic...</td>
<td align="right">26</td>
<td align="right">10</td>
<td align="right">5.80</td>
<td align="left">0.04640</td>
</tr>
<tr class="even">
<td>738</td>
<td align="left"><a href="GO:0035308" class="uri">GO:0035308</a></td>
<td align="left">negative regulation of protein dephospho...</td>
<td align="right">26</td>
<td align="right">10</td>
<td align="right">5.80</td>
<td align="left">0.04640</td>
</tr>
<tr class="odd">
<td>739</td>
<td align="left"><a href="GO:0000910" class="uri">GO:0000910</a></td>
<td align="left">cytokinesis</td>
<td align="right">85</td>
<td align="right">26</td>
<td align="right">18.97</td>
<td align="left">0.04711</td>
</tr>
<tr class="even">
<td>740</td>
<td align="left"><a href="GO:0051146" class="uri">GO:0051146</a></td>
<td align="left">striated muscle cell differentiation</td>
<td align="right">85</td>
<td align="right">26</td>
<td align="right">18.97</td>
<td align="left">0.04711</td>
</tr>
<tr class="odd">
<td>741</td>
<td align="left"><a href="GO:1901293" class="uri">GO:1901293</a></td>
<td align="left">nucleoside phosphate biosynthetic proces...</td>
<td align="right">124</td>
<td align="right">36</td>
<td align="right">27.68</td>
<td align="left">0.04730</td>
</tr>
<tr class="even">
<td>742</td>
<td align="left"><a href="GO:0008285" class="uri">GO:0008285</a></td>
<td align="left">negative regulation of cell proliferatio...</td>
<td align="right">148</td>
<td align="right">42</td>
<td align="right">33.04</td>
<td align="left">0.04827</td>
</tr>
<tr class="odd">
<td>743</td>
<td align="left"><a href="GO:0001867" class="uri">GO:0001867</a></td>
<td align="left">complement activation, lectin pathway</td>
<td align="right">7</td>
<td align="right">4</td>
<td align="right">1.56</td>
<td align="left">0.04836</td>
</tr>
<tr class="even">
<td>744</td>
<td align="left"><a href="GO:0007195" class="uri">GO:0007195</a></td>
<td align="left">adenylate cyclase-inhibiting dopamine re...</td>
<td align="right">7</td>
<td align="right">4</td>
<td align="right">1.56</td>
<td align="left">0.04836</td>
</tr>
<tr class="odd">
<td>745</td>
<td align="left"><a href="GO:0007419" class="uri">GO:0007419</a></td>
<td align="left">ventral cord development</td>
<td align="right">7</td>
<td align="right">4</td>
<td align="right">1.56</td>
<td align="left">0.04836</td>
</tr>
<tr class="even">
<td>746</td>
<td align="left"><a href="GO:0010633" class="uri">GO:0010633</a></td>
<td align="left">negative regulation of epithelial cell m...</td>
<td align="right">7</td>
<td align="right">4</td>
<td align="right">1.56</td>
<td align="left">0.04836</td>
</tr>
<tr class="odd">
<td>747</td>
<td align="left"><a href="GO:0010824" class="uri">GO:0010824</a></td>
<td align="left">regulation of centrosome duplication</td>
<td align="right">7</td>
<td align="right">4</td>
<td align="right">1.56</td>
<td align="left">0.04836</td>
</tr>
<tr class="even">
<td>748</td>
<td align="left"><a href="GO:0014032" class="uri">GO:0014032</a></td>
<td align="left">neural crest cell development</td>
<td align="right">7</td>
<td align="right">4</td>
<td align="right">1.56</td>
<td align="left">0.04836</td>
</tr>
<tr class="odd">
<td>749</td>
<td align="left"><a href="GO:0016479" class="uri">GO:0016479</a></td>
<td align="left">negative regulation of transcription fro...</td>
<td align="right">7</td>
<td align="right">4</td>
<td align="right">1.56</td>
<td align="left">0.04836</td>
</tr>
<tr class="even">
<td>750</td>
<td align="left"><a href="GO:0021846" class="uri">GO:0021846</a></td>
<td align="left">cell proliferation in forebrain</td>
<td align="right">7</td>
<td align="right">4</td>
<td align="right">1.56</td>
<td align="left">0.04836</td>
</tr>
<tr class="odd">
<td>751</td>
<td align="left"><a href="GO:0034113" class="uri">GO:0034113</a></td>
<td align="left">heterotypic cell-cell adhesion</td>
<td align="right">7</td>
<td align="right">4</td>
<td align="right">1.56</td>
<td align="left">0.04836</td>
</tr>
<tr class="even">
<td>752</td>
<td align="left"><a href="GO:0046051" class="uri">GO:0046051</a></td>
<td align="left">UTP metabolic process</td>
<td align="right">7</td>
<td align="right">4</td>
<td align="right">1.56</td>
<td align="left">0.04836</td>
</tr>
<tr class="odd">
<td>753</td>
<td align="left"><a href="GO:0046415" class="uri">GO:0046415</a></td>
<td align="left">urate metabolic process</td>
<td align="right">7</td>
<td align="right">4</td>
<td align="right">1.56</td>
<td align="left">0.04836</td>
</tr>
<tr class="even">
<td>754</td>
<td align="left"><a href="GO:0046500" class="uri">GO:0046500</a></td>
<td align="left">S-adenosylmethionine metabolic process</td>
<td align="right">7</td>
<td align="right">4</td>
<td align="right">1.56</td>
<td align="left">0.04836</td>
</tr>
<tr class="odd">
<td>755</td>
<td align="left"><a href="GO:0050879" class="uri">GO:0050879</a></td>
<td align="left">multicellular organismal movement</td>
<td align="right">7</td>
<td align="right">4</td>
<td align="right">1.56</td>
<td align="left">0.04836</td>
</tr>
<tr class="even">
<td>756</td>
<td align="left"><a href="GO:0050881" class="uri">GO:0050881</a></td>
<td align="left">musculoskeletal movement</td>
<td align="right">7</td>
<td align="right">4</td>
<td align="right">1.56</td>
<td align="left">0.04836</td>
</tr>
<tr class="odd">
<td>757</td>
<td align="left"><a href="GO:0050909" class="uri">GO:0050909</a></td>
<td align="left">sensory perception of taste</td>
<td align="right">7</td>
<td align="right">4</td>
<td align="right">1.56</td>
<td align="left">0.04836</td>
</tr>
<tr class="even">
<td>758</td>
<td align="left"><a href="GO:0060632" class="uri">GO:0060632</a></td>
<td align="left">regulation of microtubule-based movement</td>
<td align="right">7</td>
<td align="right">4</td>
<td align="right">1.56</td>
<td align="left">0.04836</td>
</tr>
<tr class="odd">
<td>759</td>
<td align="left"><a href="GO:1904398" class="uri">GO:1904398</a></td>
<td align="left">positive regulation of neuromuscular jun...</td>
<td align="right">7</td>
<td align="right">4</td>
<td align="right">1.56</td>
<td align="left">0.04836</td>
</tr>
<tr class="even">
<td>760</td>
<td align="left"><a href="GO:0017144" class="uri">GO:0017144</a></td>
<td align="left">drug metabolic process</td>
<td align="right">13</td>
<td align="right">6</td>
<td align="right">2.90</td>
<td align="left">0.04911</td>
</tr>
<tr class="odd">
<td>761</td>
<td align="left"><a href="GO:0045168" class="uri">GO:0045168</a></td>
<td align="left">cell-cell signaling involved in cell fat...</td>
<td align="right">13</td>
<td align="right">6</td>
<td align="right">2.90</td>
<td align="left">0.04911</td>
</tr>
<tr class="even">
<td>762</td>
<td align="left"><a href="GO:0043542" class="uri">GO:0043542</a></td>
<td align="left">endothelial cell migration</td>
<td align="right">44</td>
<td align="right">15</td>
<td align="right">9.82</td>
<td align="left">0.04935</td>
</tr>
<tr class="odd">
<td>763</td>
<td align="left"><a href="GO:0048167" class="uri">GO:0048167</a></td>
<td align="left">regulation of synaptic plasticity</td>
<td align="right">44</td>
<td align="right">15</td>
<td align="right">9.82</td>
<td align="left">0.04935</td>
</tr>
<tr class="even">
<td>764</td>
<td align="left"><a href="GO:0032787" class="uri">GO:0032787</a></td>
<td align="left">monocarboxylic acid metabolic process</td>
<td align="right">221</td>
<td align="right">60</td>
<td align="right">49.33</td>
<td align="left">0.04958</td>
</tr>
<tr class="odd">
<td>765</td>
<td align="left"><a href="GO:0000105" class="uri">GO:0000105</a></td>
<td align="left">histidine biosynthetic process</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="even">
<td>766</td>
<td align="left"><a href="GO:0001315" class="uri">GO:0001315</a></td>
<td align="left">age-dependent response to reactive oxyge...</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="odd">
<td>767</td>
<td align="left"><a href="GO:0002551" class="uri">GO:0002551</a></td>
<td align="left">mast cell chemotaxis</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="even">
<td>768</td>
<td align="left"><a href="GO:0002924" class="uri">GO:0002924</a></td>
<td align="left">negative regulation of humoral immune re...</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="odd">
<td>769</td>
<td align="left"><a href="GO:0003174" class="uri">GO:0003174</a></td>
<td align="left">mitral valve development</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="even">
<td>770</td>
<td align="left"><a href="GO:0003183" class="uri">GO:0003183</a></td>
<td align="left">mitral valve morphogenesis</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="odd">
<td>771</td>
<td align="left"><a href="GO:0006122" class="uri">GO:0006122</a></td>
<td align="left">mitochondrial electron transport, ubiqui...</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="even">
<td>772</td>
<td align="left"><a href="GO:0006311" class="uri">GO:0006311</a></td>
<td align="left">meiotic gene conversion</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="odd">
<td>773</td>
<td align="left"><a href="GO:0006556" class="uri">GO:0006556</a></td>
<td align="left">S-adenosylmethionine biosynthetic proces...</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="even">
<td>774</td>
<td align="left"><a href="GO:0006601" class="uri">GO:0006601</a></td>
<td align="left">creatine biosynthetic process</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="odd">
<td>775</td>
<td align="left"><a href="GO:0007488" class="uri">GO:0007488</a></td>
<td align="left">histoblast morphogenesis</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="even">
<td>776</td>
<td align="left"><a href="GO:0009153" class="uri">GO:0009153</a></td>
<td align="left">purine deoxyribonucleotide biosynthetic ...</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="odd">
<td>777</td>
<td align="left"><a href="GO:0009182" class="uri">GO:0009182</a></td>
<td align="left">purine deoxyribonucleoside diphosphate m...</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="even">
<td>778</td>
<td align="left"><a href="GO:0009699" class="uri">GO:0009699</a></td>
<td align="left">phenylpropanoid biosynthetic process</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="odd">
<td>779</td>
<td align="left"><a href="GO:0009772" class="uri">GO:0009772</a></td>
<td align="left">photosynthetic electron transport in pho...</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="even">
<td>780</td>
<td align="left"><a href="GO:0010826" class="uri">GO:0010826</a></td>
<td align="left">negative regulation of centrosome duplic...</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="odd">
<td>781</td>
<td align="left"><a href="GO:0014819" class="uri">GO:0014819</a></td>
<td align="left">regulation of skeletal muscle contractio...</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="even">
<td>782</td>
<td align="left"><a href="GO:0015669" class="uri">GO:0015669</a></td>
<td align="left">gas transport</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="odd">
<td>783</td>
<td align="left"><a href="GO:0015670" class="uri">GO:0015670</a></td>
<td align="left">carbon dioxide transport</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="even">
<td>784</td>
<td align="left"><a href="GO:0016131" class="uri">GO:0016131</a></td>
<td align="left">brassinosteroid metabolic process</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="odd">
<td>785</td>
<td align="left"><a href="GO:0016132" class="uri">GO:0016132</a></td>
<td align="left">brassinosteroid biosynthetic process</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="even">
<td>786</td>
<td align="left"><a href="GO:0018094" class="uri">GO:0018094</a></td>
<td align="left">protein polyglycylation</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="odd">
<td>787</td>
<td align="left"><a href="GO:0019755" class="uri">GO:0019755</a></td>
<td align="left">one-carbon compound transport</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="even">
<td>788</td>
<td align="left"><a href="GO:0021914" class="uri">GO:0021914</a></td>
<td align="left">negative regulation of smoothened signal...</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="odd">
<td>789</td>
<td align="left"><a href="GO:0031938" class="uri">GO:0031938</a></td>
<td align="left">regulation of chromatin silencing at tel...</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="even">
<td>790</td>
<td align="left"><a href="GO:0032877" class="uri">GO:0032877</a></td>
<td align="left">positive regulation of DNA endoreduplica...</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="odd">
<td>791</td>
<td align="left"><a href="GO:0034501" class="uri">GO:0034501</a></td>
<td align="left">protein localization to kinetochore</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="even">
<td>792</td>
<td align="left"><a href="GO:0035283" class="uri">GO:0035283</a></td>
<td align="left">central nervous system segmentation</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="odd">
<td>793</td>
<td align="left"><a href="GO:0035284" class="uri">GO:0035284</a></td>
<td align="left">brain segmentation</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="even">
<td>794</td>
<td align="left"><a href="GO:0035610" class="uri">GO:0035610</a></td>
<td align="left">protein side chain deglutamylation</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="odd">
<td>795</td>
<td align="left"><a href="GO:0035720" class="uri">GO:0035720</a></td>
<td align="left">intraciliary anterograde transport</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="even">
<td>796</td>
<td align="left"><a href="GO:0040037" class="uri">GO:0040037</a></td>
<td align="left">negative regulation of fibroblast growth...</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="odd">
<td>797</td>
<td align="left"><a href="GO:0042078" class="uri">GO:0042078</a></td>
<td align="left">germ-line stem cell division</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="even">
<td>798</td>
<td align="left"><a href="GO:0042418" class="uri">GO:0042418</a></td>
<td align="left">epinephrine biosynthetic process</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="odd">
<td>799</td>
<td align="left"><a href="GO:0042759" class="uri">GO:0042759</a></td>
<td align="left">long-chain fatty acid biosynthetic proce...</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="even">
<td>800</td>
<td align="left"><a href="GO:0043589" class="uri">GO:0043589</a></td>
<td align="left">skin morphogenesis</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="odd">
<td>801</td>
<td align="left"><a href="GO:0045719" class="uri">GO:0045719</a></td>
<td align="left">negative regulation of glycogen biosynth...</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="even">
<td>802</td>
<td align="left"><a href="GO:0045725" class="uri">GO:0045725</a></td>
<td align="left">positive regulation of glycogen biosynth...</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="odd">
<td>803</td>
<td align="left"><a href="GO:0045977" class="uri">GO:0045977</a></td>
<td align="left">positive regulation of mitotic cell cycl...</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="even">
<td>804</td>
<td align="left"><a href="GO:0046056" class="uri">GO:0046056</a></td>
<td align="left">dADP metabolic process</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="odd">
<td>805</td>
<td align="left"><a href="GO:0046606" class="uri">GO:0046606</a></td>
<td align="left">negative regulation of centrosome cycle</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="even">
<td>806</td>
<td align="left"><a href="GO:0046956" class="uri">GO:0046956</a></td>
<td align="left">positive phototaxis</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="odd">
<td>807</td>
<td align="left"><a href="GO:0048014" class="uri">GO:0048014</a></td>
<td align="left">Tie signaling pathway</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="even">
<td>808</td>
<td align="left"><a href="GO:0048843" class="uri">GO:0048843</a></td>
<td align="left">negative regulation of axon extension in...</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="odd">
<td>809</td>
<td align="left"><a href="GO:0052128" class="uri">GO:0052128</a></td>
<td align="left">positive energy taxis</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="even">
<td>810</td>
<td align="left"><a href="GO:0052173" class="uri">GO:0052173</a></td>
<td align="left">response to defenses of other organism i...</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="odd">
<td>811</td>
<td align="left"><a href="GO:0052200" class="uri">GO:0052200</a></td>
<td align="left">response to host defenses</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="even">
<td>812</td>
<td align="left"><a href="GO:0052564" class="uri">GO:0052564</a></td>
<td align="left">response to immune response of other org...</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="odd">
<td>813</td>
<td align="left"><a href="GO:0052572" class="uri">GO:0052572</a></td>
<td align="left">response to host immune response</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="even">
<td>814</td>
<td align="left"><a href="GO:0055075" class="uri">GO:0055075</a></td>
<td align="left">potassium ion homeostasis</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="odd">
<td>815</td>
<td align="left"><a href="GO:0055113" class="uri">GO:0055113</a></td>
<td align="left">epiboly involved in gastrulation with mo...</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="even">
<td>816</td>
<td align="left"><a href="GO:0060753" class="uri">GO:0060753</a></td>
<td align="left">regulation of mast cell chemotaxis</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="odd">
<td>817</td>
<td align="left"><a href="GO:0060754" class="uri">GO:0060754</a></td>
<td align="left">positive regulation of mast cell chemota...</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="even">
<td>818</td>
<td align="left"><a href="GO:0070831" class="uri">GO:0070831</a></td>
<td align="left">basement membrane assembly</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="odd">
<td>819</td>
<td align="left"><a href="GO:0070874" class="uri">GO:0070874</a></td>
<td align="left">negative regulation of glycogen metaboli...</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="even">
<td>820</td>
<td align="left"><a href="GO:0070875" class="uri">GO:0070875</a></td>
<td align="left">positive regulation of glycogen metaboli...</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="odd">
<td>821</td>
<td align="left"><a href="GO:0071459" class="uri">GO:0071459</a></td>
<td align="left">protein localization to chromosome, cent...</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="even">
<td>822</td>
<td align="left"><a href="GO:0071930" class="uri">GO:0071930</a></td>
<td align="left">negative regulation of transcription inv...</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="odd">
<td>823</td>
<td align="left"><a href="GO:0071973" class="uri">GO:0071973</a></td>
<td align="left">bacterial-type flagellum-dependent cell ...</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="even">
<td>824</td>
<td align="left"><a href="GO:0075136" class="uri">GO:0075136</a></td>
<td align="left">response to host</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="odd">
<td>825</td>
<td align="left"><a href="GO:0090306" class="uri">GO:0090306</a></td>
<td align="left">spindle assembly involved in meiosis</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="even">
<td>826</td>
<td align="left"><a href="GO:0097476" class="uri">GO:0097476</a></td>
<td align="left">spinal cord motor neuron migration</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="odd">
<td>827</td>
<td align="left"><a href="GO:0097477" class="uri">GO:0097477</a></td>
<td align="left">lateral motor column neuron migration</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="even">
<td>828</td>
<td align="left"><a href="GO:0097588" class="uri">GO:0097588</a></td>
<td align="left">archaeal or bacterial-type flagellum-dep...</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="odd">
<td>829</td>
<td align="left"><a href="GO:0098535" class="uri">GO:0098535</a></td>
<td align="left">de novo centriole assembly</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="even">
<td>830</td>
<td align="left"><a href="GO:1902076" class="uri">GO:1902076</a></td>
<td align="left">regulation of lateral motor column neuro...</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="odd">
<td>831</td>
<td align="left"><a href="GO:1902078" class="uri">GO:1902078</a></td>
<td align="left">positive regulation of lateral motor col...</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="even">
<td>832</td>
<td align="left"><a href="GO:1902668" class="uri">GO:1902668</a></td>
<td align="left">negative regulation of axon guidance</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="odd">
<td>833</td>
<td align="left"><a href="GO:1902977" class="uri">GO:1902977</a></td>
<td align="left">mitotic DNA replication preinitiation co...</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="even">
<td>834</td>
<td align="left"><a href="GO:1990034" class="uri">GO:1990034</a></td>
<td align="left">calcium ion export from cell</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="odd">
<td>835</td>
<td align="left"><a href="GO:1990035" class="uri">GO:1990035</a></td>
<td align="left">calcium ion import into cell</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="even">
<td>836</td>
<td align="left"><a href="GO:2000192" class="uri">GO:2000192</a></td>
<td align="left">negative regulation of fatty acid transp...</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="odd">
<td>837</td>
<td align="left"><a href="GO:2000251" class="uri">GO:2000251</a></td>
<td align="left">positive regulation of actin cytoskeleto...</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="even">
<td>838</td>
<td align="left"><a href="GO:2000289" class="uri">GO:2000289</a></td>
<td align="left">regulation of photoreceptor cell axon gu...</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="odd">
<td>839</td>
<td align="left"><a href="GO:2000305" class="uri">GO:2000305</a></td>
<td align="left">semaphorin-plexin signaling pathway invo...</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="even">
<td>840</td>
<td align="left"><a href="GO:2000539" class="uri">GO:2000539</a></td>
<td align="left">regulation of protein geranylgeranylatio...</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="odd">
<td>841</td>
<td align="left"><a href="GO:2000541" class="uri">GO:2000541</a></td>
<td align="left">positive regulation of protein geranylge...</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
<tr class="even">
<td>842</td>
<td align="left"><a href="GO:2000821" class="uri">GO:2000821</a></td>
<td align="left">regulation of grooming behavior</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0.45</td>
<td align="left">0.04980</td>
</tr>
</tbody>
</table>

Table 10: GO-Term compartments overrepresented in the set of
underexpressed DEGs.

<table>
<thead>
<tr class="header">
<th></th>
<th align="left">GO.ID</th>
<th align="left">Term</th>
<th align="right">Annotated</th>
<th align="right">Significant</th>
<th align="right">Expected</th>
<th align="left">classic</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>1</td>
<td align="left"><a href="GO:0005929" class="uri">GO:0005929</a></td>
<td align="left">cilium</td>
<td align="right">319</td>
<td align="right">193</td>
<td align="right">72.59</td>
<td align="left">&lt; 1e-30</td>
</tr>
<tr class="even">
<td>2</td>
<td align="left"><a href="GO:0044441" class="uri">GO:0044441</a></td>
<td align="left">ciliary part</td>
<td align="right">195</td>
<td align="right">127</td>
<td align="right">44.37</td>
<td align="left">&lt; 1e-30</td>
</tr>
<tr class="odd">
<td>3</td>
<td align="left"><a href="GO:0005930" class="uri">GO:0005930</a></td>
<td align="left">axoneme</td>
<td align="right">100</td>
<td align="right">81</td>
<td align="right">22.75</td>
<td align="left">&lt; 1e-30</td>
</tr>
<tr class="even">
<td>4</td>
<td align="left"><a href="GO:0097014" class="uri">GO:0097014</a></td>
<td align="left">ciliary cytoplasm</td>
<td align="right">100</td>
<td align="right">81</td>
<td align="right">22.75</td>
<td align="left">&lt; 1e-30</td>
</tr>
<tr class="odd">
<td>44</td>
<td align="left"><a href="GO:0000780" class="uri">GO:0000780</a></td>
<td align="left">condensed nuclear chromosome, centromeri...</td>
<td align="right">6</td>
<td align="right">6</td>
<td align="right">1.37</td>
<td align="left">0.00014</td>
</tr>
<tr class="even">
<td>45</td>
<td align="left"><a href="GO:0005583" class="uri">GO:0005583</a></td>
<td align="left">fibrillar collagen trimer</td>
<td align="right">6</td>
<td align="right">6</td>
<td align="right">1.37</td>
<td align="left">0.00014</td>
</tr>
<tr class="odd">
<td>46</td>
<td align="left"><a href="GO:0097228" class="uri">GO:0097228</a></td>
<td align="left">sperm principal piece</td>
<td align="right">6</td>
<td align="right">6</td>
<td align="right">1.37</td>
<td align="left">0.00014</td>
</tr>
<tr class="even">
<td>47</td>
<td align="left"><a href="GO:0098643" class="uri">GO:0098643</a></td>
<td align="left">banded collagen fibril</td>
<td align="right">6</td>
<td align="right">6</td>
<td align="right">1.37</td>
<td align="left">0.00014</td>
</tr>
<tr class="odd">
<td>48</td>
<td align="left"><a href="GO:0000795" class="uri">GO:0000795</a></td>
<td align="left">synaptonemal complex</td>
<td align="right">14</td>
<td align="right">10</td>
<td align="right">3.19</td>
<td align="left">0.00015</td>
</tr>
<tr class="even">
<td>49</td>
<td align="left"><a href="GO:0030992" class="uri">GO:0030992</a></td>
<td align="left">intraciliary transport particle B</td>
<td align="right">12</td>
<td align="right">9</td>
<td align="right">2.73</td>
<td align="left">0.00018</td>
</tr>
<tr class="odd">
<td>50</td>
<td align="left"><a href="GO:0097542" class="uri">GO:0097542</a></td>
<td align="left">ciliary tip</td>
<td align="right">19</td>
<td align="right">12</td>
<td align="right">4.32</td>
<td align="left">0.00018</td>
</tr>
<tr class="even">
<td>51</td>
<td align="left"><a href="GO:0097232" class="uri">GO:0097232</a></td>
<td align="left">lamellar body membrane</td>
<td align="right">8</td>
<td align="right">7</td>
<td align="right">1.82</td>
<td align="left">0.00020</td>
</tr>
<tr class="odd">
<td>52</td>
<td align="left"><a href="GO:0097233" class="uri">GO:0097233</a></td>
<td align="left">alveolar lamellar body membrane</td>
<td align="right">8</td>
<td align="right">7</td>
<td align="right">1.82</td>
<td align="left">0.00020</td>
</tr>
<tr class="even">
<td>53</td>
<td align="left"><a href="GO:0031512" class="uri">GO:0031512</a></td>
<td align="left">motile primary cilium</td>
<td align="right">10</td>
<td align="right">8</td>
<td align="right">2.28</td>
<td align="left">0.00020</td>
</tr>
<tr class="odd">
<td>54</td>
<td align="left"><a href="GO:0000793" class="uri">GO:0000793</a></td>
<td align="left">condensed chromosome</td>
<td align="right">80</td>
<td align="right">32</td>
<td align="right">18.20</td>
<td align="left">0.00038</td>
</tr>
<tr class="even">
<td>55</td>
<td align="left"><a href="GO:0005615" class="uri">GO:0005615</a></td>
<td align="left">extracellular space</td>
<td align="right">301</td>
<td align="right">93</td>
<td align="right">68.49</td>
<td align="left">0.00052</td>
</tr>
<tr class="odd">
<td>56</td>
<td align="left"><a href="GO:0042599" class="uri">GO:0042599</a></td>
<td align="left">lamellar body</td>
<td align="right">11</td>
<td align="right">8</td>
<td align="right">2.50</td>
<td align="left">0.00060</td>
</tr>
<tr class="even">
<td>57</td>
<td align="left"><a href="GO:1990718" class="uri">GO:1990718</a></td>
<td align="left">axonemal central pair projection</td>
<td align="right">5</td>
<td align="right">5</td>
<td align="right">1.14</td>
<td align="left">0.00061</td>
</tr>
<tr class="odd">
<td>58</td>
<td align="left"><a href="GO:0042383" class="uri">GO:0042383</a></td>
<td align="left">sarcolemma</td>
<td align="right">56</td>
<td align="right">24</td>
<td align="right">12.74</td>
<td align="left">0.00062</td>
</tr>
<tr class="even">
<td>59</td>
<td align="left"><a href="GO:0031224" class="uri">GO:0031224</a></td>
<td align="left">intrinsic component of membrane</td>
<td align="right">1649</td>
<td align="right">423</td>
<td align="right">375.23</td>
<td align="left">0.00070</td>
</tr>
<tr class="odd">
<td>60</td>
<td align="left"><a href="GO:0051233" class="uri">GO:0051233</a></td>
<td align="left">spindle midzone</td>
<td align="right">9</td>
<td align="right">7</td>
<td align="right">2.05</td>
<td align="left">0.00072</td>
</tr>
<tr class="even">
<td>61</td>
<td align="left"><a href="GO:0030991" class="uri">GO:0030991</a></td>
<td align="left">intraciliary transport particle A</td>
<td align="right">7</td>
<td align="right">6</td>
<td align="right">1.59</td>
<td align="left">0.00078</td>
</tr>
<tr class="odd">
<td>62</td>
<td align="left"><a href="GO:0032421" class="uri">GO:0032421</a></td>
<td align="left">stereocilium bundle</td>
<td align="right">33</td>
<td align="right">16</td>
<td align="right">7.51</td>
<td align="left">0.00101</td>
</tr>
<tr class="even">
<td>63</td>
<td align="left"><a href="GO:0030496" class="uri">GO:0030496</a></td>
<td align="left">midbody</td>
<td align="right">56</td>
<td align="right">23</td>
<td align="right">12.74</td>
<td align="left">0.00161</td>
</tr>
<tr class="odd">
<td>64</td>
<td align="left"><a href="GO:0005886" class="uri">GO:0005886</a></td>
<td align="left">plasma membrane</td>
<td align="right">1423</td>
<td align="right">365</td>
<td align="right">323.80</td>
<td align="left">0.00192</td>
</tr>
<tr class="even">
<td>65</td>
<td align="left"><a href="GO:0036157" class="uri">GO:0036157</a></td>
<td align="left">outer dynein arm</td>
<td align="right">10</td>
<td align="right">7</td>
<td align="right">2.28</td>
<td align="left">0.00194</td>
</tr>
<tr class="odd">
<td>66</td>
<td align="left"><a href="GO:0097208" class="uri">GO:0097208</a></td>
<td align="left">alveolar lamellar body</td>
<td align="right">10</td>
<td align="right">7</td>
<td align="right">2.28</td>
<td align="left">0.00194</td>
</tr>
<tr class="even">
<td>67</td>
<td align="left"><a href="GO:0016021" class="uri">GO:0016021</a></td>
<td align="left">integral component of membrane</td>
<td align="right">1610</td>
<td align="right">409</td>
<td align="right">366.35</td>
<td align="left">0.00203</td>
</tr>
<tr class="odd">
<td>68</td>
<td align="left"><a href="GO:0005819" class="uri">GO:0005819</a></td>
<td align="left">spindle</td>
<td align="right">137</td>
<td align="right">46</td>
<td align="right">31.17</td>
<td align="left">0.00226</td>
</tr>
<tr class="even">
<td>69</td>
<td align="left"><a href="GO:0000778" class="uri">GO:0000778</a></td>
<td align="left">condensed nuclear chromosome kinetochore</td>
<td align="right">4</td>
<td align="right">4</td>
<td align="right">0.91</td>
<td align="left">0.00267</td>
</tr>
<tr class="odd">
<td>70</td>
<td align="left"><a href="GO:0030089" class="uri">GO:0030089</a></td>
<td align="left">phycobilisome</td>
<td align="right">4</td>
<td align="right">4</td>
<td align="right">0.91</td>
<td align="left">0.00267</td>
</tr>
<tr class="even">
<td>71</td>
<td align="left"><a href="GO:0043256" class="uri">GO:0043256</a></td>
<td align="left">laminin complex</td>
<td align="right">4</td>
<td align="right">4</td>
<td align="right">0.91</td>
<td align="left">0.00267</td>
</tr>
<tr class="odd">
<td>72</td>
<td align="left"><a href="GO:0045263" class="uri">GO:0045263</a></td>
<td align="left">proton-transporting ATP synthase complex...</td>
<td align="right">4</td>
<td align="right">4</td>
<td align="right">0.91</td>
<td align="left">0.00267</td>
</tr>
<tr class="even">
<td>73</td>
<td align="left"><a href="GO:0072687" class="uri">GO:0072687</a></td>
<td align="left">meiotic spindle</td>
<td align="right">4</td>
<td align="right">4</td>
<td align="right">0.91</td>
<td align="left">0.00267</td>
</tr>
<tr class="odd">
<td>74</td>
<td align="left"><a href="GO:0097149" class="uri">GO:0097149</a></td>
<td align="left">centralspindlin complex</td>
<td align="right">4</td>
<td align="right">4</td>
<td align="right">0.91</td>
<td align="left">0.00267</td>
</tr>
<tr class="even">
<td>75</td>
<td align="left"><a href="GO:0030667" class="uri">GO:0030667</a></td>
<td align="left">secretory granule membrane</td>
<td align="right">68</td>
<td align="right">26</td>
<td align="right">15.47</td>
<td align="left">0.00279</td>
</tr>
<tr class="odd">
<td>76</td>
<td align="left"><a href="GO:0001669" class="uri">GO:0001669</a></td>
<td align="left">acrosomal vesicle</td>
<td align="right">21</td>
<td align="right">11</td>
<td align="right">4.78</td>
<td align="left">0.00290</td>
</tr>
<tr class="even">
<td>77</td>
<td align="left"><a href="GO:0032420" class="uri">GO:0032420</a></td>
<td align="left">stereocilium</td>
<td align="right">27</td>
<td align="right">13</td>
<td align="right">6.14</td>
<td align="left">0.00323</td>
</tr>
<tr class="odd">
<td>78</td>
<td align="left"><a href="GO:0044815" class="uri">GO:0044815</a></td>
<td align="left">DNA packaging complex</td>
<td align="right">16</td>
<td align="right">9</td>
<td align="right">3.64</td>
<td align="left">0.00378</td>
</tr>
<tr class="even">
<td>79</td>
<td align="left"><a href="GO:0005911" class="uri">GO:0005911</a></td>
<td align="left">cell-cell junction</td>
<td align="right">178</td>
<td align="right">56</td>
<td align="right">40.50</td>
<td align="left">0.00420</td>
</tr>
<tr class="odd">
<td>80</td>
<td align="left"><a href="GO:0071944" class="uri">GO:0071944</a></td>
<td align="left">cell periphery</td>
<td align="right">1519</td>
<td align="right">384</td>
<td align="right">345.65</td>
<td align="left">0.00424</td>
</tr>
<tr class="even">
<td>81</td>
<td align="left"><a href="GO:0072686" class="uri">GO:0072686</a></td>
<td align="left">mitotic spindle</td>
<td align="right">22</td>
<td align="right">11</td>
<td align="right">5.01</td>
<td align="left">0.00463</td>
</tr>
<tr class="odd">
<td>82</td>
<td align="left"><a href="GO:0005871" class="uri">GO:0005871</a></td>
<td align="left">kinesin complex</td>
<td align="right">25</td>
<td align="right">12</td>
<td align="right">5.69</td>
<td align="left">0.00477</td>
</tr>
<tr class="even">
<td>83</td>
<td align="left"><a href="GO:0005876" class="uri">GO:0005876</a></td>
<td align="left">spindle microtubule</td>
<td align="right">25</td>
<td align="right">12</td>
<td align="right">5.69</td>
<td align="left">0.00477</td>
</tr>
<tr class="odd">
<td>84</td>
<td align="left"><a href="GO:0031225" class="uri">GO:0031225</a></td>
<td align="left">anchored component of membrane</td>
<td align="right">38</td>
<td align="right">16</td>
<td align="right">8.65</td>
<td align="left">0.00601</td>
</tr>
<tr class="even">
<td>85</td>
<td align="left"><a href="GO:0097539" class="uri">GO:0097539</a></td>
<td align="left">ciliary transition fiber</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1.14</td>
<td align="left">0.01094</td>
</tr>
<tr class="odd">
<td>86</td>
<td align="left"><a href="GO:0045177" class="uri">GO:0045177</a></td>
<td align="left">apical part of cell</td>
<td align="right">133</td>
<td align="right">42</td>
<td align="right">30.26</td>
<td align="left">0.01130</td>
</tr>
<tr class="even">
<td>87</td>
<td align="left"><a href="GO:0001520" class="uri">GO:0001520</a></td>
<td align="left">outer dense fiber</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.68</td>
<td align="left">0.01176</td>
</tr>
<tr class="odd">
<td>88</td>
<td align="left"><a href="GO:0005606" class="uri">GO:0005606</a></td>
<td align="left">laminin-1 complex</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.68</td>
<td align="left">0.01176</td>
</tr>
<tr class="even">
<td>89</td>
<td align="left"><a href="GO:0005608" class="uri">GO:0005608</a></td>
<td align="left">laminin-3 complex</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.68</td>
<td align="left">0.01176</td>
</tr>
<tr class="odd">
<td>90</td>
<td align="left"><a href="GO:0031262" class="uri">GO:0031262</a></td>
<td align="left">Ndc80 complex</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.68</td>
<td align="left">0.01176</td>
</tr>
<tr class="even">
<td>91</td>
<td align="left"><a href="GO:0032133" class="uri">GO:0032133</a></td>
<td align="left">chromosome passenger complex</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.68</td>
<td align="left">0.01176</td>
</tr>
<tr class="odd">
<td>92</td>
<td align="left"><a href="GO:0035003" class="uri">GO:0035003</a></td>
<td align="left">subapical complex</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.68</td>
<td align="left">0.01176</td>
</tr>
<tr class="even">
<td>93</td>
<td align="left"><a href="GO:0043230" class="uri">GO:0043230</a></td>
<td align="left">extracellular organelle</td>
<td align="right">704</td>
<td align="right">184</td>
<td align="right">160.19</td>
<td align="left">0.01406</td>
</tr>
<tr class="odd">
<td>94</td>
<td align="left"><a href="GO:1903561" class="uri">GO:1903561</a></td>
<td align="left">extracellular vesicle</td>
<td align="right">703</td>
<td align="right">183</td>
<td align="right">159.97</td>
<td align="left">0.01678</td>
</tr>
<tr class="even">
<td>95</td>
<td align="left"><a href="GO:0005903" class="uri">GO:0005903</a></td>
<td align="left">brush border</td>
<td align="right">42</td>
<td align="right">16</td>
<td align="right">9.56</td>
<td align="left">0.01783</td>
</tr>
<tr class="odd">
<td>96</td>
<td align="left"><a href="GO:0030076" class="uri">GO:0030076</a></td>
<td align="left">light-harvesting complex</td>
<td align="right">8</td>
<td align="right">5</td>
<td align="right">1.82</td>
<td align="left">0.01820</td>
</tr>
<tr class="even">
<td>97</td>
<td align="left"><a href="GO:0030054" class="uri">GO:0030054</a></td>
<td align="left">cell junction</td>
<td align="right">487</td>
<td align="right">130</td>
<td align="right">110.82</td>
<td align="left">0.01906</td>
</tr>
<tr class="odd">
<td>98</td>
<td align="left"><a href="GO:0065010" class="uri">GO:0065010</a></td>
<td align="left">extracellular membrane-bounded organelle</td>
<td align="right">702</td>
<td align="right">182</td>
<td align="right">159.74</td>
<td align="left">0.01994</td>
</tr>
<tr class="even">
<td>99</td>
<td align="left"><a href="GO:0002139" class="uri">GO:0002139</a></td>
<td align="left">stereocilia coupling link</td>
<td align="right">11</td>
<td align="right">6</td>
<td align="right">2.50</td>
<td align="left">0.02185</td>
</tr>
<tr class="odd">
<td>100</td>
<td align="left"><a href="GO:1990075" class="uri">GO:1990075</a></td>
<td align="left">periciliary membrane compartment</td>
<td align="right">11</td>
<td align="right">6</td>
<td align="right">2.50</td>
<td align="left">0.02185</td>
</tr>
<tr class="even">
<td>101</td>
<td align="left"><a href="GO:0005882" class="uri">GO:0005882</a></td>
<td align="left">intermediate filament</td>
<td align="right">14</td>
<td align="right">7</td>
<td align="right">3.19</td>
<td align="left">0.02329</td>
</tr>
<tr class="odd">
<td>102</td>
<td align="left"><a href="GO:0061618" class="uri">GO:0061618</a></td>
<td align="left">sublamina densa</td>
<td align="right">14</td>
<td align="right">7</td>
<td align="right">3.19</td>
<td align="left">0.02329</td>
</tr>
<tr class="even">
<td>103</td>
<td align="left"><a href="GO:0070062" class="uri">GO:0070062</a></td>
<td align="left">extracellular exosome</td>
<td align="right">701</td>
<td align="right">181</td>
<td align="right">159.51</td>
<td align="left">0.02361</td>
</tr>
<tr class="odd">
<td>104</td>
<td align="left"><a href="GO:0032426" class="uri">GO:0032426</a></td>
<td align="left">stereocilium bundle tip</td>
<td align="right">6</td>
<td align="right">4</td>
<td align="right">1.37</td>
<td align="left">0.02690</td>
</tr>
<tr class="even">
<td>105</td>
<td align="left"><a href="GO:0000779" class="uri">GO:0000779</a></td>
<td align="left">condensed chromosome, centromeric region</td>
<td align="right">37</td>
<td align="right">14</td>
<td align="right">8.42</td>
<td align="left">0.02745</td>
</tr>
<tr class="odd">
<td>106</td>
<td align="left"><a href="GO:0043197" class="uri">GO:0043197</a></td>
<td align="left">dendritic spine</td>
<td align="right">51</td>
<td align="right">18</td>
<td align="right">11.60</td>
<td align="left">0.02811</td>
</tr>
<tr class="even">
<td>107</td>
<td align="left"><a href="GO:0044309" class="uri">GO:0044309</a></td>
<td align="left">neuron spine</td>
<td align="right">51</td>
<td align="right">18</td>
<td align="right">11.60</td>
<td align="left">0.02811</td>
</tr>
<tr class="odd">
<td>108</td>
<td align="left"><a href="GO:0001917" class="uri">GO:0001917</a></td>
<td align="left">photoreceptor inner segment</td>
<td align="right">24</td>
<td align="right">10</td>
<td align="right">5.46</td>
<td align="left">0.02994</td>
</tr>
<tr class="even">
<td>109</td>
<td align="left"><a href="GO:0000786" class="uri">GO:0000786</a></td>
<td align="left">nucleosome</td>
<td align="right">12</td>
<td align="right">6</td>
<td align="right">2.73</td>
<td align="left">0.03544</td>
</tr>
<tr class="odd">
<td>110</td>
<td align="left"><a href="GO:0042589" class="uri">GO:0042589</a></td>
<td align="left">zymogen granule membrane</td>
<td align="right">35</td>
<td align="right">13</td>
<td align="right">7.96</td>
<td align="left">0.03839</td>
</tr>
<tr class="even">
<td>111</td>
<td align="left"><a href="GO:0000796" class="uri">GO:0000796</a></td>
<td align="left">condensin complex</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.91</td>
<td align="left">0.03904</td>
</tr>
<tr class="odd">
<td>112</td>
<td align="left"><a href="GO:0000940" class="uri">GO:0000940</a></td>
<td align="left">condensed chromosome outer kinetochore</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.91</td>
<td align="left">0.03904</td>
</tr>
<tr class="even">
<td>113</td>
<td align="left"><a href="GO:0045095" class="uri">GO:0045095</a></td>
<td align="left">keratin filament</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.91</td>
<td align="left">0.03904</td>
</tr>
<tr class="odd">
<td>114</td>
<td align="left"><a href="GO:0032432" class="uri">GO:0032432</a></td>
<td align="left">actin filament bundle</td>
<td align="right">25</td>
<td align="right">10</td>
<td align="right">5.69</td>
<td align="left">0.04001</td>
</tr>
<tr class="even">
<td>115</td>
<td align="left"><a href="GO:0016459" class="uri">GO:0016459</a></td>
<td align="left">myosin complex</td>
<td align="right">43</td>
<td align="right">15</td>
<td align="right">9.78</td>
<td align="left">0.04732</td>
</tr>
</tbody>
</table>

Pfam domains
------------

Over

<table>
<thead>
<tr class="header">
<th></th>
<th align="left">Domain</th>
<th align="right">Count_in_DEGs</th>
<th align="right">Count_in_Background</th>
<th align="right">p_value</th>
<th align="right">Corrected_p_value</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>3</td>
<td align="left">Peptidase_C1</td>
<td align="right">18</td>
<td align="right">5</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="even">
<td>24</td>
<td align="left">DDE_1</td>
<td align="right">10</td>
<td align="right">7</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="odd">
<td>26</td>
<td align="left">Fz</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="even">
<td>42</td>
<td align="left">UCH</td>
<td align="right">7</td>
<td align="right">1</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="odd">
<td>55</td>
<td align="left">Fibrinogen_C</td>
<td align="right">7</td>
<td align="right">5</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="even">
<td>61</td>
<td align="left">Tubulin</td>
<td align="right">12</td>
<td align="right">6</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="odd">
<td>64</td>
<td align="left">GTP_EFTU_D3</td>
<td align="right">5</td>
<td align="right">2</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="even">
<td>74</td>
<td align="left">VWA</td>
<td align="right">15</td>
<td align="right">4</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="odd">
<td>79</td>
<td align="left">Ig_3</td>
<td align="right">25</td>
<td align="right">24</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="even">
<td>88</td>
<td align="left">Transposase_21</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="odd">
<td>90</td>
<td align="left">GPS</td>
<td align="right">7</td>
<td align="right">5</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="even">
<td>93</td>
<td align="left">S1</td>
<td align="right">7</td>
<td align="right">3</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="odd">
<td>104</td>
<td align="left">Inhibitor_I29</td>
<td align="right">10</td>
<td align="right">1</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="even">
<td>114</td>
<td align="left">cEGF</td>
<td align="right">4</td>
<td align="right">4</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="odd">
<td>121</td>
<td align="left">Y_phosphatase</td>
<td align="right">5</td>
<td align="right">5</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="even">
<td>123</td>
<td align="left">IBR</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="odd">
<td>140</td>
<td align="left">Roc</td>
<td align="right">6</td>
<td align="right">5</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="even">
<td>142</td>
<td align="left">Actin</td>
<td align="right">5</td>
<td align="right">2</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="odd">
<td>145</td>
<td align="left">Kelch_1</td>
<td align="right">6</td>
<td align="right">6</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="even">
<td>157</td>
<td align="left">Arf</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="odd">
<td>163</td>
<td align="left">Ig_2</td>
<td align="right">5</td>
<td align="right">3</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="even">
<td>167</td>
<td align="left">I-set</td>
<td align="right">12</td>
<td align="right">5</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="odd">
<td>179</td>
<td align="left">CARD</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="even">
<td>182</td>
<td align="left">Cu2_monooxygen</td>
<td align="right">4</td>
<td align="right">2</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="odd">
<td>197</td>
<td align="left">SRCR</td>
<td align="right">25</td>
<td align="right">19</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="even">
<td>201</td>
<td align="left">7tm_2</td>
<td align="right">6</td>
<td align="right">6</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="odd">
<td>231</td>
<td align="left">HTH_psq</td>
<td align="right">7</td>
<td align="right">0</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="even">
<td>270</td>
<td align="left">Tubulin_C</td>
<td align="right">10</td>
<td align="right">3</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="odd">
<td>274</td>
<td align="left">Chlam_PMP</td>
<td align="right">6</td>
<td align="right">5</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="even">
<td>278</td>
<td align="left">Calx-beta</td>
<td align="right">28</td>
<td align="right">17</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="odd">
<td>281</td>
<td align="left">FXa_inhibition</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="even">
<td>295</td>
<td align="left">EGF_CA</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="odd">
<td>317</td>
<td align="left">Peptidase_A17</td>
<td align="right">6</td>
<td align="right">5</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="even">
<td>330</td>
<td align="left">GAIN</td>
<td align="right">6</td>
<td align="right">3</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="odd">
<td>227</td>
<td align="left">Pkinase</td>
<td align="right">39</td>
<td align="right">61</td>
<td align="right">0.0000896</td>
<td align="right">0.0001217</td>
</tr>
<tr class="even">
<td>256</td>
<td align="left">NAC</td>
<td align="right">12</td>
<td align="right">15</td>
<td align="right">0.0003831</td>
<td align="right">0.0005181</td>
</tr>
<tr class="odd">
<td>29</td>
<td align="left">rve</td>
<td align="right">20</td>
<td align="right">29</td>
<td align="right">0.0006252</td>
<td align="right">0.0008422</td>
</tr>
<tr class="even">
<td>83</td>
<td align="left">ubiquitin</td>
<td align="right">7</td>
<td align="right">8</td>
<td align="right">0.0008187</td>
<td align="right">0.0010983</td>
</tr>
<tr class="odd">
<td>30</td>
<td align="left">NACHT</td>
<td align="right">11</td>
<td align="right">15</td>
<td align="right">0.0025496</td>
<td align="right">0.0034063</td>
</tr>
<tr class="even">
<td>22</td>
<td align="left">GTP_EFTU_D2</td>
<td align="right">5</td>
<td align="right">6</td>
<td align="right">0.0048549</td>
<td align="right">0.0064601</td>
</tr>
<tr class="odd">
<td>122</td>
<td align="left">RVT_1</td>
<td align="right">21</td>
<td align="right">35</td>
<td align="right">0.0076363</td>
<td align="right">0.0101204</td>
</tr>
<tr class="even">
<td>316</td>
<td align="left">zf-C3HC4</td>
<td align="right">4</td>
<td align="right">5</td>
<td align="right">0.0118132</td>
<td align="right">0.0155934</td>
</tr>
<tr class="odd">
<td>77</td>
<td align="left">COR</td>
<td align="right">9</td>
<td align="right">14</td>
<td align="right">0.0217390</td>
<td align="right">0.0285811</td>
</tr>
<tr class="even">
<td>261</td>
<td align="left">GTP_EFTU</td>
<td align="right">5</td>
<td align="right">7</td>
<td align="right">0.0220190</td>
<td align="right">0.0288344</td>
</tr>
<tr class="odd">
<td>54</td>
<td align="left">SAM_1</td>
<td align="right">3</td>
<td align="right">4</td>
<td align="right">0.0287298</td>
<td align="right">0.0373261</td>
</tr>
<tr class="even">
<td>325</td>
<td align="left">EF-hand_7</td>
<td align="right">3</td>
<td align="right">4</td>
<td align="right">0.0287298</td>
<td align="right">0.0373261</td>
</tr>
</tbody>
</table>

Under

<table>
<thead>
<tr class="header">
<th></th>
<th align="left">Domain</th>
<th align="right">Count_in_DEGs</th>
<th align="right">Count_in_Background</th>
<th align="right">p_value</th>
<th align="right">Corrected_p_value</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>30</td>
<td align="left">CD225</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="even">
<td>56</td>
<td align="left">Cyclin_N</td>
<td align="right">6</td>
<td align="right">6</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="odd">
<td>67</td>
<td align="left">Cu-oxidase_2</td>
<td align="right">5</td>
<td align="right">1</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="even">
<td>69</td>
<td align="left">Gal-bind_lectin</td>
<td align="right">6</td>
<td align="right">0</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="odd">
<td>79</td>
<td align="left">AAA_9</td>
<td align="right">11</td>
<td align="right">0</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="even">
<td>85</td>
<td align="left">TPH</td>
<td align="right">8</td>
<td align="right">0</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="odd">
<td>96</td>
<td align="left">MyTH4</td>
<td align="right">4</td>
<td align="right">2</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="even">
<td>109</td>
<td align="left">Carn_acyltransf</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="odd">
<td>126</td>
<td align="left">Laminin_G_2</td>
<td align="right">7</td>
<td align="right">2</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="even">
<td>163</td>
<td align="left">AAA_6</td>
<td align="right">13</td>
<td align="right">0</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="odd">
<td>186</td>
<td align="left">Arm</td>
<td align="right">5</td>
<td align="right">5</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="even">
<td>196</td>
<td align="left">DHC_N1</td>
<td align="right">10</td>
<td align="right">2</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="odd">
<td>202</td>
<td align="left">ABC2_membrane</td>
<td align="right">15</td>
<td align="right">10</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="even">
<td>221</td>
<td align="left">FERM_M</td>
<td align="right">8</td>
<td align="right">5</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="odd">
<td>248</td>
<td align="left">Cyclin_C</td>
<td align="right">6</td>
<td align="right">3</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="even">
<td>252</td>
<td align="left">HlyIII</td>
<td align="right">5</td>
<td align="right">1</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="odd">
<td>274</td>
<td align="left">Sortilin_C</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="even">
<td>289</td>
<td align="left">EBP</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="odd">
<td>293</td>
<td align="left">Integrin_beta</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="even">
<td>305</td>
<td align="left">ITI_HC_C</td>
<td align="right">4</td>
<td align="right">2</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="odd">
<td>309</td>
<td align="left">polyprenyl_synt</td>
<td align="right">5</td>
<td align="right">3</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="even">
<td>323</td>
<td align="left">Tektin</td>
<td align="right">6</td>
<td align="right">0</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="odd">
<td>330</td>
<td align="left">ADK</td>
<td align="right">5</td>
<td align="right">3</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="even">
<td>344</td>
<td align="left">G-alpha</td>
<td align="right">4</td>
<td align="right">4</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="odd">
<td>347</td>
<td align="left">Histone</td>
<td align="right">5</td>
<td align="right">3</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="even">
<td>354</td>
<td align="left">MORN</td>
<td align="right">6</td>
<td align="right">5</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="odd">
<td>357</td>
<td align="left">Kazal_2</td>
<td align="right">8</td>
<td align="right">3</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="even">
<td>374</td>
<td align="left">7tm_2</td>
<td align="right">20</td>
<td align="right">18</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="odd">
<td>387</td>
<td align="left">AAA_7</td>
<td align="right">13</td>
<td align="right">1</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="even">
<td>395</td>
<td align="left">p450</td>
<td align="right">18</td>
<td align="right">9</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="odd">
<td>404</td>
<td align="left">SHIPPO-rpt</td>
<td align="right">6</td>
<td align="right">0</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="even">
<td>405</td>
<td align="left">Dynein_heavy</td>
<td align="right">14</td>
<td align="right">0</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="odd">
<td>421</td>
<td align="left">Cu-oxidase</td>
<td align="right">5</td>
<td align="right">1</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="even">
<td>432</td>
<td align="left">Dpy-30</td>
<td align="right">7</td>
<td align="right">0</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="odd">
<td>446</td>
<td align="left">AAA_8</td>
<td align="right">12</td>
<td align="right">2</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="even">
<td>456</td>
<td align="left">FERM_C</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="odd">
<td>491</td>
<td align="left">IQ</td>
<td align="right">14</td>
<td align="right">6</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="even">
<td>503</td>
<td align="left">Abhydrolase_1</td>
<td align="right">4</td>
<td align="right">2</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="odd">
<td>527</td>
<td align="left">CHGN</td>
<td align="right">6</td>
<td align="right">3</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="even">
<td>554</td>
<td align="left">FERM_N</td>
<td align="right">4</td>
<td align="right">2</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="odd">
<td>563</td>
<td align="left">Kinesin</td>
<td align="right">13</td>
<td align="right">9</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="even">
<td>589</td>
<td align="left">ELO</td>
<td align="right">6</td>
<td align="right">0</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="odd">
<td>595</td>
<td align="left">RYDR_ITPR</td>
<td align="right">5</td>
<td align="right">5</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="even">
<td>604</td>
<td align="left">Forkhead</td>
<td align="right">7</td>
<td align="right">4</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="odd">
<td>605</td>
<td align="left">zf-MYND</td>
<td align="right">6</td>
<td align="right">5</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="even">
<td>664</td>
<td align="left">DHC_N2</td>
<td align="right">15</td>
<td align="right">2</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="odd">
<td>683</td>
<td align="left">Ins145_P3_rec</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="even">
<td>696</td>
<td align="left">Laminin_G_1</td>
<td align="right">8</td>
<td align="right">0</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="odd">
<td>700</td>
<td align="left">Y_phosphatase</td>
<td align="right">17</td>
<td align="right">15</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="even">
<td>716</td>
<td align="left">Collagen</td>
<td align="right">20</td>
<td align="right">19</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="odd">
<td>718</td>
<td align="left">TTL</td>
<td align="right">7</td>
<td align="right">3</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="even">
<td>735</td>
<td align="left">EGF_CA</td>
<td align="right">10</td>
<td align="right">9</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="odd">
<td>748</td>
<td align="left">Laminin_EGF</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="even">
<td>755</td>
<td align="left">Peptidase_M14</td>
<td align="right">4</td>
<td align="right">2</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="odd">
<td>759</td>
<td align="left">FA_desaturase</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="even">
<td>402</td>
<td align="left">Calx-beta</td>
<td align="right">50</td>
<td align="right">53</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="odd">
<td>376</td>
<td align="left">EGF</td>
<td align="right">34</td>
<td align="right">40</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="even">
<td>254</td>
<td align="left">ig</td>
<td align="right">22</td>
<td align="right">25</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="odd">
<td>122</td>
<td align="left">Ig_3</td>
<td align="right">48</td>
<td align="right">82</td>
<td align="right">0.0000000</td>
<td align="right">0.0000000</td>
</tr>
<tr class="even">
<td>443</td>
<td align="left">AAA</td>
<td align="right">33</td>
<td align="right">57</td>
<td align="right">0.0000017</td>
<td align="right">0.0000028</td>
</tr>
<tr class="odd">
<td>435</td>
<td align="left">MT</td>
<td align="right">14</td>
<td align="right">18</td>
<td align="right">0.0000031</td>
<td align="right">0.0000051</td>
</tr>
<tr class="even">
<td>334</td>
<td align="left">TPR_8</td>
<td align="right">11</td>
<td align="right">16</td>
<td align="right">0.0002040</td>
<td align="right">0.0003396</td>
</tr>
<tr class="odd">
<td>506</td>
<td align="left">Sushi</td>
<td align="right">18</td>
<td align="right">32</td>
<td align="right">0.0003727</td>
<td align="right">0.0006190</td>
</tr>
<tr class="even">
<td>680</td>
<td align="left">F-box-like</td>
<td align="right">10</td>
<td align="right">15</td>
<td align="right">0.0005295</td>
<td align="right">0.0008776</td>
</tr>
<tr class="odd">
<td>294</td>
<td align="left">Cadherin</td>
<td align="right">18</td>
<td align="right">33</td>
<td align="right">0.0006387</td>
<td align="right">0.0010563</td>
</tr>
<tr class="even">
<td>635</td>
<td align="left">fn3</td>
<td align="right">24</td>
<td align="right">48</td>
<td align="right">0.0007438</td>
<td align="right">0.0012276</td>
</tr>
<tr class="odd">
<td>703</td>
<td align="left">PID</td>
<td align="right">6</td>
<td align="right">8</td>
<td align="right">0.0010951</td>
<td align="right">0.0018035</td>
</tr>
<tr class="even">
<td>361</td>
<td align="left">AAA_5</td>
<td align="right">4</td>
<td align="right">5</td>
<td align="right">0.0021503</td>
<td align="right">0.0035191</td>
</tr>
<tr class="odd">
<td>393</td>
<td align="left">ABC2_membrane_3</td>
<td align="right">4</td>
<td align="right">5</td>
<td align="right">0.0021503</td>
<td align="right">0.0035191</td>
</tr>
<tr class="even">
<td>441</td>
<td align="left">cNMP_binding</td>
<td align="right">4</td>
<td align="right">5</td>
<td align="right">0.0021503</td>
<td align="right">0.0035191</td>
</tr>
<tr class="odd">
<td>148</td>
<td align="left">GPS</td>
<td align="right">10</td>
<td align="right">17</td>
<td align="right">0.0025973</td>
<td align="right">0.0042417</td>
</tr>
<tr class="even">
<td>103</td>
<td align="left">AIG1</td>
<td align="right">5</td>
<td align="right">7</td>
<td align="right">0.0032993</td>
<td align="right">0.0053769</td>
</tr>
<tr class="odd">
<td>290</td>
<td align="left">F-box</td>
<td align="right">12</td>
<td align="right">22</td>
<td align="right">0.0033969</td>
<td align="right">0.0055244</td>
</tr>
<tr class="even">
<td>679</td>
<td align="left">hEGF</td>
<td align="right">7</td>
<td align="right">11</td>
<td align="right">0.0036150</td>
<td align="right">0.0058546</td>
</tr>
<tr class="odd">
<td>730</td>
<td align="left">DSPc</td>
<td align="right">7</td>
<td align="right">11</td>
<td align="right">0.0036150</td>
<td align="right">0.0058546</td>
</tr>
<tr class="even">
<td>10</td>
<td align="left">SRCR</td>
<td align="right">23</td>
<td align="right">51</td>
<td align="right">0.0053630</td>
<td align="right">0.0086677</td>
</tr>
<tr class="odd">
<td>286</td>
<td align="left">Fibrinogen_C</td>
<td align="right">12</td>
<td align="right">23</td>
<td align="right">0.0057367</td>
<td align="right">0.0092524</td>
</tr>
<tr class="even">
<td>609</td>
<td align="left">Cyt-b5</td>
<td align="right">3</td>
<td align="right">4</td>
<td align="right">0.0073492</td>
<td align="right">0.0118286</td>
</tr>
<tr class="odd">
<td>500</td>
<td align="left">Myosin_head</td>
<td align="right">7</td>
<td align="right">12</td>
<td align="right">0.0080657</td>
<td align="right">0.0129550</td>
</tr>
<tr class="even">
<td>776</td>
<td align="left">DOMON</td>
<td align="right">4</td>
<td align="right">6</td>
<td align="right">0.0097570</td>
<td align="right">0.0156393</td>
</tr>
<tr class="odd">
<td>350</td>
<td align="left">LRR_9</td>
<td align="right">5</td>
<td align="right">8</td>
<td align="right">0.0099121</td>
<td align="right">0.0158553</td>
</tr>
<tr class="even">
<td>682</td>
<td align="left">Ig_2</td>
<td align="right">9</td>
<td align="right">17</td>
<td align="right">0.0105448</td>
<td align="right">0.0168328</td>
</tr>
<tr class="odd">
<td>342</td>
<td align="left">Ank_4</td>
<td align="right">8</td>
<td align="right">15</td>
<td align="right">0.0128642</td>
<td align="right">0.0204515</td>
</tr>
<tr class="even">
<td>545</td>
<td align="left">SAM_1</td>
<td align="right">8</td>
<td align="right">15</td>
<td align="right">0.0128642</td>
<td align="right">0.0204515</td>
</tr>
<tr class="odd">
<td>452</td>
<td align="left">I-set</td>
<td align="right">20</td>
<td align="right">47</td>
<td align="right">0.0177136</td>
<td align="right">0.0281037</td>
</tr>
<tr class="even">
<td>409</td>
<td align="left">HRM</td>
<td align="right">5</td>
<td align="right">9</td>
<td align="right">0.0223882</td>
<td align="right">0.0354481</td>
</tr>
<tr class="odd">
<td>168</td>
<td align="left">TPR_12</td>
<td align="right">4</td>
<td align="right">7</td>
<td align="right">0.0259011</td>
<td align="right">0.0395627</td>
</tr>
<tr class="even">
<td>371</td>
<td align="left">SAM_2</td>
<td align="right">4</td>
<td align="right">7</td>
<td align="right">0.0259011</td>
<td align="right">0.0395627</td>
</tr>
<tr class="odd">
<td>383</td>
<td align="left">DEP</td>
<td align="right">4</td>
<td align="right">7</td>
<td align="right">0.0259011</td>
<td align="right">0.0395627</td>
</tr>
<tr class="even">
<td>510</td>
<td align="left">HMG_box</td>
<td align="right">3</td>
<td align="right">5</td>
<td align="right">0.0281449</td>
<td align="right">0.0428221</td>
</tr>
<tr class="odd">
<td>535</td>
<td align="left">SET</td>
<td align="right">3</td>
<td align="right">5</td>
<td align="right">0.0281449</td>
<td align="right">0.0428221</td>
</tr>
<tr class="even">
<td>470</td>
<td align="left">ABC_tran</td>
<td align="right">10</td>
<td align="right">22</td>
<td align="right">0.0325409</td>
<td align="right">0.0494139</td>
</tr>
</tbody>
</table>

Discussion
==========

Acknowledgements
================

### Conflict of interest

None declared.

### Author contributions

Sergio Vargas designed the study, conducted the experiments and analyzed
the data. Gert Wörheide contributed reagents. Sergio Vargas and Gert
Wörheide wrote the manuscript. Sergio Vargas manages the [project
repository](https://github.com/sevragorgia/CBAS/).

### Disclaimer

None of the products or companies listed are endorsed or recommended in
anyway by the author(s) of this text.

### Scripts and Data availability

All scripts used to analyze the data, as well as some miscellaneous
scripts used to, for instance, prepare the annotation table or generate
GO-term/Pfam input files for the enrichment analyses can be found in the
[project repository](https://github.com/sevragorgia/CBAS/).

References and Notes
====================

[1] I have (repeatedly) failed to prepare libraries for Next Generation
Sequencing from bleached tissue, probably due to the presence of
secondary metabolites in the extracts.

[2] The 10 mm hole puncher has been replaced by a 6 mm version. [Click
here to see the
product.](http://www.sigmaaldrich.com/catalog/product/sigma/z708909?lang=de&region=DE)

[3] We, in fact, keep a permanent culture of "bleached" explants that
are allowed to bleach for at least 12 weeks, are fixed in liquid
nitrogen and stored at -80 °C.

[4] Zymo Research. [Click here to see the product
website.](https://www.zymoresearch.de/rna/dna-rna-co-purification/cells-tissue-rna/zr-duet-dna-rna-miniprep)

[5] EMBL GeneCore.
<http://genecore3.genecore.embl.de/genecore3/index.cfm>

[6] FastQC. <http://www.bioinformatics.babraham.ac.uk/projects/fastqc/>

[7] BioLite: <https://bitbucket.org/caseywdunn/biolite/overview>

[8] Downloaded early 2015. [Uniprot
download.](http://www.uniprot.org/downloads)

[9] Aqu2 isoforms (proteins) [*Amphimedon queenslandica* transcriptome
resource.](http://amphimedon.qcloud.qcif.edu.au/downloads.html)

[10] Blast was run from [Galaxy](https://galaxyproject.org/) and the XML
file was converted to a 25 column table by the
[ncbi\_blast\_plus](http://toolshed.g2.bx.psu.edu/view/devteam/ncbi_blast_plus)
tool.

[11] Gene Ontology Consortium. <http://geneontology.org/>

[12] The UNIPROT accession code of the blast best match associated to
each transcript in the CBAS transcriptome was used to query the GO
annotations associated with the UNIPROT protein. These GO annotations
were used to annotated the CBAS transcripts; note that these are
annotations based on annotations and should be taken with precaution.

[13] TopGO
[(https://bioconductor.org/packages/release/bioc/html/topGO.html)](https://bioconductor.org/packages/release/bioc/html/topGO.html)
expects GO terms in tab separated {Trascript, GoAnnotation} pairs. The
script produces output with a leading column indicating the UNIPROT
accession code used for each transcript, which comes handy when
controling the correct annotation was used for each transcript and can
be easily removed for further analyses.

[14] Transdecoder's Github Repository can be found
[here.](https://transdecoder.github.io/)

[15] I honestly do not know what the purpose of this blast run was, I
did it for sake of completeness and because the reads were not filtered
before assembling the conting. Yet, after seeing the results, it seems
as if not too many contigs matched the bacterial database...

[16] Parra et al. 2007. CEGMA: a pipeline to accurately annotate core
genes in eukaryotic genomes. [Bioinformatics 23 (9):
1061-1067.](http://dx.doi.org/10.1093/bioinformatics/btm071)

[17] Francis et al. 2013. A comparison across non-model animals suggests
an optimal sequencing depth for *de novo* transcriptome assembly. [BMC
Genomics 14:167](http://dx.doi.org/10.1186/1471-2164-14-167)

[18] the script is also provided as part of this repository. The
Copyright (c) of this script belongs trinityrnaseq (2014), who reserves
all rights upon the code. [See the full LICENSE
here.](https://github.com/trinityrnaseq/trinityrnaseq/blob/master/LICENSE)

[19] With these graphs were produced for transcripts that could be
successfully translated using **Transdecoder**. In the case of the
**Uniprot** and **Aqu2** annotations, the transcripts were directly
blasted using blastn. Thus could be possible to obtain **Uniprot** and
**Aqu2** annotations for transcripts that could not be translated with
**Transdecoder**. This is in fact the case, however the number of
transcripts that could not *transdecoded* but were annotated is
relatively small. In the case of **Uniprot**, the number of transcripts
that yielded annotations but were not *transdecoded* was **7535**. For
**Aqu2**, only **11915** not *transdecoded* transcripts could be
annotated.

[20] Here it is worth noting that the evalue cutoff used for blast
against the Uniprot and Aqu2 databases was 0.001. So the annotations
obtained had an evalue at least two orders of magnitud smaller than the
set cutoff evalue.
