# Description of the scripts provided

**Note** that as far as possible I have tried to provide a detailed help/usage section that will be printed if you call the scripts without options (i.e. perl name_of_script.pl). Below I try to summarize what each script does.

### Create_Transcriptome_Annotation_Table.pl

Perl script to generate an annotation metatable from a number of separate annotation files.

### Get_Pfam_Domain_List.pl

Perl script to get the list of Pfam domains contained in a file of Pfam domain by Transcript (e.g. ../Results/Pfam/all_degs_p001_2LFC.pfam). This script just returns a list of the domains present in the input file (e.g. ../Results/Pfam/all_degs_p001_2LFC.pfam.list)

### Get_DEG_Annotations.pl

Perl script that scans the metatable and, given the user provides a list of genes of interest, returns the annotations for these genes. The user can modify the output to only get the desired annotations. The default behaviour is to return all annotations for a transcript.

### Get_DEG_GO_Annotations.pl

Perl script with which the GO term annotations for the DEGs and the transcripts in the background set can be extracted from a set of GO annotation files (e.g. as is provided in Annotations/GOs)

### Get_GO_Annotaitons.pl

Perl script to scan swissprot and retrieve the GO annotations associated with a UNIPROT accession. Using this script one can use the results of a blast annotation table against UNIPROT to get the GO terms associated with the UNIPROT match of a transcript and annotate the transcript this way. **NOTE**: you are annotating against annotations! So, be careful interpreting the results.

## Other miscellaneus scripts provided

### Extract_TPM_FPKM_From_RSEM_Counts.pl

Perl script to produce a TPM/FPKM matrix from RSEM output files.

### Create_Transcript_WikiPage.pl

Perl script to generate github wikipages for a set of transcripts based on the annotations provided as a metatable. The metatable can contain DeSeq2 results appended after the last metatable columns. The pages produced can be pushed to the github wiki and will be functional (e.g. see this repo's wiki).


## Support:

**Please use the Issues tab to report a problem or ask for help.**

