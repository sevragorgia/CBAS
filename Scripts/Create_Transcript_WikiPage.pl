#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;

#ADD GNU license here!


########################
#
#
#perl Create_Transcript_WikiPage.pl --metaable ~/Desktop/CBAS/CBAS_Annotation_Metatable_20160407.csv 
#   
    #the metatable fields are as follow:
    #
    # 0.Transcript_name -> this names are including isoforms
    # 1.Length	
    # 2.Sequence
    # 3.Gene name
    # 4.Uniprot_match	
    # 5.evalue
    # 6.Uniprot_match_annotation	
    # 7.Aqu2_match	
    # 8.evalue	
    # 9.Aqu2_match_annotation	
    # 10.ORF_Type	
    # 11.Protein
    # 12.GO_Component	
    # 13.GO_Function	
    # 14.GO_Process	
    # 15.Pfam_Domains
    # 16.KEGG
    # 17.Bacteria_match	
    # 18.evalue	
    # 19.Bacteria_match_annotation	
#
###################################################################################################3



#command line options, the variable names should be self explanatory;
my $metatable;
my $is_deg_table;
my $verbose;
my $debug;

my $help = "The following options have to be provided:
    --metatable = the metatable containing all annotations for a transcriptome
    --is_deg_table = boolean indicating whether DeSeq2 results are included in the table.
A normal command line would look like:
	perl Create_Transcript_WikiPage.pl --metatable cbas_transcriptome_metatable.tsv --is_deg_table

This will generate a wiki page for each transcript in the metatable.

";

#get options from command line
GetOptions("metatable=s" => \$metatable,
	   "is_deg_table" => \$is_deg_table,
	   "verbose" => \$verbose,
       	   "debug" => \$debug,
	   "help" => \$help,
 ) or die("Error in the command line arguments\n$help");


#table reader, expects file path, a hash ref, the fields to be extracted and the field separator
sub main{

  die "\n$help\n" unless(defined $metatable);
  
  #open filehandles
  open my $metatable_handle, '<', $metatable or die "Cannot open the table file\n$!";

  chomp(my $table_header = <$metatable_handle>);

  while(<$metatable_handle>){
    chomp;
    my @table_fields = split("\t");
    #
    #the metatable fields are as follow:
    #
    # 0.Transcript_name -> this names are including isoforms
    # 1.Length	
    # 2.Sequence
    # 3.Gene name
    # 4.Uniprot_match	
    # 5.evalue
    # 6.Uniprot_match_annotation	
    # 7.Aqu2_match	
    # 8.evalue	
    # 9.Aqu2_match_annotation	
    # 10.ORF_Type	
    # 11.Protein
    # 12.GO_Component	
    # 13.GO_Function	
    # 14.GO_Process	
    # 15.Pfam_Domains
    # 16.KEGG
    # 17.Bacteria_match	
    # 18.evalue	
    # 19.Bacteria_match_annotation	
    # 20.Remarks from Uniprot (not in metatable)
    # 21.uniprot annotation notes
    # 22.baseMean
    # 23.log2FoldChange
    # 24.lfcSE
    # 25.stat
    # 26.pvalue
    # 27.padj

    #
    #build the md page
    #
    
    open my $transcript_wikipage, '>', "$table_fields[0].md" or die "Cannot open outfile\n$1";

    my $print_line = "# Transcript name: $table_fields[0]
## Gene name: $table_fields[3]
## Sequence ($table_fields[1] bp):
$table_fields[2]
## Uniprot best match:
$table_fields[4]	$table_fields[6]	evalue = $table_fields[5]
## AQU2 best match:
$table_fields[7]	$table_fields[9]		evalue = $table_fields[8]
## ORF Type = $table_fields[10]
## Protein sequence:
$table_fields[11]
## GO Terms:
### Component
$table_fields[12]
### Function
$table_fields[13]
### Process
$table_fields[14]
## PFam domains
$table_fields[15]
## KEGG
$table_fields[16]
## Bacterial match:
$table_fields[17]	$table_fields[19]	evalue = $table_fields[18]\n";

    if($is_deg_table){
    
    	$print_line .= "## DeSeq2 Results:
### Basemean = $table_fields[22]
### Log2 Fold Change = $table_fields[23] (lfcSE = $table_fields[24])
### Stat = $table_fields[25]
### p-value (BJ adjusted) = $table_fields[27] (unadjusted = $table_fields[26])";
    
    }

    print $transcript_wikipage $print_line;
    close $transcript_wikipage;

  }
  
  close $metatable_handle;
  }




#run the main
exit &main;

