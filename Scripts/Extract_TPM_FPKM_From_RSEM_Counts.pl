#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use File::Basename;

#ADD GNU license here!

################################################################################################################################################################
#blast table Output format

#Because Galaxy focuses on processing tabular data, the default output of this tool is tabular. The standard BLAST+ tabular output contains 12 columns:
#Column 	NCBI name 	Description
#1 	qseqid 	Query Seq-id (ID of your sequence)
#2 	sseqid 	Subject Seq-id (ID of the database hit)
#3 	pident 	Percentage of identical matches
#4 	length 	Alignment length
#5 	mismatch 	Number of mismatches
#6 	gapopen 	Number of gap openings
#7 	qstart 	Start of alignment in query
#8 	qend 	End of alignment in query
#9 	sstart 	Start of alignment in subject (database hit)
#10 	send 	End of alignment in subject (database hit)
#11 	evalue 	Expectation value (E-value)
#12 	bitscore 	Bit score

#The BLAST+ tools can optionally output additional columns of information, but this takes longer to calculate. Many commonly used extra columns are included by selecting the extended tabular output. The extra columns are included after the standard 12 columns. This is so that you can write workflow filtering steps that accept either the 12 or 25 column tabular BLAST output. Galaxy now uses this extended 25 column output by default.
#Column 	NCBI name 	Description
#13 	sallseqid 	All subject Seq-id(s), separated by a ';'
#14 	score 	Raw score
#15 	nident 	Number of identical matches
#16 	positive 	Number of positive-scoring matches
#17 	gaps 	Total number of gaps
#18 	ppos 	Percentage of positive-scoring matches
#19 	qframe 	Query frame
#20 	sframe 	Subject frame
#21 	qseq 	Aligned part of query sequence
#22 	sseq 	Aligned part of subject sequence
#23 	qlen 	Query sequence length
#24 	slen 	Subject sequence length
#25 	salltitles 	All subject title(s), separated by a '<>'
################################################################################################################################################################

########################
#
#
#perl Create_Transcriptome_Annotation_Table.pl --transcripts ~/Desktop/CBAS/Galaxy85-CBAS_Trinity_Reference_Transcriptome.fasta --uniprot ~/Desktop/CBAS/Annotations/Galaxy136-CBAS_Trinity_Reference_Transcriptome_Uniprot_Annotations_SingleMatch.csv --other ~/Desktop/CBAS/Annotations/Galaxy135-CBAS_Trinity_Reference_Transcriptome_Aqu2_Annotations_SingleMatch.csv --transpeps ~/Desktop/CBAS/Transdecoder/Galaxy85-CBAS.fasta.transdecoder.pep --component ~/Desktop/CBAS/GOs/component_CBAS_Trinity_Reference_Transcriptome_Uniprot_Annotations_SingleMatch_GOs.csv --process ~/Desktop/CBAS/GOs/process_CBAS_Trinity_Reference_Transcriptome_Uniprot_Annotations_SingleMatch_GOs.csv --function ~/Desktop/CBAS/GOs/function_CBAS_Trinity_Reference_Transcriptome_Uniprot_Annotations_SingleMatch_GOs.csv --bacteria ~/Desktop/CBAS/Annotations/Galaxy140-CBAS_Trinity_Reference_Transcriptome_AllBacts_Annotations.csv --degs ~/Desktop/CBAS/DeSeq2/CBAS_DeSeq_NoZeros_DESeq_Results_005_SigOnly.csv --pfam ~/Desktop/CBAS/Pfam/CBAS_pfam_scan.csv --kegg ~/Desktop/CBAS/KEGG/20160503_KEGG_Annotations.ko >CBAS_Annotation_Metatable_20160504.csv
#
#
###################################################################################################3



#command line options, the variable names should be self explanatory;
my $counts;
my $verbose;
my $debug;
my $tpm;
my $fpkm;

my $help = "The following options have to be provided:
		--counts = base filename of the counts file
		--tpm = produce a tpm matrix
		--fpkm = produce a fpkm matrix
		--verbose
 		--debug

A normal command line would look like:
	perl Extract_TPM_FPKM_From_RSEM_Counts.pl --counts CountFileBasename

Note:

This script expects a RSEM header containing the following columns:

gene_id	transcript_id(s)	length	effective_length	expected_count	TPM	FPKM

and will produce a TPM and FPKM matrix.
";

#get options from command line
GetOptions("counts=s" => \$counts,
	  "tpm" => \$tpm,
	  "fpkm" => \$fpkm,
	  "verbose" => \$verbose,
	  "debug" => \$debug,
	  "help" => \$help
) or die("Error in the command line arguments\n$help");


#table reader, expects file path, a hash ref, the fields to be extracted and the field separator
sub counts_table_reader{
  my $table_file = shift;
  my $hash_ref = shift;
  my $field = shift;
  my $sep = shift;
  
	
  #open filehandles
  open my $table, '<', $table_file or die "Cannot open the table file\n$!";

  while(<$table>){
    chomp;
    my @table_fields = split($sep);

#    print @table_fields, "\n" if($debug);
#    print $#table_fields, "\n" if($debug);
    
    my $value = $table_fields[$field];
    
    if(exists ${$hash_ref}{$table_fields[0]}){
    
      push @{${$hash_ref}{$table_fields[0]}}, $value;
    
    }else{
    
      ${$hash_ref}{$table_fields[0]}=[$value];
    
    }
    
  }
  close $table;
}

sub print_matrix{
  my $array_ref = shift;#with file names = column names
  my $hash_ref = shift;#the actual data by transcript
  my $ext = shift;#extension of the file.

  my $output_file_name = $counts . "." . $ext;
  
  #open filehandles
  open my $output_file, '>', $output_file_name or die "Cannot open file for output\n$!";
  
  #print header of outfile
  print $output_file "Transcript_ID\t";
  
  foreach my $column_name (@{$array_ref}){
    print $output_file basename($column_name), "\t";
  }
  
  print $output_file "\n";  
  
  foreach my $id (sort keys %{$hash_ref}){
    
    #print transcript id
    print $output_file $id, "\t";
    #then print the associated data
    print $output_file join("\t", @{${$hash_ref}{$id}});
    print $output_file "\n";
  }
  
}

#main loop
sub main {

	#check that some required options are defined;
	die "\n$help\n" unless defined $counts;

#	&print_call if $verbose;

	print $counts if($debug);
	
	#glob all files matching the base filename
	my @counts_files = <$counts*>;
	
	print @counts_files if($debug);
	
	my %tpm_table;
	my %fpkm_table;
		
	foreach my $counts_file (@counts_files){
	
	    counts_table_reader($counts_file, \%tpm_table, 5, "\t") if($tpm);

	    counts_table_reader($counts_file, \%fpkm_table, 6, "\t") if($fpkm);

	}
	
	print_matrix(\@counts_files, \%tpm_table, "tpm") if($tpm);
	print_matrix(\@counts_files, \%fpkm_table, "fpkm") if($fpkm);

	return 0;
}

exit &main;

