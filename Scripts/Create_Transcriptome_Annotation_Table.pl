#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;

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
my $transcripts;
my $uniprot;
my @uniprot_fields;
my $other_annotations;
my @other_fields;
my $bacteria;
my @bacteria_fields;
my $component;
my $function;
my $process;
my $transdecoder_mrna;
my $transdecoder_peptides;
my $degs;
my $verbose = "";
my $debug = "";
my $output_file = "./out.csv"; #where should I store the genomes?
my $pfam_file;
my $kegg_file;


my $help = "The following options have to be provided:
    --transcripts = the fasta file with the transcripts
		--uniprot = uniprot blast results in extended table format (i.e. 25 columns)
		--uniprot_fields = fields to kept in the annotation table, defaults to field 2, 11 and 25 (Subject Seq-id, e-value and All subject title(s))
		--other_blast = other blast results in extended table format (i.e. 25 columns), e.g. AQU2
		--other_fields = fields to kept in the annotation table, defaults to field 2, 11 and 25 (Subject Seq-id, e-value and All subject title(s))
		--bacteria = blast against bacterial database results in extended table format (i.e. 25 columns)
		--bacteria_fields = fields to kept in the annotation table, defaults to field 2, 11 and 25 (Subject Seq-id, e-value and All subject title(s))
		--out = ./out.csv
		--component = component GO annotations for the uniprot matches
		--process = component GO annotations for the uniprot matches
		--function = component GO annotations for the uniprot matches
		--transpeps = the transdecoder peptides
		--transmrna = the transdecoder mrna
		--pfam = pfam annotations
		--kegg = kegg annotations
		--degs = list of differentially expressed genes
 		--verbose F
 		--debug F

A normal command line would look like:
	perl Create_Transcriptome_Annotation_Table.pl --transcripts cbas_transcriptome.fasta --uniprot ./uniprot.csv --other ./AQU2.csv --go ./go.csv --transpeps cbas_transdec.pep --transmRNA cbas_transdec.mrna --out cbas_annotation.table

Note:

this script expects trinity like headers in the fasta files:
	>TR1|c0_g1_i1 len=396 path=[374:0-395] [-1, 374, -2]

it extracts the transcripts and uses their names (i.e. TR1|c0_g1.i1) to look for annotations. It produces a csv table joining all the annotations for the transcripts.
";

#get options from command line
GetOptions("transcripts=s" => \$transcripts,
					 "out=s" => \$output_file,
					 "uniprot=s" => \$uniprot,
					 "uniprot_fields=s{,}" => \@uniprot_fields,
					 "other=s" => \$other_annotations,
					 "other_fields=s{,}" => \@other_fields,
 					 "bacteria=s" => \$bacteria,
					 "bacteria_fields=s{,}" => \@bacteria_fields,
					 "component=s" => \$component,
					 "function=s" => \$function,
					 "process=s" => \$process,
					 "transpeps=s" => \$transdecoder_peptides,
					 "transmrna=s" => \$transdecoder_mrna,
					 "pfam=s" => \$pfam_file,
					 "kegg=s" => \$kegg_file,
					 "degs=s" => \$degs,
					 "verbose=s" => \$verbose,
					 "debug=s" => \$debug,
					 "help" => \$help
					 ) or die("Error in the command line arguments\n$help");


#fasta reader, expects file path and hash ref
sub transcript_reader{
  my $fasta_file = shift;
  my $hash_ref = shift;
  
  #initialize some counters
	my $transcript_counter = 1;
	
  #open filehandles
  open my $transcript_file, '<', $fasta_file or die "Cannot open the transcript file\n$!";
  
  my $transcript_name="";
  my $length="";
  my $sequence="";

  #read the transcript file an populate hash
  while(<$transcript_file>){
     chomp;
#    print if($debug);

     if(m/^>/){#in fasta header get line and split at spaces. replace the values in transcript_name and length
     
      #get rid of > at the begining of transcript name
      s/>//g;
    
      #first, is sequence is not empty is because we just finish reading one, store it using the previously stored transcript_name
      if($sequence){#note that this will only be true after reading the first header. Therefore the next push is safe.
      
        print $sequence, "\n" if($debug);
        push @{${$hash_ref}{$transcript_name}}, $sequence;
        
        #now clear the sequence
        $sequence = "";
      
      }
      
      ($transcript_name, $length) = split(' ');
      
      print $transcript_name, "->", $length, "\n" if($debug);
 
      #store the new sequences in hash
      ${$hash_ref}{$transcript_name}=[$length];
      
      $transcript_counter++;
    }else{#in sequence, concatenate the sequence
    
      $sequence .= $_;   
    
    }
  }
  push @{${$hash_ref}{$transcript_name}}, $sequence;#we need this push to save the last record in the file
  close $transcript_file;
}


#table reader, expects file path, a hash ref, the fields to be extracted and the field separator
sub blast_table_reader{
  my $table_file = shift;
  my $hash_ref = shift;
  my $array_ref = shift;#fields
  my $sep = shift;
	
  #open filehandles
  open my $table, '<', $table_file or die "Cannot open the table file\n$!";

  print @{$array_ref};

  while(<$table>){
    chomp;
    my @table_fields = split($sep);
    
    #at this point, the hash should be filled and all transcripts should be there.
  
    if(@{$array_ref}){
    
      my $selected_fields = "";

      foreach my $field (@{$array_ref}){
      
        #create line for debug
        $selected_fields .= $table_fields[$field] . "\t";
      
      }
      
      print $selected_fields, "\n" if($debug);
      
      ${$hash_ref}{$table_fields[0]}=$selected_fields;
      
    }else{
      
      my $default_fields = $table_fields[1] . "\t" . $table_fields[10] . "\t" . $table_fields[24]; 
      
      print $table_fields[1], "\t", $table_fields[10], "\t", $table_fields[24], "\n" if($debug);
      
      ${$hash_ref}{$table_fields[0]}=$default_fields;
      
    }
  }
  close $table;
}

#fasta reader, expects file path and hash ref
sub transdecoder_reader{
  my $fasta_file = shift;
  my $hash_ref = shift;
  
  #open filehandles
  open my $transdecoder_file, '<', $fasta_file or die "Cannot open the transdecoder file\n$!";
  
  my @transcript_header;
  my $transcript_name;
  my $sequence="";

  #read the transcript file an populate hash
  while(<$transdecoder_file>){
     chomp;
#    print if($debug);

     if(m/^>/){#in fasta header get line and split at spaces. replace the values in transcript_name and length
     
      #get rid of > at the begining of transcript name
      s/>//g;
    
      #first, if sequence is not empty is because we just finish reading one, store it using the previously stored transcript_name
      if($sequence){#note that this will only be true after reading the first header. Therefore the next push is safe.
      
        print $sequence, "\n" if($debug);
        push @{${$hash_ref}{$transcript_name}}, $sequence;
        
        #now clear the sequence
        $sequence = "";
      
      }
      
      @transcript_header = split(/\|/);
      $transcript_name = $transcript_header[0] . "|" . $transcript_header[1];

      my @type = split(' ', $transcript_header[8]);

     #store the new sequences in hash
      ${$hash_ref}{$transcript_name}=[$type[1]];
      
    }else{#in sequence, concatenate the sequence
    
      $sequence .= $_;   
    
    }
  }
  push @{${$hash_ref}{$transcript_name}}, $sequence;#we need this push to capture the last sequence in the file.
  close $transdecoder_file;
}


sub kegg_reader{
  my $fasta_file = shift;
  my $hash_ref = shift;
  
  #open filehandles
  open my $keggs, '<', $fasta_file or die "Cannot open the kegg annotation file\n$!";
  
  my @kegg_annotations;
  my $transcript_name;
  my $sequence="";

  #read the transcript file an populate hash
  while(<$keggs>){
     chomp;
#    print if($debug);
 
      @kegg_annotations = split(/\|/);
      $transcript_name = $kegg_annotations[0] . "|" . $kegg_annotations[1];
      
      my ($peptide_name, $kegg) = split(" ", $kegg_annotations[2]);

      #store the annotation in hash
      ${$hash_ref}{$transcript_name}=$kegg;
  }
  close $keggs;
}


sub go_reader{

  my $go_file = shift;
  my $hash_ref = shift;
  
  #open filehandles
  open my $go_annotations, '<', $go_file or die "Cannot open the GOs file\n$!";
  
  while(<$go_annotations>){
    chomp;
    
    my ($uniprot_accession, $transcript_id, $go_terms) = split("\t");
    
    ${$hash_ref}{$transcript_id} = $go_terms;
  
  }

  close $go_annotations;
}

sub deg_reader{
  
  my $deg_file = shift;
  my $hash_ref = shift;
  my $separator = shift;

  #open filehandles
  open my $degs_handle, '<', $deg_file or die "Cannot open the DEGs file\n$!";
  
  while(<$degs_handle>){
    chomp;
    
    my @deg_data = split($separator);
    
    ${$hash_ref}{$deg_data[0]} = "$deg_data[2]\t$deg_data[6]";
    
  }

  close $deg_file;
}


sub pfam_reader{

  my $pfam_file = shift;
  my $hash_ref = shift;
  
  
  #open filehandles
  open my $pfam_filehandle, '<', $pfam_file or die "Cannot open the pfam file\n$!";
  
  #read the transcript file an populate hash
  while(<$pfam_filehandle>){
     chomp;

     my @pfam_annotations;
     my $transcript_name;
  
     if(m/^#/){
     
     }#do nothing
     else{
      #get rid of > at the begining of transcript name
      s/\s+/ /g;
      @pfam_annotations = split(/\|/);
      $transcript_name = $pfam_annotations[0] . "|" . $pfam_annotations[1];
      
      my @fields = split(" ", $pfam_annotations[2]);
      my $pfam_phrase = $fields[5] . "," . $fields[6];

      if(exists ${$hash_ref}{$transcript_name}){

        push @{${$hash_ref}{$transcript_name}}, $pfam_phrase;
                
      }else{
      
        ${$hash_ref}{$transcript_name} = [$pfam_phrase];
        
      }
     }
   }

  close $pfam_filehandle;
}

#main loop
sub main {

	#check that some required options are defined;
	die "\n$help\n" unless defined $transcripts;

	&print_call if $verbose;

  #data structure to store the annotations. It will use the Transcript names as keys and store the values in an array
  #column order in the array will be: sequence, length, transdec_mRNA, transdec_protein, Completeness, uniprot matching accession, e-value of the match, name of uniprot match, organism, Go_Func, Go_Comp, Go_Proc, Other matching accesion, e-value, 
  my %transcript_table;
  transcript_reader($transcripts, \%transcript_table);
  
  my %uniprot_table;
  blast_table_reader($uniprot, \%uniprot_table, \@uniprot_fields, "\t") if(defined $uniprot);
  
  my %other_table;
  blast_table_reader($other_annotations, \%other_table, \@other_fields, "\t") if(defined $other_annotations);
  
  my %transdecoder_table;
  transdecoder_reader($transdecoder_peptides, \%transdecoder_table) if(defined $transdecoder_peptides);
  
  my %component_go_table;
  go_reader($component, \%component_go_table) if(defined $component);
  
  my %process_go_table;
  go_reader($process, \%process_go_table) if(defined $process);

  my %function_go_table;
  go_reader($function, \%function_go_table) if(defined $function);
  
  my %bacteria_table;
  blast_table_reader($bacteria, \%bacteria_table, \@bacteria_fields, "\t") if(defined $bacteria);
  
  my %deg_table;
  deg_reader($degs, \%deg_table, ",") if(defined $degs);
  
  my %pfam_table;
  pfam_reader($pfam_file, \%pfam_table) if(defined $pfam_file);
  
  my %kegg_table;
  kegg_reader($kegg_file, \%kegg_table) if(defined $kegg_file);

  
#join all annotations

  print "Transcript_name\tLength\tSequence\tUniprot_match\tevalue\tUniprot_match_annotation\tAqu2_match\tevalue\tAqu2_match_annotation\tORF_Type\tProtein\tGO_Component\tGO_Function\tGO_Process\tPfam_Domains\tKEGG\tBacteria_match\tevalue\tBacteria_match_annotation\tLog_Fold_Change\tAdjusted_p_value\n";
  
  my @keys = sort(keys(%transcript_table));

  foreach my $key (@keys){
  
    my $print_line = $key . "\t" . join("\t", @{$transcript_table{$key}});
  
    if($uniprot_table{$key}){#the contig has a uniprot annotation
    
      $print_line .= "\t" . $uniprot_table{$key};
    
    }else{#insert blank columns

      my $times = $#uniprot_fields > 0 ? $#uniprot_fields : 3;
      $print_line .= "\tNA" x $times;
    
    }
    
    #other blast annotations    
  
    if($other_table{$key}){#the contig has a uniprot annotation
    
      $print_line .= "\t" . $other_table{$key};
    
    }else{#insert blank columns
    
      my $times = $#uniprot_fields > 0 ? $#uniprot_fields : 3;
      $print_line .= "\tNA" x $times;
    
    }
  
    #transdecoder peptides
    if($transdecoder_table{$key}){

      $print_line .= "\t" . join("\t", @{$transdecoder_table{$key}});
      
    }else{
    
      $print_line .= "\tNA" x 2; #there are two transdecoder annotations kept
      
    }
    
    #component gos
    if($component_go_table{$key}){
    
      $print_line .= "\t" . $component_go_table{$key};
    
    }else{
    
      $print_line .= "\tNA";
    
    }
  
    #function gos
    if($function_go_table{$key}){
    
      $print_line .= "\t" . $function_go_table{$key};
    
    }else{
    
      $print_line .= "\tNA";
    
    }  
  
    #process gos
    if($process_go_table{$key}){
    
      $print_line .= "\t" . $process_go_table{$key};
    
    }else{
    
      $print_line .= "\tNA";
    
    }  
    
    #pfam annotations
    if($pfam_table{$key}){

      $print_line .= "\t" . join(";", @{$pfam_table{$key}});
    
    }else{
      
      $print_line .= "\tNA";
    
    }
    
    #kegg_annotations
    if($kegg_table{$key}){
    
      $print_line .= "\t" . $kegg_table{$key};
    
    }else{
    
      $print_line .= "\tNA";
    
    }
    
    #bacteria annotations
    if($bacteria_table{$key}){#the contig has a uniprot annotation
    
      $print_line .= "\t" . $bacteria_table{$key};
    
    }else{#insert blank columns

      my $times = $#bacteria_fields > 0 ? $#bacteria_fields : 3;
      $print_line .= "\tNA" x $times;
    
    }
    
    #degs
    $key =~ s/(.+)_i\d+/$1/g;
    if($deg_table{$key}){
    
      $print_line .= "\t" . $deg_table{$key};
    
    }else{
    
      $print_line .= "\tNA" x 2;
    
    }
    
  
    print $print_line, "\n";
  
  }
  
 	return 0;
}

exit &main;

