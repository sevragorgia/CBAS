#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;

#ADD GNU license here!


########################
#
#
#perl Get_DEG_Annotations.pl --metatable ~/Desktop/CBAS/CBAS_Annotation_Metatable_20160407.csv --transpeps --degs ~/Desktop/CBAS/DeSeq2/over_degs_p001_2LFC.csv 
#
#
###################################################################################################3



#command line options, the variable names should be self explanatory;
my $metatable;
my $sequence;
my $uniprot;
my $uniprot_annotation;
my $transdecoder_orf_type;
my $transdecoder_peptide;
my $component;
my $function;
my $process;
my $pfam;
my $bacteria;
my $bacteria_annotation;
my $kegg;
my $degs;
my $verbose;
my $debug;
my $isoforms;   
my $output_file = "./out.csv"; #where should I store the genomes?



my $help = "The following options have to be provided:
    --metatable = the metatable containing all annotations for a transcriptome

the following flags may be activated and will change the behaviour of the program and its output
    
    --sequence = the contig sequence
		--uniprot = matching uniprot accession
		--uniprot_annotation = annotation derived from the uniprot accession
		--transorf = the type of transdecoder orf
		--transpeps = the transdecoder peptide
		--component = component GO annotations for the uniprot match
		--process = component GO annotations for the uniprot match
		--function = component GO annotations for the uniprot match
		--pfam = pfam annotations
		--kegg = kegg annotations
		--bacteria = matching bacterial accession 
		--bacteria_annotation = annotation derived from the bacterial accession
		--isoforms = deactivate merging of isoforms 
		--degs = list of differentially expressed genes. If not provided the table will not be filtered.
 		--verbose
 		--debug

A normal command line would look like:
	perl Get_Annotations.pl --metatable cbas_transcriptome_metatable.tsv --uniprot --transpeps --pfam --degs cbas_degs.txt

This will generate three output files for the uniprot accession numbers, the transdecoder peptides (the longest will be chosen unless --isoforms is active) and the pfam domains. These files can be further use in other programs.


";

#get options from command line
GetOptions("metatable=s" => \$metatable,
					 "uniprot" => \$uniprot,
					 "uniprot_annotation" => \$uniprot_annotation,
					 "component" => \$component,
					 "function" => \$function,
					 "process" => \$process,
					 "transpeps" => \$transdecoder_peptide,
					 "transorf" => \$transdecoder_orf_type,
					 "pfam" => \$pfam,
 					 "bacteria" => \$bacteria,
					 "bacteria_annotation" => \$bacteria_annotation,
					 "degs=s" => \$degs,
					 "kegg" => \$kegg,
					 "isoforms" => \$isoforms,
					 "verbose" => \$verbose,
					 "debug" => \$debug,
					 "help" => \$help
					 ) or die("Error in the command line arguments\n$help");


sub deg_reader{
  
  my $deg_file = shift;
  my $array_ref = shift;

  #open filehandles
  open my $degs_handle, '<', $deg_file or die "Cannot open the DEGs file\n$!";
  
  while(<$degs_handle>){

    chomp;
    push(@{$array_ref}, $_);
    
  }

  close $deg_file;
}

#table reader, expects file path, a hash ref, the fields to be extracted and the field separator
sub main{

  die "\n$help\n" unless(defined $metatable);
  
  my @list_of_degs;
  &deg_reader($degs, \@list_of_degs);
  
  my @extracted_degs;
  
  #open filehandles
  open my $metatable_handle, '<', $metatable or die "Cannot open the table file\n$!";

  while(<$metatable_handle>){
    chomp;
    my @table_fields = split("\t");
    #
    #the metatable fields are as follow:
    #
    # 0.Transcript_name -> this names are including isoforms
    # 1.Length	
    # 2.Sequence
    # 3.Uniprot_match	
    # 4.evalue
    # 5.Uniprot_match_annotation	
    # 6.Aqu2_match	
    # 7.evalue	
    # 8.Aqu2_match_annotation	
    # 9.ORF_Type	
    # 10.Protein
    # 11.GO_Component	
    # 12.GO_Function	
    # 13.GO_Process	
    # 14.Pfam_Domains
    # 15.KEGG
    # 16.Bacteria_match	
    # 17.evalue	
    # 18.Bacteria_match_annotation	
    #
    
    my $in_list = 0;
    #at this point, the hash should be filled and all transcripts should be there.
    if(@list_of_degs){
      $table_fields[0] =~ s/_i\d+//g unless($isoforms);#get rid of the isoform part of the name.
      
      $in_list = grep {$table_fields[0] eq $_ } @list_of_degs;
      
    }

    if($in_list){
    
      unless(grep {$table_fields[0] eq $_ } @extracted_degs){
    
        print $table_fields[0];
      
        print "\t $table_fields[2]" if($sequence);
      
        print "\t $table_fields[3]" if($uniprot); 

        print "\t $table_fields[5]" if($uniprot_annotation);

        print "\t $table_fields[9]" if($transdecoder_orf_type);

        print "\t $table_fields[10]" if($transdecoder_peptide);

        print "\t $table_fields[11]" if($component);
      
        print "\t $table_fields[12]" if($function);
      
        print "\t $table_fields[13]" if($process);
      
        print "\t $table_fields[14]" if($pfam);

        print "\t $table_fields[15]" if($kegg);

        print "\t $table_fields[16]" if($bacteria);

        print "\t $table_fields[18]" if($bacteria_annotation);
      
        print "\n";
        
        push(@extracted_degs, $table_fields[0]) unless($isoforms);#add the name of the gene to the @extracted_degs array to avoid printing the same info for each isoform unless isoforms are wanted
      }
    }
  }
  
  close $metatable_handle;
}




#run the main
exit &main;

