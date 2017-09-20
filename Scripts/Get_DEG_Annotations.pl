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
my $sequence;
my $length;
my $gene_name;
my $uniprot;
my $uniprot_evalue;
my $uniprot_annotation;
my $aqu2;
my $aqu2_evalue;
my $aqu2_annotation;
my $transdecoder_orf_type;
my $transdecoder_peptide;
my $component;
my $function;
my $process;
my $pfam;
my $bacteria;
my $bacteria_evalue;
my $bacteria_annotation;
my $kegg;
my $degs;
my $verbose;
my $debug;
my $isoforms;  
my $all;
my $output_file = "./out.csv"; #where should I store the genomes?



my $help = "The following options have to be provided:
    --metatable = the metatable containing all annotations for a transcriptome
    --degs = list of differentially expressed genes or more generally genes of interest to look up in the table

the following flags may be activated and will change the behaviour of the program and its output
    
		--sequence = the contig sequence
		--length = contig length in bp
		--gene_name = name of the trinity gene
		--uniprot = matching uniprot accession
		--uniprot_evalue = evalue of the uniprot best match
		--uniprot_annotation = annotation derived from the uniprot accession
		--aqu2 = matching Amphimedon queenslandica gene
		--aqu2_evalue = evalue of best A. queenslandica match
		--aqu2_annotation = annotation derived from the A. queenslandica match 
		--transorf_type = type of transdecoder orf
		--transpeps = the transdecoder peptide
		--component = component GO annotations for the uniprot match
		--process = component GO annotations for the uniprot match
		--function = component GO annotations for the uniprot match
		--pfam = pfam annotations
		--kegg = kegg annotations
		--bacteria = matching bacterial accession 
		--bacteria_evalue = evalue of the best matching
		--bacteria_annotation = annotation derived from the bacterial accession
		--isoforms = deactivate merging of isoforms 
 		--all = print all available annotations for DEGs
		--verbose
 		--debug

A normal command line would look like:
	perl Get_Annotations.pl --metatable cbas_transcriptome_metatable.tsv --uniprot --transpeps --pfam --degs cbas_degs.txt

This will generate three output files for the uniprot accession numbers, the transdecoder peptides (the longest will be chosen unless --isoforms is active) and the pfam domains. These files can be further use in other programs.


";

#get options from command line
GetOptions("metatable=s" => \$metatable,
<<<<<<< Updated upstream
					 "length" => \$length,
=======
>>>>>>> Stashed changes
					 "sequence" => \$sequence,
					 "uniprot" => \$uniprot,
					 "uniprot_annotation" => \$uniprot_annotation,
					 "uniprot_evalue" => \$uniprot_evalue,
					 "aqu2" => \$aqu2,
					 "aqu2_evalue" => \$aqu2_evalue,
					 "aqu2_annotation" => \$aqu2_annotation,
					 "gene_name" => \$gene_name,
					 "component" => \$component,
					 "function" => \$function,
					 "process" => \$process,
					 "transpeps" => \$transdecoder_peptide,
					 "transorf_type" => \$transdecoder_orf_type,
					 "pfam" => \$pfam,
 					 "bacteria" => \$bacteria,
					 "bacteria_annotation" => \$bacteria_annotation,
					 "bacteria_evalue" => \$bacteria_evalue,
					 "degs=s" => \$degs,
					 "kegg" => \$kegg,
					 "isoforms" => \$isoforms,
					 "verbose" => \$verbose,
					 "debug" => \$debug,
					 "help" => \$help,
					 "all" => \$all
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

  die "\n$help\n" unless(defined $metatable && defined $degs);
  
  my @list_of_degs;
  &deg_reader($degs, \@list_of_degs);
  
  my @extracted_degs;
  
  #open filehandles
  open my $metatable_handle, '<', $metatable or die "Cannot open the table file\n$!";

  chomp(my $table_header = <$metatable_handle>);

  my @table_header_fields = split("\t", $table_header);

  #yes, I know this is very bad programming.
  
  print $table_header_fields[0];

  print "\t $table_header_fields[1]" if($length || $all);
  print "\t $table_header_fields[2]" if($sequence || $all);
  print "\t $table_header_fields[3]" if($gene_name || $all);
  print "\t $table_header_fields[4]" if($uniprot || $all);
  print "\t $table_header_fields[5]" if($uniprot_evalue || $all);
  print "\t $table_header_fields[6]" if($uniprot_annotation || $all);
  print "\t $table_header_fields[7]" if($aqu2 || $all);
  print "\t $table_header_fields[8]" if($aqu2_evalue || $all);
  print "\t $table_header_fields[9]" if($aqu2_annotation || $all);
  print "\t $table_header_fields[10]" if($transdecoder_orf_type || $all);
  print "\t $table_header_fields[11]" if($transdecoder_peptide || $all);
  print "\t $table_header_fields[12]" if($component || $all);
  print "\t $table_header_fields[13]" if($function || $all);
  print "\t $table_header_fields[14]" if($process || $all);
  print "\t $table_header_fields[15]" if($pfam || $all);
  print "\t $table_header_fields[16]" if($kegg || $all);
  print "\t $table_header_fields[17]" if($bacteria || $all);
  print "\t $table_header_fields[18]" if($bacteria_evalue || $all);
  print "\t $table_header_fields[19]" if($bacteria_annotation || $all);

  print "\n"; 

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
      
	print "\t $table_fields[1]" if($length || $all);

	print "\t $table_fields[2]" if($sequence || $all);
     
	print "\t $table_fields[3]" if($gene_name || $all);

        print "\t $table_fields[4]" if($uniprot || $all); 

	print "\t $table_fields[5]" if($uniprot_evalue || $all);

        print "\t $table_fields[6]" if($uniprot_annotation || $all);

	print "\t $table_fields[7]" if($aqu2 || $all);

	print "\t $table_fields[8]" if($aqu2_evalue || $all);

	print "\t $table_fields[9]" if($aqu2_annotation || $all);

        print "\t $table_fields[10]" if($transdecoder_orf_type || $all);

        print "\t $table_fields[11]" if($transdecoder_peptide || $all);

        print "\t $table_fields[12]" if($component || $all);
      
        print "\t $table_fields[13]" if($function || $all);
      
        print "\t $table_fields[14]" if($process || $all);
      
        print "\t $table_fields[15]" if($pfam || $all);

        print "\t $table_fields[16]" if($kegg || $all);

        print "\t $table_fields[17]" if($bacteria || $all);

	print "\t $table_fields[18]" if($bacteria_evalue || $all);

        print "\t $table_fields[19]" if($bacteria_annotation || $all);
      
        print "\n";
        
        push(@extracted_degs, $table_fields[0]) unless($isoforms);#add the name of the gene to the @extracted_degs array to avoid printing the same info for each isoform unless isoforms are wanted
      }
    }
  }
  
  close $metatable_handle;
}




#run the main
exit &main;

