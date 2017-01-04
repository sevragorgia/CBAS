#!/usr/bin/perl

#Sergio Vargas R. 2015. Munich, Germany

use strict;
use warnings;
use Getopt::Long;
use File::Basename;


#ADD GNU license here!

#command line options, the variable names should be self explanatory;
my $verbose;
my $debug;
my $save_tmps = "";
my $back_file = undef;
my $deg_file = undef;
my $component_file;
my $function_file;
my $process_file;
my $fuse_isoforms;

my $help = "Get_GO_Annotations.pl

The following options have to be provided:
		--background; a file listing the names of the background transcripts
    --degs; a file listing the names of the differentially expressed transcripts
    --component; a file with the GO component annotations 
    --function; a file with the GO function annotations 
    --process; a file with the GO process annotations 
    --fuse_isoforms; fuse trinity isoforms in 'genes'
 		--verbose F
 		--debug F
 		
A normal command line would look like:
	perl Get_GO_Annotations.pl --background back.csv --degs up.csv --component component_go_annotations.csv

The options --background and --degs requires you to pass a list of transcripts of interest which GOs will be fetched.

The list is simply a list of transcripts with no headers, footers or blank lines, e.g.:

TR4|c9_g1_i1
TR41535|c9_g1_i1
TR41536|c9_g1_i1

Note:

this script expects trinity transcripts annotated using blast agains SWISSPROT (with output 6). This looks like:

B8GNX2	TR3|c0_g1_i1	GO:0005524,GO:0008270,GO:0031072,GO:0046872,GO:0051082,GO:0051082

or 

TR3|c0_g1_i1	GO:0005524,GO:0008270,GO:0031072,GO:0046872,GO:0051082,GO:0051082

it extracts the transcript name and the associated annotations and saves them if the transcript is list of background, up or down regulated genes.

The script produces output files ending wiht .function.gos, .process.gos and component.gos with the go terms for each transcript as follows:

TR41535|c9_g1_i1	GO:XXXX,GO:YYYY

these files should be then usable for enrichment analyses pipelines like top_go.

if the flag --fuse_isoforms is used. The script will fuse trinity isoforms in genes, for instance, TR41535|c9_g1_i1 TR41535|c9_g1_i2 TR41535|c9_g1_i3 will be collapsed to TR41535|c9_g1

";

#get options from command line
GetOptions(
           "background=s" => \$back_file,
           "degs=s" => \$deg_file,
           "component=s" => \$component_file,
           "function=s" => \$function_file,
           "process=s" => \$process_file,
           "fuse_isoforms" => \$fuse_isoforms,
					 "verbose" => \$verbose,
					 "debug" => \$debug,
					 "help" => \&help
					 ) or die("Error in the command line arguments\n$help");

sub help{
	print $help;
}

#taken and modified from http://search.cpan.org/~szabgab/Array-Unique-0.08/lib/Array/Unique.pm
sub Unique {
  my $ref_array = shift;
	my %seen = ();
	my @unique = ();
	foreach my $item (@{$ref_array}) {
	  print $item, "\n" if($debug);
		unless ($seen{$item}){
       # if we get here we have not seen it before
       $seen{$item} = 1;
       push (@unique, $item);
  	}
  }
  return @unique; 
}

sub go_reader{

  my $go_file = shift;
  my $hash_ref = shift;
  
  #open filehandles
  open my $go_annotations, '<', $go_file or die "Cannot open the GOs file\n$!";
  
  while(<$go_annotations>){
    chomp;
    
    my ($uniprot_accession, $transcript_id, $go_terms) = split("\t");
    
    $transcript_id =~ s/_i\d+//g if($fuse_isoforms);
    
    print "transcript_id collapsed to ", $transcript_id, "\n" if($debug && $fuse_isoforms);
    
    my @go_terms = split(",", $go_terms);
    
    if(exists ${$hash_ref}{$transcript_id}){

      print "pushing @go_terms in array\n" if($debug);
      push @{${$hash_ref}{$transcript_id}}, @go_terms;
    
    }else{
      
      print "creating hash key with value ", $transcript_id, "\n" if($debug);
      ${$hash_ref}{$transcript_id} = \@go_terms;
     
    }
  }

  close $go_annotations;

}

sub Get_List_of_Queries{#expects a file and a ref to an array
	my $list_file = shift;
  my $ref_array = shift;

  open my $list_file_filehandle, '<' ,$list_file or die "Cannot open list file: ", $!;
	
	my $targets=0;
	
	while(<$list_file_filehandle>){
		chomp;
#		s/\|//g;
#		my ($id, $transcript) = split("\t");
#		print $transcript, "\n" if($debug);
#		push(@target_transcripts, $transcript);

		print $_, "\n" if($debug);
		push(@{$ref_array}, $_);

		$targets++;
		}

	print "$targets in list file\n" if($verbose);

  close $list_file_filehandle;
}

sub Extract_GOs_for_list{#receives a list of genes of interest, and three hash refs where to look for annotations

  my $outname = shift;
  my $list_ref = shift;
  my $functions_ref = shift;
  my $components_ref = shift;  
  my $processes_ref = shift;

  open my $functions_outfile, '>' , ($outname . ".functions.gos") or die "Cannot open functions outfile: ", $!;
  open my $components_outfile, '>' , ($outname . ".components.gos") or die "Cannot open components outfile: ", $!;
  open my $processes_outfile, '>' , ($outname . ".processes.gos") or die "Cannot open processes outfile: ", $!;

  foreach my $item (@{$list_ref}){
    if(exists ${$functions_ref}{$item}){
      if (@{${$functions_ref}{$item}}){
        
        print $item, "\t", Unique(${$functions_ref}{$item}), "\n" if($debug);
        print $functions_outfile $item, "\t", join(",", Unique(${$functions_ref}{$item})), "\n"; 
      
      }else{
      
        print $item, "\thas no annotations\n"  if($debug);
        print $functions_outfile $item, "\tNA\n";
        
      }
      

    }else{
    
      print $item, "\tnot found\n" if($debug);
    
    }


 ######################################
 
    if(exists ${$components_ref}{$item}){
      if (@{${$components_ref}{$item}}){
        
        print $item, "\t", Unique(${$components_ref}{$item}), "\n" if($debug);
        print $components_outfile $item, "\t", join(",", Unique(${$components_ref}{$item})), "\n"; 
      
      }else{
      
        print $item, "\thas no annotations\n" if($debug);
        print $components_outfile $item, "\tNA\n";
        
      }
      

    }else{
    
      print $item, "\tnot found\n" if($debug);
    
    } 
 
 ################################################################
 
    if(exists ${$processes_ref}{$item}){
      if (@{${$processes_ref}{$item}}){
        
        print $item, "\t", Unique(${$processes_ref}{$item}), "\n" if($debug);
        print $processes_outfile $item, "\t", join(",", Unique(${$processes_ref}{$item})), "\n"; 
      
      }else{
      
        print $item, "\thas no annotations\n" if($debug);
        print $processes_outfile $item, "\tNA\n";
                
      }
      

    }else{
    
      print $item, "\tnot found\n" if($debug);
    
    } 
 
  }
  close $functions_outfile; 
  close $components_outfile; 
  close $processes_outfile; 
} 


sub main{
  
  die "\n",$help,"\n" unless($component_file or $function_file or $process_file);

  my %components;
  if(defined $component_file){
    print "reading component GOs\n" if($verbose);
    go_reader($component_file, \%components);  
  }
  
  my %functions;
  if(defined $function_file){
    print "reading function GOs\n" if($verbose);
    go_reader($function_file, \%functions);
  }
  
  my %processes;
  if(defined $process_file){
    print "reading process GOs\n" if($verbose);
    go_reader($process_file, \%processes);  
  }
  
	my @background_list;

  if(defined $back_file){
		
		print "opening $back_file\n" if $verbose;
		Get_List_of_Queries($back_file, \@background_list);
	
	}

	my @deg_list;
	
	if(defined $deg_file){
	
		print "opening $deg_file\n" if $verbose;	
		Get_List_of_Queries($deg_file, \@deg_list);
	
	}

  my ($back_name, $back_dir, $back_suffix) = fileparse($back_file, (".csv"));
  Extract_GOs_for_list(($back_dir . $back_name), \@background_list, \%functions, \%components, \%processes);

  my ($degs_name, $degs_dir, $degs_suffix) = fileparse($deg_file, (".csv"));
  Extract_GOs_for_list(($degs_dir . $degs_name), \@deg_list, \%functions, \%components, \%processes);

	return 0;
}

exit &main();
