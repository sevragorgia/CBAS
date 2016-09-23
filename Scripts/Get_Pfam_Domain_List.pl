#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;

#ADD GNU license here!


########################
#
#
#perl Get_Pfam_Domain_List.pl --pfam
#
#
###################################################################################################3



#command line options, the variable names should be self explanatory;
my $pfam;
my $verbose;
my $debug;


my $help = "The following options have to be provided:
    --pfam = the pfam annotations per gene/transcript

A normal command line would look like:
	perl Get_Pfam_Domain_List.pl --pfam pfam_domain_per_transcript_file.txt >pfam.list

This will generate the list of pfam domains contained in the given list.


";

#get options from command line
GetOptions("pfam=s" => \$pfam,
 					 "verbose" => \$verbose,
					 "debug" => \$debug,
					 "help" => \$help
					 ) or die("Error in the command line arguments\n$help");


sub pfam_reader{
  
  my $pfam_file = shift;
  my $array_ref = shift;

  #open filehandles
  open my $pfam_handle, '<', $pfam_file or die "Cannot open the Pfam file\n$!";
  
  while(<$pfam_handle>){

    chomp;
    my ($transcript, $domains) = split("\t");
    
    push(@{$array_ref}, $domains);
    
  }

  close $pfam_file;
}

#table reader, expects file path, a hash ref, the fields to be extracted and the field separator
sub main{

  die "\n$help\n" unless(defined $pfam);
  
  my %domain_list;
  my @pfam_domains;
  &pfam_reader($pfam, \@pfam_domains);
  
  foreach my $domain_annotation (@pfam_domains){
  
    $domain_annotation =~ s/\s+//g;
  
    my @domains = split(";", $domain_annotation);#this is a list of PFcodes,domain names
    
    foreach my $domain (@domains){
    
      my ($pf_code, $domain_name) = split(",", $domain);
      
      $domain_list{$pf_code}=$domain_name unless(exists $domain_list{$pf_code});
    
    }
  }  
  
  #print the hash,
  
  foreach my $key (keys %domain_list){
    
    print $key, "\t", $domain_list{$key}, "\n" unless($key eq "NA");
  
  }
    
}




#run the main
exit &main;

