#! /usr/bin/perl -w
# Matthew Bashton 9th of May 2013
# Takes the tab delimited output of cuffdiff and adds ensembl annotation needs a biomart dump of the required annotation

use strict;

# Get command line arguments

my $anno_file = $ARGV[0];
my $exp_file = $ARGV[1];

unless (defined $anno_file && $exp_file) {
    die "\n*** You need to supply both a annotation file from ensembl and the results of cuffdiff e.g. gene_exp.diff ***\n\n";
}

my %annotation_data;
my %alt_annotation_data;

# Read in all the ensemble annotation data NOTE: you'll need to hack this bit if you export differing amounts of cols from ensemble, 13 cols in total
# Also you'll need to get rid of the empty fields between successive characters using perl -p -e "s/ (?:^|(?<=\t)) (?=\t|$) /null/xg" CPannotation.tab
# Should really put that code into this script.

open (INPUT1, "$anno_file");
while (<INPUT1>) {
    if (/^(ENSCPOG\S+)\t(\S+)\t(.+)\t(MT|scaffold\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)/) {
    #      1              2     3       4                5      6      7      8      9      10     11     12     13
	# $1 Ensembl Gene ID Key
	# $2 Ensembl Transcript ID 0
	# $3 Description 1
	# $4 Chromosome Name 2
	# $5 Gene Start (bp) 3
	# $6 Gene End (bp) 4
	# $7 Strand 5
	# $8 Transcript Start (bp) 6
	# $9 Transcript End (bp) 7
	# $10 Associated Gene Name 8
	# $11 Associated Transcript Name 9
	# $12 UniProt/SwissProt Accession, might not be present 10
	# $13 UniProt/TrEMBL Accession, might not be present 11

	#print "$1\t$2\t$3\t$4\t$5\t$6\t$7\t$8\t$9\t$10\t$11\t$12\t$13\n";
	$annotation_data{$1} = [$2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13];
	$alt_annotation_data{$10} = [$1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13];
    }
}
close INPUT1;


# Next open up the cuffdiff output whilst working on this will integrate the annotation

open (INPUT2, "$exp_file");
while (<INPUT2>) {
    if (/^(XLOC\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)/) {
	#  1           2      3      4      5      6      7      8      9      10     11     12     13     14
        # 1 test_id
	# 2 gene_id
	# 3 gene
	# 4 locus   
	# 5 sample_1
	# 6 sample_2
	# 7 status
	# 8 value_1
	# 9 value_2
	# 10 log2(fold_change)
	# 11 test_stat
	# 12 p_value
	# 13 q_value
	# 14 significant

	# Print the cuffdiff stuff
	print "$1\t$2\t$3\t$4\t$5\t$6\t$7\t$8\t$9\t$10\t$11\t$12\t$13\t$14";

	my $CD1 = $1;
	my $CD2 = $2;
	my $CD3 = $3;
	my $CD4 = $4;
	my $CD5 = $5;
	my $CD6 = $6;
	my $CD7 = $7;
	my $CD8 = $8;
	my $CD9 = $9;
	my $CD10 = $10;
	my $CD11 = $11;
	my $CD12 = $12;
	my $CD13 = $13;
	my $CD14 = $14;

	#print "$1:  ";
	# Print the extra annotation based on matching $1 to the key of %annotation_data

	# Flags
	my $anno = 1;
	my $alt_anno = 1;

	# Alter gene to single entry allow a tab or a comma
	$CD3 =~ /(\w+)/;
	my $gene = $1;

	# Check to see if any annotation is around for the transcript id
	if (!defined $annotation_data{$gene}[0]) {
	    $anno = 0;
	}
	
	if (!defined $alt_annotation_data{$gene}[0]) {
	    $alt_anno = 0;
	}
    

	if ($anno == 0 && $alt_anno == 0) {
	    print "\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\n";
	    #        1  2  3  4  5  6  7  8  9  10 11 
	    next;
	}
	
	if ($anno == 1) {
	    print "\t$gene\t$annotation_data{$gene}[0]\t$annotation_data{$gene}[1]\t$annotation_data{$gene}[2]\t$annotation_data{$gene}[3]\t$annotation_data{$gene}[4]\t$annotation_data{$gene}[5]\t$annotation_data{$gene}[6]\t$annotation_data{$gene}[7]\t$annotation_data{$gene}[8]\t$annotation_data{$gene}[9]\t$annotation_data{$gene}[10]\t$annotation_data{$gene}[11]\n";
	}
	
	if ($anno == 0 && $alt_anno == 1) {
	    print "\t$alt_annotation_data{$gene}[0]\t$alt_annotation_data{$gene}[1]\t$alt_annotation_data{$gene}[2]\t$alt_annotation_data{$gene}[3]\t$alt_annotation_data{$gene}[4]\t$alt_annotation_data{$gene}[5]\t$alt_annotation_data{$gene}[6]\t$alt_annotation_data{$gene}[7]\t$alt_annotation_data{$gene}[8]\t$alt_annotation_data{$gene}[9]\t$alt_annotation_data{$gene}[10]\t$alt_annotation_data{$gene}[11]\t$alt_annotation_data{$gene}[12]\n";
	}
    }
}
close INPUT2;
