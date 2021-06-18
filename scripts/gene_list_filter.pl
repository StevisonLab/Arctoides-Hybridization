#!/usr/bin/perl

#use strict;
#use warnings;

my $input="/home/aublys/Macaque_New/ensGene.filtered.bed";
my $output="gene_list.bed"; 

open(INPUT,"$input");
open(OUTPUT,">$output");

#initialize variables
my $last_gene="";
my $last_txstart="1";
my $last_txend="1";
my $last_chr="";
my $last_gene2="";

$lc=0;
$gene_1=1;

while (<INPUT>){
    chomp;
    my @input_array=split("\t",$_);

    my $chr=$input_array[0];
    my $txstart=$input_array[1];
    my $txend=$input_array[2];
    my $gene=$input_array[4];
    my $gene2=$input_array[5];
    my $ori=$input_array[6];

    if ($lc <= 10) {
	print STDERR "Gene: $gene ; Last gene: $last_gene ; Orientation: $ori\n";
	$lc++;
    } 

    if ( $gene eq $last_gene) { 
	if ( $txstart < $last_txstart)  {
	    $max_txstart=$txstart;
	}
	
	if ( $txend > $last_txend) {
	    $max_txend=$txend;
	}

	if ($gene_1==1) {
	    $gene_1=0;
	    $max_txstart=$txstart;
	    $max_txend=$txend;
	}

#    } elsif ( $gene eq $last_gene && $ori eq "-" ) { 
#	if ( $txstart > $last_txstart) {
#	    $max_txstart=$txstart;
#	}
#	
#	if ( $txend < $last_txend) {
#	    $max_txend=$txend;
#	}
#
#	if ($gene_1==1) {
#	    $gene_1=0;
#	    $max_txstart=$txstart;
#	    $max_txend=$txend;
#	}

    } elsif ( $gene ne $last_gene && $lc > 1) {
	#add flanking
	$start_adj=$max_txstart - 5000;
	$end_adj=$max_txend + 5000;

	if ($start_adj <0) {
	    $start_adj=0;
	}

	#print last gene
#	print OUTPUT "$last_chr\t$max_txstart\t$max_txend\t$last_gene\t$last_gene2\t$last_ori\n";

	#print adjusted gene coordinates
	print OUTPUT "$last_chr\t$start_adj\t$end_adj\t$last_gene\t$last_gene2\t$last_ori\n";
	$gene_1=1;
	$max_txstart=$txstart;
	$max_txend=$txend;
    } else {
#	print STDERR "LC: $lc\n";
	$max_txstart=$txstart;
	$max_txend=$txend;
    } 
    
    $last_ori=$ori;
    $last_gene=$gene;
    $last_txstart=$txstart;
    $last_txend=$txend;
    $last_chr=$chr;
    $last_gene2=$gene2;
	
}
