#! /usr/bin/perl

#open axt file and convert to VCF output
#last updated July 14, 2019


$input=$ARGV[0];
$output=$ARGV[1];

unless ($#ARGV==1) {
    print STDERR "Please provide input and output filenames on command line.\n\n";
    die;
} #end unless 

open(INPUT, $input);
open(OUTPUT1, ">$output");

#print OUTPUT1 "##fileformat=VCFv4.0\n";
#print OUTPUT1 "#CHROM\tPOS\tID\trheMac8\tpapAnu4\tQUAL\tFILTER\tINFO\n";

$new_alignment=0;
$score=0;
$new_identity=0;
$identity=0;

while (<INPUT>) {
    chomp;

    # identity=255/258 (98.8%)
    # coverage=258/53671032 (0.0%)
    # continuity=258/258 (100.0%)
    # cigar=258m
    #a score=23884

    if($_=~/^# identity/) { #capture percent identity for new alignment
	@identity_array=split(/\s+/,$_);
	$identity_array[2]=~s/\(//;
	$identity_array[2]=~s/%\)//;
	$new_identity=$identity_array[2];
    } #end if

    if ($_=~/^a/) {  #new alignment

	$new_alignment=1;
	@sequence_one=split("",$sequence1);
	@sequence_two=split("",$sequence2);

	for ($x=0; $x<=$#sequence_one;$x++) {
	    unless ($sequence_one[$x]=~/-/) {
		print OUTPUT1 "$ref\t$sequence_start\t.\t$sequence_one[$x]\t$sequence_two[$x]\t$score\tCHROM.POS:%ID:ORI:LEN\t$alt.$alt_seq_start:$identity:$ali_ori:$ali_length\n";
		$sequence_start++;
		if($ori==1) {
		    $alt_seq_start++;
		} else {
		    $alt_seq_start--;
		} #end else
	    } #end unless
	} #end for

	#store score of new alignment
	@score_array=split(/\s+/,$_);
	$score_array[1]=~s/score=//;
	$score=$score_array[1];

	#store new identity as current identity
	$identity=$new_identity;

	next;
    } #end if

    if ($_=~/^s/ && $new_alignment==1) {  #species line #1
	$new_alignment=0;
	
	#s chr19 33779994 211 +  53671032 

	@input=split(/\s+/,$_);
#	@chr_name=split(/\./,$input[1]);
	$ref=$input[1];
	$sequence1=$input[6];
	$sequence_start=$input[2];
	$ali_length=$input[3];

#	print STDERR "Species 1: $input[1]\n";
    } elsif ($_=~/^s/ && $new_alignment==0) {  #species line #2
	@input=split(/\s+/,$_);
#	@chr_name=split(/\./,$input[1]);
	$alt=$input[1];
        $sequence2=$input[6];
	$alt_seq_start=$input[2];

	$ali_ori=$input[4];

	if($input[4]=~/-/) {
	    $ori=0;
	} else {
	    $ori=1;
	} #end if
#	print STDERR "Species 2: $input[1]\n";
    } #end if

} #end while

@sequence_one=split("",$sequence1);
@sequence_two=split("",$sequence2);

for ($x=0; $x<=$#sequence_one;$x++) {
    unless ($sequence_one[$x]=~/-/) {
	print OUTPUT1 "$ref\t$sequence_start\t.\t$sequence_one[$x]\t$sequence_two[$x]\t$score\tCHROM.POS:%ID:ORI:LEN\t$alt.$alt_seq_start:$identity:$ali_ori:$ali_length\n";
	$sequence_start++;
	if ($ori==1) {
	    $alt_seq_start++;
	} else {
	    $alt_seq--;
	} #end else
    } #end unless
} #end for
