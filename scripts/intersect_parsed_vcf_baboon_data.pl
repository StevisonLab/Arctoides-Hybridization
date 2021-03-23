#!/usr/bin/perl

#Part of Step4 (old Step5) - see Step4 shell script in main directory!!!

#Program reads in intersection VCF file from Step3
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  BGI-CR-5        M_Arctoides     M_Assamensis    M_Fascicularis  M_Thibetana
#chr20   2       .       A       .       41.31   .       AN=10;DP=169    GT:AD:DP:RGQ    0/0:61,0:61:120 0/0:21,0:21:63  ./.     0/0:26,0:26:78  0/0:42,0:42:120
#chr20   4       .       A       .       43.31   .       AN=10;DP=176    GT:AD:DP:RGQ    0/0:61,0:61:120 0/0:25,0:25:75  ./.     0/0:28,0:28:84  0/0:42,0:42:120
#chr20   6       .       G       .       52.31   .       AN=10;DP=189    GT:AD:DP:RGQ    0/0:61,0:61:120 0/0:27,0:27:81  ./.     0/0:37,0:37:111 0/0:42,0:42:120
#chr20   8       .       C       .       49.31   .       AN=10;DP=198    GT:AD:DP:RGQ    0/0:61,0:61:120 0/0:27,0:27:81  ./.     0/0:37,0:37:111 0/0:50,0:50:120
#chr20   9       .       T       .       50.31   .       AN=10;DP=198    GT:AD:DP:RGQ    0/0:61,0:61:120 0/0:27,0:27:81  ./.     0/0:37,0:37:111 0/0:50,0:50:120
#chr20   12      .       C       .       32.31   .       AN=10;DP=200    GT:AD:DP:RGQ    0/0:61,0:61:120 0/0:27,0:27:81  ./.     0/0:37,0:37:111 0/0:50,0:50:120
#chr20   13      .       A       .       32.31   .       AN=10;DP=200    GT:AD:DP:RGQ    0/0:61,0:61:120 0/0:27,0:27:81  ./.     0/0:37,0:37:111 0/0:50,0:50:120


#Program reads in baboon vcf from lastz alignments between baboon and rhesus
##fileformat=VCFv4.0
#CHROM  POS     ID      rheMac2 papAnu2 QUAL    FILTER  INFO
#chr19   429     .       T       T       40080   CHROM.POS:%ID:ORI:LEN   chr4.43936:91.9:-:511
#chr19   430     .       A       A       40080   CHROM.POS:%ID:ORI:LEN   chr4.43935:91.9:-:511
#chr19   431     .       C       C       40080   CHROM.POS:%ID:ORI:LEN   chr4.43934:91.9:-:511


#input1: Parsed VCF file
$vcf=$ARGV[0];

#input2: Baboon VCF file
$baboon=$ARGV[1];

$output=$ARGV[2];

#usage
unless ($#ARGV==2) {
    print STDERR "Error, please provide input parsed VCF file and baboon VCF file.\n\n";
    die;
} #end unless

#open inputs
open(VCF,$vcf);
open(BABOON,$baboon);

#open output
#$output=$vcf;
#$output=~s/3_first_intersected_VCFs/4_parsed_VCFs/;
#$output=~s/.output.vcf/.wBaboon.vcf/;

$started_baboon = 0;

open(OUTPUT,">$output");

print STDERR "Opened two inputs $vcf and $baboon.\n";

while (<VCF>) {
    chomp;
    
    if ($_=~/^##/) {
	print OUTPUT "$_\n"; #print header lines to output
    } elsif ($_=~/^#/) {
	print OUTPUT "$_\tpapAnu4\n";  #print main header line to output adding Baboon sample
    } else {
	@input_array=split(/\s+/,$_);


	$POSITION_BED="$input_array[0]\t$input_array[1]\t$input_array[2]";

	$GENOTYPES="$input_array[9]";

	for ($g=10; $g<=$#input_array;$g++){
	    $GENOTYPES=$GENOTYPES . "\t$input_array[$g]";
	} #end for

	#test last baboon position against current line of VCF
	if ($input_array[1]==$baboon_position && $input_array[3]=~/$baboon_array[3]/ && $started_baboon==1) {
	    if ($input_array[3]=~/$baboon_array[4]/) {  #baboon base matches reference allele
		print OUTPUT "$POSITION_BED\t$input_array[3]\t$input_array[4]\t$input_array[5]\t$input_array[6]\t$input_array[7]\t$input_array[8]\t$GENOTYPES\t0/0:$baboon_geno_info\n";
	    } elsif ($input_array[4]=~/$baboon_array[4]/ && $input_array[4]!~/[\.,]/) { #baboon base matches alt allele, alt not multiples or empty
		print OUTPUT "$POSITION_BED\t$input_array[3]\t$input_array[4]\t$input_array[5]\t$input_array[6]\t$input_array[7]\t$input_array[8]\t$GENOTYPES\t1/1:$baboon_geno_info\n";
	    } elsif ($input_array[4]=~/$baboon_array[4]/ && $input_array[4]=~/,/) { #baboon base matches alt allele, but multiples exist
		@allele_array=split(",",$input_array[4]);
		for ($i=0; $i<=$#allele_array; $i++) {
		    if ($allele_array[$i]=~/$baboon_array[4]/) {
			$k=$i+1;
			print OUTPUT "$POSITION_BED\t$input_array[3]\t$input_array[4]\t$input_array[5]\t$input_array[6]\t$input_array[7]\t$input_array[8]\t$GENOTYPES\t$k/$k:$baboon_geno_info\n";
		    } #end if
		} #end for
	    } elsif ($input_array[4]!~/$baboon_array[4]/ &&  $input_array[4]=~/\./) { #no match, make new alt allele
		print OUTPUT "$POSITION_BED\t$input_array[3]\t$baboon_array[4]\t$input_array[5]\t$input_array[6]\t$input_array[7]\t$input_array[8]\t$GENOTYPES\t1/1:$baboon_geno_info\n";
	    } elsif ($input_array[4]!~/$baboon_array[4]/ &&  $input_array[4]!~/[\.,]/) { #no match, make 2nd alt allele
		print OUTPUT "$POSITION_BED\t$input_array[3]\t$input_array[4],$baboon_array[4]\t$input_array[5]\t$input_array[6]\t$input_array[7]\t$input_array[8]\t$GENOTYPES\t2/2:$baboon_geno_info\n";
	    } else {
		print STDERR "Error: site $input_array[1] and $baboon_position match, but nothing written to output:$baboon_geno_info\n";
	    } #end if
	    next;  #next line of VCF, skip going into baboon loop
	} elsif ($input_array[1]<$baboon_position) {  #baboon position still greater than VCF position
	    print OUTPUT "$POSITION_BED\t$input_array[3]\t$input_array[4]\t$input_array[5]\t$input_array[6]\t$input_array[7]\t$input_array[8]\t$GENOTYPES\t./.:.:.:.:.\n";
	    next;  #next line of VCF, skip going into baboon loop
	} #end if

	while (<BABOON>) {
		chomp;
		$started_baboon=1;
		@baboon_array=split(/\s+/,$_);
		$baboon_position = ($baboon_array[1]/1)+1; #add 1 to baboon positions as lastz coordinates are zero based, not 1-based

		@baboon_info=split(/:/,$baboon_array[7]);

		$baboon_geno_info=".:$baboon_array[5]:$baboon_info[1]:$baboon_info[2]";

		if ($input_array[1]==$baboon_position && $input_array[3]=~/$baboon_array[3]/) {
		    if ($input_array[3]=~/$baboon_array[4]/) {  #baboon base matches reference allele
		    print OUTPUT "$POSITION_BED\t$input_array[3]\t$input_array[4]\t$input_array[5]\t$input_array[6]\t$input_array[7]\t$input_array[8]\t$GENOTYPES\t0/0:$baboon_geno_info\n";
		    } elsif ($input_array[4]=~/$baboon_array[4]/ && $input_array[4]!~/[\.,]/) { #baboon base matches alt allele, alt not multiples or empty
		    print OUTPUT "$POSITION_BED\t$input_array[3]\t$input_array[4]\t$input_array[5]\t$input_array[6]\t$input_array[7]\t$input_array[8]\t$GENOTYPES\t1/1:$baboon_geno_info\n";
		    } elsif ($input_array[4]=~/$baboon_array[4]/ && $input_array[4]=~/,/) { #baboon base matches alt allele, but multiples exist
			@allele_array=split(",",$input_array[4]);
			for ($i=0; $i<=$#allele_array; $i++) {
			    if ($allele_array[$i]=~/$baboon_array[4]/) {
				$k=$i+1;
				print OUTPUT "$POSITION_BED\t$input_array[3]\t$input_array[4]\t$input_array[5]\t$input_array[6]\t$input_array[7]\t$input_array[8]\t$GENOTYPES\t$k/$k:$baboon_geno_info\n";
			    } #end if
			} #end for
		    } elsif ($input_array[4]!~/$baboon_array[4]/ &&  $input_array[4]=~/\./) { #no match, make new alt allele
		    print OUTPUT "$POSITION_BED\t$input_array[3]\t$baboon_array[4]\t$input_array[5]\t$input_array[6]\t$input_array[7]\t$input_array[8]\t$GENOTYPES\t1/1:$baboon_geno_info\n";
		    } elsif ($input_array[4]!~/$baboon_array[4]/ &&  $input_array[4]!~/\./) { #no match, make 2nd alt allele
		    print OUTPUT "$POSITION_BED\t$input_array[3]\t$input_array[4],$baboon_array[4]\t$input_array[5]\t$input_array[6]\t$input_array[7]\t$input_array[8]\t$GENOTYPES\t2/2:$baboon_geno_info\n";
		    } else {
			print STDERR "Error: site $input_array[1] and $baboon_position match, but nothing written to output:$baboon_geno_info\n";
		    } #end if
		    last;
		} elsif ($input_array[1]<$baboon_position) {
		    print OUTPUT "$POSITION_BED\t$input_array[3]\t$input_array[4]\t$input_array[5]\t$input_array[6]\t$input_array[7]\t$input_array[8]\t$GENOTYPES\t./.:.:.:.:.\n";
		    last;
#		} elsif ($input_array[1]>$baboon_position) { #skipped site of baboon position, now past it!
#		    print OUTPUT "$POSITION_BED\t$input_array[3]\t$input_array[4]\t$input_array[5]\t$input_array[6]\t$input_array[7]\t$input_array[8]\t$GENOTYPES\t./.:$baboon_geno_info\n";
#		    last;
		} else {
		    next;
		} #end if
		
	} #end else
	    
    } #end BABOON VCF while
    
} #end MAIN VCF while
