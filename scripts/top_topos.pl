#! /usr/bin/perl

use List::Util qw(max);

#program reads in output from Twisst and looks for consective windows with high support for one topology

$input="merged.twisst_final.tsv";

open(INPUT,$input);
open(OUTPUT, ">haplo_summary3.tsv");

#scaffold start end mid sites lnL topo1 topo2 topo3
#    chr16735570284750-242.456990384192
#    chr163818645765150-231.750530288288
#    chr18651119901014150-308.00204200132244
#    chr112024148421281750-242.79882048096

$header=1;
$firstline=1;

while (<INPUT>) {
    chomp;

    if ($header==1) {
	$header=0;
	next;
    } #end if

    @input_array=split(/\s+/,$_);
    @topo=($input_array[6],$input_array[7],$input_array[8]);

    $max=max @topo;

    @topo_short = grep {!/$max/} @topo;

    $second_max=max @topo_short;

    $highest_topo_weight = $input_array[6] + $input_array[7] + $input_array[8];

    $GThalf=1+($highest_topo_weight/2);

    $per_diff=$max/$highest_topo_weight;

    $cutoff = 1;  #1 is for extreme; otherwise use 2/3

#    $per_diff=($max-$second_max)/576;


    if ($firstline==1) {
	print STDERR "Topo: $topo[0],$topo[1],$topo[2]\n";
	print STDERR "Max: $max; Second highest: $second_max; Cutoff = $cutoff\n";
	$firstline=0;
    } #end if

    if ($topo[0]>$topo[1] && $topo[0]>$topo[2] && $per_diff>=$cutoff) { #topo1 highest
	print OUTPUT "$input_array[0]\t$input_array[1]\t$input_array[2]\ttopo1\t$max\t$second_max\n";
    } elsif ($topo[1]>$topo[0] && $topo[1]>$topo[2] && $per_diff>=$cutoff) { # topo2 highest
	print OUTPUT "$input_array[0]\t$input_array[1]\t$input_array[2]\ttopo2\t$max\t$second_max\n";
    } elsif ($topo[2]>$topo[1] && $topo[2]>$topo[0] && $per_diff>=$cutoff) { #topo3 highest
	print OUTPUT "$input_array[0]\t$input_array[1]\t$input_array[2]\ttopo3\t$max\t$second_max\n";
    } elsif ($topo[0]==$topo[1] || $topo[1]==$topo[2] || $topo[2]==$topo[0] || $per_diff<$cutoff) { #two values match or max is below 300
	print OUTPUT "$input_array[0]\t$input_array[1]\t$input_array[2]\tunresolved\t$max\t$second_max\n";
    } else { #there shouldn't be an else
	print STDERR "Error: made it to the else\n";
    } #end if

} #end while

close OUTPUT;

open(NEW,"haplo_summary3.tsv");
open(FINAL,">top_topos_extreme.tsv");

$top=1;

while (<NEW>) {
    chomp;
    @line=split(/\s+/,$_);

    if ($top==1) {
	$last_topo=$line[3];
	$last_chr=$line[0];
	$open=$line[1];
	$close=$line[2];
	$top=0;
	next;
    } #end if

    if ($last_topo=~/$line[3]/ && $line[0]=~/$last_chr/) { #same topo and same chr
	$close=$line[2];
    } else { #new topo
	print FINAL "$last_chr\t$open\t$close\t$last_topo\n";
	$last_topo=$line[3];
	$last_chr=$line[0];
	$open=$line[1];
	$close=$line[2];
    } #end if

} #end while

print FINAL "$last_chr\t$open\t$close\t$last_topo\n";

