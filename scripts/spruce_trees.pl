#! /usr/bin/perl

#reads in geno files and only keeps sites with AT LEAST 1 PER GROUP based on pop file

#$pop_file="../split_by_chr/pop_list_new.txt";
$pop_file=$ARGV[0];

#$geno_file="../split_by_chr/geno_files/Macaque_merged.chrM.geno";
$geno_file=$ARGV[1];

#unzip geno file input
system("gunzip ${geno_file}.gz");

open(POP,$pop_file);
%pops=();

while (<POP>) {
    chomp;
    @line=split(/\s+/,$_);
    $pops{$line[0]} = $line[1];
} #end while

$output="${geno_file}.spruced";
open(OUTPUT,">$output");
open(INPUT,$geno_file);

%cols=();

while (<INPUT>) {
    chomp;
    @input_line=split(/\s+/,$_);
    $skip_site=0;

    if ($_=~/^#/) { #header line
	for ($i=2;$i<=$#input_line;$i++) {
	    $group=$pops{$input_line[$i]};  #get column number for sample and extract group from pops hash
	    push @{$cols{$group}},$i;       #add column number to cols hash as an array for each group
	} #end for
	print OUTPUT "$_\n";  #print header to output file
	next;
    } #end if

    foreach $var (keys %cols) {  #loop through individuals within each group
	$missing=0;  #set missing to zero for this group

	for ($k=0;$k<=$#{$cols{$var}};$k++) {  #loop through the group
	    $column=$cols{$var}[$k];  #extract column number 

	    if ($input_line[$column]=~/N\/N/) {  #check corresponding columns for missingness
		$missing+=1; #increment if missing
	    } #end if

	} #end for

	if ($missing>$#{$cols{$var}}) { #all individuals within ANY group are missing, so skip this line!
	    $skip_site=1;
	} #end if

    } #end foreach

    if ($skip_site==0) {
	print OUTPUT "$_\n";  #print line that passes missing check to output file
    } #end if

} #end while

#zip geno file input
system("gzip ${geno_file}");
