#!/usr/bin/perl
use lib '/home/dp/ipsa';
use Perl::utils;

if(@ARGV==0) {
    print STDERR "This utility computes (a) the global exon inclusion and processing rates for a given set of annotated exons and (b) inclusion and processing rates of SJs\n";
}

parse_command_line(ssj      => {description=>'the input ssj (BED) file', ifunreadable=>'input not specified'},
                   ssc      => {description=>'the input ssc (BED) file'},
		   mincount => {default=>10,  description=>'the min value of the denominator of the fraction'}
                  );

#######################################################################################################################

# Read SJs from ssj bed file and memorize splice sites
print STDERR "[<$ssj";
open FILE, $ssj || die();
while($line=<FILE>) {
    chomp $line;
    ($id, $count) = split /\t/, $line;
    ($chr, $beg, $end, $str) = split /\_/, $id;
    ($D, $A) = $str eq "+" ? (D, A) : (A, D);
    $site{$chr}{$str}{$beg}{$D} = 1;
    $site{$chr}{$str}{$end}{$A} = 1;
}
close FILE;
print STDERR ", indexing";
foreach $chr(sort keys(%site)) {
    foreach $str(sort keys(%{$site{$chr}})) {
	foreach $pos(sort {$a<=>$b} keys(%{$site{$chr}{$str}})) {
	    $index{$chr}{$str}{$pos} = ++$n;
	    @{$xendi[$n]}= ($chr, $str, $pos);
	}
    }
}
print STDERR ", n=$n]\n";

#######################################################################################################################

# Read ssj table again and gather inclusion and exclusion counts
print STDERR "[<$ssj";
open FILE, $ssj || die("Can't open $ssj\n");
while($line=<FILE>) {
    chomp $line;
    ($id, $count) = split /\t/, $line;
    ($chr, $beg, $end, $str) = split /\_/, $id;
    $start = $index{$chr}{$str}{$beg};
    $stop  = $index{$chr}{$str}{$end};
    next unless($start && $stop && $start < $stop);
    ($D, $A) = $str eq "+" ? (D, A) : (A, D);
    $inc{$chr}{$str}{$beg}{$D}+=$count;
    $inc{$chr}{$str}{$end}{$A}+=$count;
    for($k = $start + 1; $k<$stop; $k++) {
	($chr, $str, $pos) = @{$xendi[$k]}; 
	$exc{$chr}{$str}{$pos}+=$count;
    }
}
close FILE;
print STDERR "]\n";

if(-e $ssc) {
    print STDERR "[<$ssc";
    open FILE, $ssc || die();
    while($line=<FILE>) {
    	chomp $line;
	($id, $count) = split /\t/, $line;
    	($chr, $pos, $str) = split /\_/, $id;
	$ret{$chr}{$str}{$pos}+=$count;
    }
    print STDERR "]\n";
    close FILE;
}

#######################################################################################################################

print STDERR "[>stdout";
for($i=1; $i<=$n; $i++) {
    ($chr, $str, $pos) = @{$xendi[$i]};
    foreach $type(keys(%{$inc{$chr}{$str}{$pos}})) {
        print join("\t", join("_", $chr, $pos, $str, $type), 0+$inc{$chr}{$str}{$pos}{$type}, 0+$exc{$chr}{$str}{$pos}, 0+$ret{$chr}{$str}{$pos}), "\n";
    }
}

