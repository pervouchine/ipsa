#!/usr/bin/perl
use Perl::utils;

if(@ARGV==0) {
    print STDERR "This utility computes (a) the global exon inclusion and processing rates for a given set of annotated exons and (b) inclusion and processing rates of SJs\n";
}

parse_command_line(
                   ssj      => {description=>'the input ssj (BED) file', ifunreadable=>'input not specified'},
		   vip	    => {description=>'splicing graph (BED) file', ifunreadable=>'splicing graph not specified'},
		   mincount => {default=>10,  description=>'the min value of the denominator of the fraction'}
                  );

#######################################################################################################################

# Read ssj table again and gather inclusion and exclusion counts
print STDERR "[<$ssj";
open FILE, $ssj || die("Can't open $ssj\n");
while($line=<FILE>) {
    chomp $line;
    ($id, $count) = split /\t/, $line;
    $data{$id} = $count;
}
close FILE;
print STDERR "]\n";

print STDERR "[<$vip";
open FILE, $vip || die("Can't open $vip\n");
while($line=<FILE>) {
    chomp $line;
    ($chr, $beg, $end, $name, $score, $str, $paths) = split /\t/, $line;
    ($sj1, $sj2) = split /\:/, $paths;
    @a = split /\,/, $sj1;
    for($x=0, $i=0; $i<@a; $i+=2) {
	$x+= $data{join("_", $chr, $a[$i], $a[$i+1], $str)}
    }
    @b = split /\,/, $sj2;
    for($y=0, $i=0; $i<@b; $i+=2) {
        $y+= $data{join("_", $chr, $b[$i], $b[$i+1], $str)}
    }
    $x/=(@a/2);
    $y/=(@b/2);
    print "$chr:$paths:$str\t",frac($x,$y),"\n";
}
close FILE;
print STDERR "]\n";

