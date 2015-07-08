#!/usr/bin/perl
use Perl::utils;

if(@ARGV==0) {
}

parse_command_line( ssj => {description=>'the splice junction file', ifunreadable=>'bed not specified'},
                    minstaggered=>{description=>'he minimum umber of staggered reads', default=>2},
                    nucleotides =>{description=>'the splice site nucleotides', default=>GTAG});

open FILE, $ssj || die();
while($line=<FILE>) {
    chomp $line;
    ($id, $total, $staggered, $entropy, $annot, $nuc) = split /\t/, $line;
    $jnc{$id}++ if($nuc eq $nucleotides && $staggered>=$minstaggered);
}
close FILE;

while($line=<STDIN>) {
    chomp $line;
    ($id, $deg, $offset, $count) = split /\t/, $line;
    ($chr, $x, $y, $z, $t, $strand) = split /\_/, $id;
    foreach $str("+", "-") {
        if($strand eq $str || $strand eq '.') {
    	    print "$line\n" if($jnc{join("_", $chr, $x, $y, $str)} && $jnc{join("_", $chr, $z, $t, $str)});
	}
    }
}
