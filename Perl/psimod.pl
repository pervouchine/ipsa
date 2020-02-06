#!/usr/bin/perl
use Perl::utils;


if(@ARGV==0) {
    print STDERR "This utility computes (a) the global exon inclusion and processing rates for a given set of annotated exons and (b) inclusion and processing rates of SJs\n";
}

parse_command_line(annot    => {description=>'the annotation (GTF) file', ifunreadable=>'annotation not specified'},
                   ssj      => {description=>'the input ssj  (tsv) file', ifunreadable=>'input not specified'},
                   mincount => {default=>10,  description=>'the min value of the denominator of the fraction'},
                  );

open FILE, $annot || die();
print STDERR "[<$annot";
while($line=<FILE>) {
    chomp $line;
    ($chr, $trash, $element, $beg, $end, $trash, $str, $frame, $attr) = split /\t/, $line;
    $f{$element}{$chr}{$str}{$beg}{$end}++; # forward
    $b{$element}{$chr}{$str}{$end}{$beg}++; # backward
}
close FILE;
print STDERR "]\n";

print STDERR "[<$ssj";
open FILE, $ssj || die();
while($line=<FILE>) {
    chomp $line;
    ($chr, $beg, $end, $str, $count) = split /[\t\_]/, $line;
    $data{$chr}{$str}{$beg}{$end} = $count;
    $site{$chr}{$str}{$beg} += $count;
    $site{$chr}{$str}{$end} += $count;
    $f{intron}{$chr}{$str}{$beg}{$end}++;
    $b{intron}{$chr}{$str}{$end}{$beg}++;
}
close FILE;
print STDERR "]\n";

foreach $chr(sort keys(%{$f{exon}})) {
    foreach $str(sort keys(%{$f{exon}{$chr}})) {
	foreach $x(sort{$a<=>$b} keys(%{$f{exon}{$chr}{$str}})) {
	    foreach $y(sort{$a<=>$b} keys(%{$f{exon}{$chr}{$str}{$x}})) {
		$inc{$chr}{$str}{$x}{$y} = $site{$chr}{$str}{$x} + $site{$chr}{$str}{$y};
		foreach $a(sort{$a<=>$b} keys(%{$b{intron}{$chr}{$str}{$x}})) {
		    foreach $b(sort {$a<=>$b} keys(%{$f{intron}{$chr}{$str}{$y}})) {
			$exc{$chr}{$str}{$x}{$y} += $data{$chr}{$str}{$a}{$b};
		    }
		}
	    }
	}
    }
}

foreach $chr(sort keys(%{$f{exon}})) {
    foreach $str(sort keys(%{$f{exon}{$chr}})) {
        foreach $x(sort{$a<=>$b} keys(%{$f{exon}{$chr}{$str}})) {
            foreach $y(sort{$a<=>$b} keys(%{$f{exon}{$chr}{$str}{$x}})) {
		$inc = 0+$inc{$chr}{$str}{$x}{$y};
		$exc = 0+$exc{$chr}{$str}{$x}{$y};
		$psi  = frac($inc, 2*$exc);
		next unless($psi =~ /\d/);
        	print join("\t", $chr, 'SJPIPE', 'exon', $x, $y, int(1000*$psi), $str, '.', set_attributes(psi=>$psi,inc=>$inc,exc=>$exc)), "\n";
	    }
	}
    }
}



