#!/usr/bin/perl
use Perl::utils;

if(@ARGV==0) {
}

parse_command_line(annot    => {description=>'the annotation (GFX) file', ifunreadable=>'annotation not specified'},
                   counts   => {description=>'the input ssj counts master file', ifunreadable=>'input not specified'},
                   mincount => {default=>10,  description=>'the min value of the denominator of the fraction'},
                  );

open FILE, $annot || die();
print STDERR "[<$annot";
while($line=<FILE>) {
    chomp $line;
    ($chr, $trash, $element, $beg, $end, $trash, $str, $frame, $attr) = split /\t/, $line;
    $elem{$element}{$chr}{$str}{$beg}{$end}++;
    $mele{$element}{$chr}{$str}{$end}{$beg}++;
}
close FILE;
print STDERR "]\n";


print STDERR "[<$counts";
open FILE, $counts || die();
@header = split /\t/, <FILE>;
while($line=<FILE>) {
    chomp $line;
    ($chr, $beg, $end, $str, $rest) = split /[\t\_]/, $line, 5;
    $data{$chr}{$str}{$beg}{$end} = [split /\t/, $rest];
}
close FILE;
print STDERR "]\n";


print STDERR "[>";
print join("\t", @header);
foreach $chr(sort keys(%{$elem{exon}})) {
    foreach $str(keys(%{$elem{exon}{$chr}})) {
	foreach $x(sort {$a<=>$b} keys(%{$elem{exon}{$chr}{$str}})) {
	    foreach $y(sort {$a<=>$b} keys(%{$elem{exon}{$chr}{$str}{$x}})) {
		foreach $a(keys(%{$mele{intron}{$chr}{$str}{$x}})) {
		    foreach $b(keys(%{$elem{intron}{$chr}{$str}{$y}})) {
			next unless($elem{intron}{$chr}{$str}{$a}{$b});
			print join("_", $chr, $a, $x, $y, $b, $str);
			for($i=0;$i<@header;$i++) {
			    $inc = $data{$chr}{$str}{$a}{$x}->[$i] + $data{$chr}{$str}{$y}{$b}->[$i];
			    $exc = $data{$chr}{$str}{$a}{$b}->[$i];
                	    print "\t", frac($inc, 2*$exc);
			}
			print "\n";
		    }
		}
	    }
	}
    }
}
print STDERR "]\n";




