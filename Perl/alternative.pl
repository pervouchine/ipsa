#!/usr/bin/perl
use Perl::utils;

if(@ARGV==0) {
    print STDERR "This utility computes splicing graph\n";
}

parse_command_line(annot   => {description=>'the annotation (gtf)', ifunreadable=>'annotation not specified'});

print STDERR "[<=$annot";
open FILE, $annot || die;
while($line=<FILE>) {
    chomp $line;
    ($chr, $source, $feature, $beg, $end, $score, $str, $frame, $attr) = split /\t/, $line;
    next if($f{$feature}{$chr}{$str}{$beg}{$end});
    $f{$feature}{$chr}{$str}{$beg}{$end}++;
    $g{$feature}{$chr}{$str}{$beg}++;
    $g{$feature}{$chr}{$str}{$end}++;
}
close FILE;
print STDERR "]\n";

foreach $chr(keys(%{$f{intron}})) {
    foreach $str(keys(%{$f{intron}{$chr}})) {
	foreach $beg(keys(%{$f{intron}{$chr}{$str}})) {
	    foreach $end(keys(%{$f{intron}{$chr}{$str}{$beg}})) {
		print join("_",$chr,$beg, $end,$str), "\n" if($g{intron}{$chr}{$str}{$beg}>1 || $g{intron}{$chr}{$str}{$end}>1);
	    } 
	}
    }
}

foreach $chr(keys(%{$f{exon}})) {
    foreach $str(keys(%{$f{exon}{$chr}})) {
	foreach $beg(keys(%{$f{exon}{$chr}{$str}})) {
	    foreach $end(keys(%{$f{exon}{$chr}{$str}{$beg}})) {
		next unless($g{intron}{$chr}{$str}{$beg}>=1 && $g{intron}{$chr}{$str}{$end}>=1);
		print join("_",$chr,$beg, $end,$str), "\n" if($g{intron}{$chr}{$str}{$beg}>1 || $g{intron}{$chr}{$str}{$end}>1);
	    }
	}
    }
}
    





