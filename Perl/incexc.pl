#!/usr/bin/perl
use Perl::utils;

if(@ARGV==0) {
    print STDERR "This utility computes (a) the global exon inclusion and processing rates for a given set of annotated exons and (b) inclusion and processing rates of SJs\n";
}

parse_command_line(annot    => {description=>'the annotation (gfx) file', ifunreadable=>'annotation not specified'},
		   exons    => {description=>'exons file with exons to quantify'},
                   ssj      => {description=>'the input ssj (BED) file', ifunreadable=>'input not specified'},
                   ssc      => {description=>'the input ssc (BED) file'},
		   mincount => {default=>10,  description=>'the min value of the denominator of the fraction'},
                  );

#######################################################################################################################

# Read annotation and splice junctions and memorize them

# Annotation
print STDERR "[<$annot";
open FILE, $annot || die();
while($line=<FILE>) {
    chomp $line;
    ($chr, $source, $feature, $beg, $end, $score, $str) = split /\t/, $line;
    ($beg, $end) = reverse ($beg, $end) if(strand_c2i($str)<0);
    $exon{$chr}{$str}{$beg}{$end}{$source}++ if($feature eq "exon");
    $count53{$chr}{$str}{$beg}{$end} = $count35{$chr}{$str}{$end}{$beg} = 0 if($feature eq "intron");
}
close FILE;
print STDERR "]\n";

# Junctions
print STDERR "[<$ssj";
open FILE, $ssj || die();
while($line=<FILE>) {
    chomp $line;
    ($id, $count) = split /\t/, $line;
    ($chr, $beg, $end, $str) = split /\_/, $id;
    ($beg, $end) = reverse ($beg, $end) if(strand_c2i($str)<0);
    $count53{$chr}{$str}{$beg}{$end} = $count35{$chr}{$str}{$end}{$beg} = $count;
    $count5X{$chr}{$str}{$beg} += $count;
    $countX3{$chr}{$str}{$end} += $count;
}
close FILE;
print STDERR "]\n";

#Mini exons
if(-r $exons) {
    print STDERR "[<$exons";
    open FILE, $exons || die();
    while($line=<FILE>) {
        chomp $line;
        ($id) = split /\t/, $line;
	($chr, $beg, $end, $str) = split /\_/, $id;
	($beg, $end) = reverse ($beg, $end) if(strand_c2i($str)<0);
	$exon{$chr}{$str}{$beg}{$end}{MINIEX}++;
    }
    close FILE;
    print STDERR "]\n";
}

#Splice sites
if(-r $ssc) {
    print STDERR "[<$ssc";
    open FILE, $ssc || die();
    while($line=<FILE>) {
        chomp $line;
        ($id, $count) = split /\t/, $line;
        ($chr, $pos, $str) = split /\_/, $id;
        $count00{$chr}{$str}{$pos} = $count;
    }
    print STDERR "]\n";
    close FILE;
}

#######################################################################################################################

print STDERR "[>exons";
foreach $chr(keys(%exon)) {
    foreach $str(keys(%{$exon{$chr}})) {
	foreach $beg(sort {$a<=>$b} keys(%{$exon{$chr}{$str}})) {
	    foreach $end(sort {$a<=>$b} keys(%{$exon{$chr}{$str}{$beg}})) {
    		$inc = $exc = 0;
    		foreach $y(keys(%{$count53{$chr}{$str}{$end}})) {
		    foreach $x(keys(%{$count35{$chr}{$str}{$beg}})) {
	    		if($count53{$chr}{$str}{$x}{$y} ne undef) {
	    		    $inc += $count35{$chr}{$str}{$beg}{$x} + $count53{$chr}{$str}{$end}{$y};
			    $exc += $count53{$chr}{$str}{$x}{$y};
	    		}
		    }
    		}
    		$ret = $count00{$chr}{$str}{$beg} + $count00{$chr}{$str}{$end} + 0;
    		$psi  = frac($inc, 2*$exc);
    		$cosi = frac($inc + 2*$exc, $ret);
    		next unless(join(undef, $psi, $cosi) =~ /\d/);
		($x, $y) = sort {$a<=>$b} ($beg, $end);
    		print join("\t", $chr, join(",",keys(%{$exon{$chr}{$str}{$beg}{$end}})), 'exon', $x, $y, int(1000*$psi), $str, '.', set_attributes(psi=>$psi, cosi=>$cosi, inc=>$inc, exc=>$exc, ret=>$ret)), "\n";
	    }
	}
    }
}

print STDERR ",introns";
foreach $chr(keys(%count53)) {
    foreach $str(keys(%{$count53{$chr}})) {
        foreach $beg(keys(%{$count53{$chr}{$str}})) {
            foreach $end(keys(%{$count53{$chr}{$str}{$beg}})) {
                $nDA = $count53{$chr}{$str}{$beg}{$end};
                $nDX = $count5X{$chr}{$str}{$beg} - $nDA;
                $nXA = $countX3{$chr}{$str}{$end} - $nDA;
                $nD  = $count00{$chr}{$str}{$beg} + 0;
                $nA  = $count00{$chr}{$str}{$end} + 0;
                $psi5  = frac($nDA, $nDX);
                $psi3  = frac($nDA, $nXA);
                $cosi5 = frac($nDA + $nDX, $nD);
                $cosi3 = frac($nDA + $nXA, $nA);
		next unless(join(undef, $psi5,$psi3,$cosi5,$cosi3) =~ /\d/);
		($x, $y) = sort {$a<=>$b} ($beg, $end);
                print join("\t", $chr, 'SJPIPE', 'intron', $x, $y, $psi5=~/\d/ ? int(500*($psi5 + $psi3)) : '.', $str, '.',
                        set_attributes(psi5=>$psi5, psi3=>$psi3, cosi5=>$cosi5, cosi3=>$cosi3, nDA=>$nDA, nDX=>$nDX, nXA=>$nXA, nD=>$nD, nA=>$nA)), "\n";
            }
        }
    }
}
print STDERR "]\n";



