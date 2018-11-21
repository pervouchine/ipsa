#!/usr/bin/perl
use Perl::utils;

if(@ARGV==0) {
    print STDERR "This utility computes inclusion and processing rates of SJs (psi_5,3 and cosi_5,3)\n";
}

parse_command_line(ssj      => {description=>'the input ssj (BED) file', ifunreadable=>'input not specified'},
                   ssc      => {description=>'the input ssc (BED) file'},
		   mincount => {default=>10,  description=>'the min value of the denominator of the fraction'},
                   stranded => {default=>1,   description=>'1(yes) or 0(no)'}
                  );

#######################################################################################################################

print STDERR "[<$ssj";
open FILE, $ssj || die("Can't open $ssj\n");
while($line=<FILE>) {
    chomp $line;
    ($chr, $beg, $end, $s, $count) = split /\t/, $line;
    $str = strand_c2i($s) * $stranded;
    ($beg, $end) = reverse ($beg, $end) if($str<0);
    $count53{$chr}{$str}{$beg}{$end}+=$count;
    $count5X{$chr}{$str}{$beg}+=$count;
    $countX3{$chr}{$str}{$end}+=$count;
}
close FILE;
print STDERR "]\n";

if(-e $ssc) {
    print STDERR "[<$ssc";
    open FILE, $ssc || die();
    while($line=<FILE>) {
    	chomp $line;
    	($chr, $pos, $trash, $s, $count) = split /\t/, $line;
	$str = strand_c2i($s) * $stranded;
        $count00{$chr}{$str}{$pos}+=$count;
    }
    print STDERR "]\n";
    close FILE;
}

#######################################################################################################################

print STDERR "[>stdout";
foreach $chr(sort keys(%count53)) {
    foreach $str(sort {$a<=>$b} keys(%{$count53{$chr}})) {
        foreach $beg(sort {$a<=>$b} keys(%{$count53{$chr}{$str}})) {
            foreach $end(sort {$a<=>$b} keys(%{$count53{$chr}{$str}{$beg}})) {
                $nDA = $count53{$chr}{$str}{$beg}{$end};
                $nDX = $count5X{$chr}{$str}{$beg} - $count53{$chr}{$str}{$beg}{$end};
                $nXA = $countX3{$chr}{$str}{$end} - $count53{$chr}{$str}{$beg}{$end};
                $nD = $count00{$chr}{$str}{$beg} + 0;
                $nA = $count00{$chr}{$str}{$end} + 0;
                $psi5  = frac($nDA, $nDX);
                $psi3  = frac($nDA, $nXA);
                $cosi5 = frac($nDA + $nDX, $nD);
                $cosi3 = frac($nDA + $nXA, $nA);
		next unless(join(undef, $psi5,$psi3,$cosi5,$cosi3) =~ /\d/);
		($x, $y) = sort {$a<=>$b} ($beg, $end);
                print join("\t", $chr, 'SJPIPE', 'intron', $x, $y, $psi5=~/\d/ ? int(500*($psi5 + $psi3)) : '.', strand_i2c($str), '.',
                        set_attributes(psi5=>$psi5, psi3=>$psi3, cosi5=>$cosi5, cosi3=>$cosi3, nDA=>$nDA, nDX=>$nDX, nXA=>$nXA, nD=>$nD, nA=>$nA)), "\n";
            }
        }
    }
}
print STDERR "]\n";



