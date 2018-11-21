#!/usr/bin/perl
use lib '/home/dp/ipsa';
use Perl::utils;

if(@ARGV==0) {
    print STDERR "This utility takes an aggregated TSV file, the genomic annotation, and the genome, and outputs a TSV with two more columns: ";
    print STDERR "(8) annotation status and (9) splice sites\n";
    print STDERR "If the input was strandless ('.' in column 4) then each line will be reported twice, one for '+' and and for '-' strand\n";
}

parse_command_line(	in	=> {description=>'the input A06 junction', ifunreadable=>'input not specified'},
			annot	=> {description=>'the annotation (gfx)', ifunreadable=>'annotation not specified'},
			w	=> {description=>'window', default=>100});

print STDERR "[$annot";
open FILE, $annot;
while($line=<FILE>) {
    chomp $line;
    ($chr, $source, $feature, $beg, $end, $score, $str, $frame, $group) = split /\t/, $line;
    next unless($feature eq "intron");
    ($beg, $end) = reverse ($beg, $end) if($str eq "-");
    $f{"$chr;$str;D"}{$beg}{1}++;
    $f{"$chr;$str;A"}{$end}{1}++;
}
close FILE;
print STDERR "]\n";

print STDERR "[$in";
open FILE, $in;
while($line=<FILE>) {
    chomp $line;
    ($jnc, $count, $stag, $ent, $ann, $nuc) = split /\t/, $line;
    next unless($nuc eq "GTAG");
    ($chr, $beg, $end, $str) = split /\_/, $jnc;
    ($beg, $end) = reverse ($beg, $end) if($str eq "-");
    $f{"$chr;$str;D"}{$beg}{0}+=$count;
    $f{"$chr;$str;A"}{$end}{0}+=$count;
}
close FILE;
print STDERR "]\n";

foreach $id(keys(%f)) {
    ($chr, $str, $type) = split /;/, $id;
    @coord = sort {$a<=>$b} keys(%{$f{$id}});
    for($i=0;$i<@coord;$i++) {
        $pos = $coord[$i];
        if($f{$id}{$pos}{0}) {
            @b = ();
            for($j=$i;$j>=0;$j--) {
		last if(abs($coord[$j]-$pos)>$w);
                if($f{$id}{$coord[$j]}{1}) {
                    push @b, $coord[$j]-$pos;
                }
            }
            for($j=$i;$j<@coord;$j++) {
		last if(abs($coord[$j]-$pos)>$w);
                if($f{$id}{$coord[$j]}{1}) {
		    push @b, $coord[$j]-$pos;
                }
            }
            @c = sort {abs($a)<=>abs($b)} @b;
            $shift = $c[0];
	    next if(@c==0 || $shift==0);
            print join("\t", $chr, $pos, $str, $type, $shift, 0+$f{$id}{$pos}{0}, 0+$f{$id}{$pos+$shift}{0}), "\n";
        }
    }
}


