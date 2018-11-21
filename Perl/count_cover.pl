#!/usr/bin/perl
use Perl::utils;

#Read gtf exons (or CDS) and save (exon,transcript) relation
#Also save splice sites for future 
while($line=<STDIN>) {
    chomp $line;
    ($chr, $source, $feature, $beg, $end, $score, $str, $frame, $attr) = split /\t/, $line;
    %attr = get_attributes($attr);
    next unless($id = $attr{'transcript_id'});
    push @{$exons{$id}}, [$beg, $end];
    $chr{$id}{$chr} = $str{$id}{$str} = 1;
    $site{$chr}{$str}{$beg} = $site{$chr}{$str}{$end} = 0;
    $freq{$chr}{$str}{$beg}{$end}++;
    $m++;
}

# indexing enumerates all splice sites in the increasing order on each chr-str
# index array is the inverse function
print STDERR "[Indexing ";
foreach $chr(sort keys(%site)) {
    foreach $str(sort keys(%{$site{$chr}})) {
	@pos = sort {$a<=>$b} keys(%{$site{$chr}{$str}});
	foreach $pos(@pos) {
	    $site{$chr}{$str}{$pos} = ++$n;
	    $index[$n] = $pos;
	}
    }
}
print STDERR "$n]\n";

#for each transcript which is on one strand and chromosome (no trans-splicing)
#run over exons that belong to the transcript and increase corresponding count
#run over exons that are covered by the given transcript and increase corresponding count
foreach $id(keys(%exons)) {
    if(keys(%{$chr{$id}})==1 && keys(%{$str{$id}})==1) {
	($chr, $str) = (keys(%{$chr{$id}}), keys(%{$str{$id}}));
	@array = sort {$a->[0]<=>$b->[0]} @{$exons{$id}};
	for($j=$site{$chr}{$str}{$array[0]->[0]};$j<$site{$chr}{$str}{$array[-1]->[1]};$j++) {
	    foreach $end(keys(%{$freq{$chr}{$str}{$index[$j]}})) {
		$cover{$chr}{$str}{$index[$j]}{$end}++;
	    }
	}
    }
    else {
	$trans++;
    }
}
print STDERR "[WARNING: $trans genes trans-spliced]\n" if($trans>0);

foreach $chr(keys(%cover)) {
    foreach $str(keys(%{$cover{$chr}})) {
	foreach $beg (keys(%{$cover{$chr}{$str}})) {
	    foreach $end (keys(%{$cover{$chr}{$str}{$beg}})) {
    		print join("\t", $chr, $beg, $end, $str, 0 + $freq{$chr}{$str}{$beg}{$end}, 0 + $cover{$chr}{$str}{$beg}{$end}), "\n";
	    }
	}
    }
}
