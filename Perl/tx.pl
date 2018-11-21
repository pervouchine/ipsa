#!/usr/bin/perl
use Perl::utils;

if(@ARGV==0) {
    print STDERR "This utility computes exon inclusion rates based on the quantifiaction data\n"; 
}

parse_command_line(annot=>{description=>'the genomic annotation', ifunreadable=>'the genomic annotation is missing'},
		   quant=>{description=>'the transcript quantification file (gtf/gff)', ifunreadable=>'the quantification file missing'},
		   field=>{description=>'comma separated list of fields to be averaged', default=>"RPKM"},
		   minc =>{description=>'the minimum sum of rpkms to consider the inclusion ratio as reliable', default=>0.5});
		
#######################################################################################################################

@field = split /\,/, $field;

print STDERR "[<$annot";
open FILE, $annot || die();
while($line=<FILE>) {
    chomp $line;
    ($chr, $trash, $element, $beg, $end, $trash, $str, $frame, $attr) = split /\t/, $line;
    next unless($element eq "exon");
    $eid = join("_", $chr, $beg, $end, $str);
    %attr = get_attributes($attr);
    foreach $tid(split /\,/, $attr{'transcript_id'}) {
	$TE{$tid}{$eid}++;
    }
    $exon{$chr}{$str}{$beg}{$end}++;
}
close FILE;
print STDERR ", indexing";

foreach $chr(keys(%exon)) {
    foreach $str(keys(%{$exon{$chr}})) {
    	foreach $beg(sort {$a<=>$b} keys(%{$exon{$chr}{$str}})) {
	    foreach $end(sort {$a<=>$b} keys(%{$exon{$chr}{$str}{$beg}})) {
	    	push @{$loe{$chr}{$str}}, [$beg, $end];
	    }
	}
    }
}
print STDERR "]\n";

print STDERR "[<$quant";
open FILE, "sort -k1,1 -k4,4n $quant |" || die();
while($line=<FILE>) {
    chomp $line;
    ($chr, $source, $element, $beg, $end, $name, $str, $trash, $attr) = split /\t/, $line;
    next unless($element eq "transcript");
    %attr = get_attributes($attr);

    $tid = $attr{'transcript_id'};
    $abundance = avg(@attr{@field});

    while($curr{$chr}{$str} < @{$loe{$chr}{$str}} && $loe{$chr}{$str}->[$curr{$chr}{$str}]->[0]<$beg) {$curr{$chr}{$str}++;}
    for($i=$curr{$chr}{$str}; $i<@{$loe{$chr}{$str}} && $loe{$chr}{$str}->[$i]->[0]<=$end; $i++) {
	next unless($loe{$chr}{$str}->[$i]->[1]<=$end);
	$eid = join("_", $chr, @{$loe{$chr}{$str}->[$i]}, $str);
	$tot{$eid}+=$abundance;
    }
    foreach $eid(keys(%{$TE{$tid}})) {
	$inc{$eid} += $abundance;
    }
}
close FILE;
print STDERR "]\n";


print STDERR "[>stdout";
foreach $chr(keys(%exon)) {
    foreach $str(keys(%{$exon{$chr}})) {
	for($i=0; $i<@{$loe{$chr}{$str}}; $i++) {
	    $eid = join("_", $chr, @{$loe{$chr}{$str}->[$i]}, $str);
	    $psi  = frac($inc{$eid}, $tot{$eid} - $inc{$eid});
	    next unless($psi =~ /\d/);
	    print join("\t", $chr, 'TXPIPE', 'exon', @{$loe{$chr}{$str}->[$i]}, int(1000*$psi), $str, '.', set_attributes(psitx=>$psi, inc=>$inc{$eid}, tot=>$tot{$eid})), "\n";
	}	
    }
}
print STDERR "]\n";



