#!/usr/bin/perl
use Perl::utils;

if(@ARGV==0) {
    print STDERR "This is a line-based utility for extracting pA sites from SAM files\n";
}

parse_command_line( read1 => {description=>'flip read1 yes/no (1/0)', default=>1},
		    read2 => {description=>'flip read2 yes/no (1/0)', default=>0},
		    readLength => {description=>'read length', ifabsent=>'not specified'},
		    minmatch => {description=>'min number of nucleoties for match part', default=>31},
		    mintail  => {description=>'min number of nucleoties for polyA part', default=>4},
		    minpercenta => {description=>'min percent A in polyA tail', default=>80},
		    unstranded => {description=>'unstranded flag',store=>T},
		    bam   => {description=>'input bam file', ifabsent=>'not specified'},
		    lim   => {description=>'stop after this number of lines (for debug)',default=>0});
 
$BAM_FREAD1 = 0x40;
$BAM_FREAD2 = 0x80;
$BAM_FREVERSE = 0x10;
@STRAND = ("+", "-");

@read = ($read1, $read2);
for($s=0; $s<2; $s++) {
    print STDERR "[Warning: will take reverse complement of read ", $s+1, "]\n" if($read[$s] % 2);
}

open FILE, "sjcount-3.1/samtools-0.1.18/samtools view $bam |" || die();
while(<FILE>){
    ($id, $flag, $ref, $pos, $qual, $cigar, $z, $z, $z, $seq) = split /\t/;
    $s = (($flag & $BAM_FREVERSE)>0);
    $strand = ($flag & $BAM_FREAD1) ? ($s + $read[0]) & 1 : ($s + $read[1]) & 1;

    $n++;
    last if($n>$lim && $lim>0);

    if($cigar=~/^(\d+)M(\d+)S$/) {
	$x = $1;
	$y = $2;
	if($x >= $minmatch && $y >= $mintail) {
	    if(perc(substr($seq,-$y,$y),'A') >= $minpercenta) {
		print "3p\t$strand\n";
	    }
	}
    }
    if($cigar=~/^(\d+)S(\d+)M$/) {
	$x = $2;
	$y = $1;
	if($x >= $minmatch && $y >= $mintail) {
            if(perc(substr($seq,0,$y),'T') >= $minpercenta) {
		print "5p\t$strand\n";
	    }
	}
    }   
}
close FILE;


foreach $key(sort keys(%count)) {
    print "$key\t$count{$key}\n";
}


sub perc {
    my $t = @_[0];
    $t =~ s/@_[1]//g;
    return(int(100*(1-length($t)/length(@_[0]))));
}

