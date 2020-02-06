#!/usr/bin/perl
use Perl::utils;

if(@ARGV==0) {
    print STDERR "This utility creates a makefile for the sj pipeline, taking the index file from STDIN and printing the makefile to STDOUT\n";
}

parse_command_line(     in      => {description=>'the index file'},
     			dir     => {description=>'the output directory', ifabsent=>'output directory not specified'},
                        repository => {description=>'the repository subdirectory for bam files'},
                        param   => {description=>'parameters passed to sjcount'},
                        id      => {description=>'sample id field', default=>''},
                        idsep   => {description=>'id field separator', default=>'_'},
                        margin  => {description=>'margin for aggregate', default=>5},
                        entropy => {description=>'entropy lower threshold', default=>1.5},
                        deltaSS => {description=>'distance threshold for splice sites', default=>10},
                        minstatus=>{description=>'annotation status lower threshold', default=>0},
                        mincount=> {description=>'min number of counts for the denominator', default=>10},
                        annot   => {description=>'the annotation (gfx)', ifabsent=>'annotation not specified'},
                        genome  => {description=>'the genome (without .dbx or .idx)', ifabsent=>'genome not specified'},
                        merge   => {description=>'the name of the output to merge in case if blocks are missing', default=>"all"},
			prefix  => {description=>'chromosome name prefix'},
                        SJCOUNTDIR =>{variable=>T,ifabsent=>'SJCOUNTDIR not specified'});

$program = $SJCOUNTDIR."sjcount";

@id = split /\,/, $id;
unless(@id>0) {
    print STDERR "[WARNING: id field not specified, will use file name instead]\n";
}

if(-e -r $in) {
    open $f, $in || die("Can't read $in exiting\n");
} else  {
    $f=*STDIN;
}

#make(script=>"vip1.pl", input=>{'<'=>$annot}, output=>{'>'=>$dir."sg.vip"}, mkdir=>T);

while($line=<$f>) {
    chomp $line;
    ($file, $attr) = split /\t/, $line;
    %attr = get_features($attr);

    @array = split /\//, $file;
    $target = pop(@array);
    $target =~ s/\.\w*$//;
    $key = join($idsep, @attr{@id});
    $key = $target unless($key);
    die("Key is not unique: $key\n") if($hasbeen{$attr{'type'}}{$key}++);

    if(($attr{'type'} eq "bam" || $attr{'file_format'} eq "bam") && ($attr{'view'} eq "Alignments" || $attr{'output_type'} eq "alignments")) {
	if($attr{'readType'}=~/(\d+)(D*)$/) {
	    $readLength = $1;
	    $stranded   = ($2 eq "D" ? "" : "-unstranded");
	}
	else {
	    die("Read length not specified");
	}
	die("Incorrect read length") unless($readLength>0);

        if($file=~/^http/ || $file=~/^ftp/) {
            make(script=>"wget", before=>$file, output=>{-O=>"$repository$target.bam"}, mkdir=>T, endpoint=>'download');
            $file = "$repository$target.bam";
        }

##
	make(script=>"pA.pl", input=>{-bam=>$file}, after=>"-readLength $readLength $stranded", output=>{'>'=>fn($key,P01,tsv)}, endpoint=>P01);
	make(script=>"aggregate.awk", input=>{''=>fn($key,P01,tsv)}, output=>{'>'=>fn($key,P02,tsv)}, before=>"-v degree=0", endpoint=>P02);
##

	$param = "-read1 1 -read2 0" if($attr{'sense'} eq 'MATE2_SENSE');
	$param = "-read1 0 -read2 1" if($attr{'sense'} eq 'MATE1_SENSE');
	make(script=>$program, input=>{-bam=>$file}, output=>{-ssj=>fn($key,A01,ssj,tsv), -ssc=>fn($key,A01,ssc,tsv), -log=>fn($key,A01,ssj,'log')},
             after=>"-nbins $readLength $param $stranded -quiet", mkdir=>T, endpoint=>A01);

	$prm = "-v readLength=$readLength -v margin=$margin -v prefix=$prefix";
	make(script=>"aggregate.awk", input=>{''=>fn($key,A01,ssj,tsv)}, output=>{'>'=>fn($key,A02,ssj,tsv)}, before=>"-v degree=1 $prm", endpoint=>A02);

	make(script=>"aggregate.awk", input=>{''=>fn($key,A01,ssc,tsv)}, output=>{'>'=>fn($key,A02,ssc,tsv)}, before=>"-v degree=0 $prm", endpoint=>A02);

	make(script=>"aggregate.awk", input=>{''=>fn($key,A01,ssj,tsv)}, output=>{'>'=>fn($key,D01,tsv)},     before=>"-v degree=2",      endpoint=>D01);
	make(script=>"annotate.pl", input=>{-in=>fn($key,A02,ssj,tsv), -annot=>$annot, -dbx=>"$genome.dbx", -idx=>"$genome.idx"}, output=>{'>'=>fn($key,A03,ssj,tsv)}, after=>"-deltaSS $deltaSS", endpoint=>A03);
	make(script=>"choose_strand.awk", input=>{''=>fn($key,A03,ssj,tsv)}, output=>{'>'=>fn($key,A04,ssj,tsv)}, endpoint=>A04);
	make(script=>"constrain_ssc.awk", input=>{''=>fn($key,A02,ssc,tsv), jncfile=>fn($key,A04,ssj,tsv)}, output=>{'>'=>fn($key,A04,ssc,tsv)}, endpoint=>A04);
	make(script=>"constrain_mex.awk", input=>{''=>fn($key,D01,tsv), jncfile=>fn($key,A04,ssj,tsv)}, output=>{'>'=>fn($key,D02,tsv)}, endpoint=>D02);
	make(script=>"extract_mex.awk", input=>{''=>fn($key,D02,tsv)}, output=>{'>'=>fn($key,D06,tsv)}, endpoint=>D06);

	make(script=>"awk", before=>"'\$\$4>=$entropy && \$\$5>=$minstatus'", input=>{''=>fn($key,A04,ssj,tsv)}, output=>{'>'=>fn($key,A06,ssj,tsv)}, endpoint=>A06);
        make(script=>"awk", before=>"'\$\$4>=$entropy'", input=>{''=>fn($key,A04,ssc,tsv)}, output=>{'>'=>fn($key,A06,ssc,tsv)}, endpoint=>A06);
	make(script=>"shifts.pl", input=>{-in=>fn($key,A06,ssj,tsv), -annot=>$annot}, output=>{'>'=>fn($key,G06,ssj,tsv)}, endpoint=>shifts);

        make(script=>'tsv2bed.pl', input=>{'<'=>fn($key,A06,ssj,tsv)}, output=>{'>'=>fn($key,E06,ssj,bed)}, between=>"-extra 2,3,4,5,6,7", endpoint=>E06);
        make(script=>'tsv2bed.pl', input=>{'<'=>fn($key,A06,ssc,tsv)}, output=>{'>'=>fn($key,E06,ssc,bed)}, between=>"-extra 2 -ssc", endpoint=>E06);
        make(script=>'tsv2gff.pl', input=>{'<'=>fn($key,A06,ssj,tsv)}, output=>{'>'=>fn($key,E06,ssj,gff)}, between=>"-o count 2 -o stagg 3 -o entr 4 -o annot 5 -o nucl 6 -o IDR 7", endpoint=>E06);

	$prm = "-mincount $mincount";

    	make(script=>"zeta.pl", input=>{-ssj=>fn($key,A06,ssj,tsv), -ssc=>fn($key,A06,ssc,tsv), -annot=>$annot, -exons=>fn($key,D06,tsv)}, output=>{'>'=>fn($key,A07,gff)}, between=>$prm, endpoint=>A07);

	make(script=>"alef.pl", input=>{-ssj=>fn($key,A06,ssj,tsv),-vip=>$dir."sg.vip"}, output=>{'>'=>fn($key,A08,tsv)}, endpoint=>A08);

    	make(script=>"psimod.pl", input=>{-ssj=>fn($key,A06,ssj,tsv), -annot=>$annot}, output=>{'>'=>fn($key,B07,gff)}, endpoint=>B07);

        $merge_tsv{A}{ssj}{fn($key,A06,ssj,tsv)} 	= $key;
        $merge_tsv{A}{ssc}{fn($key,A06,ssc,tsv)} 	= $key;

        $merge_gff{A}{'psi,cosi'}{fn($key,A07,gff)} 	= $key;
        $merge_gff{A}{'inc,exc,ret'}{fn($key,A07,gff)} 	= $key;
        $merge_gff{A}{'psi5,psi3'}{fn($key,A07,gff)} 	= $key;
        $merge_gff{A}{'cosi5,cosi3,cosit'}{fn($key,A07,gff)} 	= $key;
        #$merge_gff{B}{'psicas,psiloc'}{fn($key,B07,gff)} = $key;
    }

    if(($attr{'type'} eq "gff" || $attr{'type'} eq "gtf") && $attr{'view'} =~ /^Transcript/) { 
        make(script=>"tx.pl", input=>{-quant=>$file, -annot=>$annot}, output=>{'>'=>fn($key,C07,gff)}, between=>"-field RPKM1,RPKM2", endpoint=>C07);
	$merge_gff{C}{psitx}{fn($key,C07,gff)} = $key;
    }
}

####################################################################################################################################################

foreach $endpoint(keys(%merge_tsv)) {
    foreach $arm(keys(%{$merge_tsv{$endpoint}})) {
	make2(script=>"merge_tsv.pl", inputs=>{-i=>\%{$merge_tsv{$endpoint}{$arm}}}, outputs=>{''=>{'>'=>fn($merge,counts,$arm,tsv)}}, endpoint=>$endpoint);
    }
}

foreach $endpoint(keys(%merge_gff)) {
    foreach $arms(keys(%{$merge_gff{$endpoint}})) {
	%outputs=();
	foreach $arm(split /\,/, $arms) {
	    $outputs{$arm} = fn($merge,$endpoint,$arm,tsv);
	}
        make2(script=>"merge_gff.pl", inputs=>{-i=>\%{$merge_gff{$endpoint}{$arms}}}, outputs=>{-o=>\%outputs}, endpoint=>$endpoint);
    }
}

####################################################################################################################################################

print "all :: A \n";

sub fn {
    return(@_[1]=~/^[A-Z]\d+$/ ? join(undef, $dir, @_[1], "/", join('.', @_)) : join(undef, $dir, join('.', @_)));
}
