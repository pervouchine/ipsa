#!/usr/bin/perl
use Perl::utils;

@suffixes = ('ssj','ssc');

if(@ARGV==0) {
    print STDERR "This utility creates a makefile for the sj pipeline, taking the index file from STDIN and printing the makefile to STDOUT\n";
}

parse_command_line(     in      => {description=>'the index file', ifabsent=>'index file not specified'},
                        dir     => {description=>'the output directory', ifabsent=>'output directory not specified'},
                        repository => {description=>'the repository subdirectory for bam files'},
                        param   => {description=>'parameters passed to sjcount'},
                        id   	=> {description=>'sample id field', default=>''},
                        idsep   => {description=>'id field separator', default=>'_'},
                        margin  => {description=>'margin for aggregate', default=>5},
                        entropy => {description=>'entropy lower threshold', default=>1.5},
                        deltaSS => {description=>'distance threshold for splice sites', default=>10},
                        minstatus=>{description=>'annotation status lower threshold', default=>0},
                        mincount=> {description=>'min number of counts for the denominator', default=>10},
                        annot   => {description=>'the annotation (gfx)', ifabsent=>'annotation not specified'},
                        genome  => {description=>'the genome (without .dbx or .idx)', ifabsent=>'genome not specified'},
                        merge   => {description=>'the name of the output to merge in case if blocks are missing', default=>"all"},
			other	=> {description=>'subset of exons and introns'},
			prefix  => {description=>'chromosome name prefix'},
                        SJCOUNTDIR =>{variable=>T,ifabsent=>'SJCOUNTDIR not specified'}
		  );

$program = $SJCOUNTDIR."sjcount";

@id = split /\,/, $id;
unless(@id>0) {
    print STDERR "[WARNING: id field not specified, will use file name instead]\n";
}

if($in) {
    open $f, $in || die("Can't read $in exiting\n");
} else  {
    $f=*STDIN;
}

while($line=<$f>) {
    chomp $line;
    ($file, $attr) = split /\t/, $line;
    %attr = get_features($attr);

    @array = split /\//, $file;
    ($target) = split /\./, pop(@array);
    $key = join($idsep, @attr{@id});
    $key = $target unless($key);
    die("Key is not unique: $key\n") if($hasbeen{$key}++);

    if($attr{'type'} eq "bam" && $attr{'view'}=~/^Alignments/) {
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

        $param = "-read1 1 -read2 0" if($attr{'sense'} eq 'MATE2_SENSE');
        $param = "-read1 0 -read2 1" if($attr{'sense'} eq 'MATE1_SENSE');

	make(script=>$program, input=>{-bam=>$file}, output=>{-ssj=>fn($key,A01,ssj,tsv), -ssc=>fn($key,A01,ssc,tsv), -log=>fn($key,A01,ssj,'log')},
             after=>"-nbins $readLength $param $stranded -quiet", mkdir=>T, endpoint=>A01);

        $prm = "-v readLength=$readLength -v margin=$margin -v prefix=$prefix";
        make(script=>"aggregate.awk", input=>{''=>fn($key,A01,ssj,tsv)}, output=>{'>'=>fn($key,A02,ssj,tsv), logfile=>fn($key,A02,ssj,'log')}, before=>"-v degree=1 $prm", endpoint=>A02);
        make(script=>"aggregate.awk", input=>{''=>fn($key,A01,ssc,tsv)}, output=>{'>'=>fn($key,A02,ssc,tsv), logfile=>fn($key,A02,ssc,'log')}, before=>"-v degree=0 $prm", endpoint=>A02);
        make(script=>"aggregate.awk", input=>{''=>fn($key,A01,ssj,tsv)}, output=>{'>'=>fn($key,D01,tsv)}, before=>"-v degree=2", endpoint=>D01);
        make(script=>"annotate.pl", input=>{-in=>fn($key,A02,ssj,tsv), -annot=>$annot, -dbx=>"$genome.dbx", -idx=>"$genome.idx"}, output=>{'>'=>fn($key,A03,ssj,tsv)}, after=>"-deltaSS $deltaSS", endpoint=>A03);
        make(script=>"choose_strand.awk", input=>{''=>fn($key,A03,ssj,tsv)}, output=>{'>'=>fn($key,A04,ssj,tsv)}, endpoint=>A04);
        make(script=>"constrain_ssc.awk", input=>{''=>fn($key,A02,ssc,tsv), jncfile=>fn($key,A04,ssj,tsv)}, output=>{'>'=>fn($key,A04,ssc,tsv)}, endpoint=>A04);
        make(script=>"constrain_mex.awk", input=>{''=>fn($key,D01,tsv), jncfile=>fn($key,A04,ssj,tsv)}, output=>{'>'=>fn($key,D02,tsv)}, endpoint=>D02);
        make(script=>"extract_mex.awk", input=>{''=>fn($key,D02,tsv)}, output=>{'>'=>fn($key,D06,tsv)}, endpoint=>D06);

        make(script=>"awk", before=>"'\$\$4>=$entropy && \$\$5>=$minstatus'", input=>{''=>fn($key,A04,ssj,tsv)}, output=>{'>'=>fn($key,A06,ssj,tsv)}, endpoint=>A06);
        make(script=>"awk", before=>"'\$\$4>=$entropy'", input=>{''=>fn($key,A04,ssc,tsv)}, output=>{'>'=>fn($key,A06,ssc,tsv)}, endpoint=>A06);

        make(script=>'tsv2bed.pl', input=>{'<'=>fn($key,A06,ssj,tsv)}, output=>{'>'=>fn($key,E06,ssj,bed)}, between=>"-extra 2,3,4,5,6,7", endpoint=>E06);
        make(script=>'tsv2bed.pl', input=>{'<'=>fn($key,A06,ssc,tsv)}, output=>{'>'=>fn($key,E06,ssc,bed)}, between=>"-extra 2 -ssc", endpoint=>E06);
        make(script=>'tsv2gff.pl', input=>{'<'=>fn($key,A06,ssj,tsv)}, output=>{'>'=>fn($key,E06,ssj,gff)}, between=>"-o count 2 -o stagg 3 -o entr 4 -o annot 5 -o nucl 6 -o IDR 7", endpoint=>E06);

        $prm = "-mincount $mincount";

        make(script=>"zeta.pl", input=>{-ssj=>fn($key,A06,ssj,tsv), -ssc=>fn($key,A06,ssc,tsv), -annot=>$annot, -exons=>fn($key,D06,tsv)}, output=>{'>'=>fn($key,A07,gff)}, between=>$prm, endpoint=>A07);
        #make(script=>"psicas.pl", input=>{-ssj=>fn($key,A06,ssj,tsv), -annot=>$annot}, output=>{'>'=>fn($key,B07,gff)}, endpoint=>B07);
    }

    if($attr{'type'} eq "gff" || $attr{'type'} eq "gtf") { 
        make(script=>"tx.pl", input=>{-quant=>$file, -annot=>$annot}, output=>{'>'=>fn($grp,C07,gff)}, endpoint=>C07);
    }
}
close FILE;


#######################################################################################################################################################################

print "all :: A07\n";

$z = $id=~/\w/ ? "-id $id" : undef;

foreach $x([psi, cosi, A07, exon], [inc, exc, A07, exon], [psi5, psi3, A07, intron], [cosi5, cosi3, A07, intron]) {
    ($one, $two, $branch, $f) = @{$x};
    make2(script=>'merge_large_gff.pl', inputs=>{'<'=>{''=>$in}}, before=>"-ext .$branch.gff $z -dir $dir$branch/ -f $f ". ($f eq "exon" ? $other : "-percent 0.25"), 
	 outputs=>{-out=>{$one=>fn($merge, $one, $branch, tsv), $two=>fn($merge, $two, $branch, tsv)}}, endpoint=>master);
}

make(script=>'merge_large_tsv.pl', input=>{'<'=>$in}, before=>"-ext .A06.ssj.tsv $z -dir $dir"."A06/ ", output=>{'>'=>fn($merge,counts,ssj,tsv)}, endpoint=>master);
make(script=>'merge_large_tsv.pl', input=>{'<'=>$in}, before=>"-ext .A06.ssc.tsv $z -dir $dir"."A06/ ", output=>{'>'=>fn($merge,counts,ssc,tsv)}, endpoint=>master);

sub fn {
    return(@_[1]=~/^[A-Z]\d+$/ ? join(undef, $dir, @_[1], "/", join('.', @_)) : join(undef, $dir, join('.', @_)));
    
}


