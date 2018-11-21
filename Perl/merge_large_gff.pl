#!/usr/bin/perl
use Perl::utils;

if(@ARGV==0) {
}


parse_command_line(dir    => {description=>'directory name', ifabsent=>'need directory name'},
                   ext    => {description=>'file extention', ifabsent=>'need file extention'},
		   id     => {description=>'sample id field', default=>''},
                   idsep  => {description=>'id field separator', default=>'_'},
		   feature=> {description=>'feature'},
		   subset => {description=>'file to subset id on'},
		   lim	  => {description=>'limit number of files for debug', default=>0},
                   out    => {description=>'endpoint and output tsv file name', array=>hash},
                   percent=> {description=>'do not report rows that have more than <value> NAs', default=>0});

%output = @out;

@id = split /\,/, $id;
unless(@id>0) {
    print STDERR "[WARNING: id field not specified, will use file name instead]\n";
}

if($subset) {
    open FILE, $subset || die("Can't read $subset\n");
    while($line=<FILE>) {
    	chomp $line;
    	$subset{$line}++;
    }
    close FILE;
}

while($line=<STDIN>) {
    chomp $line;
    ($file, $attr) = split /\t/, $line;
    %attr = get_features($attr);

    @array = split /\//, $file;
    ($target) = split /\./, pop(@array);
    $key = join($idsep, @attr{@id});
    $key = $target unless($key);
    die("Key is not unique: $key\n") if($hasbeen{$key}++);

    next unless($key);
    if($attr{'type'} eq "bam" && $attr{'view'}=~/^Alignments/) {
        $file = "$dir$key$ext";
	$input{$file} = $key;
    }
}

foreach $file(keys(%input)) {
    $name = $input{$file};
    print STDERR "[",++$N," $file $name";
    open FILE, $file;
    while($line = <FILE>) {
	($chr, $source, $f, $beg, $end, $score, $str, $frame, $attr) = split /\t/, $line;
	next if($feature && $feature ne $f);
        $id = join("_", $chr, $beg, $end, $str);
	next if($subset && $subset{$id} eq undef);
        %attr = get_attributes($attr);
        foreach $key(keys(%output)) {
            next unless($attr{$key}=~/\d/);
            $data{$key}{$name}{$id} = $attr{$key};
            $rows{$key}{$id}++;
            $cols{$key}{$name}++;
        }
    }
    print STDERR "]\n";
    close FILE; 
    last if($N>$lim && $lim>0);
}

foreach $key(sort keys(%output)) {
    @c = sort keys(%{$cols{$key}});
    @r = sort keys(%{$rows{$key}});
    print STDERR "[$key=>$output{$key}";
    open FILE, ">$output{$key}" || die("can't write to this file");
    print FILE join("\t", @c),"\n";
    foreach $id(@r) {
	@arr = ($id);
	$num = 0;
    	foreach $name(@c) {
	    $value = $data{$key}{$name}{$id};
            push @arr, $value =~/\d/ ? $value : "NA";
            $num++ if($value =~/\d/);
        }
        print FILE join("\t", @arr), "\n" if($num>=$percent*@c);
    }
    close FILE;
    print STDERR "]\n";
}


