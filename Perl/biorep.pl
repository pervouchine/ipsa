#!/usr/bin/perl
use Perl::utils;

$rep = $ARGV[0];

open FILE, $rep || die;
while($line=<FILE>) {
    chomp $line;
    @array = split /\t/, $line;
    $label = join("_", @array[0..1]);
    foreach $elem(@array) {
	$arm{$elem} = $label;
    }
    $arm{join(",", @array[0..1])} = $label;
}
close FILE;

while($line=<STDIN>) {
    chomp $line;
    ($file, $attr) = split /\t/, $line;
    %attr = get_features($attr);
    next unless($arm{$attr{'labExpId'}});
    $line.= " idrGroup=$arm{$attr{'labExpId'}};";
    print "$line\n";
}





