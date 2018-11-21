use Perl::utils;

if(@ARGV==0) {
    print STDERR "This utility computes inclusion and processing rates of SJs (psi_5,3 and cosi_5,3)\n";
}

parse_command_line(ssj      => {description=>'the input ssj (BED) file', ifunreadable=>'input not specified'},
                   ssc      => {description=>'the input ssc (BED) file'},
		   psi5	    => {description=>'psi5 output', ifabsent=>"output 5 not specified"},
		   psi3     => {description=>'psi3 output', ifabsent=>"output 3 not specified"},
                   mincount => {default=>10,  description=>'the min value of the denominator of the fraction'}
                  );

#######################################################################################################################

print STDERR "[<$ssj";
open FILE, $ssj || die("Can't open $ssj\n");
$line=<FILE>;
while($line=<FILE>) {
    chomp $line;
    @array = split /\t/, $line;
    ($chr, $beg, $end, $str) = split /\_/, (shift @array);
    #($beg, $end) = reverse ($beg, $end) if($str eq "-");
    @{$data{$id}} = @array;
next;
    for($i=0;$i<@array;$i++) {
	$DA{join("_",$chr, $beg, $end, $str)}[$i]+=$array[$i];
	$DX{join("_",$chr, $beg, $str)}[$i]+=$array[$i];
	$XA{join("_",$chr, $end, $str)}[$i]+=$array[$i];
    }
    print STDERR ".";
}
print STDERR "]";

exit(0);

open PSI5, ">$psi5";
open PSI3, ">$psi3";
foreach $id(sort keys(%DA)) {
    ($chr, $beg, $end, $str) = split /\_/, $id;
    ($beg1, $end1) = $str eq "+" ? ($beg, $end) : ($end, $beg);
    print PSI5 join("_", $chr, $beg1, $end1, $str);
    print PSI3 join("_", $chr, $beg1, $end1, $str);
    for($i=0;$i<@array;$i++) {
	print PSI5 "\t", frac($DA{join("_",$chr, $beg, $end, $str)}[$i], $DX{join("_",$chr, $beg, $str)}[$i]-$DA{join("_",$chr, $beg, $end, $str)}[$i]);
	print PSI3 "\t", frac($DA{join("_",$chr, $beg, $end, $str)}[$i], $XA{join("_",$chr, $end, $str)}[$i]-$DA{join("_",$chr, $beg, $end, $str)}[$i]);
    }
    print PSI5 "\n";
    print PSI3 "\n";
}
close PSI5;
close PSI3;
