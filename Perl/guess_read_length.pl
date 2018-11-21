use Perl::utils;

while($line=<STDIN>) {
    chomp $line;
    ($file, $attr) = split /\t/, $line;
    %attr = get_features($attr);
    $attr{'readType'} = $ARGV[0] if($ARGV[0] =~/\w/);

    if($attr{'readType'}=~/(\d+)(D*)$/) {
	#print join("\t", $file, $attr), "\n";
	next;
    }
    if($attr{'readType'}=~/(2x)*NA(D*)$/) {
	$PE = $1;
	$ST = $2;

	$s = `sjcount-3.1/samtools-0.1.18/samtools view $file | awk '\$6~/^[0-9]+M\$/' | head -n 100 | cut -f6`;
    	%f=();
    	foreach $x(split /\n/, $s) {
	    $f{$x}++;
    	}
    	$z = (sort {$f{$b}<=>$f{$a}} keys(%f))[0];
    	substr($z,-1,1)=undef;
    	next unless($z>0);
    	$attr{'readType'} = "$PE$z$ST";
print STDERR $attr{'readType'};

	$attr = set_features(%attr);
    	print join("\t", $file, $attr), "\n";
	next; 
    }
}
