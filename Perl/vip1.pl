#!/usr/bin/perl

# Reading transcript structure from GFX
while($line=<STDIN>) {
    chomp $line;
    ($chr, $source, $feature, $bed, $end, $name, $str, $frame, $attr)  = split /\t/, $line;
    next unless($feature eq "exon");
    if($attr=~/transcript_id \"(.*)\"/) {
	foreach $id(split /\,/, $1) {
	    push @{$exons{$id}}, ($bed, $end);
	    $loc{$id}{"$chr,$str"}++;	# this maps id->chr,str
	    $csi{"$chr,$str"}{$id}++;   # this maps chr,str->list of ids
	}
    }
}

foreach $chrstr(sort keys(%csi)) {
    ($chr,$str) = split /\,/, $chrstr;
    print STDERR "$chr $str\n";
    @site = %site = %paths = %from  = %to = ();

    foreach $id(keys(%{$csi{$chrstr}})) {
    	if(keys(%{$loc{$id}})>1) {
	    $trans_splicing++;
	    next;
    	}
        @list = sort {$a <=> $b} @{$exons{$id}};
        @list = reverse @list if($str eq "-");

        for($i=0;$i<@list;$i++) {
	    $t = ($i%2 == 0) ? (($i==0) ? "S" : 5) : (($i==@list-1) ? "E" : 3); 
	    $key = "$list[$i];$t";
	    $site{$key} = $list[$i];
	    $type{$key} = $t;
	    $list[$i] = $key;
    	}

    	for($i=0;$i<@list-1;$i++) {    
	    $from{$list[$i]}{$list[$i+1]} = $to{$list[$i+1]}{$list[$i]} = TRUE;
	    for($j=$i+1;$j<@list;$j++) {
	    	$p = join(",", @list[$i..$j]);
	    	$paths{$list[$i]}{$list[$j]}{$p} = TRUE;
	    }
    	}

    }

    foreach $x(keys(%paths)) {
	next unless(keys(%{$from{$x}})>1);
	foreach $y(keys(%{$paths{$x}})) {
	    next unless(keys(%{$to{$y}})>1);
	    @arr = sort {lex($a,$b)} keys(%{$paths{$x}{$y}});
	    for($i=0;$i<@arr;$i++) {
		for($j=$i+1;$j<@arr;$j++) {
		    %f = ();
		    @a = split /\,/, $arr[$i];
		    @b = split /\,/, $arr[$j];
		    foreach $z(@a,@b) {$f{$z}++;}
		    if(keys(%f)+2==@a+@b) {
			$t = join(":", join(undef, @type{@a}),join(undef, @type{@b}));
			$p1 = join(",", sort {$a<=>$b} @site{@a});
			$p2 = join(",", sort {$a<=>$b} @site{@b});
			print join("\t", $chr, (sort {$a<=>$b} ($site{$x}, $site{$y})), $t, 2, $str, join(":", $p1, $p2)),"\n";
		    }
		}
	    }  
	}
    }
}

print STDERR "Trans-spliced transcripts : $trans_splicing\n" if($trans_splicing);



sub lex {                                                       # the input is two pointers to arrays
    my $flag = -1;
    my @x = split /\,/, @_[0];
    my @y = split /\,/, @_[1];

    return(-1) if(@x < @y);                                     # if the first array is shorter return -1
    return(1)  if(@x > @y);                                     # if the second is shorter, return 1

    for(my $i=0;$i<@x;$i++) {                                   # if both arrays are of the same length
        $flag = 1 unless($x[$i] le $y[$i]);                     # flag comes out -1 if all p[i]<=q[i]
    }

    return($str eq "+" ? $flag : -$flag);                       # this changes if the strand in -1
}



