#!/usr/bin/perl
use Perl::utils;

if(@ARGV==0) {
    print STDERR "This utility computes splicing graph\n";
}

parse_command_line(annot   => {description=>'the annotation (gtf)', ifunreadable=>'annotation not specified'},
		   lim     => {description=>'number of lines from gtf to read (for debug)', default=>0},
		   maxpath => {description=>'max path length', default=>10},
		   mindim  => {description=>'min bubble dimension to be reported', default=>2},
		   fraction=> {description=>'fractionate exonic segments', store=>TRUE},
		   only    => {description=>'consider only annotated paths', store=>TRUE});


#################################################################################################################

print STDERR "[<=$annot";
open FILE, $annot || die;
while($line=<FILE>) {
    chomp $line;
    ($chr, $source, $feature, $beg, $end, $score, $str, $frame, $attr) = split /\t/, $line;
    next unless($feature eq "exon");
    $str = strand_c2i($str);
    ($beg, $end) = reverse ($beg, $end) if($str<0);
    %attr = get_attributes($attr);
    foreach $transcript(split /\,/, $attr{'transcript_id'}) {
	push @{$exons{$chr}{$str}{$transcript}}, [$beg, $end];
    }
    $n++;
    last if($lim>0 && $n>$lim);
}
print STDERR "]\n";

foreach $chr(sort keys(%exons)) {
    foreach $str(sort keys(%{$exons{$chr}})) {
	%boundary = %site = %index = ();
	%paths = %from = %to = ();
	
	foreach $transcript(keys(%{$exons{$chr}{$str}})) {							# memorize exon boundaries
	    @segm = sort {$a->[0] <=> $b->[0]} @{$exons{$chr}{$str}{$transcript}};
	    @segm = reverse @segm if($str<0);
	    for($i=0;$i<@segm;$i++) {
		$boundary{$segm[$i]->[0]}{5} = 1;								# 5 = 5' exon boundary = A = acceptor site 
		$boundary{$segm[$i]->[1]}{3} = 1;								# 3 = 3' exon boundary = D = donor site
	    }
	    $boundary{$segm[0]->[0]}{B} = $boundary{$segm[-1]->[0]}{E} = 1;					# B and E = specisl labels for transcript start and end
        }

        @pos = ();                                                      					# @pos keeps positions
        foreach $x(sort{$a<=>$b} keys(%boundary)) {
            push @pos, $x unless($boundary{$x}{5} && $boundary{$x}{3});						# except ones that are double boundaries
        }
        @pos = reverse @pos if($str<0);
        for($i=0;$i<@pos;$i++) {										# index is created here to know
            $index{$pos[$i]} = $i;										# for each boundary which is the next boundary
	    $site{$pos[$i]} = $boundary{$pos[$i]}{5} ? 5 : 3;							# %site keeps the information about 5/3 exon boundary
        }

	$m=0;
        %stub = ();
	foreach $transcript(keys(%{$exons{$chr}{$str}})) {
	    progressbar(++$m, 0+keys(%{$exons{$chr}{$str}}), join("\t",$chr,strand_i2c($str), undef));	
	    %trace = ();
	    @segm = sort {$a->[0] <=> $b->[0]} @{$exons{$chr}{$str}{$transcript}};
	    @segm = reverse @segm if($str<0);
	    for($i=0;$i<@segm;$i++) {
		if($i>0 && $index{$segm[$i]->[0]} && $index{$segm[$i-1]->[1]}) {
		    $from{$segm[$i-1]->[1]}{$segm[$i]->[0]} = $to{$segm[$i]->[0]}{$segm[$i-1]->[1]} = 1;	# 1 means intron
		}
		next unless($index{$segm[$i]->[0]} && $index{$segm[$i]->[1]});
                if($fraction) {                                                                                 # if graph type == 5533
                    for($j = $index{$segm[$i]->[0]}; $j<$index{$segm[$i]->[1]}; $j++) {                         #   divide each exon into atomic segments
                        $from{$pos[$j]}{$pos[$j+1]} = $to{$pos[$j+1]}{$pos[$j]} = 0;                            #   0 means exon
                        $trace{$pos[$j]} = $trace{$pos[$j+1]} = 1;                                              #   memorize all vertices for the path
                    }
                }
                else {                                                                                          # if graph type == 5353
                    $from{$segm[$i]->[0]}{$segm[$i]->[1]} = $to{$segm[$i]->[1]}{$segm[$i]->[0]} = 0;            #   0 means exon
                    $trace{$segm[$i]->[0]} = $trace{$segm[$i]->[1]} = 1;                                        #   memorize all vertices for the path
                }
            }
	    next unless($only);
            @arr = sort {$a<=>$b} keys(%trace);                                                             	# trace keeps all vertices of the given transcript
            @arr = reverse @arr if($str<0);                                                                 	#
            for($i = 0; $i < @arr; $i++) {                                                                  	# for each
                for($j=$i + 1; $j < @arr; $j++) {                                                           	#   pair of vertices
                    $stub{join("\t", @arr[$i .. $j])}++ unless($j-$i>$maxpath);                             	#     keep path from i to j in stubs
                }
            }
    	}

	$m=0;
        if($only) {                                                                             		#Compute splicing graph from stubs
            foreach $s(keys(%stub)) {										# for each stub
		progressbar(++$m, 0+keys(%stub), "Paths from stub\t");
                @arr = split /\t/, $s;										#  take the path
                push @{$paths{$arr[0]}{$arr[-1]}}, [@arr];							#    set path{a}{b} = the path
            }
        }
        else {                                                                                  		#Compute splicing graph from scratch
            for($i=@pos-1;$i>=0;$i--) {                                                         		# going 3' to 5'
		progressbar(++$m, 0+@pos, "All paths\t");
                $x = $pos[$i];                                                                  		# for each site x 
                foreach $y(keys(%{$from{$x}})) {                                                		#   for each site y: x~y (~ denotes intron or exon)
                    push @{$paths{$x}{$y}}, [$x, $y];                                           		#     add path {x,y}
                    foreach $z(keys(%{$paths{$y}})) {                                           		#     for each z for which a path from y exists
                        foreach $p(@{$paths{$y}{$z}}) {                                         		#       for each path p={y,...,z}   
                            push @{$paths{$x}{$z}}, [$x, @{$p}] unless(@{$p}>=$maxpath);        		#       add path {x,y,...,z} unless p is too long
                        }
                    }
                }
            }
        }

	#Analyze splicing graph
	$k = 0;
	foreach $x(keys(%paths)) {                                                      # for each site x
	    progressbar(++$k, 0+keys(%paths), "Splice graph\t");			#
	    next unless(keys(%{$from{$x}}) >= $mindim);					#   drop x if less than mindim edges start from it.
            foreach $y(keys(%{$paths{$x}})) {                                           #   for each site y ; x->...-> y
		next unless(keys(%{$to{$y}}) >= $mindim);				#     drop y if less than mindim edges start from it
		%belongs = ();								#     initialize the internal nodes relationship
		for(my $i=0;$i<@{$paths{$x}{$y}};$i++) {				#     for each i
		    my $p = $paths{$x}{$y}->[$i];                                       #       and path p[i]
		    for(my $j=1;$j<@{$p}-1;$j++) {					#	and for each its INTERNAL node
			$belongs{$p->[$j]}{$i} = 1;					#	  make them incident %belongs{node}{path}
		    }
		}
		%equiv = ();								# an equivqlence relationship over paths %equiv{path}{path}
		for(my $i=0;$i<@{$paths{$x}{$y}};$i++) {				# for each path p[i] 
		    for(my $j=$i+1;$j<@{$paths{$x}{$y}};$j++) {				#   for each path q[j], j>i
			foreach $node(keys(%belongs)) {					#     run over all nodes
			    $equiv{$j}{$i} = $equiv{$i}{$j} = 1 if($belongs{$node}{$i} && $belongs{$node}{$j});
			}								#     make p ~ q if p and q have a common node
		    }									#
		    $equiv{$i}{$i} = 1;							#   p ~ p
		}

                my @clusters = cliques(\%equiv);					# @clusters are equivalence classes of paths from x to y factorized by ~

                if(@clusters >= $mindim) {						# if there are at last mindim such classes
                    my @representative = ();                                            #   array of representative paths, one from each cluster
                    my %model = ();                                                     #   model array for drawing
                    foreach $cluster(@clusters) {                                       #   foreach equivalence class
                        @arr = sort {lex($paths{$x}{$y}->[$a],$paths{$x}{$y}->[$b])} @{$cluster}; # sort paths in the class LEXICOGRAPHICALLY 
                        $p = $paths{$x}{$y}->[shift(@arr)];                             #     take the shortest in the cluster as a representative
                        push @representative, $p;
                        @model{@{$p}} = @site{@{$p}};
                    }
                    @positions = sort {$a<=>$b} keys(%model);                           # positions are common from the %model hash 
                    @positions = reverse @positions if($str<0);                         # sort positions according to the strand
		    @arr = @crd = ();
		    @representative = sort{lex($a,$b)} @representative;
                    foreach $p(@representative) {
                        push @arr, join("", @model{@{$p}});
			push @crd, join(",", @{$p});
                    }
		    print join("\t", $chr, $x, $y, strand_i2c($str), 0+@clusters, join(":", @arr), join(":", @crd)), "\n";
                }
	    }
	}
    }
}




sub lex {										# the input is two pointers to arrays
    my $flag = -1;

    return(-1) if(@{@_[0]} < @{@_[1]});							# if the first array is shorter return -1
    return(1)  if(@{@_[0]} > @{@_[1]});							# if the second is shorter, return 1

    for(my $i=0;$i<@{@_[0]};$i++) {							# if both arrays are of the same length
        $flag = 1 unless(@_[0]->[$i] le @_[1]->[$i]);					# flag comes out -1 if all p[i]<=q[i]
    }

    return(@_[0]->[1] > @_[0]->[0] ? $flag : -$flag);					# this changes if the strand in -1 (indirectly checking by array order)
}


##############################################################################################################
# I/O functions
##############################################################################################################

sub print_relation {
    my $x;
    my $y;
    foreach $x (keys(%{@_[0]})) {
	foreach $y (keys(%{${@_[0]}{$x}})) {
	    print "$x\t$y\n";
	}
    }
    print "\n";
}

##############################################################################################################
# Algebraic functions
##############################################################################################################

sub symmetrify {
    # input: pointer to a relation
    foreach my $x (keys(%{@_[0]})) {
	next unless($x);
	foreach my $y (keys(%{${@_[0]}{$x}})) {
	    next unless($y);
	    ${@_[0]}{$y}{$x}=${@_[0]}{$x}{$y};
	}
	${@_[0]}{$x}{$x}=1
    }
}

sub product {
    # input: pointer to the first term, pointer to the second term, pointer to the output
    foreach my $x (keys(%{@_[0]})) {
	next unless($x);
	foreach my $y (keys(%{${@_[0]}{$x}})) {
	    next unless($y);
	    foreach my $z (keys(%{${@_[1]}{$y}})) {
		next unless($z);
		${@_[2]}{$x}{$z}=1;
	    }
	}
    }
}

sub adjoint {
    foreach my $x (keys(%{@_[0]})) {
	foreach my $y (keys(%{${@_[0]}{$x}})) {
	    ${@_[1]}{$x}{$y}=1;
	}
    }
}

##############################################################################################################
# Transitive closure 
##############################################################################################################

sub getall {
    my $x = @_[1];
    my @res = ($x);
    $BINREL_FLAG{$x}=1;
    foreach my $y (keys(%{${@_[0]}{$x}})) {
	next if($y eq undef);
	push @res, getall(@_[0],$y) unless($BINREL_FLAG{$y});
    }
    return(@res);
}

sub transitive_closure {
    %BINREL_FLAG=();
    my $x;
    my %res=();
    foreach $x (keys(%{@_[0]})) {
	next if($x eq undef);
	unless($BINREL_FLAG{$x}) {
	    foreach $y(getall(@_[0],$x)) {
		next if($y eq undef);
		$res{$x}{$y}=1;
	    }
	}
    }
    return(%res);
}

sub cliques {
    %BINREL_FLAG=();
    my @res = ();
    foreach my $x (keys(%{@_[0]})) {
	next if($x eq undef);
	unless($BINREL_FLAG{$x}) {
	    push @res, [getall(@_[0],$x)];
	}
    }
    return(@res);
}
