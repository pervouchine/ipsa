1) aggregate.pl

This utility aggregates tsv output of sjcount by the 5th column (offset) and outputs BED6+3 with three extra columns being (7) total count, (8) staggered read count, (9) entropy

	-margin ..., the margin for offset, default=0
	-maxintron ..., max intron length, default=0
	-minintron ..., min intron length, default=0
	-readLength ..., the read length, default=0
	-tsv ..., the input file, obligatory

2) annotate.pl

This utility takes an aggregated BED6+3 file, the genomic annotation, and the genome, and outputs BED6+3+2 with two more columns: (10) annotation status and (11) splice sites
If BED6+3 was strandless ('.' in column 6) then each line will be reported twice, one for '+' and and for '-' strand

	-annot ..., the annotation (gtf), obligatory
	-bed ..., the input file, obligatory
	-dbx ..., the genome (dbx), obligatory
	-idx ..., the genome (idx), obligatory

3) biorep.pl
*** this utility is deprecated *** 

4) choose_strand.pl

This utility takes a BED6+3+2 file, where (10) is the annotation status and (11) is splice sites, and selects for each pair of beg/end the strand based on these two columns

	-annot ..., annotation column, default=10
	-bed ..., input file, obligatory
	-sites ..., splice site column, default=11

5) constrain_ssc.pl

This utility takes a BED6 ssc file and constraints its content to splice sites which are present in BED6 ssj file

	-annot ..., annotation column, default=10
	-sites ..., splice site column, default=11
	-ssc ..., input ssc (bed) file, obligatory
	-ssj ..., input ssj (bed) file, obligatory

6) make.pl

This utility creates a makefile for the sj pipeline, taking the index file from STDIN and printing the makefile to STDOUT

	-annot ..., the annotation (gtf), obligatory
	-block ..., the blocking field for merge
	-dir ..., the output directory, obligatory
	-entropy ..., entropy lower threshold, default=3
	-genome ..., the genome (without .dbx or .idx), obligatory
	-group ..., the grouping field for IDR, obligatory
	-idr ..., IDR upper threshold, default=0.1
	-margin ..., margin for aggregate, default=5
	-merge ..., the name of the output to merge in case if blocks are missing
	-mincount ..., min number of counts for the denominator, default=20
	-param ..., parameters passed to sjcount
	-repository ..., the repository subdirectory for bam files

7) merge_gff.pl

This utility merges gtf/gff files into rectangular tables with colnames

	Example: Perl/merge_gff.pl -i <gtf_1> <label_1> -i <gtf_2> <label_2> ... -o feature_1 <tsv1> -o feature_2 <tsv2>
	-i ..., input gtf file name and label, array=hash
	-o ..., feature and output tsv file name, array=hash
	-percent ..., do not report rows that have more than <value> NAs

8) sam2ssj.pl

This is a line-based utility for extracting splice junctions from SAM files

	-lim ..., stop after this number of lines (for debug), default=0
	-read1 ..., flip read1 yes/no (1/0), default=1
	-read2 ..., flip read2 yes/no (1/0), default=0


9) transcript_elements.pl

This utility takes a gtf annotation and reformats it into a more compact, quickly readable file. Only exons are taken into account.

	-gtf ..., the input file, obligatory
	-source ..., the content of the source field, default=SJPIPE

10) tx.pl

This utility computes exon inclusion rates based on the quantifiaction data

	-annot ..., the genomic annotation, obligatory
	-minc ..., the minimum sum of rpkms to consider the inclusion ratio as reliable, default=0
	-quant ..., the transcript quantification file (gtf/gff), obligatory

11) zeta.pl

This utility computes (a) the global exon inclusion and processing rates for a given set of annotated exons and (b) inclusion and processing rates of SJs

	-annot ..., the annotation (GTF) file, obligatory
	-mincount ..., the min value of the denominator of the fraction, default=20
	-ssc ..., the input ssc (BED) file
	-ssj ..., the input ssj (BED) file, obligatory
	-stranded ..., 1(yes) or 0(no), default=1

