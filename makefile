DIR=example/
MAPTOOLSDIR=maptools-2.0/
SJCOUNTDIR=sjcount-3.1/

###############################################################################################

all :: ${SJCOUNTDIR}sjcount ${MAPTOOLSDIR}bin/transf ${MAPTOOLSDIR}bin/getsegm

${SJCOUNTDIR}sjcount : 
	wget https://github.com/pervouchine/sjcount/archive/v3.1.tar.gz -O v3.1.tar.gz
	tar -xf v3.1.tar.gz 
	rm -f v3.1.tar.gz
	make -C ${SJCOUNTDIR} all

${MAPTOOLSDIR}bin/transf ${MAPTOOLSDIR}bin/getsegm:
	wget https://github.com/pervouchine/maptools/archive/v2.0.tar.gz -O v2.0.tar.gz
	tar -xf v2.0.tar.gz 
	rm -f v2.0.tar.gz
	mkdir -p ${MAPTOOLSDIR}bin/
	make -C ${MAPTOOLSDIR} all

clean ::
	rm -f -r ${SJCOUNTDIR} ${MAPTOOLSDIR} example.dat example.mk

###############################################################################################

all :: example.mk
	# Installation complete. Type 'make run' to execute a test run.

run ::	example.mk ${DIR}hg19.idx ${DIR}hg19.dbx ${DIR}hg19v18.gfx
	make -f example.mk all 
	# Test run complete. The output is in the output/ directory.
	ls -l ${DIR}output/

###############################################################################################

${DIR}gencode.v18.annotation.gtf :
	mkdir -p ${DIR}
	wget ftp://ftp.sanger.ac.uk/pub/gencode/release_18/gencode.v18.annotation.gtf.gz -O ${DIR}gencode.v18.annotation.gtf.gz
	gunzip -f ${DIR}gencode.v18.annotation.gtf.gz

${DIR}chr22Y/chr22.fa :
	mkdir -p ${DIR}
	wget http://genome.crg.eu/~dmitri/export/ipsa/chr22Y.tar.gz -O ${DIR}chr22Y.tar.gz
	tar -xf ${DIR}chr22Y.tar.gz -C ${DIR}
	rm -f ${DIR}chr22Y.tar.gz

${DIR}hg19.idx ${DIR}hg19.dbx : ${DIR}chr22Y/chr22.fa
	mkdir -p ${DIR}
	${MAPTOOLSDIR}bin/transf -dir ${DIR}chr22Y/chr22.fa -idx ${DIR}hg19.idx -dbx ${DIR}hg19.dbx

${DIR}hg19v18.gfx : ${DIR}gencode.v18.annotation.gtf
	perl Perl/transcript_elements.pl - < ${DIR}gencode.v18.annotation.gtf > ${DIR}hg19v18.gfx

example.dat : 
	wget http://genome.crg.eu/~dmitri/export/ipsa/example_ipsa.dat -O example.dat

example.mk : example.dat Perl/make.pl
	perl Perl/make.pl -repository ${DIR}input/ -dir ${DIR}output/ -group idrGroup -param '-read1 0' -annot ${DIR}hg19v18.gfx -genome ${DIR}hg19 -merge pooled -in  example.dat > example.mk

###############################################################################################

