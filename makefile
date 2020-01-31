DIR=example/
MAPTOOLSDIR=maptools-3.2/
SJCOUNTDIR=sjcount-full/
WEB=http://cb.skoltech.ru/dp/ipsa/

###############################################################################################

all :: ${SJCOUNTDIR}sjcount ${MAPTOOLSDIR}bin/transf ${MAPTOOLSDIR}bin/getsegm

${SJCOUNTDIR}sjcount : 
	wget https://github.com/pervouchine/sjcount/archive/v3.1.tar.gz -O v3.1.tar.gz
	tar -xf v3.1.tar.gz 
	rm -f v3.1.tar.gz
	make -C ${SJCOUNTDIR} all

${MAPTOOLSDIR}bin/transf ${MAPTOOLSDIR}bin/getsegm:
	wget https://github.com/pervouchine/maptools/archive/v3.2.tar.gz -O v3.2.tar.gz
	tar -xf v3.2.tar.gz 
	rm -f v3.2.tar.gz
	mkdir -p ${MAPTOOLSDIR}bin/
	make -C ${MAPTOOLSDIR} all

clean ::
	rm -f -r ${SJCOUNTDIR} ${MAPTOOLSDIR} example.dat example.mk

###############################################################################################

all :: example.mk
	# Installation complete. Type 'make run' to execute a test run.

run ::	example.mk ${DIR}hg19.idx ${DIR}hg19.dbx ${DIR}hg19v19.gfx
	make -f example.mk all 
	# Test run complete. The output is in the output/ directory.
	ls -l ${DIR}output/

###############################################################################################

${DIR}gencode.v19.annotation.gtf :
	mkdir -p ${DIR}
	wget ${WEB}gencode.v19.annotation.gtf.gz -O ${DIR}gencode.v19.annotation.gtf.gz
	gunzip -f ${DIR}gencode.v19.annotation.gtf.gz

${DIR}chr22Y.fa :
	mkdir -p ${DIR}
	wget ${WEB}chr22Y.fa.gz -O ${DIR}chr22Y.fa.gz
	gunzip ${DIR}chr22Y.fa.gz

${DIR}hg19.idx ${DIR}hg19.dbx : ${DIR}chr22Y.fa
	mkdir -p ${DIR}
	${MAPTOOLSDIR}bin/transf -exactdir -dir ${DIR}chr22Y.fa -idx ${DIR}hg19.idx -dbx ${DIR}hg19.dbx

${DIR}hg19v19.gfx : ${DIR}gencode.v19.annotation.gtf
	perl Perl/transcript_elements.pl - < ${DIR}gencode.v19.annotation.gtf > ${DIR}hg19v19.gfx

example.dat : 
	wget ${WEB}example.dat -O example.dat

example.mk : example.dat Perl/make.pl
	perl Perl/make.pl -repository ${DIR}input/ -dir ${DIR}output/ -id id -annot ${DIR}hg19v19.gfx -genome ${DIR}hg19 -in  example.dat > example.mk

###############################################################################################

