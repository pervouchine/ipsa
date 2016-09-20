suppressPackageStartupMessages(library("optparse"))
option_list <- list(make_option(c("-t", "--tsv"),  help="tsv input"),
                    make_option(c("-p", "--pdf"), help="pdf output"))

opt <- parse_args(OptionParser(option_list=option_list))

print(opt$tsv)
data = read.delim(opt$tsv, header=F)
colnames(data) = c('file','degree','offset','count')
with(data,  aggregate(list(count=count), by=list(offset=offset), sum)) -> df
pdf(opt$pdf)
with(df, barplot(count,names.arg=offset,main=opt$tsv))
dev.off()

