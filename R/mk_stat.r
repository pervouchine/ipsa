suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(ggplot2))
option_list <- list(make_option(c("-i", "--tsv"),  help="tsv input"),
		    make_option(c("-f", "--fct"),  help="every f-th row", default=1, type="double"),
                    make_option(c("-o", "--pdf"), help="pdf output"))

opt <- parse_args(OptionParser(option_list=option_list))

print(opt$tsv)
data = read.delim(opt$tsv, header=F)
colnames(data) = c('ann','log2count','nsp','nsj')

print(subset(data, log2count==0 & nsp==1))

data = subset(data, (nsp-1) %% opt$fct == 0)

#data$ann = factor(data$ann, levels=c(0,1,2,3), labels=c('Both unknown','One known','Both known','Intron known'))

pdf(opt$pdf, width=4, height=4)
the_base_size = 10
theme_set(theme_bw(base_size = the_base_size))
p = ggplot(data, aes(x=log2count,y=nsp, fill=nsj)) + geom_tile() + facet_wrap(~ann) + xlab(expression(log[2](count))) + ylab("# of samples")

p +  scale_fill_gradient("# of SJ", low='white', high='red')
p + scale_fill_gradientn("# of SJ", colours = rev(rainbow(7)[1:6]))
p +  scale_fill_gradient("# of SJ", low='white', high='red', trans="log10")
p + scale_fill_gradientn("# of SJ", colours = rev(rainbow(7)[1:6]), trans="log10")

dev.off()

