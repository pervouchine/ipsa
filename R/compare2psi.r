suppressPackageStartupMessages(library("optparse"))
option_list <- list(make_option(c("-x", "--file1"), help="first gff"),
                    make_option(c("-y", "--file2"), help="second gff"),
		    make_option(c("-a", "--name1"), help="first name"),
                    make_option(c("-b", "--name2"), help="second name"),
		    make_option(c("-o", "--out"),   help="pdf output")
)

opt <- parse_args(OptionParser(option_list=option_list))

#file1 = "data/human_sv2/A07/LID16627.A07.gff"
#file2 = "data/human_sv2/C07/LID16627_TranscriptGencV19.C07.gff" 
#name1 = "psi" 
#name2 = "psitx"

file1 = opt$file1
file2 = opt$file2 
name1 = opt$name1 
name2 = opt$name2

pdf(opt$out)
require(VennDiagram)
require(ggplot2)

A = na.omit(read.delim(pipe(paste("awk '$3==\"exon\"' ", file1, " | print_gff_attributes.pl INT", name1)), header=F))
B = na.omit(read.delim(pipe(paste("awk '$3==\"exon\"' ", file2, " | print_gff_attributes.pl INT", name2)), header=F))

X = merge(A,B,by=1)

grid.newpage()
venn.plot <- draw.pairwise.venn(dim(A)[1], dim(B)[1], dim(X)[1], c(name1, name2),fill= c("blue", "red"),alpha=0.3)
grid.draw(venn.plot)
ggplot(X, aes(x=V2.x,y=V2.y)) + geom_point(size=1) + xlab(name1) + ylab(name2)
ggplot(X, aes(x=V2.x,y=V2.y)) + stat_bin2d() + scale_fill_gradient(low="white", high="blue") + xlab(name1) + ylab(name2)
ggplot(X, aes(x=V2.x,y=V2.y)) + stat_bin2d() + scale_fill_gradient(low="white", high="blue", trans="log10") + xlab(name1) + ylab(name2)

apply(X[,c('V2.x','V2.y')],1,function(x){abs(x[2]-x[1])}) -> diff
plot(ecdf(diff),xlab=paste("|",name1,"-",name2,"|"))

dev.off()
