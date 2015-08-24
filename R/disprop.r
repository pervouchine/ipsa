suppressPackageStartupMessages(library("ggplot2"))

cmd_args = commandArgs();
pdf = cmd_args[length(cmd_args)]

df = data.frame();
for(i in 6:(length(cmd_args)-1)) {
    print(cmd_args[i])
    a = read.delim(cmd_args[i], header=F)
    a$id  = unlist(lapply(strsplit(as.character(a$V1),"_"),function(x){paste(x[1:3],collapse="_")}))
    a$str = unlist(lapply(strsplit(as.character(a$V1),"_"),function(x){x[4]})) 
    with(a, aggregate(list(count=V2),by=list(id=id),FUN=function(x){log10(x[2])-log10(x[1])})) -> z
    hist(z$count, breaks=100) -> h
    b = data.frame(file=cmd_args[i], mids = h$mids, density = h$density)
    df = rbind(df, b)
}

unlist(lapply(strsplit(as.character(df$file),"/"),function(x){x[length(x)]})) -> x
gsub(".A03.ssj.tsv","", x) -> df$sample

pdf(pdf)
theme_set(theme_bw(base_size = 12))
ggplot(df,aes(x=mids,y=density,colour=sample)) + geom_line() + geom_point() + xlab("log10(+) - log10(-)")
dev.off()

