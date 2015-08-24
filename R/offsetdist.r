suppressPackageStartupMessages(library("ggplot2"))

cmd_args = commandArgs();
pdf = cmd_args[length(cmd_args)]

df = data.frame();
for(i in 6:(length(cmd_args)-1)) {
    print(cmd_args[i])
    a = read.delim(cmd_args[i], header=F)
    a$V4 = a$V3*100/sum(a$V3)
    df = rbind(df, a)
}

colnames(df) = c('file','offset','count','freq')

unlist(lapply(strsplit(as.character(df$file),"/"),function(x){x[length(x)]})) -> x
gsub(".A02.ss[jc].log","", x) -> df$sample

pdf(pdf)
theme_set(theme_bw(base_size = 12))
ggplot(df,aes(x=offset,y=count,colour=sample)) + geom_line() + geom_point() + ylab("Number of reads")
ggplot(df,aes(x=offset,y=freq,colour=sample)) + geom_line() + geom_point()  + ylab("Proportion of reads, %")
dev.off()

