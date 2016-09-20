library(abind)

opt=c()

opt$inc = "data/human_sv2/fourcelllines.A.inc.tsv"
opt$exc = "data/human_sv2/fourcelllines.A.exc.tsv"
opt$meta = "meta.dat"
opt$formula = "cell"

inc = read.delim(opt$inc)
exc = read.delim(opt$exc)
meta = read.delim(opt$meta)
form = opt$formula

inc[is.na(inc)] <- 0
exc[is.na(exc)] <- 0

apply(inc*exc,1,function(x){length(x[x>0])}) -> nz

data = abind(inc[nz>0.75*ncol(inc),], exc[nz>0.75*ncol(exc),], along = 3)

diffsplice <- function(x) {
    w = data.frame(meta[rownames(x),])
    colnames(w) = colnames(meta)
    glm(as.formula(paste("x", form, sep = "~")), data = w, family=binomial()) -> model
    melt(summary(model)$coefficients[,c('Estimate','Pr(>|z|)')])
}

apply(data, 1, diffsplice) -> df
do.call(rbind, df) -> df1
as.character(lapply(strsplit(rownames(df1),"\\."),function(x){x[1]})) -> df1$exon


