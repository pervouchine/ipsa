A = read.delim("data/shin/hs.A.psi.tsv")
B = read.delim("data/shin/mm.A.psi.tsv")

map = read.delim("map1.tsv", header=F)
C = subset(map, V2==1)[,c(1,8)]

f<-function(X){apply(X,1,function(x){length(x[!is.na(x)])})->n;X[n>=0.9*ncol(X),]}
g<-function(X){apply(X,1,function(x){max(x,na.rm=T)-min(x,na.rm=T)})->d;X[d>=0.05,]}

A1 = g(f(A))
B1 = g(f(B))

merge(merge(C, A1, by.x=1, by.y=0), B1, by.x=2, by.y=0) -> X
print(dim(A1))
print(dim(B1))
print(dim(X))
X[is.na(X)]<-0.5
write.table(as.matrix(dist(t(X[,c(-1,-2)]))), file="hsmm.tsv", col.names=T, row.names=T, quote=F, sep="\t")

