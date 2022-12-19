require(RColorBrewer)

args <- commandArgs(trailing=TRUE)

## args <- c("THG_pruned.fam", "selected_indi_pop.list", "THG_pruned.2.Q_1")

pops <- as.character(read.table(args[1])[,1])
q <- read.table(args[2])
outpng <- args[3]



## t <- read.table(args[1])
## n <- as.character(read.table(args[2])[,1])
## q <- read.table(args[3])
## outpng <- args[4]


## indiPop <- as.character(t[,2])
## names(indiPop) <- t[,1]
## pops <- as.character(indiPop[n])

## write.table(cbind(n, pops, q), file=paste0(outpng, ".data"), quote=F, row.names=F, col.names=F, sep="\t") 

ncolors <- ncol(q)
if(ncolors<3){
    colors <- brewer.pal(3, "Set1")[1:2]
}else {
    colors <- brewer.pal(ncolors, "Set1")
}


reorderQ <- order(q[,1])
pops <- pops[reorderQ]
q <- q[reorderQ,]

reorderPops <- order(pops)
pops <- pops[reorderPops]
q <- q[reorderPops,]

poporder <- unique(pops)
popN <- table(pops)[poporder]

vbar <- cumsum(popN)
labelidx <- round(vbar - popN/2)
vbar <- vbar[-length(vbar)]

firstrow = 1:length(poporder) %% 2 == 0

print(popN)

## temp
spaces=rep(0,sum(popN))
spaces[vbar+1] = 1
## , ylab = "Ancestral Proportions"
bitmap(outpng, res=400, heigh=2, width=6)
par(mar=c(2,2,1,0))
bp <- barplot(t(as.matrix(q)), col=colors, border=NA, xaxt="n", yaxt="n", space=spaces)
axis(1, at=bp[labelidx][firstrow], labels=poporder[firstrow], tick=F, line=-1)
axis(1, at=bp[labelidx][!firstrow], labels=poporder[!firstrow], tick=F, line=0)
axis(2, pos=0)
mtext("Ancestral Proportions", side=2, line=1)
dev.off()
