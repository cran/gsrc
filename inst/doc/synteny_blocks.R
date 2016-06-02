## ---- results="hide"-----------------------------------------------------
library(openxlsx)
library(dbscan)
library(gsrc)

## ------------------------------------------------------------------------
tf1 <- tempfile()
utils::download.file(url="http://www.sciencedirect.com/science/MiamiMultiMediaURL/1-s2.0-S2352340915000062/1-s2.0-S2352340915000062-mmc6.xlsx/311593/html/S2352340915000062/5528a88e468e01866744f28fc42139c6/mmc6.xlsx",
              destfile = tf1, method = "internal")
synA <- openxlsx::read.xlsx(tf1)
utils::download.file(url = "http://www.sciencedirect.com/science/MiamiMultiMediaURL/1-s2.0-S2352340915000062/1-s2.0-S2352340915000062-mmc7.xlsx/311593/html/S2352340915000062/d616ded5ee9399d4fe29df72435780c9/mmc7.xlsx",
              destfile = tf1, method = "internal")
synC <- openxlsx::read.xlsx(tf1)
unlink(tf1)

## ------------------------------------------------------------------------
synA <- synA[synA$unigene %in% synC$unigene, ]
synC <- synC[synC$unigene %in% synA$unigene, ]

## ------------------------------------------------------------------------
synA <- synA[!is.na(synA$unigene),]
synC <- synC[!is.na(synC$unigene),]

## ------------------------------------------------------------------------
synA <- synA[, c(2,5,6,7)]
synC <- synC[, c(2,5,6,7)]

## ------------------------------------------------------------------------
syn_uni <- merge(synA, synC, by = "unigene")
syn_uni <- syn_uni[order(syn_uni$A.Chr, syn_uni$A.start),]

## ------------------------------------------------------------------------
syn_uni$Alev <- as.factor(syn_uni$A.Chr)
syn_uni$Clev <- as.factor(syn_uni$C.Chr)
syn_uni$AGlo <- syn_uni$A.start
syn_uni$CGlo <- syn_uni$C.start

## ------------------------------------------------------------------------
amaxs <- sapply(levels(syn_uni$Alev), function(x) max(syn_uni$A.start[syn_uni$A.Chr==x]))
for(i in 2:length(levels(syn_uni$Alev))){
  syn_uni$AGlo[syn_uni$A.Chr == levels(syn_uni$Alev)[i]] <- syn_uni$AGlo[syn_uni$A.Chr == levels(syn_uni$Alev)[i]] + sum(amaxs[1:(i-1)])
}

## ------------------------------------------------------------------------
cmaxs <- sapply(levels(syn_uni$Clev), function(x) max(syn_uni$C.start[syn_uni$C.Chr==x]))
for(i in 2:length(levels(syn_uni$Clev))){
  syn_uni$CGlo[syn_uni$C.Chr == levels(syn_uni$Clev)[i]] <- syn_uni$CGlo[syn_uni$C.Chr == levels(syn_uni$Clev)[i]] + sum(cmaxs[1:(i-1)])
}

## ------------------------------------------------------------------------
csamax <- cumsum(amaxs)
cscmax <- cumsum(cmaxs)

## ---- fig.show = "hold", fig.width = 10, fig.height = 10-----------------
plot(syn_uni$AGlo, syn_uni$CGlo, pch = 19,cex=0.2, col = rgb(0,0,0,alpha = 0.05))
abline(v = c(0, csamax), h=c(0, cscmax))

## ------------------------------------------------------------------------
syn_uni2 <- syn_uni[, c(2, 3, 5, 6)]
names(syn_uni2) <- c("chr1", "pos1", "chr2", "pos2")
synteny_blocks <- find_blocks(syn_uni2, eps = 2000000, minPts = 50,
                              minLength = 1000000, maxLength = 10000000)

## ------------------------------------------------------------------------
synteny_blocks$blocks$off1 <- c(0,cumsum(amaxs))[match(synteny_blocks$blocks$chr1, names(amaxs))]
synteny_blocks$blocks$off2 <- c(0,cumsum(cmaxs))[match(synteny_blocks$blocks$chr2, names(cmaxs))]

## ---- fig.show = "hold", fig.width = 10, fig.height = 5------------------
plot.new()
cols <- rainbow(n = 10, start = 0, end = 1, alpha = 0.2)
max1 <- max(synteny_blocks$blocks$end1) 
max2 <- max(synteny_blocks$blocks$end2) 
axis(3, at = (csamax - min(amaxs) / 2) / max(csamax), 
     labels = unique(synA$A.Chr), tick = FALSE, cex.axis = 0.8, las = 1)
axis(1, at = (cscmax - min(cmaxs) / 2) / max(cscmax), labels = unique(synC$C.Chr), tick = FALSE, cex.axis = 0.8, las = 1)
y <- c(1, 1, 0, 0)
for(i in 1:nrow(synteny_blocks$blocks)){
  x <- c(synteny_blocks$blocks$start1[i], 
         synteny_blocks$blocks$end1[i], 
         synteny_blocks$blocks$end2[i], 
         synteny_blocks$blocks$start2[i])
  x <- c(x[1:2] / max1, x[3:4] / max2)
  polygon(x, y, col = cols[unique(synteny_blocks$blocks$chr1) %in% synteny_blocks$blocks$chr1[i]])
}


## ---- cache=FALSE--------------------------------------------------------
sessionInfo()

