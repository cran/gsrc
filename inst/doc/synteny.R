## ---- results="hide"-----------------------------------------------------
library(openxlsx)
library(dbscan)
library(gsrc)

## ------------------------------------------------------------------------
tf1 <- tempfile()
url1 <- paste("http://www.sciencedirect.com/science/MiamiMultiMediaURL/",
              "1-s2.0-S2352340915000062/1-s2.0-S2352340915000062-mmc6.xlsx/",
              "311593/html/S2352340915000062/5528a88e468e01866744f28fc42139c6/",
              "mmc6.xlsx", sep = "")
utils::download.file(url = url1, destfile = tf1, method = "internal", mode = "wb")
synA <- openxlsx::read.xlsx(tf1)
url2 <- paste("http://www.sciencedirect.com/science/MiamiMultiMediaURL/",
              "1-s2.0-S2352340915000062/1-s2.0-S2352340915000062-mmc7.xlsx/",
              "311593/html/S2352340915000062/d616ded5ee9399d4fe29df72435780c9/",
              "mmc7.xlsx", sep = "")
utils::download.file(url = url2, destfile = tf1, method = "internal", mode = "wb")
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


## ---- results="hide", eval=FALSE-----------------------------------------
#  library(DECIPHER)
#  library(dbscan)
#  library(gsrc)

## ---- eval = FALSE-------------------------------------------------------
#  fas <- c(Genome1="PATHTOSUBGENOMEA",
#           Genome2="PATHTOSUBGENOMEB")
#  db <- "PATHTODB"
#  for (i in seq_along(fas)) {
#    Seqs2DB(fas[i], "FASTA", db, names(fas[i]))
#  }
#  synteny <- FindSynteny(db, processors = NULL, storage = 5)
#  save(synteny, file = "PATHTOSYNTENYFILE")

## ---- results = "hide"---------------------------------------------------
load("brassica_synteny_filtered.rda")

## ---- eval = FALSE-------------------------------------------------------
#  syn_uni <- as.data.frame(synteny[[2]])
#  syn_uni <- syn_uni[syn_uni$score>1500,]
#  syn_uni <- syn_uni[, c(1,5,7,2,6,8)]
#  syn_uni <- cbind(syn_uni, rowMeans(syn_uni[, 2:3]), rowMeans(syn_uni[, 5:6]))
#  syn_uni <- syn_uni[, c(1,7,4,8)]
#  colnames(syn_uni) <- c("index1", "start1", "index2", "start2")
#  syn_uni <- syn_uni[order(syn_uni$index1, syn_uni$start1),]
#  syn_uni$Alev <- as.factor(syn_uni$index1)
#  syn_uni$Clev <- as.factor(syn_uni$index2)
#  syn_uni$AGlo <- syn_uni$start1
#  syn_uni$CGlo <- syn_uni$start2

## ------------------------------------------------------------------------
amaxs <- sapply(levels(syn_uni$Alev), function(x) max(syn_uni$start1[syn_uni$index1==x]))
for(i in 2:length(levels(syn_uni$Alev))){
  syn_uni$AGlo[syn_uni$index1 == levels(syn_uni$Alev)[i]] <- syn_uni$AGlo[syn_uni$index1 == levels(syn_uni$Alev)[i]] + sum(amaxs[1:(i-1)])
}

cmaxs <- sapply(levels(syn_uni$Clev), function(x) max(syn_uni$start2[syn_uni$index2==x]))
for(i in 2:length(levels(syn_uni$Clev))){
  syn_uni$CGlo[syn_uni$index2 == levels(syn_uni$Clev)[i]] <- syn_uni$CGlo[syn_uni$index2 == levels(syn_uni$Clev)[i]] + sum(cmaxs[1:(i-1)])
}
csamax <- cumsum(amaxs)
cscmax <- cumsum(cmaxs)


## ---- fig.show = "hold", fig.width = 8, fig.height = 8-------------------
plot(syn_uni$AGlo, syn_uni$CGlo, pch = 19, cex = 0.2, col = rgb(0, 0, 0, alpha = 0.05))
abline(v = c(0, csamax), h=c(0, cscmax))

## ------------------------------------------------------------------------
syn_uni2 <- syn_uni[, c(1, 2, 3, 4)]
names(syn_uni2) <- c("chr1", "pos1", "chr2", "pos2")

## ------------------------------------------------------------------------
blocks <- gsrc::find_blocks(syn_uni2, eps = 1500000, minPts = 20, minLength = 1000000, maxLength = 10000000)

## ---- fig.show = "hold", fig.width = 9, fig.height = 5-------------------
plot.new()
cols <- rainbow(
  n = 10, start = 0, end = 1, alpha = 0.2
)
max1 <- max(blocks$blocks$end1)
max2 <- max(blocks$blocks$end2)
axis(
  3, at = (csamax - min(amaxs) / 2) / max(csamax),
  labels = sort(unique(blocks$blocks$chr1)), tick = FALSE, cex.axis = 0.8, las = 1
)
axis(
  1, at = (cscmax - min(cmaxs) / 2) / max(cscmax),
  labels = sort(unique(blocks$blocks$chr2)), tick = FALSE, cex.axis = 0.8, las = 1
)
y <- c(1, 1, 0, 0)
for (i in 1:nrow(blocks$blocks)) {
  x <- c(
    blocks$blocks$start1[i],
    blocks$blocks$end1[i],
    blocks$blocks$end2[i],
    blocks$blocks$start2[i]
  )
  x <- c(x[1:2] / max1, x[3:4] / max2)
  polygon(x, y, col = cols[unique(blocks$blocks$chr1) %in% blocks$blocks$chr1[i]])
}

## ---- results="hide", eval=FALSE-----------------------------------------
#  library(DECIPHER)
#  library(dbscan)
#  library(gsrc)

## ---- eval = FALSE-------------------------------------------------------
#  fas <- c(Genome1="PATHTOGENOME1",
#           Genome2="PATHTOGENOME2")
#  db <- "PATHTODB"
#  for (i in seq_along(fas)) {
#    Seqs2DB(fas[i], "FASTA", db, names(fas[i]))
#  }
#  synteny <- FindSynteny(db, processors = NULL, storage = 5)
#  save(synteny, file = "PATHTOSYNTENYFILE")

## ---- results = "hide"---------------------------------------------------
load("cotton_synteny_filtered.rda")

## ---- eval = FALSE-------------------------------------------------------
#  syn_uniC <- as.data.frame(synteny[[2]])
#  syn_uniC <- syn_uniC[syn_uniC$score>1500,]
#  syn_uniC <- syn_uniC[, c(1,5,7,2,6,8)]
#  syn_uniC <- cbind(syn_uniC, rowMeans(syn_uniC[, 2:3]), rowMeans(syn_uniC[, 5:6]))
#  syn_uniC <- syn_uniC[, c(1,7,4,8)]
#  colnames(syn_uniC) <- c("index1", "start1", "index2", "start2")
#  syn_uniC <- syn_uniC[order(syn_uniC$index1, syn_uniC$start1),]
#  syn_uniC$Alev <- as.factor(syn_uniC$index1)
#  syn_uniC$Clev <- as.factor(syn_uniC$index2)
#  syn_uniC$AGlo <- syn_uniC$start1
#  syn_uniC$CGlo <- syn_uniC$start2

## ------------------------------------------------------------------------
amaxs <- sapply(levels(syn_uniC$Alev), function(x) max(syn_uniC$start1[syn_uniC$index1==x]))
for(i in 2:length(levels(syn_uniC$Alev))){
  syn_uniC$AGlo[syn_uniC$index1 == levels(syn_uniC$Alev)[i]] <- syn_uniC$AGlo[syn_uniC$index1 == levels(syn_uniC$Alev)[i]] + sum(amaxs[1:(i-1)])
}

cmaxs <- sapply(levels(syn_uniC$Clev), function(x) max(syn_uniC$start2[syn_uniC$index2==x]))
for(i in 2:length(levels(syn_uniC$Clev))){
  syn_uniC$CGlo[syn_uniC$index2 == levels(syn_uniC$Clev)[i]] <- syn_uniC$CGlo[syn_uniC$index2 == levels(syn_uniC$Clev)[i]] + sum(cmaxs[1:(i-1)])
}
csamax <- cumsum(amaxs)
cscmax <- cumsum(cmaxs)

## ---- fig.show = "hold", fig.width = 8, fig.height = 8-------------------
plot(syn_uniC$AGlo, syn_uniC$CGlo, pch = 19, cex = 0.2, col = rgb(0, 0, 0, alpha = 0.05))
abline(v = c(0, csamax), h=c(0, cscmax))

## ------------------------------------------------------------------------
syn_uni2 <- syn_uniC[, c(1, 2, 3, 4)]
names(syn_uni2) <- c("chr1", "pos1", "chr2", "pos2")

## ------------------------------------------------------------------------
blocks <- gsrc::find_blocks(syn_uni2,eps = 5000000, minPts = 20, minLength = 10000000, maxLength = 50000000)

## ---- fig.show = "hold", fig.width = 9, fig.height = 5-------------------
plot.new()
cols <- rainbow(
  n = 13, start = 0, end = 1, alpha = 0.2
)
max1 <- max(blocks$blocks$end1)
max2 <- max(blocks$blocks$end2)
axis(
  3, at = (csamax - min(amaxs) / 2) / max(csamax),
  labels = unique(blocks$blocks$chr1), tick = FALSE, cex.axis = 0.8, las = 1
)
axis(
  1, at = (cscmax - min(cmaxs) / 2) / max(cscmax),
  labels = unique(blocks$blocks$chr2), tick = FALSE, cex.axis = 0.8, las = 1
)
y <- c(1, 1, 0, 0)
for (i in 1:nrow(blocks$blocks)) {
  x <- c(
    blocks$blocks$start1[i],
    blocks$blocks$end1[i],
    blocks$blocks$end2[i],
    blocks$blocks$start2[i]
  )
  x <- c(x[1:2] / max1, x[3:4] / max2)
  polygon(x, y, col = cols[unique(blocks$blocks$chr1) %in% blocks$blocks$chr1[i]])
}

