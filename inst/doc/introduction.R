## ---- results="hide"-----------------------------------------------------
library(gsrc)
require(devtools)
devtools::install_github("grafab/brassicaData")

## ---- eval = FALSE-------------------------------------------------------
#  files <- list.files("/YOUR/DATA/REPOSITORY/",
#                      pattern = "idat",full.names = TRUE)

## ---- eval = FALSE-------------------------------------------------------
#  files <- list.files(system.file("extdata",
#                                  package = "brassicaData"),
#                      full.names = TRUE,
#                      pattern = "idat")

## ---- eval = FALSE-------------------------------------------------------
#  samples <- read_sample_sheets(files =
#                                  list.files(system.file("extdata",package = "brassicaData"),
#                                             full.names = TRUE,
#                                             pattern = "csv"))

## ---- eval = FALSE-------------------------------------------------------
#  controls <- grep("H2O", samples$Names)
#  if(length(controls) > 0) samples <- samples[-controls, ]
#  files <- grep(paste(samples$ID, collapse = "|"), files, value = TRUE)

## ---- eval = FALSE-------------------------------------------------------
#  column_names <- sapply(strsplit(files, split = "/"), FUN=function(x) x[length(x)])

## ------------------------------------------------------------------------
data(dictionary, package = "brassicaData", envir = environment())
head(dictionary)
data(chrPos, package = "brassicaData", envir = environment())
head(chrPos)

## ---- eval = FALSE-------------------------------------------------------
#  raw_data <- read_intensities(files = files,
#                               dict = dictionary,
#                               cnames = column_names,
#                               pos = chrPos)

## ---- eval = FALSE-------------------------------------------------------
#  str(raw_data)

## ---- eval = FALSE-------------------------------------------------------
#  raw_data <- rename_samples(raw_data,
#                             samples = samples[,2:1],
#                             suffix = c("_Grn", "_Red"))

## ------------------------------------------------------------------------
data(raw_napus, package = "brassicaData", envir = environment())

## ---- fig.show = "hold", fig.width = 10, fig.height = 10, fig.cap = "Raw Data Histogram"----
check_raw(raw_napus, thresh = 28000, breaks = 20)

## ---- eval = TRUE--------------------------------------------------------
length(raw_napus$samples)
raw_napus <- filt_samp(raw_napus, check_raw(raw = raw_napus, plot = FALSE, thresh = 28000))
length(raw_napus$samples)

## ---- eval = TRUE, fig.show = "hold", fig.width = 10, fig.height = 10, fig.cap = "Boxplot comparing green and red signal distibutions"----
boxplot(as.vector(raw_napus$raw[, seq(1, length(raw_napus$samples), 2)]),
        as.vector(raw_napus$raw[, seq(2, length(raw_napus$samples), 2)]),
        names = c("Green", "Red"))

## ---- eval = TRUE--------------------------------------------------------
norm_dat <- intens_theta(raw_napus, norm = "both", scaling = "mean", transf = "log")
str(norm_dat)

## ---- eval = TRUE--------------------------------------------------------
head(norm_dat$samples)
norm_dat <- remove_suffix(norm_dat, "_Grn")
head(norm_dat$samples)

## ---- eval = TRUE, fig.show = "hold", fig.width = 10, fig.height = 10, fig.cap="Intensity histogram"----
hist(norm_dat$intensity, breaks = 1000)

## ---- eval = TRUE, fig.show = "hold", fig.width = 10, fig.height = 10, fig.cap="Theta histogram"----
hist(norm_dat$theta, breaks = 1000)

## ---- eval = FALSE-------------------------------------------------------
#  rm(raw_napus)

## ---- eval = TRUE--------------------------------------------------------
norm_dat <- geno_baf_rratio(norm_dat, delthresh = 11)
str(norm_dat)

## ---- eval = TRUE, fig.show = "hold", fig.width = 10, fig.height = 10, fig.cap="B-allele frequency histogram"----
hist(norm_dat$baf, breaks = 1000)

## ---- eval = TRUE, fig.show = "hold", fig.width = 10, fig.height = 10, fig.cap="Barplot of genotypes"----
tmp <- table(norm_dat$geno, useNA = "ifany")
barplot(tmp, names.arg = c(names(tmp)[1:4], "NA"))

## ---- eval = TRUE, fig.show = "hold", fig.width = 10, fig.height = 10, fig.cap="Log R ratio histogram"----
hist(norm_dat$rratio, breaks = 1000)

## ---- eval = TRUE--------------------------------------------------------
norm_dat$theta <- norm_dat$intensities <- NULL

## ---- eval = TRUE--------------------------------------------------------
length(norm_dat$snps)
norm_dat <- filt_snps(norm_dat, norm_dat$snps[is.na(rowMeans(norm_dat$baf, na.rm = TRUE))])
length(norm_dat$snps)

## ---- eval = TRUE--------------------------------------------------------
norm_dat <- segm(norm_dat)
str(norm_dat)

## ---- eval = TRUE--------------------------------------------------------
norm_dat <- cnv(norm_dat, dup = 0.03, del = -0.06)
str(norm_dat)

## ---- eval = TRUE, fig.show = "hold", fig.width = 10, fig.height = 10, fig.cap="Barplot of CNVs"----
barplot(table(norm_dat$cnv))

## ---- eval = TRUE--------------------------------------------------------
data(synteny_blocks, package = "brassicaData", envir = environment())

## ---- eval = TRUE--------------------------------------------------------
norm_dat <- trans_location(norm_dat, synteny_blocks, min1 = 5, min2 = 5, maxdiff = 10)

## ---- eval = TRUE, fig.show = "hold", fig.width = 10, fig.height = 10, fig.cap="Genome structure rearrangements within one sample"----
plot_gsr(norm_dat, sb = synteny_blocks, samp = 1)

## ---- eval = TRUE, fig.show = "hold", fig.width = 10, fig.height = 10, fig.cap="Genome structure rearrangements within one sample, inluding B-allele frequencies and translocations"----
plot_gsr(norm_dat, sb = synteny_blocks, samp = 1, baf = TRUE, tl =TRUE)

## ---- eval = TRUE, fig.show = "hold", fig.width = 10, fig.height = 10, fig.cap="Mean B-alllele frequencies and Log R rations of the population"----
plot_global(norm_dat, sb = synteny_blocks)

