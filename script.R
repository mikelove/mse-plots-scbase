library(LoomExperiment)
files <- list.files(".", ".loom")
dat <- list()
nms <- sub("deng.(.*).loom","\\1",files)
for (i in 1:4) {
  dat[[nms[i]]] <- import(files[i])
}
sapply(dat, dim)
lapply(dat, assayNames)
mean(colSums(assay(dat[["sim.em"]], "matrix"))) # 147538

h5ls("deng.sim.em.loom")

get.mse <- function(x,y) {
  s <- rowSums((x-y)^2, na.rm=TRUE)
  n <- rowSums(!is.nan(x-y))
  s/n
}

counts <- assay(dat[["sim.uq"]], "matrix")
 
# EM version

idx <- mcols(dat[["sim.em"]])$Rhat_ase < 1.1 &
       mcols(dat[["sim.em"]])$Selected > 0
table(idx)
true <- (assay(dat[["full.em"]], "F") / assay(dat[["full.em"]], "matrix"))[idx,]
nopp <- (assay(dat[["sim.em"]], "F") / assay(dat[["sim.em"]], "matrix"))[idx,]
pp <- assay(dat[["sim.em"]], "p_k")[idx,]
cts <- rowMeans(counts[idx,])

# uniq version

idx <- mcols(dat[["sim.uq"]])$Rhat_ase < 1.1 &
       mcols(dat[["sim.uq"]])$Selected > 0
table(idx)
true <- ( assay(dat[["full.uq"]], "uniq_F") /
          (assay(dat[["full.uq"]], "uniq_F") +
           assay(dat[["full.uq"]], "uniq_B")) )[idx,]
nopp <- (assay(dat[["sim.uq"]], "uniq_F") / assay(dat[["sim.uq"]], "uniq_total"))[idx,]
pp <- assay(dat[["sim.uq"]], "p_k")[idx,]
cts <- rowMeans(counts[idx,])

# calculate MSE

mse.nopp <- get.mse(nopp, true)
mse.pp <- get.mse(pp, true)
keep <- rowSums(!is.nan(nopp)) >= 15 & !is.nan(mse.nopp - mse.pp)
table(keep) # 6751 / 5759

library(ggplot2)
d <- data.frame(log10ave=log10(cts + .1),
                nopp.minus.pp=(mse.nopp - mse.pp))
ggplot(d[keep,], aes(x=log10ave,y=nopp.minus.pp)) +
  geom_hex(bins=120) +
  geom_hline(yintercept=0, alpha=.6,color="grey") +
  ylim(-.4,.4) +
  scale_fill_gradient2(low="black",mid="red",high="white",midpoint=20)
