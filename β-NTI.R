setwd("/Users/MycologyLab/yanqi/DADA2/OIL_SPIL_SAMPLE/new_0410/")

#建树 已经有了asv的树 asv_tree,nwk

library(ape)
library(picante)
library(vegan)
library(reshape2)

dat <- read.csv("asv_t1o.csv", row.names = 1, check.names = FALSE)
seqs <- dat$Sequence
comm <- dat[, !(colnames(dat) %in% "Sequence")]
comm <- t(comm)

# 去掉全0 ASV
comm <- comm[, colSums(comm) > 0]

# 读树
tree <- read.tree("asv_tree.nwk")

# 匹配ASV
common_asv <- intersect(colnames(comm), tree$tip.label)

comm <- comm[, common_asv]
tree <- drop.tip(tree, setdiff(tree$tip.label, common_asv))

pd <- cophenetic(tree)

bmntd.obs <- as.matrix(
  comdistnt(
    comm,
    pd,
    abundance.weighted = TRUE
  )
)
#stage4
rand <- 999
n <- nrow(comm)

bmntd.null <- array(NA, dim = c(n, n, rand))

for (i in 1:rand) {
  
  rand.tree <- tipShuffle(tree)
  pd.rand <- cophenetic(rand.tree)
  
  bmntd.null[,,i] <- as.matrix(
    comdistnt(
      comm,
      pd.rand,
      abundance.weighted = TRUE
    )
  )
}

# 👇 一定在循环外！！
bmntd.mean <- apply(bmntd.null, c(1,2), mean)
bmntd.sd   <- apply(bmntd.null, c(1,2), sd)

betaNTI <- (bmntd.obs - bmntd.mean) / bmntd.sd

#Step 5️⃣ 计算 RCbray
# Bray-Curtis
bray.obs <- as.matrix(vegdist(comm, method = "bray"))

# null
bray.null <- array(NA, dim = c(n, n, rand))

for (i in 1:rand) {
  rand.comm <- randomizeMatrix(comm, null.model = "independentswap")
  
  bray.null[,,i] <- as.matrix(
    vegdist(rand.comm, method = "bray")
  )
}

# RCbray
RCbray <- matrix(NA, n, n)

for (i in 1:n) {
  for (j in 1:n) {
    
    null.dist <- bray.null[i,j,]
    obs <- bray.obs[i,j]
    
    RCbray[i,j] <- (sum(null.dist < obs) - sum(null.dist > obs)) / rand
  }
}

#6
process <- matrix(NA, n, n)

process[betaNTI > 2]  <- "HeS"
process[betaNTI < -2] <- "HoS"

mid <- abs(betaNTI) <= 2

process[mid & RCbray > 0.95]  <- "DL"
process[mid & RCbray < -0.95] <- "HD"
process[mid & abs(RCbray) <= 0.95] <- "DR"
#stage7
df <- melt(process)
colnames(df) <- c("sample1","sample2","process")

# 去掉对角线
df <- df[df$sample1 != df$sample2, ]

prop <- prop.table(table(df$process))
prop

#stage8
barplot(prop,
        col = c("red","blue","green","purple","grey"),
        ylab = "Proportion",
        main = "Assembly Processes")
