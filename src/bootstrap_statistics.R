setwd("~/GitHub/tango")
vals <- as.matrix(read.table('data/bootstrap_trees/H60_1kboot_data.txt',skip = 5))

#splits that define clusters
# 39: F1/F2
# 37: F1
# 14: G
# 43: H1/H2
# 38: H1
# 42: B2/B4
# 12: B2
# 22: B4
# 3: B1
# 41: I
# 17: B3
# 44: A
# 27: A1
# 37: A2/A3
# 36: A2
# 16: D1
# 19: E
# 40: C
# 9: C1

cluster_splits <- c(39,37,14,43,38,42,12,22,3,41,17,44,27,37,36,16,19,40,9)
vals.splits <- vals[cluster_splits,]
vals.nosplits <- vals[!(1:59 %in% cluster_splits), ]
vals.all <- vals[,]

par(mfrow=c(2,2))
hist(1-vals.splits[, 1])
hist(1-vals.nosplits[, 1])
hist(1-vals[, 1])


