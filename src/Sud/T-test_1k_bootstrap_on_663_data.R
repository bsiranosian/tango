#setwd("~/GitHub/tango")
#vals <- as.matrix(read.table('figures/bootstrap_trees/663_1000boot_data.txt',skip = 6))
  
#vals <- as.matrix(read.table('~/Desktop/1K Bootstrap on 663 Data Matrix',skip = 6))  
  
vals <- as.matrix(read.table("~/Desktop/1kbootstrap",skip = 6))  
  
cluster_splits <- c(661,655,639,643,617,625,653,637,650,659,646,644,640,635,633,642)
subcluster_splits <- c(630,628,638,623,626,624,620,632,608,609,631,618,614,615,636,627,616,602,601,598,603,611,595)
subcluster_separate_splits <- c(645,634,621,610,649,604)
singleton_splits <- c(641,648,657,654,647,621)

vals.splits.clusters <- vals[cluster_splits,]
vals.splits.subclusters <- vals[subcluster_splits,]
vals.splits.separates <- vals[cluster_separate_splits,]
vals.splits.singletons <- vals[singleton_splits,]

vals.nosplits.clusters <- vals[!(1:662 %in% cluster_splits), ]
vals.nosplits.subclusters <- vals[!(1:662 %in% subcluster_splits), ]
vals.nosplits.separates <- vals[!(1:662 %in% subcluster_separate_splits), ]
vals.nosplits.singletons <- vals[!(1:662 %in% singleton_splits), ]
vals.all <- vals[,]


par(mfrow=c(2,2))

hist(1-vals.splits.clusters[, 1])
hist(1-vals.nosplits.clusters[, 1])
hist(1-vals.all[, 1])

t.test(vals.splits.clusters[,1], vals.nosplits.clusters[,1])
t.test(vals.splits.clusters[,1], vals[,1])



