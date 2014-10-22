# set the working directory so filenames are relative to it
setwd("~/GitHub/tango")
# read in the data from the bootstrapping tree. This reads it in as a matrix. 
# skip=5 skips the first 5 lines of the file because they aren't relevant 
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

# this is a list of the split numbers that define clusters
cluster_splits <- c(39,37,14,43,38,42,12,22,3,41,17,44,27,37,36,16,19,40,9)
# here i'm selecting a subset of the vals matrix. This picks out the rows with numbers defined in cluster_splits
# in general you can subset any matrix with the form matrix[rows, columns]
# here there is a comma after the rows subset because I want all the columns  
vals.splits <- vals[cluster_splits,]
# select rows that aren't in cluster_splits
# here I'm subsetting on members of the sequence 1:59 that are not in cluster_splits
# i.e. the rows not in the list. 
vals.nosplits <- vals[!(1:59 %in% cluster_splits), ]

# this sets the plotting window to have 2 rows and 2 columns, meaning it can hold 4 plots.
# can change to any number of dimensions
par(mfrow=c(2,2))
# plot a histogram of the p-values from the first column (au values) 
# I'm doing [,1] because we want all rows and only the first column. 
# subtrating from 1 so things are more interpretable as p-values
hist(1-vals.splits[, 1])
# plot a histogram of the splits not defining a cluster
hist(1-vals.nosplits[, 1])
# plot a histogram of all the data
hist(1-vals[, 1])

# use the t.test function to test for differences in the means of the two datasets
# this compares p-values of cluster splits to p-values of the other splits
t.test(vals.splits[,1], vals.nosplits[,1])
# this compares cluster splits to the entire distribution (not sure if this is valid to do)
t.test(vals.splits[,1], vals[,1])

