# TUD variation tests - within cluster variation of certain signals, etc
tud <- as.matrix(read.table('~/GitHub/tango/data/all_phages_TUD.tsv'))
tudf <- as.data.frame(read.table('~/GitHub/tango/data/all_phages_TUD.tsv'))
strsplit("224(E)",split="\\(")[[1]][2]

clusters <- sapply(rownames(tud), function(i) strsplit(x=strsplit(i,split="\\(")[[1]][2], split='\\)')[[1]][1])
cf <- factor(clusters)

tudf["cluster"] <- clusters

cnums <- sapply(levels(cf), function(x) sum(tudf["cluster"]==x))
b1var <- apply(tudf[tudf["cluster"]=="B1",1:256], 2, var)
b1mean <- apply(tudf[tudf["cluster"]=="B1",1:256], 2, mean)
b1median <- apply(tudf[tudf["cluster"]=="B1",1:256], 2, median)

a1var <- apply(tudf[tudf["cluster"]=="A1",1:256], 2, var)
a1mean <- apply(tudf[tudf["cluster"]=="A1",1:256], 2, mean)
a1median <- apply(tudf[tudf["cluster"]=="A1",1:256], 2, median)

allvar <- apply(tudf[,1:256], 2, var)

tcor <- cor(tud)
samplecor <- cor(t(tud))
pc <- prcomp(tud)
predict1 <- predict(pc)[,1]
