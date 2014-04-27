# repeat some of the analysis on the data with reverse complements added
tud <- as.data.frame(read.table('~/GitHub/tango/data/all_phages_TUD.tsv'))
tudrc <- as.data.frame(read.table('~/GitHub/tango/data/with_reverse_complement//all_phages_TUD_4_RC.tsv'))
d<- dist(tudrc)
plot(hist(d))
#holy shit it looks even better this way. 

clusters <- sapply(rownames(tudrc), function(i) strsplit(x=strsplit(i,split="\\(")[[1]][2], split='\\)')[[1]][1])
cf <- factor(clusters)

tudrc["cluster"] <- clusters

tudVars <- sort(apply(tudrc[1:256],2,var), decreasing=T)
#TCGA has the most variance and is also very bimodal. Try split on this
TCGA1 <- tudrc[tudrc[,'TCGA']<2.2,"cluster"]
TCGA2 <- tudrc[tudrc[,'TCGA']>2.2,"cluster"]
# Literally a perfect split between clusters.
table(TCGA1)
table(TCGA2)
# only singletons are in both tables
sum(names(table(TCGA1)) %in% names(table(TCGA2))) + sum(names(table(TCGA2)) %in% names(table(TCGA1)))

#B3 cluster can be differentiated just on the basis of TUD in GATC
length(tudrc[tudrc[,'GATC']>3.5,"cluster"])
sum(tudrc[,"cluster"]=="B3")

#Try GAAG
plot(hist(tudrc[,"GAAG"],breaks=50))
# GAAG is characteristic of B1 phage
plot(hist(tudrc[tudrc[,"cluster"]=="B1","GAAG"],breaks=50))
GAAG1 <- tudrc[tudrc[,'GAAG']<0.49,"cluster"]
GAAG2 <- tudrc[tudrc[,'GAAG']>0.49,"cluster"]
#only the Ps contaminate this split, otherwise perfect
table(GAAG1)
table(GAAG2)

# are tetras normally distributed within a given cluster?
# check 
qqnorm(tudrc[tudrc[,"cluster"]=="B1",1])
plot(hist(tudrc[tudrc[,"cluster"]=="B1",1]))

