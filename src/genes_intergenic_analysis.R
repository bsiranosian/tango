setwd("~/GitHub/tango")
gd.3 <- as.matrix(read.table('data/gene_analysis/wlidcat_gene_dev3.txt', sep=','))
gd.4 <- as.matrix(read.table('data/gene_analysis/wlidcat_gene_dev4.txt', sep=','))
gd.3[gd.3==0] <- NA
gd.4[gd.4==0] <- NA

c3 <- cor(t(gd.3))
heatmap(c3, Rowv = NA, Colv=NA)

c4 <- cor(t(gd.4))
  heatmap(c4, Rowv = NA, Colv=NA)
