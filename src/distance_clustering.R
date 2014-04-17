# clustering on TUD distance 
tud <- as.matrix(read.table('~/GitHub/tango/src/hatful60TUD.tsv'))
d <- dist(tud, method='euclidean',)
# try different clustering methods
fit1 <- hclust(d,method='ward')
fit2 <- hclust(d*d,method='ward')
fit3 <- hclust(d,method='complete')
fit4 <- hclust(d,method='single')
par(mfrow=c(2,2))
plot(fit1)
plot(fit2)
plot(fit3)
plot(fit4)

