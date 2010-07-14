###################################################
### chunk number 1: 
###################################################
options(width=60)
ps.options(family="sans")


###################################################
### chunk number 2: 
###################################################
library(BioNet)
library(DLBCL)
data(dataLym)
data(interactome)


###################################################
### chunk number 3: 
###################################################
pvals <- cbind(t=dataLym$t.pval, s=dataLym$s.pval)
rownames(pvals) <- dataLym$label
pval <- aggrPvals(pvals, order=2, plot=FALSE)


###################################################
### chunk number 4: 
###################################################
subnet <- subNetwork(dataLym$label, interactome)
subnet <- rmSelfLoops(subnet)
subnet


###################################################
### chunk number 5: 
###################################################
fb <- fitBumModel(pval, plot=FALSE)
scores <- scoreNodes(subnet, fb, fdr=0.001)


###################################################
### chunk number 6: 
###################################################
module <- runFastHeinz(subnet, scores)
logFC <- dataLym$diff
names(logFC) <- dataLym$label


###################################################
### chunk number 7: 
###################################################
plotModule(module, scores=scores, diff.expr=logFC)


###################################################
### chunk number 8: 
###################################################
library(BioNet)
library(DLBCL)
data(exprLym)
data(interactome)


###################################################
### chunk number 9: 
###################################################
exprLym


###################################################
### chunk number 10: 
###################################################
interactome


###################################################
### chunk number 11: 
###################################################
network <- subNetwork(featureNames(exprLym), interactome)
network


###################################################
### chunk number 12: 
###################################################
network <- largestComp(network)
network


###################################################
### chunk number 13: 
###################################################
library(genefilter)
library(impute)
expressions <- impute.knn(exprs(exprLym))$data
t.test <- rowttests(expressions, fac=exprLym$Subgroup)    


###################################################
### chunk number 14: 
###################################################
t.test[1:10, ]


###################################################
### chunk number 15: 
###################################################
library(xtable)
top.table <- xtable(t.test[1:10,], display=c("s", "f", "f", "f"))
print(top.table, floating=FALSE)                                    


###################################################
### chunk number 16: 
###################################################
data(dataLym)
ttest.pval <- t.test[, "p.value"]
surv.pval <- dataLym$s.pval
names(surv.pval) <- dataLym$label
pvals <- cbind(ttest.pval, surv.pval)


###################################################
### chunk number 17: 
###################################################
pval <- aggrPvals(pvals, order=2, plot=FALSE)


###################################################
### chunk number 18: 
###################################################
fb <- fitBumModel(pval, plot=FALSE)
fb   


###################################################
### chunk number 19: 
###################################################
dev.new(width=13, height=7)
par(mfrow=c(1,2))
hist(fb)
plot(fb)

dev.off()


###################################################
### chunk number 20: 
###################################################
plotLLSurface(pval, fb)


###################################################
### chunk number 21: 
###################################################
scores <- scoreNodes(network=network, fb=fb, fdr=0.001)


###################################################
### chunk number 22: 
###################################################
network <- rmSelfLoops(network)                                                                                                     
writeHeinzEdges(network=network, file="lymphoma_edges_001", use.score=FALSE)
writeHeinzNodes(network=network, file="lymphoma_nodes_001", node.scores = scores)


###################################################
### chunk number 23: 
###################################################
datadir <- file.path(.path.package("BioNet"), "extdata") 
dir(datadir)


###################################################
### chunk number 24: 
###################################################
module <- readHeinzGraph(node.file=file.path(datadir, "lymphoma_nodes_001.txt.0.hnz"), network=network)
diff <- t.test[, "dm"]
names(diff) <- rownames(t.test) 


###################################################
### chunk number 25: 
###################################################
plotModule(module, diff.expr=diff, scores=scores)


###################################################
### chunk number 26: 
###################################################
plot3dModule(module, windowSize = c(100,100,1500,1000), diff.or.score=scores)


###################################################
### chunk number 27: 
###################################################
sum(scores[nodes(module)])
sum(scores[nodes(module)]>0)
sum(scores[nodes(module)]<0)


###################################################
### chunk number 28: 
###################################################
library(BioNet)
library(DLBCL)
library(ALL)
data(ALL)
data(interactome)


###################################################
### chunk number 29: 
###################################################
ALL


###################################################
### chunk number 30: 
###################################################
interactome


###################################################
### chunk number 31: 
###################################################
mapped.eset = mapByVar(ALL, network=interactome, attr="geneID")
mapped.eset[1:5,1:5]


###################################################
### chunk number 32: 
###################################################
length(intersect(rownames(mapped.eset), nodes(interactome)))


###################################################
### chunk number 33: 
###################################################
network = subNetwork(rownames(mapped.eset), interactome)
network
network <- largestComp(network)
network <- rmSelfLoops(network)
network


###################################################
### chunk number 34: 
###################################################
library(limma)
design= model.matrix(~ -1+ factor(c(substr(unlist(ALL$BT), 0, 1))))
colnames(design)= c("B", "T")
contrast.matrix <- makeContrasts(B-T, levels=design)
contrast.matrix 
fit <- lmFit(mapped.eset, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)


###################################################
### chunk number 35: 
###################################################
pval = fit2$p.value[,1]                                


###################################################
### chunk number 36: 
###################################################
fb <- fitBumModel(pval, plot=FALSE)
fb


###################################################
### chunk number 37: 
###################################################
dev.new(width=13, height=7)
par(mfrow=c(1,2))
hist(fb)
plot(fb)


###################################################
### chunk number 38: 
###################################################
scores <- scoreNodes(network=network, fb=fb, fdr=1e-14)


###################################################
### chunk number 39: 
###################################################
writeHeinzEdges(network=network, file="ALL_edges_001", use.score=FALSE)
writeHeinzNodes(network=network, file="ALL_nodes_001", node.scores = scores)


###################################################
### chunk number 40: 
###################################################
datadir <- file.path(.path.package("BioNet"), "extdata") 
module <- readHeinzGraph(node.file=file.path(datadir, "ALL_nodes_001.txt.0.hnz"), network=network)


###################################################
### chunk number 41: 
###################################################
nodeDataDefaults(module, attr="diff") <- ""
nodeData(module, n=nodes(module), attr="diff") <- fit2$coefficients[nodes(module),1]
nodeDataDefaults(module, attr="score") <- ""
nodeData(module, n=nodes(module), attr="score") <- scores[nodes(module)]
nodeData(module)[1]


###################################################
### chunk number 42: 
###################################################
saveNetwork(module, file="ALL_module", type="XGMML")


