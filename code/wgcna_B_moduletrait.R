# 打印全部模块-性状相关（复用 wgcna_B.R 的参数重算）
suppressMessages(library(WGCNA))
options(stringsAsFactors = FALSE)
allowWGCNAThreads()
dd <- "C:/Users/98757/Desktop/生信分析/data/r_crosscheck"
expr <- read.csv(file.path(dd, "gse83453_expr.csv"), row.names = 1, check.names = FALSE)
grp <- read.csv(file.path(dd, "gse83453_groups.csv"), row.names = 1)
datExpr <- as.data.frame(t(expr))[rownames(grp), ]
trait <- as.integer(grp$group == "Calcified")
madv <- apply(datExpr, 2, mad)
sel <- union(names(sort(madv, decreasing = TRUE))[1:5000], c("CD74", "MMP9"))
datExpr <- datExpr[, sel]
net <- blockwiseModules(datExpr, power = 12, TOMType = "unsigned", minModuleSize = 30,
                        mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, verbose = 0)
modColors <- labels2colors(net$colors)
MEs <- orderMEs(moduleEigengenes(datExpr, colors = modColors)$eigengenes)
mtc <- cor(MEs, trait, use = "p")
mtp <- corPvalueStudent(mtc, nrow(datExpr))
res <- data.frame(module = gsub("ME", "", rownames(mtc)), cor = round(mtc[, 1], 3), P = signif(mtp[, 1], 2))
res <- res[order(-abs(res$cor)), ]
print(res, row.names = FALSE)
write.csv(res, "C:/Users/98757/Desktop/生信分析/论文_选题B/TableS_WGCNA_module_trait.csv", row.names = FALSE)
cat("done\n")
