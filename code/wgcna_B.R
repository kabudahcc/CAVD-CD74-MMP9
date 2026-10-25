# 选题B：WGCNA 共表达网络分析（GSE83453，19 钙化 vs 8 正常）
# 输出：Fig4 面板素材（模块-性状热图、关键模块 GS-MM 散点）+ 模块归属表 + JSON 摘要
suppressMessages(library(WGCNA))
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

dd  <- "C:/Users/98757/Desktop/生信分析/data/r_crosscheck"
out <- "C:/Users/98757/Desktop/生信分析/论文_选题B"

expr <- read.csv(file.path(dd, "gse83453_expr.csv"), row.names = 1, check.names = FALSE)
grp  <- read.csv(file.path(dd, "gse83453_groups.csv"), row.names = 1)
datExpr <- as.data.frame(t(expr))           # 样本 × 基因
datExpr <- datExpr[rownames(grp), ]
trait <- data.frame(Calcified = as.integer(grp$group == "Calcified"))
rownames(trait) <- rownames(grp)

# 取 MAD 前 5000 基因；强制保留 CD74 / MMP9
madv <- apply(datExpr, 2, mad)
top  <- names(sort(madv, decreasing = TRUE))[1:min(5000, ncol(datExpr))]
sel  <- union(top, intersect(c("CD74", "MMP9"), colnames(datExpr)))
datExpr <- datExpr[, sel]

# 软阈值
powers <- c(1:10, seq(12, 20, 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 0)
png(file.path(out, "Fig4_WGCNA_softthreshold.png"), width = 1600, height = 800, res = 200)
par(mfrow = c(1, 2))
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft threshold (power)", ylab = "Scale-free topology fit (signed R^2)", type = "n",
     main = "Scale independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], labels = powers, col = "#C0504D")
abline(h = 0.80, col = "grey", lty = 2)
plot(sft$fitIndices[, 1], sft$fitIndices[, 5], xlab = "Soft threshold (power)",
     ylab = "Mean connectivity", type = "n", main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, col = "#C0504D")
dev.off()

softPower <- sft$powerEstimate
if (is.na(softPower) || softPower < 3) {
  cand <- sft$fitIndices[, 1][-sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2] >= 0.8]
  softPower <- if (length(cand) > 0) min(cand) else 6
}
cat("softPower =", softPower, "\n")

net <- blockwiseModules(datExpr, power = softPower, TOMType = "unsigned",
                        minModuleSize = 30, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE, verbose = 0)
modColors <- labels2colors(net$colors)
names(modColors) <- colnames(datExpr)
MEs <- orderMEs(moduleEigengenes(datExpr, colors = modColors)$eigengenes)

# 模块-性状相关
modTraitCor    <- cor(MEs, trait, use = "p")
modTraitPvalue <- corPvalueStudent(modTraitCor, nrow(datExpr))
txt <- matrix(paste0(sprintf("%.2f", modTraitCor), "\n(", signif(modTraitPvalue, 1), ")"),
              nrow = nrow(modTraitCor))
png(file.path(out, "Fig4_WGCNA_module_trait.png"), width = 1200, height = 1800, res = 200)
par(mar = c(6, 9, 3, 1))
labeledHeatmap(Matrix = modTraitCor, xLabels = "Calcified", yLabels = names(MEs),
               ySymbols = names(MEs), colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = txt, setStdMargins = FALSE, cex.text = 0.8, cex.lab.y = 0.75,
               zlim = c(-1, 1), main = "Module–trait relationships")
dev.off()

# 关键模块 = 与钙化相关性最高且 P<0.05 的模块
ord <- order(-abs(modTraitCor[, 1]))
keyIdx <- ord[modTraitPvalue[ord, 1] < 0.05][1]
keyME  <- names(MEs)[keyIdx]
keyCol <- gsub("ME", "", keyME)
cat("key module:", keyCol, "| cor =", round(modTraitCor[keyIdx, 1], 3),
    "| P =", signif(modTraitPvalue[keyIdx, 1], 3), "\n")

# GS / MM
GS <- cor(datExpr, trait$Calcified, use = "p")[, 1]
MM <- cor(datExpr, MEs, use = "p")
inKey <- modColors == keyCol
png(file.path(out, "Fig4_WGCNA_GS_MM.png"), width = 1100, height = 1100, res = 200)
plot(abs(MM[inKey, keyME]), abs(GS[inKey]), pch = 16, col = "#78909C",
     xlab = paste0("Module membership (", keyCol, ")"), ylab = "Gene significance for calcification",
     main = paste0(keyCol, " module: GS vs MM (r = ",
                   sprintf("%.2f", cor(abs(MM[inKey, keyME]), abs(GS[inKey]))), ")"))
for (g in c("CD74", "MMP9")) {
  if (g %in% colnames(datExpr) && modColors[colnames(datExpr) == g] == keyCol) {
    points(abs(MM[g, keyME]), abs(GS[g]), pch = 17, col = "#C0504D", cex = 1.6)
    text(abs(MM[g, keyME]), abs(GS[g]), g, pos = 3, col = "#C0504D", font = 2)
  }
}
dev.off()

# 模块归属表
memb <- data.frame(gene = colnames(datExpr), module = modColors,
                   MM_key = round(as.numeric(MM[, keyME]), 3),
                   GS = round(as.numeric(GS), 3))
memb <- memb[order(memb$module, -abs(memb$MM_key)), ]
write.csv(memb, file.path(out, "TableS_WGCNA_module_membership.csv"), row.names = FALSE)

res <- list(softPower = softPower, n_genes = ncol(datExpr), n_modules = length(unique(modColors)),
            key_module = keyCol,
            key_module_cor = round(modTraitCor[keyIdx, 1], 3),
            key_module_P = signif(modTraitPvalue[keyIdx, 1], 3),
            module_sizes = as.list(table(modColors)),
            CD74_module = if ("CD74" %in% colnames(datExpr)) modColors["CD74"] else NA,
            MMP9_module = if ("MMP9" %in% colnames(datExpr)) modColors["MMP9"] else NA,
            CD74_GS = if ("CD74" %in% colnames(datExpr)) round(GS["CD74"], 3) else NA,
            MMP9_GS = if ("MMP9" %in% colnames(datExpr)) round(GS["MMP9"], 3) else NA,
            CD74_MM_key = if ("CD74" %in% colnames(datExpr)) round(MM["CD74", keyME], 3) else NA,
            MMP9_MM_key = if ("MMP9" %in% colnames(datExpr)) round(MM["MMP9", keyME], 3) else NA)
jsonlite::write_json(res, file.path(out, "results_wgcna.json"), auto_unbox = TRUE, pretty = TRUE)
print(res)
cat("WGCNA done\n")
