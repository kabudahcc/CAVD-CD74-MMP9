# 选题B 升级9c-1：multiMiR 查询 CD74/MMP9 的 miRNA（miRTarBase 验证 + TargetScan/miRDB 预测）
suppressMessages(library(multiMiR))
genes <- c("CD74", "MMP9")

val <- get_multimir(org = "hsa", target = genes, table = "validated", summary = TRUE)
vd <- val@data
cat("validated rows:", nrow(vd), "\n")

pre <- get_multimir(org = "hsa", target = genes, table = "predicted",
                    summary = TRUE, predicted.cutoff = 20, predicted.cutoff.type = "p",
                    predicted.site = "conserved")
pd <- pre@data
cat("predicted rows:", nrow(pd), "\n")

# 每个基因：验证型 miRNA + 双库(TargetScan+miRDB)共同预测的 miRNA
pick <- function(g) {
  v <- unique(vd$ mature_mirna_id[vd$target_symbol == g])
  p <- pd[pd$target_symbol == g, ]
  cnt <- table(p$mature_mirna_id)
  both_db <- names(cnt)[cnt >= 2]  # targetscan 与 mirdb 同时预测
  list(validated = v, predicted_both = both_db)
}
cd <- pick("CD74"); mm <- pick("MMP9")
cat("CD74 validated:", length(cd$validated), " predicted_both:", length(cd$predicted_both), "\n")
cat("MMP9 validated:", length(mm$validated), " predicted_both:", length(mm$predicted_both), "\n")

shared_val <- intersect(cd$validated, mm$validated)
cd_set <- union(cd$validated, cd$predicted_both)
mm_set <- union(mm$validated, mm$predicted_both)
shared_all <- intersect(cd_set, mm_set)
cat("shared validated:", paste(shared_val, collapse = ", "), "\n")
cat("shared all:", length(shared_all), "\n")

out <- data.frame(
  miRNA = shared_all,
  CD74_validated = shared_all %in% cd$validated,
  MMP9_validated = shared_all %in% mm$validated,
  CD74_predicted = shared_all %in% cd$predicted_both,
  MMP9_predicted = shared_all %in% mm$predicted_both)
write.csv(out, file.path(Sys.getenv("WS", "."), "mirna_shared.csv"), row.names = FALSE)
cat("saved mirna_shared.csv:", nrow(out), "miRNAs\n")
