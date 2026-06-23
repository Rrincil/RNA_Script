library(ggvenn)
library(ggplot2)
library(dplyr)
library(data.table)
rm(list = ls())
setwd("/media/nfs/nfs02/wangyi/AD_Human_Canine/01人犬相似亚型交集基因")
# 1. 读取数据
AD_Human_Group1 <- fread("/media/nfs/nfs02/wangyi/AD_Human_logFpkmAnd001_New/10差异基因与maker基因取交集/AD_Group1_Intersection_1.2FC2.txt") %>% as.data.frame()
AD_Canine_Group1 <- fread("/media/nfs/nfs02/wangyi/AD_Cains_logFpkmAnd001_New/10差异基因与maker基因取交集/AD_Group1_Intersection_1.2FC2.txt") %>% as.data.frame()

# 2. 取基因
AD_Human_gene <- AD_Human_Group1$Gene
AD_Canine_gene <- AD_Canine_Group1$Gene

# 3. 基因列表去重
AD_Human_gene <- unique(as.character(AD_Human_gene))
AD_Canine_gene <- unique(as.character(AD_Canine_gene))
intersection <- intersect(AD_Human_gene, AD_Canine_gene)

# 4. 输出统计
cat(paste0("AD_Human_vs_", "差异基因数：", length(AD_Human_gene), "\n"))
cat(paste0("AD_Canine_vs_", "差异基因数：", length(AD_Canine_gene), "\n"))
cat(paste0("交集基因数：", length(intersection), "\n"))

# 5. 保存结果
write.table(
  data.frame(Gene = intersection, stringsAsFactors = FALSE),
  paste0("AD_Human_Group1_vs_AD_Canine_Group1_Intersection_.txt"),
  row.names = F, col.names = T, quote = F, sep = "\t"
)
cat(paste0("✅ 结果保存：AD_Human_Group1_vs_AD_Canine_Group1_Intersection_.txt\n"))

# 6. 韦恩图（无语法错误）
venn_data <- list(
  "AD_Human_Group1_DEGs",
  "AD_Canine_Group1_DEGs"
)
names(venn_data) <- venn_data
venn_data[[1]] <- AD_Human_gene
venn_data[[2]] <- AD_Canine_gene

p <- ggvenn(
  venn_data,
  columns = names(venn_data),
  show_elements = F,
  show_percentage = T,
  fill_color = c("#E41A1C", "#1E90FF"),
  fill_alpha = 0.5,
  set_name_size = 4,
  text_size = 4
)
ggsave(
  paste0("AD_Human_Group1_vs_AD_Canine_Group1_Intersection_.pdf"),
  plot = p, width = 7, height = 7, dpi = 300
)
cat(paste0("✅ 韦恩图保存：AD_Human_Group1_vs_AD_Canine_Group1_Intersection_.pdf\n")) 

# 预览交集基因
cat("\n交集基因前20个：\n")
print(if (length(intersection) > 0) head(intersection, 20) else "无")


cat("\n=======================================\n")
cat("所有存在的组都处理完成！\n")
cat("=======================================\n")