# ==============================================
# 完整WGCNA分析代码（4个AD亚型特异性模块与Hub基因识别）
# ==============================================
rm(list = ls())
options(stringsAsFactors = FALSE)

# 1. 加载所需包
library(WGCNA)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)  # 人类基因注释库（犬用org.Cf.eg.db）

# 2. 设置工作目录（替换为你的路径）
setwd("/media/nfs/nfs02/wangyi/AD_Human_logFpkmAnd001_New/10WGCNA")

# ==============================================
# 步骤1：数据预处理
# ==============================================
# 2.1 读取表达矩阵
expr_matrix <- fread("/media/nfs/nfs02/wangyi/AD_Human_logFpkmAnd001_New/07AD分组差异分析/AllFPKM/AD_Human_Groups_rownams.xls")
expr_matrix <- as.data.frame(expr_matrix)
expr_matrix <- expr_matrix[,!grepl("^Healthy", colnames(expr_matrix))]  # 剔除健康组（如需保留可注释）

# 2.2 格式调整（基因名设为行名）
rownames(expr_matrix) <- expr_matrix$gene
expr_matrix <- expr_matrix[,-1]

# 2.3 样本分组标注
sample_names <- colnames(expr_matrix)
sample_group <- case_when(
  grepl("AD_Human_Group1", sample_names) ~ "AD_Group1",
  grepl("AD_Human_Group2", sample_names) ~ "AD_Group2",
  grepl("AD_Human_Group3", sample_names) ~ "AD_Group3",
  grepl("AD_Human_Group4", sample_names) ~ "AD_Group4"
)
table(sample_group)  # 检查分组分布

# 2.4 过滤低表达基因（可选，已注释，如需启用取消注释）
# keep_gene <- rowSums(expr_matrix > 0.1) >= ncol(expr_matrix)*0.5
# expr_matrix <- expr_matrix[keep_gene, ]
# cat("过滤后基因数：", nrow(expr_matrix), "\n")

# 2.5 转置矩阵（WGCNA要求：行=样本，列=基因）
expr_matrix_t <- t(expr_matrix)
cat("矩阵维度：", nrow(expr_matrix_t), "样本 ×", ncol(expr_matrix_t), "基因\n")

# ==============================================
# 步骤2：软阈值选择（构建无标度网络）
# ==============================================
enableWGCNAThreads(6)  # 开启多线程（根据服务器配置调整）

# 计算软阈值
powers <- c(1:10, seq(12,20,2))
sft <- pickSoftThreshold(
  data = expr_matrix_t,
  powerVector = powers,
  verbose = 5,
  networkType = "unsigned"
)

# 可视化软阈值（修复β编码问题，用beta）
pdf("1_SoftThreshold_Selection.pdf", width=10, height=5)
par(mfrow=c(1,2))
# 无标度拟合指数图
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (beta)", ylab="Scale Free Topology Fit (R²)",
     type="n", main="Scale Free Topology Fit")
points(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       pch=21, bg="red", cex=2)
abline(h=0.8, col="blue", lty=2)

# 平均连通性图
plot(sft$fitIndices[,1], sft$fitIndices[,4],
     xlab="Soft Threshold (beta)", ylab="Mean Connectivity",
     type="n", main="Mean Connectivity")
points(sft$fitIndices[,1], sft$fitIndices[,4], pch=21, bg="green", cex=2)
dev.off()

# 选择最优软阈值
optimal_beta <- ifelse(max(-sign(sft$fitIndices[,3])*sft$fitIndices[,2])>=0.8,
                       sft$powerEstimate, max(powers))
cat("最优软阈值beta：", optimal_beta, "\n")

# ==============================================
# 步骤3：构建共表达模块
# ==============================================
net <- blockwiseModules(
  datExpr = expr_matrix_t,          # 修复参数名（新版本用datExpr）
  power = optimal_beta,
  minModuleSize = 30,              # 最小模块基因数
  mergeCutHeight = 0.25,           # 模块合并阈值
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = "AD_TOM",
  verbose = 3
)

# 模块颜色转换
module_colors <- labels2colors(net$colors)
module_gene_count <- table(module_colors)
cat("模块统计：\n"); print(module_gene_count)

# 可视化模块聚类树
pdf("2_Module_Dendrogram.pdf", width=12, height=6)
plotDendroAndColors(
  net$dendrograms[[1]],
  colors = module_colors,
  groupLabels = "Modules",
  addGuide = TRUE,
  main = "Gene Dendrogram & Module Colors"
)
dev.off()

# ==============================================
# 步骤4：模块与亚型关联分析
# ==============================================
# 4.1 构建亚型二分类表型矩阵
sample_info <- data.frame(
  SampleID = rownames(expr_matrix_t),
  Group = sample_group,
  AD_Group1 = ifelse(sample_group=="AD_Group1",1,0),
  AD_Group2 = ifelse(sample_group=="AD_Group2",1,0),
  AD_Group3 = ifelse(sample_group=="AD_Group3",1,0),
  AD_Group4 = ifelse(sample_group=="AD_Group4",1,0),
  row.names = rownames(expr_matrix_t)
)

# 4.2 计算模块特征值（ME）
MEs <- moduleEigengenes(expr_matrix_t, module_colors)$eigengenes
MEs <- orderMEs(MEs)

# 4.3 模块-亚型相关性计算
module_subtype_cor <- cor(MEs, sample_info[,4:7], use="p")
module_subtype_p <- corPvalueStudent(module_subtype_cor, nrow(expr_matrix_t))

# 可视化相关性热图
pdf("3_Module_Subtype_Correlation.pdf", width=10, height=8)
labeledHeatmap(
  Matrix = module_subtype_cor,
  xLabels = colnames(module_subtype_cor),
  yLabels = rownames(module_subtype_cor),
  colors = greenWhiteRed(50),
  textMatrix = paste0(signif(module_subtype_cor,2),"\n(",signif(module_subtype_p,1),")"),
  margins = c(10,12),
  main = "Module-Subtype Correlation (r/P)"
)
dev.off()

# 4.4 筛选亚型特异性模块（|r|>0.7, P<0.05）
subtype_modules <- list()
for (subtype in colnames(module_subtype_cor)) {
  sig_mods <- rownames(module_subtype_cor)[abs(module_subtype_cor[,subtype])>0.7 & module_subtype_p[,subtype]<0.05]
  subtype_modules[[subtype]] <- sig_mods
  cat(subtype, "特异性模块：", paste(sig_mods, collapse=", "), "\n")
}

# ==============================================
# 步骤5：筛选Hub基因
# ==============================================
hub_genes <- list()
for (subtype in names(subtype_modules)) {
  sig_mods <- subtype_modules[[subtype]]
  if (length(sig_mods)==0) next
  
  # 基因显著性（GS）：基因与亚型的相关性
  gs <- as.numeric(cor(expr_matrix_t, sample_info[,subtype], use="p"))
  names(gs) <- colnames(expr_matrix_t)
  
  # 提取模块基因
  mod_genes <- lapply(sig_mods, function(m) colnames(expr_matrix_t)[module_colors==m])
  mod_genes <- unique(unlist(mod_genes))
  
  # 模块隶属度（MM）：基因与模块ME的相关性
  mm <- cor(expr_matrix_t[,mod_genes], MEs[,sig_mods])
  
  # 筛选Hub基因（MM>0.8, GS>0.6）
  hub <- c()
  for (m in sig_mods) {
    genes <- colnames(expr_matrix_t)[module_colors==m]
    m_mm <- mm[genes, m]
    m_gs <- gs[genes]
    hub <- c(hub, names(m_mm)[m_mm>0.8 & m_gs>0.6])
  }
  
  hub_genes[[subtype]] <- unique(hub)
  cat(subtype, "Hub基因数：", length(hub_genes[[subtype]]), "\n")
}

# 保存Hub基因
saveRDS(hub_genes, "4_Subtype_HubGenes.rds")
write.table(
  do.call(rbind, lapply(names(hub_genes), function(x) data.frame(Subtype=x, HubGene=hub_genes[[x]]))),
  "4_Subtype_HubGenes.txt", sep="\t", quote=F, row.names=F
)

# ==============================================
# 步骤6：Hub基因功能富集
# ==============================================
enrich_func <- function(genes) {
  # 基因名转换为ENTREZID
  entrez <- mapIds(org.Hs.eg.db, keys=genes, keytype="SYMBOL", column="ENTREZID", multiVals="first")
  entrez <- na.omit(entrez)
  
  if (length(entrez)<5) return(NULL)
  
  # KEGG富集
  kegg <- enrichKEGG(gene=entrez, organism="hsa", pvalueCutoff=0.05)
  # GO-BP富集
  go <- enrichGO(gene=entrez, OrgDb=org.Hs.eg.db, ont="BP", pvalueCutoff=0.05, readable=T)
  
  return(list(KEGG=kegg, GO=go))
}

# 批量富集分析
enrich_res <- lapply(hub_genes, enrich_func)
saveRDS(enrich_res, "5_HubGene_Enrichment.rds")

# 可视化AD_Group1富集结果（示例）
if (!is.null(enrich_res[["AD_Group1"]]$KEGG)) {
  pdf("5_AD_Group1_KEGG_Enrichment.pdf", width=12, height=8)
  dotplot(enrich_res[["AD_Group1"]]$KEGG, showCategory=10, title="AD_Group1 Hub Gene KEGG Enrichment")
  dev.off()
}

# ==============================================
# 步骤7：结果汇总输出
# ==============================================
# 模块-亚型关联汇总表
module_summary <- data.frame(
  Module = rownames(module_subtype_cor),
  AD_Group1_r = module_subtype_cor[,"AD_Group1"],
  AD_Group1_p = module_subtype_p[,"AD_Group1"],
  AD_Group2_r = module_subtype_cor[,"AD_Group2"],
  AD_Group2_p = module_subtype_p[,"AD_Group2"],
  AD_Group3_r = module_subtype_cor[,"AD_Group3"],
  AD_Group3_p = module_subtype_p[,"AD_Group3"],
  AD_Group4_r = module_subtype_cor[,"AD_Group4"],
  AD_Group4_p = module_subtype_p[,"AD_Group4"],
  Gene_Count = module_gene_count[rownames(module_subtype_cor)]
)
write.table(module_summary, "6_Module_Summary.txt", sep="\t", quote=F, row.names=F)

# 完成提示
cat("\n========== WGCNA分析完成！核心输出文件：==========\n")
cat("1. 1_SoftThreshold_Selection.pdf：软阈值选择图\n")
cat("2. 2_Module_Dendrogram.pdf：模块聚类树\n")
cat("3. 3_Module_Subtype_Correlation.pdf：模块-亚型相关性热图\n")
cat("4. 4_Subtype_HubGenes.txt：各亚型Hub基因列表\n")
cat("5. 5_HubGene_Enrichment.rds：Hub基因功能富集结果\n")
cat("6. 6_Module_Summary.txt：模块关联汇总表\n")