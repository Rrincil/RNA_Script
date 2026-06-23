setwd("/media/nfs/nfs02/wangyi/AD_Human_Canine/03人犬cellchat/Canine")
rm(list = ls())
# ======================
# 加载包+版本检查
# ======================
library(CellChat)
library(patchwork)
library(dplyr)
library(circlize)
library(RColorBrewer)


cat("=== 版本信息 ===\n")
cellchat_version <- as.character(packageVersion("CellChat"))
cat("CellChat版本：", cellchat_version, "\n")
cat("R版本：", R.version.string, "\n\n")

# ======================
# 全局工具函数（核心修复：解决文件名非法字符问题）
# ======================
# 1. 清理文件名字符（兼容Windows/Linux/Mac）
clean_filename_char <- function(name) {
  # 替换所有非法字符：/ \ : * ? " < > | + 空格 等，统一改为下划线
  clean_name <- gsub("[\\/\\\\:\\*\\?\\\"\\<\\>\\|\\+\\s]", "_", name)
  # 替换连续下划线为单个，移除首尾下划线
  clean_name <- gsub("_+", "_", clean_name)
  clean_name <- gsub("^_|_$", "", clean_name)
  return(clean_name)
}

# 2. 修复高度计算
calc_safe_height <- function(n_sources, multiplier = 0.5) {
  height <- n_sources * 1
  if (height < 3) height <- 3
  return(height)
}

# 3. 通路名校验+匹配
validate_pathways <- function(input_paths, cellchat_pathways) {
  input_paths_clean <- gsub("-", "_", input_paths) %>% tolower()
  cellchat_paths_clean <- gsub("-", "_", cellchat_pathways) %>% tolower()
  matched_idx <- match(input_paths_clean, cellchat_paths_clean)
  matched_paths <- cellchat_pathways[!is.na(matched_idx)]
  
  cat("=== Cellchat中实际的通路名列表 ===\n")
  print(head(cellchat_pathways, 20)) # 只打印前20个避免刷屏
  cat("\n=== 通路名匹配结果 ===\n")
  cat("输入的通路名：", paste(input_paths, collapse = ", "), "\n")
  cat("匹配到的通路名：", paste(matched_paths, collapse = ", "), "\n")
  cat("未匹配到的通路名：", paste(setdiff(input_paths, matched_paths), collapse = ", "), "\n")
  
  return(matched_paths)
}

# ======================
# 关键参数设置
# ======================
group_specific_cells <- list(
  Group1 = c("ACTG2+MC/ACR_SMC","APOD+OxPhos_MC","Bas_Immune","Bas_Quiescent","CASQ2+MD/MC_SMC",
              "CLEC14A+ImmuReg_Endo","COL14A1+EcmMus_Mast","COL6A3+EcmCma_SMC","CTLac","DC3","DUSP+Treg_NK",
              "ENO1+EnergyProc_Mast","HF","IGFBP2+MesenchymalDiff_Fib","KRT14+EpiKrt_Endo","Mac3",
              "MRC1+HighEnergy_Endo","PRG4+ECMReg_Fib","SELE+EndoRegAdhesion_Fib","SEPT11+AutoVesT_MC",
              "Spi_AntiMicro","Tcm","Th17","Th2"),
  # Group1 = c( "CLEC14A+ImmuReg_Endo", "COL6A3+EcmCma_SMC","PRG4+ECMReg_Fib","Tcm", "Bas_Quiescent"),              
  
  Group2 = c("C7+ImmuResp_Fib","CD163_Mac","CD16a+Cyto_NK","COL11A1+CollagenECM_Fib","H2AZ2+ImmuRge-NK",
             "IL1B+ECMReInflam_Fib","KRAS+ImmuSig_Mast","LC1","MKI67+HighProlif_Fib","Mac1",
             "S100A8+MetSt_MC","Th22","Treg","cDC1","cDC2","LC2","B cells","Tex"),
  # Group2 = c("KRAS+ImmuSig_Mast", "LC1", "LC2","MKI67+HighProlif_Fib","S100A8+MetSt_MC","Th22","Tex"),

  Group3 = c("MKI67+CellDiv_Endo","COL23A1+ECMBuild_Fib","GAS6+MelSyn_MC","MHCII+AP_NK",
            "ACTA2+HighEnergyActin_Fib","CTLex","DARS1+RiboTrans_Mast","GFRA2+NeuroAdh_MC","Gra",
            "HDC+ImmuReg_Mast","KRT5+EpiProt_MC","LCE3D+EpiDiffCornif_Fib","Mac2","Mac4",
            "Mac_inf","S100A8+BarrierProtec_NK","ncTRM","cTRM"),

  Group4 = c("NKT","STING1+TransReg_Fib","ASXL1+SplicePig_MC","Bas_Prolif","IRF8+Inflam_NK",
              "LC3","MAG+NeuroMyel_MC","POLR1G+RiboRNAProc_Endo","REL+StressProtReg_Mast",
              "SFN+RSR/HM_SMC","SFRP5+WntRegAdhesion_Fib","SOX10+GlialDiff_Fib",
              "STING1+CSR/AR_SMC","Spi_Early_mid","Spi_mig","Th1","Tmm")
)

immune_cell_types <- c(  
  "COL14A1+EcmMus_Mast", "DARS1+RiboTrans_Mast", "ENO1+EnergyProc_Mast", 
  "HDC+ImmuReg_Mast", "KRAS+ImmuSig_Mast", "REL+StressProtReg_Mast",
  "CD163_Mac","Mac_inf","Mac1","Mac2","Mac3","Mac4",
  "cDC1","cDC2","DC3","LC1","LC2","LC3",
  "B cells", "Th1", "Th17","Th2","Th22","Treg", 
  "CTLac","CTLex","cTRM","ncTRM","NKT", "Tex", "Tcm","Tmm"
)
AD_Canine1_Human1_cell_types <- c( "CLEC14A+ImmuReg_Endo", "COL6A3+EcmCma_SMC","PRG4+ECMReg_Fib","Tcm", "Bas_Quiescent")
AD_Canine2_Human3_cell_types <- c("KRAS+ImmuSig_Mast", "LC1", "LC2","MKI67+HighProlif_Fib","S100A8+MetSt_MC","Th22","Tex")

min_val <- 1e-6  
main_result_dir <- "Immune_Cell_Communication_Results"
if (!dir.exists(main_result_dir)) dir.create(main_result_dir, recursive = TRUE)

# ======================
# 1. 数据准备+LR基因预提取
# ======================
cat("=== 数据准备 ===\n")
raw_expression <- read.csv("Deconv_Predicted_CellType_Expression.csv", row.names = 1, check.names = FALSE)
if (grepl("cell|Cell|TYPE|Type", toString(rownames(raw_expression)[1])) || 
    grepl("gene|Gene", toString(colnames(raw_expression)[1]))) {
  raw_expression <- t(raw_expression)
  cat("✓ 已转置表达矩阵（基因在行，细胞在列）\n")
}
expr_avg <- as.matrix(raw_expression)

# 处理负值+NA+Inf
expr_avg[expr_avg < 0] <- min_val
expr_avg[is.na(expr_avg) | is.infinite(expr_avg)] <- min_val

all_cell_types <- colnames(expr_avg)
gene_names <- rownames(expr_avg)
cat(sprintf("✓ 原始表达矩阵（清理后）：%d基因 × %d细胞\n", nrow(expr_avg), ncol(expr_avg)))

# 提取LR基因
lr_db <- CellChatDB.human$interaction
lr_genes_all <- unique(c(lr_db$ligand, lr_db$receptor))
lr_genes_matched <- intersect(lr_genes_all, gene_names)
cat(sprintf("✓ CellChatDB中匹配的LR基因数：%d/%d\n", length(lr_genes_matched), length(lr_genes_all)))

# 保留LR基因+填充极小值
if (length(lr_genes_matched) > 0) {
  expr_avg <- expr_avg[lr_genes_matched, ]
  expr_avg[expr_avg == 0] <- min_val
  cat(sprintf("✓ 保留LR基因并填充极小值后矩阵：%d基因 × %d细胞\n", nrow(expr_avg), ncol(expr_avg)))
} else {
  stop("❌ 无LR基因匹配，请检查基因名格式")
}

# 读取比例矩阵
ad_proportions <- read.csv("AD_Sample_CellType_Proportions.csv", row.names = 1, check.names = FALSE)
ad_proportions[is.na(ad_proportions) | ad_proportions == 0] <- min_val
cat(sprintf("✓ 比例矩阵（清理后）：%d样本 × %d细胞\n", nrow(ad_proportions), ncol(ad_proportions)))

# 样本分组信息
sample_groups <- c(
  # 第1组
  "AD-24"=1,"AD-30"=1,"AD-29"=1,"AD-12"=1,"AD-11"=1,"AD-7"=1,"AD-10"=1,"AD-9"=1,"AD-27"=1,
  # 第2组
  "AD-5"=2,"AD-18"=2,"AD-28"=2,"AD-21"=2,"AD-14"=2,"AD-16"=2,"AD-20"=2,"AD-23"=2,
  # 第3组
  "AD-1"=3,"AD-26"=3,"AD-19"=3,"AD-17"=3,"AD-6"=3,"AD-2"=3,"AD-4"=3,
  # 第4组
  "AD-8"=4,"AD-25"=4,"AD-22"=4,"AD-15"=4,"AD-3"=4,"AD-13"=4
)

# ======================
# 2. 按组计算加权表达矩阵
# ======================
cat("\n=== 计算加权表达矩阵（仅LR基因+组专属细胞）===\n")
group_weighted_expr <- list()
group_avg_props <- list()
group_matched_cells <- list()

for (g in 1:4) {
  group_name <- paste0("Group", g)
  target_cells <- group_specific_cells[[group_name]]
  group_samples <- names(sample_groups)[sample_groups == g]
  group_samples <- intersect(group_samples, rownames(ad_proportions))
  cat(sprintf("[%s] 样本数：%d\n", group_name, length(group_samples)))
  
  # 计算平均比例
  group_prop_all <- colMeans(ad_proportions[group_samples, , drop = FALSE], na.rm = TRUE)
  group_prop_all[is.na(group_prop_all)] <- min_val
  
  # 筛选专属细胞
  matched_cells <- intersect(target_cells, names(group_prop_all))
  unmatched_cells <- setdiff(target_cells, names(group_prop_all))
  group_matched_cells[[group_name]] <- matched_cells
  cat(sprintf("[%s] 预设细胞数：%d，实际匹配数：%d\n", group_name, length(target_cells), length(matched_cells)))
  
  if (length(unmatched_cells) > 0) {
    warning(sprintf("[%s] 未匹配的细胞：%s", group_name, paste(unmatched_cells[1:5], collapse = ", "))) # 只打印前5个
  }
  if (length(matched_cells) < 2) {
    warning(sprintf("[%s] 有效细胞数<2，无法计算通讯，跳过该组！", group_name))
    next
  }
  
  # 筛选比例+表达矩阵
  group_prop <- group_prop_all[matched_cells]
  expr_group <- expr_avg[, matched_cells, drop = FALSE]
  
  # 批量计算加权表达
  weighted_expr <- expr_group * matrix(rep(group_prop, nrow(expr_group)), nrow = nrow(expr_group), byrow = TRUE)
  # 清理加权矩阵
  weighted_expr[weighted_expr < 0] <- min_val
  weighted_expr[is.na(weighted_expr) | is.infinite(weighted_expr)] <- min_val
  
  # 存储结果
  group_weighted_expr[[group_name]] <- weighted_expr
  group_avg_props[[group_name]] <- group_prop
  
  # 打印预览
  cat(sprintf("[%s] 加权矩阵：%d LR基因×%d专属细胞\n", group_name, nrow(weighted_expr), ncol(weighted_expr)))
  if (nrow(weighted_expr) >= 3) {
    for (gene in rownames(weighted_expr)[1:3]) {
      cat(sprintf("[%s] %s：加权后均值=%.6f\n", group_name, gene, mean(weighted_expr[gene, ])))
    }
  }
}

# 过滤无效组
valid_groups <- names(group_weighted_expr)[sapply(group_weighted_expr, function(x) !is.null(x) && ncol(x)>=2)]
if (length(valid_groups) == 0) stop("❌ 所有组的有效细胞数均<2，无法进行后续分析！")
cat("\n✅ 有效分析组：", paste(valid_groups, collapse = ", "), "\n")

# ======================
# 3. 免疫细胞匹配检查
# ======================
cat("\n=== 免疫细胞匹配检查（按组）===\n")
group_immune_matched <- list()
for (group_name in valid_groups) {
  matched_cells <- group_matched_cells[[group_name]]
  immune_matched <- intersect(matched_cells, immune_cell_types)
  group_immune_matched[[group_name]] <- immune_matched
  cat(sprintf("[%s] 组专属细胞中免疫细胞数：%d/%d\n", group_name, length(immune_matched), length(matched_cells)))
  if (length(immune_matched) == 0) {
    warning(sprintf("[%s] 无匹配的免疫细胞，后续bubble图可能无法生成", group_name))
  }
}

# ======================
# 4. 核心分析函数
# ======================
analyze_group_full_cellchat <- function(group_name, weighted_expr, matched_cells, group_prop, immune_matched) {
  # cat(sprintf("\n=== 处理%s ===\n", group_name))
  result_dir <- "Immune_Cell_Communication_Results"
  group_name <- "Group2"
  weighted_expr <- group_weighted_expr[[group_name]]
  group_prop <- group_avg_props[[group_name]]
  cat(sprintf("\n=== 处理%s ===\n", group_name))
  
  cat(sprintf("\n=== 处理%s ===\n", group_name))
  # 子结果目录（按分组创建，避免文件混乱）
  group_subdir <- file.path(main_result_dir, group_name)
  if (!dir.exists(group_subdir)) dir.create(group_subdir, recursive = TRUE)
  
  # 清理加权矩阵
  weighted_expr_clean <- weighted_expr
  weighted_expr_clean[weighted_expr_clean < 0] <- min_val
  weighted_expr_clean[is.na(weighted_expr_clean) | is.infinite(weighted_expr_clean)] <- min_val
  weighted_expr_clean[weighted_expr_clean == 0] <- min_val
  cat("✓ 加权矩阵清理完成，关键统计：\n")
  cat(sprintf("  最小值：%.6f，最大值：%.6f，是否有负值：%s\n", 
              min(weighted_expr_clean), max(weighted_expr_clean), any(weighted_expr_clean < 0)))

  # 未过滤的细胞列表
  current_cells <- colnames(weighted_expr_clean)
  cat(sprintf("✓ 原始矩阵（未过滤细胞）：%d基因 × %d专属细胞\n", nrow(weighted_expr_clean), ncol(weighted_expr_clean)))

  # 创建CellChat对象
  cat("✓ 创建CellChat对象...\n")
  meta <- data.frame(
    cellType = current_cells,
    size = group_prop[current_cells],
    is_immune = ifelse(current_cells %in% immune_cell_types, "Immune", "Non-Immune"),
    is_AD_Canine1_Human1_cell_types = ifelse(current_cells %in% AD_Canine1_Human1_cell_types, "AD_Canine1_Human1", "Non-AD_Canine1_Human1"),
    is_AD_Canine2_Human3_cell_types = ifelse(current_cells %in% AD_Canine2_Human3_cell_types, "AD_Canine2_Human3", "Non-AD_Canine2_Human3"),
    row.names = current_cells,
    check.names = FALSE
  )
  cellchat <- createCellChat(
    object = weighted_expr_clean, 
    meta = meta, 
    group.by = "cellType",
    do.sparse = FALSE
  )
  
  # 手动设置CellChatDB
  cellchat@DB <- CellChatDB.human
  cellchat <- subsetData(cellchat, features = rownames(weighted_expr_clean))
  cat(sprintf("✓ CellChat对象：%d专属细胞，%d LR基因\n", length(cellchat@idents), nrow(cellchat@data.signaling)))
  
  # 筛选高表达基因
  cat("✓ 筛选高表达基因...\n")
  cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = NULL)
  cellchat@var.features$features <- unique(c(cellchat@var.features$features, rownames(weighted_expr_clean)))
  cat(sprintf("✓ 强制保留LR基因后特征数：%d\n", length(cellchat@var.features$features)))
  
  # 识别LR互作
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cat(sprintf("✓ 识别到LR互作数：%d\n", nrow(cellchat@LR$LR)))
  cellchat@LR$LRsig <- cellchat@LR$LR
  
  # 计算通讯概率
  cat("✓ 计算通讯概率...\n")
  cellchat <- computeCommunProb(
    cellchat,
    population.size = FALSE,
    type = "truncatedMean",
    trim = 0.01,
    nboot = 5
  )
  
  # 通路分析和网络聚合
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  total_edges <- sum(cellchat@net$weight > 0, na.rm = TRUE)
  cat(sprintf("✓ 通路数：%d，总通讯边数：%d\n", length(cellchat@netP$pathways), total_edges))
  
  if (total_edges == 0) {
    cat("✗ 无有效通讯边，跳过可视化\n")
    return(cellchat)
  }
  
  # ========== 1. 网络图绘制 ==========
  # 全部细胞网络图
  all_cells <- rownames(cellchat@net$weight)
  if (length(all_cells) >= 2) {
    net_all <- cellchat@net$weight[all_cells, all_cells, drop = FALSE]
    net_all[is.na(net_all)] <- 0
    if (sum(net_all > 0) > 0) {
      pdf(file.path(group_subdir, paste0(group_name, "_All_Cells_Network.pdf")), width = 8, height = 7)
      netVisual_circle(
        net_all,
        vertex.weight = group_prop[all_cells],
        weight.scale = TRUE,
        label.edge = FALSE,
        vertex.label.cex = 0.55,
        edge.width.max = 3,
        title.name = paste0(group_name, " 全部细胞通讯网络")
      )
      dev.off()
      cat("✓ 全部细胞网络图已保存\n")
    }
  }

  # 组专属细胞网络图
  pdf(file.path(group_subdir, paste0(group_name, "_Specific_Cells_Network.pdf")), width = 16, height = 16)
  netVisual_circle(
    cellchat@net$weight,
    vertex.weight = group_prop[current_cells],
    weight.scale = TRUE,
    label.edge = FALSE,
    vertex.label.cex = 0.8,
    title.name = paste0(group_name, " 专属细胞通讯网络")
  )
  dev.off()
  cat("✓ 专属细胞网络图已保存\n")
  
  # 主要细胞网络图
  cell_comm_strength <- rowSums(cellchat@net$weight, na.rm = TRUE) + colSums(cellchat@net$weight, na.rm = TRUE)
  main_cells <- names(sort(cell_comm_strength, decreasing = TRUE)[1:min(15, length(cell_comm_strength))])
  if (length(main_cells) >= 2) {
    net_main <- cellchat@net$weight[main_cells, main_cells, drop = FALSE]
    net_main[is.na(net_main)] <- 0
    pdf(file.path(group_subdir, paste0(group_name, "_Main_Specific_Cells_Network.pdf")), width = 14, height = 14)
    netVisual_circle(
      net_main,
      vertex.weight = group_prop[main_cells],
      weight.scale = TRUE,
      label.edge = FALSE,
      vertex.label.cex = 1.0,
      title.name = paste0(group_name, " 主要专属细胞通讯网络")
    )
    dev.off()
    cat("✓ 主要专属细胞网络图已保存\n")
  }
  
  # 免疫细胞网络图
  if (length(immune_matched) >= 2) {
    existing_immune <- intersect(immune_matched, rownames(cellchat@net$weight))
    if (length(existing_immune) >= 2) {
      net_immune <- cellchat@net$weight[existing_immune, existing_immune, drop = FALSE]
      net_immune[is.na(net_immune)] <- 0
      pdf(file.path(group_subdir, paste0(group_name, "_Immune_Specific_Cells_Network.pdf")), width = 14, height = 14)
      netVisual_circle(
        net_immune,
        vertex.weight = group_prop[existing_immune],
        weight.scale = TRUE,
        label.edge = FALSE,
        vertex.label.cex = 0.9,
        title.name = paste0(group_name, " 专属免疫细胞通讯网络")
      )
      dev.off()
      cat("✓ 专属免疫细胞网络图已保存\n")
    } else {
      cat("✗ 专属细胞中无足够免疫细胞，跳过免疫网络图\n")
    }
  }
  # AD_Canine1_Human1 网络图
  if (length(intersect(current_cells, AD_Canine1_Human1_cell_types)) >= 2) {
    ad1_cells <- intersect(AD_Canine1_Human1_cell_types, current_cells)
    net_ad1 <- cellchat@net$weight[ad1_cells, ad1_cells, drop = FALSE]
    net_ad1[is.na(net_ad1)] <- 0
    if (sum(net_ad1 > 0) > 0) {
      pdf(file.path(group_subdir, paste0(group_name, "_AD_Canine1_Human1_Network.pdf")), width = 5, height = 5)
      netVisual_circle(
        net_ad1,
        # top = 10,
        # edge.width.max = 3,
        vertex.weight = group_prop[ad1_cells],
        weight.scale = TRUE,
        label.edge = TRUE,
        vertex.label.cex = 0.7,
        title.name = paste0(group_name, " AD_Canine1_Human1 细胞通讯网络")
      )
      dev.off()
      cat("✓ AD_Canine1_Human1 网络图已保存\n")
    }
  }
  # AD_Canine2_Human3 网络图
  if (length(intersect(current_cells, AD_Canine2_Human3_cell_types)) >= 2) {
    ad2_cells <- intersect(AD_Canine2_Human3_cell_types, current_cells)
    net_ad2 <- cellchat@net$weight[ad2_cells, ad2_cells, drop = FALSE]
    net_ad2[is.na(net_ad2)] <- 0
    if (sum(net_ad2 > 0) > 0) {
      pdf(file.path(group_subdir, paste0(group_name, "_AD_Canine2_Human3_Network.pdf")), width = 6, height = 6)
      netVisual_circle(
        net_ad2,
        vertex.weight = group_prop[ad2_cells],
        weight.scale = TRUE,
        label.edge = TRUE,
        vertex.label.cex = 0.7,
        title.name = paste0(group_name, " AD_Canine2_Human3 细胞通讯网络")
      )
      dev.off()
      cat("✓ AD_Canine2_Human3 网络图已保存\n")
    }
  }
  

  # 保存通讯数据
  df.net <- subsetCommunication(cellchat)
  write.csv(df.net, file.path(group_subdir, paste0(group_name, "_net_lr.csv")), quote = FALSE)
  df.netp <- subsetCommunication(cellchat, slot.name = "netP")
  write.csv(df.netp, file.path(group_subdir, paste0(group_name, "_net_pathway.csv")), quote = FALSE)

  # ========== 从df.net中筛选AD_Canine1_Human1相关互作 ==========
  # ========== 修正：严格筛选AD_Canine1_Human1细胞间的互作（仅source+target都在列表内） ==========
  if(group_name == "Group1"){
    cat(sprintf("\n=== 从df.net中筛选%s AD_Canine1_Human1细胞间的互作 ===\n", group_name))
    print(AD_Canine1_Human1_cell_types)
    ad1_cells <- AD_Canine1_Human1_cell_types
    if (length(ad1_cells) == 0) {
      cat("✗ 当前组无AD_Canine1_Human1细胞，跳过筛选\n")
    } else {
      # 关键修改：将 | 改为 &，仅保留source和target都在列表内的互作
      df.net_ad1 <- df.net %>%
        dplyr::filter(source %in% ad1_cells & target %in% ad1_cells)

      cat(sprintf("✓ 共筛选出 %d 条AD_Canine1_Human1细胞间的互作\n", nrow(df.net_ad1)))
      write.csv(df.net_ad1,
                file.path(group_subdir, paste0(group_name, "_AD_Canine1_Human1_net_lr.csv")),
                row.names = FALSE, quote = FALSE)

      # 打印筛选结果（仅列表内细胞）
      cat("筛选后source唯一值（仅AD列表内）：\n")
      print(unique(df.net_ad1$source))
      cat("筛选后target唯一值（仅AD列表内）：\n")
      print(unique(df.net_ad1$target))
    }
    df.net <- df.net_ad1
  }else if(group_name == "Group2"){
    cat(sprintf("\n=== 从df.net中筛选%s AD_Canine2_Human3细胞间的互作 ===\n", group_name))
    print(AD_Canine2_Human3_cell_types)
    ad2_cells <- AD_Canine2_Human3_cell_types
    if (length(ad2_cells) == 0) {
      cat("✗ 当前组无AD_Canine2_Human3细胞，跳过筛选\n")
    } else {
      # 关键修改：将 | 改为 &，仅保留source和target都在列表内的互作
      df.net_ad2 <- df.net %>%
        dplyr::filter(source %in% ad2_cells & target %in% ad2_cells)

      cat(sprintf("✓ 共筛选出 %d 条AD_Canine2_Human3细胞间的互作\n", nrow(df.net_ad2)))
      write.csv(df.net_ad2,
                file.path(group_subdir, paste0(group_name, "_AD_Canine2_Human3_net_lr.csv")),
                row.names = FALSE, quote = FALSE)

      # 打印筛选结果（仅列表内细胞）
      cat("筛选后source唯一值（仅AD列表内）：\n")
      print(unique(df.net_ad2$source))
      cat("筛选后target唯一值（仅AD列表内）：\n")
      print(unique(df.net_ad2$target))  
    }
    df.net <- df.net_ad2
  }

  # ========== 筛选Top5细胞（通讯强度） ==========
  cat("✓ 筛选通讯强度最强的前5个细胞类型...\n")
  source_prob <- df.net %>% 
    dplyr::group_by(cell = source) %>% 
    dplyr::summarise(source_prob = sum(prob), .groups = "drop")

  target_prob <- df.net %>% 
    dplyr::group_by(cell = target) %>% 
    dplyr::summarise(target_prob = sum(prob), .groups = "drop")

  cell_comm_total <- dplyr::full_join(source_prob, target_prob, by = "cell") %>%
    dplyr::mutate(
      source_prob = ifelse(is.na(source_prob), 0, source_prob),
      target_prob = ifelse(is.na(target_prob), 0, target_prob),
      total_prob = source_prob + target_prob
    ) %>%
    dplyr::select(cell, total_prob)

  cell_comm_weighted <- cell_comm_total %>%
    dplyr::left_join(
      data.frame(cell = names(group_prop), cell_prop = group_prop),
      by = "cell"
    ) %>%
    dplyr::mutate(
      cell_prop = ifelse(is.na(cell_prop), min_val, cell_prop),
      weighted_score = total_prob * cell_prop
    )

  if(group_name == "Group2"){
    top5_cells <- cell_comm_weighted %>%
      dplyr::arrange(dplyr::desc(weighted_score)) %>%
      dplyr::slice_head(n = 7) %>%
      dplyr::pull(cell)
  }else{
    top5_cells <- cell_comm_weighted %>%
      dplyr::arrange(dplyr::desc(weighted_score)) %>%
      dplyr::slice_head(n = 5) %>%
      dplyr::pull(cell)
  }



  # 打印Top5结果
  cat("✓ 【", group_name, "】通讯强度最强的前5个细胞类型：\n")
  for (i in 1:length(top5_cells)) {
    cell_score <- cell_comm_weighted %>% dplyr::filter(cell == top5_cells[i])
    cat(sprintf("  %d. %s（总通讯强度：%.6f | 细胞占比：%.6f | 加权得分：%.6f）\n",
                i, top5_cells[i], cell_score$total_prob[1], cell_score$cell_prop[1], cell_score$weighted_score[1]))
  }


  if(group_name == "Group2"){
    # 保存Top5细胞基础信息
    top5_detail <- cell_comm_weighted %>%
      dplyr::filter(cell %in% top5_cells) %>%
      dplyr::arrange(dplyr::desc(weighted_score)) %>%
      dplyr::mutate(rank = 1:7) %>%
      dplyr::select(rank, cell_type = cell, total_comm_strength = total_prob, 
                    cell_proportion = cell_prop, weighted_score)
  }else{
    # 保存Top5细胞基础信息
    top5_detail <- cell_comm_weighted %>%
      dplyr::filter(cell %in% top5_cells) %>%
      dplyr::arrange(dplyr::desc(weighted_score)) %>%
      dplyr::mutate(rank = 1:5) %>%
      dplyr::select(rank, cell_type = cell, total_comm_strength = total_prob, 
                    cell_proportion = cell_prop, weighted_score)
  }

  write.csv(
    top5_detail,
    file.path(group_subdir, paste0(group_name, "_Top5_Comm_CellTypes.csv")),
    row.names = FALSE, quote = FALSE
  )

  # ========== 筛选Top5细胞分角色的前3通路 ==========
  cat("\n✓ 开始筛选Top5细胞分角色的前3核心通路...\n")
  top5_cell_role_pathway <- data.frame()

  for (cell_name in top5_cells) {
    cat(sprintf("\n--- 分析细胞：%s ---\n", cell_name))

    # 发送者前3通路
    source_interactions <- df.net %>% dplyr::filter(source == cell_name)
    if (nrow(source_interactions) == 0) {
      cat(sprintf("  【发送者】%s无作为发送者的通讯互作，跳过\n", cell_name))
      source_pathway_record <- data.frame(
        cell_name = cell_name,
        role = "sender",
        pathway_rank = 1:3,
        pathway_name = "No_Interaction",
        interaction_count = 0,
        total_prob = 0,
        stringsAsFactors = FALSE
      )
    } else {
      source_pathway <- source_interactions %>%
        dplyr::group_by(pathway_name) %>%
        dplyr::summarise(
          interaction_count = dplyr::n(),
          total_prob = sum(prob),
          .groups = "drop"
        ) %>%
        dplyr::arrange(dplyr::desc(total_prob), dplyr::desc(interaction_count)) %>%
        dplyr::mutate(pathway_rank = dplyr::row_number()) %>%
        dplyr::slice_head(n = 3)

      if (nrow(source_pathway) < 3) {
        fill_rows <- data.frame(
          pathway_name = rep("No_Interaction", 3 - nrow(source_pathway)),
          interaction_count = 0,
          total_prob = 0,
          pathway_rank = (nrow(source_pathway)+1):3
        )
        source_pathway <- dplyr::bind_rows(source_pathway, fill_rows)
      }
      source_pathway_record <- source_pathway %>%
        dplyr::mutate(
          cell_name = cell_name,
          role = "sender"
        ) %>%
        dplyr::select(cell_name, role, pathway_rank, pathway_name, interaction_count, total_prob)
      
      cat("  【发送者】前3通路：\n")
      for (j in 1:3) {
        rec <- source_pathway_record %>% dplyr::filter(pathway_rank == j)
        cat(sprintf("    %d. %s（互作数：%d | 总通讯强度：%.6f）\n",
                    j, rec$pathway_name[1], rec$interaction_count[1], rec$total_prob[1]))
      }
    }

    # 接收者前3通路
    target_interactions <- df.net %>% dplyr::filter(target == cell_name)
    if (nrow(target_interactions) == 0) {
      cat(sprintf("  【接收者】%s无作为接收者的通讯互作，跳过\n", cell_name))
      target_pathway_record <- data.frame(
        cell_name = cell_name,
        role = "receiver",
        pathway_rank = 1:3,
        pathway_name = "No_Interaction",
        interaction_count = 0,
        total_prob = 0,
        stringsAsFactors = FALSE
      )
    } else {
      target_pathway <- target_interactions %>%
        dplyr::group_by(pathway_name) %>%
        dplyr::summarise(
          interaction_count = dplyr::n(),
          total_prob = sum(prob),
          .groups = "drop"
        ) %>%
        dplyr::arrange(dplyr::desc(total_prob), dplyr::desc(interaction_count)) %>%
        dplyr::mutate(pathway_rank = dplyr::row_number()) %>%
        dplyr::slice_head(n = 3)

      if (nrow(target_pathway) < 3) {
        fill_rows <- data.frame(
          pathway_name = rep("No_Interaction", 3 - nrow(target_pathway)),
          interaction_count = 0,
          total_prob = 0,
          pathway_rank = (nrow(target_pathway)+1):3
        )
        target_pathway <- dplyr::bind_rows(target_pathway, fill_rows)
      }
      target_pathway_record <- target_pathway %>%
        dplyr::mutate(
          cell_name = cell_name,
          role = "receiver"
        ) %>%
        dplyr::select(cell_name, role, pathway_rank, pathway_name, interaction_count, total_prob)
      
      cat("  【接收者】前3通路：\n")
      for (j in 1:3) {
        rec <- target_pathway_record %>% dplyr::filter(pathway_rank == j)
        cat(sprintf("    %d. %s（互作数：%d | 总通讯强度：%.6f）\n",
                    j, rec$pathway_name[1], rec$interaction_count[1], rec$total_prob[1]))
      }
    }
    
    top5_cell_role_pathway <- dplyr::bind_rows(top5_cell_role_pathway, source_pathway_record, target_pathway_record)
  }

  # 保存分角色通路结果
  write.csv(
    top5_cell_role_pathway,
    file.path(group_subdir, paste0(group_name, "_Top5_Cells_Role_Specific_Pathways.csv")),
    row.names = FALSE, quote = FALSE
  )

  top5_cell_role_pathway_valid <- top5_cell_role_pathway %>%
    dplyr::filter(pathway_name != "No_Interaction")
  write.csv(
    top5_cell_role_pathway_valid,
    file.path(group_subdir, paste0(group_name, "_Top5_Cells_Role_Specific_Pathways_Valid.csv")),
    row.names = FALSE, quote = FALSE
  )

  cat(sprintf("\n✓ Top5细胞分角色通路结果已保存至：%s\n", group_subdir))

  # ========== 绘制Top5细胞分角色通路图（核心修复文件名） ==========
  # 创建通路图子目录
  pathway_plot_dir <- file.path(group_subdir, "Pathway_Plots")
  if (!dir.exists(pathway_plot_dir)) dir.create(pathway_plot_dir, recursive = TRUE)
  head(cellchat@netP)
  head(cellchat@netP$pathways)
  # cellchat@idents 结构：名称=细胞barcode，值=cellType（即current_cells）


  extract_cellchat_subset_debug <- function(cellchat, target_celltypes) {
    cat(sprintf("\n=== 【终极调试版】子集化 + 容错修复 ===\n"))
    cat("目标细胞类型：", paste(target_celltypes, collapse = ", "), "\n")

    # 1. 精准匹配目标细胞类型
    all_celltypes <- colnames(cellchat@data)  
    target_celltypes_exist <- intersect(target_celltypes, all_celltypes)

    if (length(target_celltypes_exist) == 0) {
      stop("❌ 无匹配的目标细胞类型！请检查名称拼写")
    }
    cat(sprintf("✅ 匹配到%d个有效细胞类型：%s\n", length(target_celltypes_exist), paste(target_celltypes_exist, collapse = ", ")))

    # 2. 初始化子集对象（深度拷贝避免修改原对象）
    cellchat_subset <- cellchat
    cellchat_subset@data <- cellchat@data[, target_celltypes_exist, drop = FALSE]
    cellchat_subset@data.signaling <- cellchat@data.signaling[, target_celltypes_exist, drop = FALSE]

    # 核心修复1：重新设置idents（保持因子水平仅为目标细胞）
    cellchat_subset@idents <- factor(target_celltypes_exist, levels = target_celltypes_exist)
    names(cellchat_subset@idents) <- target_celltypes_exist
    cellchat_subset@meta <- cellchat@meta[target_celltypes_exist, , drop = FALSE]

    # 打印基础信息
    cat("\n【调试1】筛选前表达矩阵维度：", nrow(cellchat@data), "×", ncol(cellchat@data), "\n")
    cat("【调试1】筛选后表达矩阵维度：", nrow(cellchat_subset@data), "×", ncol(cellchat_subset@data), "\n")
    cat("【调试1】@idents 最终长度：", length(cellchat_subset@idents), "| 水平数：", nlevels(cellchat_subset@idents), "\n")

    # ---- 2.2 @net 完整筛选（核心修复：强制设置dimnames + 命名属性）----
    if (!is.null(cellchat@net)) {
      cellchat_subset@net <- list()
      # 处理weight
      if (!is.null(cellchat@net$weight)) {
        orig_weight <- cellchat@net$weight[target_celltypes_exist, target_celltypes_exist, drop = FALSE]
        # 修复1：填充空值 + 强制设置dimnames
        orig_weight[is.na(orig_weight)] <- 1e-6
        orig_weight[orig_weight == 0] <- 1e-6
        dimnames(orig_weight) <- list(target_celltypes_exist, target_celltypes_exist) # 关键：命名必须匹配
        cellchat_subset@net$weight <- orig_weight
        cat("【调试2】@net$weight 筛选后维度：", nrow(orig_weight), "×", ncol(orig_weight), "\n")
        cat("【调试2】@net$weight 非0值数量：", sum(orig_weight > 1e-6, na.rm = TRUE), "\n")
      }
      # 处理count
      if (!is.null(cellchat@net$count)) {
        orig_count <- cellchat@net$count[target_celltypes_exist, target_celltypes_exist, drop = FALSE]
        orig_count[is.na(orig_count)] <- 0
        dimnames(orig_count) <- list(target_celltypes_exist, target_celltypes_exist)
        cellchat_subset@net$count <- orig_count
        cat("【调试2】@net$count 筛选后维度：", nrow(orig_count), "×", ncol(orig_count), "\n")
      }
      # 处理prob（核心修复：维度+命名双校验）
      if (!is.null(cellchat@net$prob)) {
        net_prob_dim <- dim(cellchat@net$prob)
        cat("【调试2】原始@net$prob 维度：", paste(net_prob_dim, collapse = "×"), "\n")

        if (length(net_prob_dim) == 2) {
          orig_prob <- cellchat@net$prob[target_celltypes_exist, target_celltypes_exist, drop = FALSE]
          dimnames(orig_prob) <- list(target_celltypes_exist, target_celltypes_exist)
        } else if (length(net_prob_dim) == 3) {
          orig_prob <- cellchat@net$prob[target_celltypes_exist, target_celltypes_exist, , drop = FALSE]
          dimnames(orig_prob) <- list(
            target_celltypes_exist, 
            target_celltypes_exist, 
            dimnames(cellchat@net$prob)[[3]] # 保留原第三维命名
          )
        } else {
          orig_prob <- cellchat@net$prob
          warning("⚠️ @net$prob 维度异常（", paste(net_prob_dim, collapse = "×"), "），跳过筛选")
        }
        orig_prob[is.na(orig_prob)] <- 1e-6
        orig_prob[orig_prob == 0] <- 1e-6
        cellchat_subset@net$prob <- orig_prob
        cat("【调试2】@net$prob 筛选后维度：", paste(dim(orig_prob), collapse = "×"), "\n")
      } else {
        cat("【调试2】@net$prob 不存在，跳过筛选\n")
      }
    }

    # ---- 2.3 @netP 完整修复（核心：同步所有子槽位 + 命名强制匹配）----
    cat("\n【调试3】原始@netP$prob 维度：", paste(dim(cellchat@netP$prob), collapse = "×"), "\n")
    if (!is.null(cellchat@netP$prob)) {
    # 1. 筛选netP$prob + 强制设置dimnames
    orig_netP_prob <- cellchat@netP$prob[target_celltypes_exist, target_celltypes_exist, , drop = FALSE]
    orig_netP_prob[is.na(orig_netP_prob)] <- 1e-6
    orig_netP_prob[orig_netP_prob == 0] <- 1e-6
    # 关键：手动设置三维矩阵的dimnames
    dimnames(orig_netP_prob) <- list(
      source = target_celltypes_exist,
      target = target_celltypes_exist,
      pathway = cellchat@netP$pathways # 保留原通路名
    )
    cellchat_subset@netP$prob <- orig_netP_prob

    # 2. 打印筛选后维度信息
    netP_prob_dim <- dim(cellchat_subset@netP$prob)
    cat("【调试3】筛选后@netP$prob 维度：", paste(netP_prob_dim, collapse = "×"), "\n")
    cat("【调试3】@netP$prob source数：", netP_prob_dim[1], "\n")
    cat("【调试3】@netP$prob target数：", netP_prob_dim[2], "\n")
    cat("【调试3】@netP$prob 通路数：", netP_prob_dim[3], "\n")

    # 3. 检查每个通路的非0值数量
    pathway_non0 <- sapply(1:netP_prob_dim[3], function(p) {
      sum(cellchat_subset@netP$prob[, , p] > 1e-6, na.rm = TRUE)
    })
    names(pathway_non0) <- cellchat@netP$pathways
    cat("【调试3】各通路非0通讯数（前10）：\n")
    print(head(sort(pathway_non0, decreasing = TRUE), 10))

    # 4. 过滤全0通路 + 同步更新所有netP子槽位
    valid_pathway_idx <- which(pathway_non0 > 0)
    if (length(valid_pathway_idx) == 0) {
      warning("⚠️ 所有通路的通讯概率全为0，创建占位通路避免报错")
      # 修复：创建有命名的占位通路（核心：命名长度必须匹配数值维度）
      cellchat_subset@netP$prob <- array(1e-6, dim = c(length(target_celltypes_exist), length(target_celltypes_exist), 1))
      dimnames(cellchat_subset@netP$prob) <- list(
        source = target_celltypes_exist,
        target = target_celltypes_exist,
        pathway = "Placeholder_Pathway"
      )
      cellchat_subset@netP$pathways <- "Placeholder_Pathway"
      # 同步设置其他netP子槽位（避免绘图时读取空值）
      cellchat_subset@netP$count <- array(0, dim = dim(cellchat_subset@netP$prob))
      dimnames(cellchat_subset@netP$count) <- dimnames(cellchat_subset@netP$prob)
      cellchat_subset@netP$pairwise <- data.frame(
        source = target_celltypes_exist[1],
        target = target_celltypes_exist[1],
        pathway = "Placeholder_Pathway",
        prob = 1e-6,
        count = 0,
        stringsAsFactors = FALSE
      )
    } else {
      # 5. 保留有效通路 + 同步筛选所有netP子槽位
      cellchat_subset@netP$prob <- cellchat_subset@netP$prob[, , valid_pathway_idx, drop = FALSE]
      cellchat_subset@netP$pathways <- cellchat@netP$pathways[valid_pathway_idx]
      
      # 同步筛选netP$count（关键：dimnames必须匹配）
      if (!is.null(cellchat@netP$count)) {
        cellchat_subset@netP$count <- cellchat@netP$count[target_celltypes_exist, target_celltypes_exist, valid_pathway_idx, drop = FALSE]
        cellchat_subset@netP$count[is.na(cellchat_subset@netP$count)] <- 0
        dimnames(cellchat_subset@netP$count) <- list(
          source = target_celltypes_exist,
          target = target_celltypes_exist,
          pathway = cellchat_subset@netP$pathways
        )
        cat("【调试3】过滤后@netP$count 维度：", paste(dim(cellchat_subset@netP$count), collapse = "×"), "\n")
      }
      
      # 同步筛选netP$pairwise（绘图核心依赖！）
      if (!is.null(cellchat@netP$pairwise)) {
        cellchat_subset@netP$pairwise <- cellchat@netP$pairwise[
          cellchat@netP$pairwise$source %in% target_celltypes_exist &
          cellchat@netP$pairwise$target %in% target_celltypes_exist &
          cellchat@netP$pairwise$pathway %in% cellchat_subset@netP$pathways, 
          , drop = FALSE
        ]
        # 修复：空pairwise填充默认值
        if (nrow(cellchat_subset@netP$pairwise) == 0) {
          cellchat_subset@netP$pairwise <- data.frame(
            source = target_celltypes_exist[1],
            target = target_celltypes_exist[1],
            pathway = cellchat_subset@netP$pathways[1],
            prob = 1e-6,
            count = 0,
            stringsAsFactors = FALSE
          )
        }
      }
      
      # 同步筛选netP$LR/netP$gene（避免绘图时找不到基因信息）
      if (!is.null(cellchat@netP$LR)) {
        cellchat_subset@netP$LR <- cellchat@netP$LR[cellchat@netP$LR$pathway_name %in% cellchat_subset@netP$pathways, , drop = FALSE]
      }
      if (!is.null(cellchat@netP$gene)) {
        cellchat_subset@netP$gene <- cellchat@netP$gene[cellchat@netP$gene$pathway %in% cellchat_subset@netP$pathways, , drop = FALSE]
      }
    }
    
    # 6. 重新聚合网络 + 校验聚合结果
    cellchat_subset <- aggregateNet(cellchat_subset)
    # 强制修复聚合后的pairwise命名（最终兜底）
    if (!is.null(cellchat_subset@netP$pairwise)) {
      rownames(cellchat_subset@netP$pairwise) <- 1:nrow(cellchat_subset@netP$pairwise)
    }
    
    cat("【调试3】过滤后通路数：", length(cellchat_subset@netP$pathways), "\n")
    cat("【调试3】过滤后@netP$prob 维度：", paste(dim(cellchat_subset@netP$prob), collapse = "×"), "\n")
    cat("【调试3】过滤后@netP$pairwise 行数：", nrow(cellchat_subset@netP$pairwise), "\n")
    }

    # ---- 2.4 修复LR和特征基因（保证命名一致性）----
    if (!is.null(cellchat@LR)) {
    cellchat_subset@LR <- cellchat@LR
    # 筛选LRsig并保证非空
    if (!is.null(cellchat_subset@netP$pathways)) {
      cellchat_subset@LR$LRsig <- cellchat_subset@LR$LR[
        cellchat_subset@LR$LR$pathway_name %in% cellchat_subset@netP$pathways, 
        , drop = FALSE
      ]
      if (nrow(cellchat_subset@LR$LRsig) == 0) {
        # 填充占位LR信息（匹配占位通路）
        cellchat_subset@LR$LRsig <- cellchat_subset@LR$LR[1, , drop = FALSE]
        cellchat_subset@LR$LRsig$pathway_name <- cellchat_subset@netP$pathways[1]
      }
    }
    # 强制设置rownames
    rownames(cellchat_subset@LR$LRsig) <- 1:nrow(cellchat_subset@LR$LRsig)
    }

    if (!is.null(cellchat@var.features)) {
    lr_genes <- unique(c(cellchat_subset@LR$LR$ligand, cellchat_subset@LR$LR$receptor))
    cellchat_subset@var.features$features <- intersect(cellchat@var.features$features, lr_genes)
    if (length(cellchat_subset@var.features$features) == 0) {
      cellchat_subset@var.features$features <- if (length(lr_genes) > 0) lr_genes[1:min(5, length(lr_genes))] else "Placeholder_Gene"
    }
    # 关键：features必须有命名（长度匹配）
    names(cellchat_subset@var.features$features) <- cellchat_subset@var.features$features
    cat("【调试4】筛选后特征基因数：", length(cellchat_subset@var.features$features), "\n")
    }

    # ---- 2.5 重置全局参数 + 最终校验 ----
    cellchat_subset@options$cell.use <- target_celltypes_exist
    cellchat_subset@options$population.size <- FALSE 
    cellchat_subset@options$nCell <- length(target_celltypes_exist) # 补充细胞数参数

    # 最终校验：打印核心数据的命名状态
    cat("\n【最终校验】子集化后核心数据状态：\n")
    cat("  @net$weight dimnames：", paste(sapply(dimnames(cellchat_subset@net$weight), length), collapse = "×"), "\n")
    cat("  @netP$prob dimnames：", paste(sapply(dimnames(cellchat_subset@netP$prob), length), collapse = "×"), "\n")
    cat("  @netP$pathways数：", length(cellchat_subset@netP$pathways), "\n")
    cat("  @netP$pairwise行数：", nrow(cellchat_subset@netP$pairwise), "\n")
    cat("  @LR$LRsig行数：", nrow(cellchat_subset@LR$LRsig), "\n")

    return(cellchat_subset)
  }

  if(group_name == "Group1"){
    cellchat_AD1_subset <- extract_cellchat_subset_debug(
      cellchat = cellchat,
      target_celltypes = AD_Canine1_Human1_cell_types
    )
    head(cellchat_AD1_subset@idents)
  }else if(group_name == "Group2"){
    cellchat_AD1_subset <- extract_cellchat_subset_debug(
      cellchat = cellchat,
      target_celltypes = AD_Canine2_Human3_cell_types
    )
    head(cellchat_AD1_subset@idents)
  }else{
    cellchat_AD1_subset <- cellchat
  }
  # cellchat <- cellchat_AD1_subset
  # head(cellchat@idents)
  # if(group_name == "Group1"){
  #   cellchat <- subsetData(cellchat, idents.use = AD_Canine1_Human1_cell_types)
  # }else if(group_name == "Grou2"){
  #   cellchat <- subsetData(cellchat, idents.use =AD_Canine2_Human3_cell_types)
  # }
  # 目标细胞列表
  # add_cells <- c("Th17", "Th22", "Th1", "Th2", "Treg")
  # target_cell_list <- unique(c(top5_cells, add_cells))
  target_cell_list <- c(top5_cells)
  target_cell_list

  for (target_cell in target_cell_list) {
    # target_cell <- "Bas_Quiescent"
    if (!(target_cell %in% current_cells)) {
      cat(sprintf("\n✗ 目标细胞【%s】不在当前cellchat中，跳过\n", target_cell))
      next
    }
    
    # 清理细胞名（核心修复！！！）
    target_cell_clean <- clean_filename_char(target_cell)
    
    # 提取前3通路
    sender_pathways <- top5_cell_role_pathway %>%
      dplyr::filter(cell_name == target_cell, role == "sender", pathway_name != "No_Interaction") %>%
      dplyr::arrange(pathway_rank) %>%
      dplyr::pull(pathway_name) %>%
      head(3)

    receiver_pathways <- top5_cell_role_pathway %>%
      dplyr::filter(cell_name == target_cell, role == "receiver", pathway_name != "No_Interaction") %>%
      dplyr::arrange(pathway_rank) %>%
      dplyr::pull(pathway_name) %>%
      head(3)

    # ========== 绘制发送者通路图 ==========
    if (length(sender_pathways) > 0) {
      cat(sprintf("\n--- 绘制细胞【%s（发送者）】的前3通路图 ---\n", target_cell))
      for (pathway in sender_pathways) {
        # pathway <- "LAMININ"
        if (!(pathway %in% cellchat_AD1_subset@netP$pathways)) {
          cat(sprintf("✗ 通路【%s】不存在，跳过\n", pathway))
          next
        }

        # 生成安全的文件名前缀
        file_prefix <- paste0(group_name, "_", target_cell_clean, "_Sender_", clean_filename_char(pathway))
        
        # LR贡献图
        tryCatch({
          Sys.sleep(1) 
          gc(verbose = FALSE)
          pdf(file.path(pathway_plot_dir, paste0(file_prefix,"_LR_contribution.pdf")), width = 5, height = 2.5)
          p1 <- netAnalysis_contribution(cellchat_AD1_subset, signaling = pathway)
          print(p1)
          dev.off()
          Sys.sleep(1)
          cat(sprintf("✓ 已保存%s通路分析结果\n", pathway))
        }, error = function(e) {
          cat(sprintf("✗ 绘制【%s】LR贡献图失败：%s\n", pathway, e$message))
        })

        # Heatmap图
        tryCatch({
          Sys.sleep(1) 
          gc(verbose = FALSE)
          save_path <- file.path(pathway_plot_dir, paste0(file_prefix, "_Heatmap.pdf"))
          pdf(
            save_path,
            width = 3.5, height =3,
            onefile = TRUE,
            useDingbats = FALSE,
            family = "Helvetica"
          )
          p <- netVisual_heatmap(cellchat_AD1_subset, signaling = pathway, color.heatmap = "Reds")
          print(p)
          dev.off()
          Sys.sleep(1)
          cat(sprintf("✓ 已保存%s Heatmap图\n", pathway))
        }, error = function(e) {
          cat(sprintf("✗ 绘制【%s】Heatmap图失败：%s\n", pathway, e$message))
        })

        # Hierarchy图
        tryCatch({
          cat(sprintf("正在绘制：%s_Hierarchy.pdf\n", file_prefix))
          save_path <- file.path(pathway_plot_dir, paste0(file_prefix, "_Hierarchy.pdf"))
          while (dev.cur() > 1) dev.off()
          pdf(save_path, width = 9, height = 8)
          netVisual_aggregate(
            object = cellchat_AD1_subset,
            signaling = pathway,
            layout = "hierarchy",
            vertex.weight = 20,
            vertex.size.max = 20,
            vertex.label.cex = 0.8,
            edge.width.max = 10,
            show.legend = TRUE,
            vertex.receiver = seq(1:min(3, length(current_cells))),
            legend.pos.x = 6, 
            legend.pos.y = 1
          )
          dev.off()
          cat(sprintf("✓ 已保存%s Hierarchy图\n", pathway))
        }, error = function(e) {
          cat(sprintf("✗ 绘制【%s】Hierarchy图失败：%s\n", pathway, e$message))
        })

        # Chord图
        tryCatch({
          cat(sprintf("正在绘制：%s_Chord.pdf\n", file_prefix))
          pdf(file.path(pathway_plot_dir, paste0(file_prefix, "_Chord.pdf")), width = 10, height = 8)
          netVisual_aggregate(cellchat_AD1_subset, signaling = pathway, layout = "chord")
          dev.off()
          cat(sprintf("✓ 已保存%s Chord图\n", pathway))
        }, error = function(e) {
          cat(sprintf("✗ 绘制【%s】Chord图失败：%s\n", pathway, e$message))
        })

        cat(sprintf("✓ 已保存【%s（发送者）-%s】的3类图\n", target_cell, pathway))
      }
    } else {
      cat(sprintf("\n✗ 细胞【%s】无有效发送者通路，跳过发送者绘图\n", target_cell))
    }

    # ========== 绘制接收者通路图 ==========
    if (length(receiver_pathways) > 0) {
      cat(sprintf("\n--- 绘制细胞【%s（接收者）】的前3通路图 ---\n", target_cell))
      for (pathway in receiver_pathways) {
        if (!(pathway %in% cellchat_AD1_subset@netP$pathways)) {
          cat(sprintf("✗ 通路【%s】不存在，跳过\n", pathway))
          next
        }

        # 生成安全的文件名前缀
        file_prefix <- paste0(group_name, "_", target_cell_clean, "_Receiver_", clean_filename_char(pathway))
        
        # Heatmap图
        tryCatch({
          Sys.sleep(1) 
          gc(verbose = FALSE)
          save_path <- file.path(pathway_plot_dir, paste0(file_prefix, "_Heatmap.pdf"))
          pdf(
            save_path,
            width = 3.5, height = 3,
            onefile = TRUE,
            useDingbats = FALSE,
            family = "Helvetica"
          )
          p <- netVisual_heatmap(cellchat_AD1_subset, signaling = pathway, color.heatmap = "Reds")
          print(p)
          dev.off()
          Sys.sleep(1)
          cat(sprintf("✓ 已保存%s Heatmap图\n", pathway))
        }, error = function(e) {
          cat(sprintf("✗ 绘制【%s】Heatmap图失败：%s\n", pathway, e$message))
        })

        # Hierarchy图
        tryCatch({
          cat(sprintf("正在绘制：%s_Hierarchy.pdf\n", file_prefix))
          save_path <- file.path(pathway_plot_dir, paste0(file_prefix, "_Hierarchy.pdf"))
          while (dev.cur() > 1) dev.off()
          pdf(save_path, width = 9, height = 8)
          netVisual_aggregate(
            object = cellchat_AD1_subset,
            signaling = pathway,
            layout = "hierarchy",
            vertex.weight = 20,
            vertex.size.max = 20,
            vertex.label.cex = 0.8,
            edge.width.max = 10,
            show.legend = TRUE,
            vertex.receiver = seq(1:min(3, length(current_cells))),
            legend.pos.x = 6, 
            legend.pos.y = 1
          )
          dev.off()
          cat(sprintf("✓ 已保存%s Hierarchy图\n", pathway))
        }, error = function(e) {
          cat(sprintf("✗ 绘制【%s】Hierarchy图失败：%s\n", pathway, e$message))
        })

        # Chord图
        tryCatch({
          cat(sprintf("正在绘制：%s_Chord.pdf\n", file_prefix))
          pdf(file.path(pathway_plot_dir, paste0(file_prefix, "_Chord.pdf")), width = 10, height = 8)
          netVisual_aggregate(cellchat_AD1_subset, signaling = pathway, layout = "chord")
          dev.off()
          cat(sprintf("✓ 已保存%s Chord图\n", pathway))
        }, error = function(e) {
          cat(sprintf("✗ 绘制【%s】Chord图失败：%s\n", pathway, e$message))
        })

        cat(sprintf("✓ 已保存【%s（接收者）-%s】的3类图\n", target_cell, pathway))
      }
    } else {
      cat(sprintf("\n✗ 细胞【%s】无有效接收者通路，跳过接收者绘图\n", target_cell))
    }
  }
  # cat("\n✓ 所有目标细胞的分角色通路图绘制完成！结果保存至：", pathway_plot_dir, "\n")

  # ========== 绘制Bubble图 ==========
  # 创建Bubble图子目录
  bubble_plot_dir <- file.path(group_subdir, "Bubble_Plots")
  if (!dir.exists(bubble_plot_dir)) dir.create(bubble_plot_dir, recursive = TRUE)
  
  target_th_list <- c(top5_cells, "Th17", "Th22", "Th1", "Th2", "Treg")
  immune_indices <- which(current_cells %in% target_th_list)
  
  if (length(immune_indices) == 0) {
    cat("✗ 无有效免疫细胞，跳过Bubble图\n")
  } else {
      # 处理Top5细胞
      cat("=== 开始绘制Top5细胞的发送者/接收者Bubble图 ===\n")
      for (target_cell in top5_cells) {
        # target_cell <- "Bas_Quiescent"
        if (!(target_cell %in% current_cells)) {
          cat(sprintf("\n✗ Top5细胞【%s】不在cellchat中，跳过\n", target_cell))
          next
        }
        clean_cell <- clean_filename_char(target_cell)
        # 接收者Bubble图
        cat(sprintf("\n--- 绘制【%s（接收者）】的Bubble图 ---\n", target_cell))
        sources_use_receiver <- setdiff(cellchat_AD1_subset@idents, target_cell)
        if (group_name == "Group1") {
          sources_use_receiver <- AD_Canine1_Human1_cell_types
        }else if(group_name == "Group2"){
          sources_use_receiver <- AD_Canine2_Human3_cell_types
        }
        if (length(sources_use_receiver) > 0) {
          # 全通路
          tryCatch({
            cat(sprintf("正在绘制：%s→%s（接收者）全通路Bubble图...\n", group_name, target_cell))
            gg_receiver_all <- netVisual_bubble(
              cellchat,
              sources.use = sources_use_receiver,
              targets.use = target_cell,
              remove.isolate = FALSE,
              title.name = paste0(group_name, " - Others→", target_cell, " (Receiver - All Pathways)")
            )
            ggsave(
              file.path(bubble_plot_dir, paste0(group_name, "_Receiver_", clean_cell, "_All_bubble.pdf")),
              gg_receiver_all, 
              width = 4, 
              height = 40,
              limitsize = FALSE
            )
          }, error = function(e) {
            cat(sprintf("✗ 绘制【%s】接收者全通路Bubble图失败：%s\n", target_cell, e$message))
          })

          # 前3通路
          receiver_top3_pathways <- top5_cell_role_pathway %>%
            dplyr::filter(cell_name == target_cell, role == "receiver", pathway_name != "No_Interaction") %>%
            dplyr::arrange(pathway_rank) %>%
            dplyr::pull(pathway_name) %>%
            head(3)

          if (length(receiver_top3_pathways) > 0) {
          tryCatch({
            cat(sprintf("正在绘制：%s→%s（接收者）前3通路（%s）Bubble图...\n", 
                        group_name, target_cell, paste(receiver_top3_pathways, collapse = "/")))
            gg_receiver_top3 <- netVisual_bubble(
              cellchat,
              sources.use = sources_use_receiver,
              targets.use = target_cell,
              signaling = receiver_top3_pathways,
              remove.isolate = FALSE,
              title.name = paste0(group_name, " - Others→", target_cell, "Top3 Pathways"),
              font.size.title = 9                 
            )
            ggsave(
              file.path(bubble_plot_dir, paste0(group_name, "_Receiver_", clean_cell, "_Top3_bubble.pdf")),
              gg_receiver_top3, 
              width = 3, 
              height = calc_safe_height(length(sources_use_receiver), 0.5)*2,
              # height = calc_safe_height(length(sources_use_receiver), 0.5)*3,
              limitsize = FALSE
            )
            cat(sprintf("✓ 已保存【%s（接收者）】前3通路Bubble图\n", target_cell))
          }, error = function(e) {
            cat(sprintf("✗ 绘制【%s】接收者前3通路Bubble图失败：%s\n", target_cell, e$message))
          })
          } else {
            cat(sprintf("✗ 【%s（接收者）】无有效前3通路，跳过\n", target_cell))
          }
        }

        # 发送者Bubble图
        cat(sprintf("\n--- 绘制【%s（发送者）】的Bubble图 ---\n", target_cell))
        targets_use_sender <- setdiff(cellchat@idents, target_cell)
        if (group_name == "Group1") {
          targets_use_sender <- AD_Canine1_Human1_cell_types
        }else if(group_name == "Group2"){
          targets_use_sender <- AD_Canine2_Human3_cell_types
        }
        if (length(targets_use_sender) > 0) {
          # 全通路
          tryCatch({
            cat(sprintf("正在绘制：%s（发送者）→%s 全通路Bubble图...\n", target_cell, group_name))
            gg_sender_all <- netVisual_bubble(
              cellchat,
              sources.use = target_cell,
              targets.use = targets_use_sender,
              remove.isolate = FALSE,
              title.name = paste0(group_name, " - ", target_cell, "→Others (Sender - All Pathways)")
            )
            ggsave(
              file.path(bubble_plot_dir, paste0(group_name, "_Sender_", clean_cell, "_All_bubble.pdf")),
              gg_sender_all, 
              width = 4, 
              height = 40,
              limitsize = FALSE
            )
          }, error = function(e) {
            cat(sprintf("✗ 绘制【%s】发送者全通路Bubble图失败：%s\n", target_cell, e$message))
          })

          # 前3通路
          sender_top3_pathways <- top5_cell_role_pathway %>%
            dplyr::filter(cell_name == target_cell, role == "sender", pathway_name != "No_Interaction") %>%
            dplyr::arrange(pathway_rank) %>%
            dplyr::pull(pathway_name) %>%
            head(3)

          if (length(sender_top3_pathways) > 0) {
            tryCatch({
              cat(sprintf("正在绘制：%s（发送者）→%s 前3通路（%s）Bubble图...\n", 
                          target_cell, group_name, paste(sender_top3_pathways, collapse = "/")))
              gg_sender_top3 <- netVisual_bubble(
                cellchat,
                sources.use = target_cell,
                targets.use = targets_use_sender,
                signaling = sender_top3_pathways,
                remove.isolate = FALSE,
                title.name = paste0(group_name, " - ", target_cell, "→Others Top3 Pathways"),
                font.size.title = 9                 
              )
              ggsave(
                file.path(bubble_plot_dir, paste0(group_name, "_Sender_", clean_cell, "_Top3_bubble.pdf")),
                gg_sender_top3, 
                width = 3, 
                height = calc_safe_height(length(targets_use_sender), 0.5)*2,
                # height = calc_safe_height(length(targets_use_sender), 0.5)*3,
                limitsize = FALSE
              )
              cat(sprintf("✓ 已保存【%s（发送者）】前3通路Bubble图\n", target_cell))
            }, error = function(e) {
              cat(sprintf("✗ 绘制【%s】发送者前3通路Bubble图失败：%s\n", target_cell, e$message))
            })
          } else {
            cat(sprintf("✗ 【%s（发送者）】无有效前3通路，跳过\n", target_cell))
          }
        }

        cat(sprintf("✓ 完成【%s】的发送者/接收者Bubble图绘制\n", target_cell))
      }

      # 处理Th系列细胞
      cat("\n=== 开始绘制Th系列细胞的接收者Bubble图 ===\n")
      th_cells <- c("Th17", "Th22", "Th1", "Th2", "Treg")
      spec_pathways_input <- c("THBS","SELPLG","SELL","IL1","IL4","MHC-II","FN1","CCL","AGT","ADIPONECTIN")
      spec_pathways_valid <- validate_pathways(spec_pathways_input, cellchat@netP$pathways)

      for (target_th in th_cells) {
        if (!(target_th %in% current_cells)) {
          cat(sprintf("\n✗ Th细胞【%s】不在cellchat中，跳过\n", target_th))
          next
        }

        clean_th <- clean_filename_char(target_th)
        sources_use_th <- setdiff(cellchat@idents, target_th)
        if (length(sources_use_th) == 0) {
          cat(sprintf("✗ 除%s外无其他细胞，跳过\n", target_th))
          next
        }

        # 全通路
        tryCatch({
          cat(sprintf("\n正在绘制：%s→%s（接收者）全通路Bubble图...\n", group_name, target_th))
          gg_th_all <- netVisual_bubble(
            cellchat,
            sources.use = sources_use_th,
            targets.use = target_th,
            remove.isolate = FALSE,
            title.name = paste0(group_name, " - Others→", target_th, " (Th Cell - All Pathways)")
          )
          ggsave(
            file.path(bubble_plot_dir, paste0(group_name, "_Others_to_", clean_th, "_All_bubble.pdf")),
            gg_th_all, 
            width = 8, 
            height = 40,
            limitsize = FALSE
          )
        }, error = function(e) {
          cat(sprintf("✗ 绘制【%s】全通路Bubble图失败：%s\n", target_th, e$message))
        })

        # 特异性通路
        if (length(spec_pathways_valid) > 0) {
          tryCatch({
            cat(sprintf("正在绘制：%s→%s 特异性通路Bubble图...\n", group_name, target_th))
            gg_th_spec <- netVisual_bubble(
              cellchat,
              sources.use = sources_use_th,
              targets.use = target_th,
              signaling = spec_pathways_valid,
              remove.isolate = FALSE,
              title.name = paste0(group_name, " - Others→", target_th, " (Th Cell - Specific Pathways)")
            )
            ggsave(
              file.path(bubble_plot_dir, paste0(group_name, "_Others_to_", clean_th, "_Specific_bubble.pdf")),
              gg_th_spec, 
              width = 5, 
              height = calc_safe_height(length(sources_use_th), 0.5),
              limitsize = FALSE
            )
            cat(sprintf("✓ 已保存【%s】特异性通路Bubble图\n", target_th))
          }, error = function(e) {
            cat(sprintf("✗ 绘制【%s】特异性通路Bubble图失败：%s\n", target_th, e$message))
          })
        } else {
          cat(sprintf("✗ 【%s】无有效特异性通路，跳过\n", target_th))
        }

        cat(sprintf("✓ 完成【%s】的Bubble图绘制\n", target_th))
      }
  }
  cat("\n=== 所有Bubble图绘制完成！结果目录：", bubble_plot_dir, "===\n")

  # 保存CellChat对象
  saveRDS(cellchat, file.path(group_subdir, paste0(group_name, "_CellChat_results.rds")))
  cat("✓ 保存完整分析结果：", file.path(group_subdir, paste0(group_name, "_CellChat_results.rds")), "\n\n")
  
  return(cellchat)  
}

# ======================
# 5. 批量处理有效组
# ======================
cat("\n=== 批量处理有效组 ===\n")
cellchat_results <- list()

for (group_name in valid_groups) {
  cat(sprintf("\n====================================\n"))
  cellchat_results[[group_name]] <- analyze_group_full_cellchat(
    group_name = group_name,
    weighted_expr = group_weighted_expr[[group_name]],
    matched_cells = group_matched_cells[[group_name]],
    group_prop = group_avg_props[[group_name]],
    immune_matched = group_immune_matched[[group_name]]
  )
}

cat("\n=== 所有分组分析完成！主结果目录：", main_result_dir, "===\n")