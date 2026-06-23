  library(gmodels)
  library(ggplot2)
  library(RColorBrewer)
  library(ggsci)
  library(ggrepel)
  library(data.table)
  library(Rtsne)
  library(dendextend)
  library(ape)
  library(stats)

  # 读取命令行参数
  args <- commandArgs()

  # 读取表达数据
  expr <- fread(args[6])
  expr <- as.data.frame(expr)

  # 分组筛选（假设输入2个分组），根据命令行参数args[7]、args[8]筛选列
  if (length(args) == 8) {
    colnum1 <- grep(paste0(args[7], ""), colnames(expr))
    colnum2 <- grep(paste0(args[8], ""), colnames(expr))
    expr <- expr[, c(colnum1, colnum2)]
  }

  # 定义分组
  group <- rep(0, ncol(expr))
  sorted <- rep(0, 2)
  for (m in 7:(7 + 1)) {
    M <- m - 6
    sorted[M] <- args[m]
    colnum <- grep(args[m], colnames(expr))
    for (n in 1:length(colnum)) {
      group[colnum[n]] <- args[m]
    }
  }

  # 筛选表达量均值 > 0.5 的基因
  keep <- apply(expr, 1, mean) > 0.5
  expr <- expr[keep, ]
  data <- t(as.matrix(expr))

  # PCA 分析
  data.pca <- fast.prcomp(data, retx = TRUE, scale = FALSE, center = TRUE)
  a <- summary(data.pca)
  tmp <- a[["importance"]]
  pro1 <- as.numeric(sprintf("%.3f", tmp[2, 1])) * 100
  pro2 <- as.numeric(sprintf("%.3f", tmp[2, 2])) * 100

  pc <- as.data.frame(a[["x"]])
  pc$group <- group
  pc$names <- rownames(pc)

  xlab <- paste0("PC1 (", pro1, "%)")
  ylab <- paste0("PC2 (", pro2, "%)")
 # 确保加载依赖包
library(ggplot2)
library(ggrepel)
library(dplyr)  # 用于分组筛选

  # 绘制 PCA 图
  pca_plot <- ggplot(pc, aes(PC1, PC2)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
    stat_ellipse(geom = "polygon", alpha = 0.1, 
                aes(fill = group, color = group), 
                linewidth = 0.3) +
    geom_point(pch = 21, colour = "black", #点
              aes(fill = group)) +
    # 添加避免重叠的标签，字体大小设为4
    # 移除连线，仅显示不重叠或关键位置的标签
# 修改 geom_text_repel 部分：
    geom_text_repel(
      aes(label = names),
      size = 1.8,          # 缩小字体，减少空间占用
      max.overlaps = 50,   # 允许更多标签尝试显示（默认10，调大）
      force = 20,          # 增强排斥力（默认1，调大至20+）
      box.padding = 0.6,   # 标签与边界的间距
      point.padding = 0.6, # 标签与点的间距
      segment.colour = NA, # 移除连接线
      direction = "both",  # 允许标签向任意方向移动（上下左右）
      seed = 123,          # 固定随机种子，布局可复现
      data = pc %>%        # 仅显示组内 PC1 或 PC2 极值的点
        group_by(group) %>% 
        filter(PC1 %in% range(PC1) | PC2 %in% range(PC2))
    ) +
    labs(x = xlab, y = ylab) +
    theme_bw(base_size = 12) +  
    theme(
      panel.grid = element_blank(), 
      axis.text = element_text(size = 8, color = "black"),
      axis.title = element_text(size = 10, color = "black", face = "bold"),
      legend.title = element_text(size = 10, color = "black", face = "bold"), #图例
      legend.text = element_text(size = 8, color = "black"),            #图例    
      legend.background = element_rect(fill = "white",                      
                                      colour = "black", 
                                      linewidth = 0.3),
      legend.key.size = unit(0.6, "lines"),   #图例                                
      legend.position = "right",                                            
      legend.box = "vertical",                                              
      legend.box.margin = margin(0, 0, 0, 0, "cm"),                        
      plot.margin = margin(t = 0.3, r = 1.2, b = 0.3, l = 0.3, "cm"),
      aspect.ratio = 1
    ) +
    labs(fill = "Group", color = "Group") +
    geom_rug(aes(colour = group), linewidth = 0.3) +
    scale_color_aaas(limits = sorted) +
    scale_fill_aaas(limits = sorted)

  # 保存 PCA 图像
  out1 <- paste(sorted, collapse = "-vs-")
  out2 <- paste0("PCA.", out1)
  ggsave(#图
    filename = paste0(out2, ".pdf"), 
    plot = pca_plot, 
    width = 3,          
    height = 3,         
    units = "in", 
    dpi = 600
  )

  # t-SNE 分析
  data_unique <- unique(data)
  set.seed(12345)
  tsne_out <- Rtsne(as.matrix(data_unique), perplexity = 1)

  dm <- as.data.frame(tsne_out$Y)
  dm$group <- group[match(rownames(data_unique), rownames(data))]
  dm$names <- rownames(data_unique)
  colnames(dm)[1:2] <- c("DM1", "DM2")

  # t-SNE 绘图
  tsne <- ggplot(dm, aes(DM1, DM2)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_point(size = 2, aes(colour = group)) +
    labs(x = "Dimension 1", y = "Dimension 2") +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(size = 10, color = "black"),
      axis.title = element_text(size = 11, color = "black", face = "bold"),
      legend.title = element_text(size = 12, color = "black", face = "bold"),
      legend.text = element_text(size = 11, color = "black"),
      legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.3),
      legend.key.size = unit(0.5, "lines"),
      legend.position = "right",
      plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
      aspect.ratio = 1  
    ) +
    labs(color = "Group") +
    stat_ellipse(geom = "polygon", alpha = 0.1, linetype = "dashed", aes(fill = group, colour = group), linewidth = 0.3) +
    geom_rug(aes(color = group), linewidth = 0.3) +
    scale_color_nejm(limits = sorted, guide = "none") +
    scale_fill_nejm(limits = sorted, guide = "none")

  out3 <- paste0("t-sne.", out1)
  ggsave(paste0(out3, ".pdf"), tsne, width = 5.5, height = 5.5, units = "in", dpi = 600)

  # 层次聚类分析
  method <- c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")
  for (i in 1:length(method)) {
    hc <- hclust(dist(data), method = method[i])

    # 颜色映射
    color <- rep(0, length(labels(hc)))
    for (m in 7:(7 + 1)) {
      M <- m - 6
      colnum <- grep(args[m], labels(hc))
      color[colnum] <- pal_aaas("default")(2)[M]  # 调整为2组颜色
    }

    dend <- as.dendrogram(hc) %>%
      set("labels_col", color) %>%
      set("leaves_pch", 17) %>%
      set("labels_cex", 0.8) %>%  
      set("leaves_cex", 0.8) %>%
      set("leaves_col", color)

    out4 <- paste0("hclust_", method[i], ".", out1)
    rw <- length(labels(hc)) * 0.2  
    # pdf(file = paste0(out4, ".pdf"), width = max(rw, 2.7), height = max(rw, 4)) 
    pdf(file = paste0(out4, ".pdf"), width = 6, height = 3.5)   
    par(mar = c(12, 2.5, 0.5, 0), oma = c(0, 1, 0, 0))  

    dend_height <- attr(dend, "height")
    ylim_max <- dend_height * 1.2  
    plot(dend, ylab = "Height", cex.axis = 0.9, cex.lab = 1.0, ylim = c(0, ylim_max))
    dev.off()
  }