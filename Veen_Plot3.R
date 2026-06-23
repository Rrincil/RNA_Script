# 下载包（已注释，假设已安装）
#install.packages("readxl")
#install.packages("UpSetR")
rm(list = ls())
# 载入包
library(UpSetR)
library(openxlsx)
library(RColorBrewer)
library(readxl)
library(data.table)

# 数据集文件路径列表
file_paths <- c(
  "data/GSE120977/gsea_report_for_control_1740761010238.xls",
  # "data/GSE133566/gsea_report_for_WT_Chow_1740993812565.xls",
  "data/GSE164085/gsea_report_for_control_1740853152557.xls",
  "data/GSE159676/gsea_report_for_healthy_1741003644405.xls"
)

# 数据集名称列表
dataset_names <- c(
  "GSE120977",
  # "GSE133566",
  "GSE164085",
  "GSE159676"
)

# 读取数据并筛选 NOM p-val 小于 0.05 的行，存储在列表中
data_list <- lapply(file_paths, function(file_path) {
  data <- fread(file_path)
  return(data[data$`NOM p-val` < 0.05, ])
})
names(data_list) <- dataset_names

# 找到所有数据集中 NAME 列的最大长度
max_length <- max(sapply(data_list, function(x) length(x$NAME)))

# 创建一个空的数据框，预先分配足够的行数并设置列名
G_UpGene <- data.frame(matrix(ncol = length(dataset_names), nrow = max_length))
colnames(G_UpGene) <- dataset_names
G_UpGene[] <- lapply(G_UpGene, as.character)

# 为数据框添加列数据
for (i in seq_along(data_list)) {
  G_UpGene[[dataset_names[i]]] <- c(data_list[[i]]$NAME, rep(NA, max_length - length(data_list[[i]]$NAME)))
}

data <- G_UpGene

# 以下是集合图相关代码，保持不变
# 调整与美化后的集合图#
upset(fromList(data))#基础图

#高亮显示特定几个集合的交集
upset(fromList(data),      
      nsets=length(data),#显示数据集的所有数据,nsets = 数值调整可视化数据集数量
      nintersects=180,#显示前多少个
      sets=c("GSE120977","GSE164085","GSE159676"
      ), # 指定集合或用keep.order = TRUE保持集合按输入的顺序排序
      number.angles = 0, #交互集合柱状图的柱标倾角
      point.size=2.5, #图中点的大小
      line.size=1, #图中连接线粗细
      mainbar.y.label="Intersection size", #y轴的标签
      main.bar.color = 'black', #y轴柱状图颜色
      matrix.color="black", #x轴点的颜色
      sets.x.label="Set size",   #x轴的标签
      sets.bar.color=brewer.pal(3,"Set3"),#x轴柱状图的颜色;Set1中只有9个颜色，Set3中有12个颜色，Paired中有12个颜色
      mb.ratio = c(0.7, 0.3), #bar plot和matrix plot图形高度的占比
      # order.by = "degree",
      order.by = "freq", #y轴矩阵排序,如"freq"频率，"degree"程度
      text.scale=c(1.5,1.5,1.5,1.5,1.5,1), #6个参数intersection size title（y标题大小）,intersection size tick labels（y刻度标签大小）, set size title（set标题大小）, set size tick labels（set刻度标签大小）, set names（set 分类标签大小）, numbers above bars（柱数字大小）的设置
      shade.color="red", #图中阴影部分的颜色
      queries=list(list(query=intersects,params=list("GSE120977","GSE164085","GSE159676"),color="red",active=T)#设置自己想要展示的特定组的交集，通过queries参数进行设置，需要展示几个关注组合的颜色，就展示几个
      )
      # queries=list(list(query=intersects,params=list("GSE18224","GSE18801"),color="red",active=T),#设置自己想要展示的特定组的交集，通过queries参数进行设置，需要展示几个关注组合的颜色，就展示几个
      #              list(query=intersects,params=list("GSE18801","GSE24489"),color="blue",active=T),
      #              list(query=intersects,params=list("GSE18801","GSE37597"),color="green",active=T),
      #              list(query=intersects,params=list("GSE37597","GSE56348_Sham-vs-CH","GSE56348_Sham-vs-HF"),color="yellow",active=T),
      #              list(query=intersects,params=list("GSE84142","GSE56348_Sham-vs-CH","GSE56348_Sham-vs-HF"),color="purple",active=T),
      #              list(query=intersects,params=list("GSE201764","GSE56348_Sham-vs-CH","GSE56348_Sham-vs-HF"),color="orange",active=T) )
)

# 简化交集计算部分，使用循环处理
pathway_lists <- G_UpGene
pathways <- Reduce(intersect, pathway_lists)
# pathways <- as.data.frame(pathways)

# 其他交集计算部分，可根据需要继续优化

# 简化创建 _intersect 数据框部分
intersect_data_list <- list()
for (i in seq_along(dataset_names)) {
  subset_data <- data_list[[i]][data_list[[i]]$NAME %in% pathways,]
  intersect_data <- data.frame(
    Person = dataset_names[i],
    Day = subset_data$NAME,
    Value = subset_data$`NOM p-val`,
    NES = subset_data$NES
  )
  intersect_data_list[[i]] <- intersect_data
}

full_intersec_data <- do.call(rbind, intersect_data_list)
# full_intersec_data <- full_intersec_data[!is.na(full_intersec_data$Value),]



# 保存数据部分，保持不变
# 大于等于3条以上才有交集
pathways <- as.data.frame(pathways)
colnames(pathways)[1] <- "pathways"
write.table(pathways,"./AllPathways_intersect_pathways.xlsx",sep = '\t',quote = F,row.names = F)

write.table(full_intersec_data,"./full_intersec_data.xlsx",sep = '\t',quote = F,row.names = F)

All_DownPathways <- as.data.frame(data)
write.table(All_DownPathways,"./All_DownPathways.xlsx",sep = '\t',quote = F,row.names = F)

AllPathways_UpPathways <- as.data.frame(data)
# 去除NA
AllPathways_UpPathways[is.na(AllPathways_UpPathways)] <- ""

write.table(AllPathways_UpPathways,"./AllPathways.xlsx",sep = '\t',quote = F,row.names = F)

All_AllPathways <- as.data.frame(data)
# 去除NA
All_AllPathways2 <- All_AllPathways
All_AllPathways2[is.na(All_AllPathways2)] <- ""
write.table(All_AllPathways2,"./AllPathways_AllPathways.xlsx",sep = '\t',quote = F,row.names = F)