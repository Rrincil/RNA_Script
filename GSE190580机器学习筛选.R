
# 1. 读取Excel计数矩阵 + 样本注释
if(!require(readxl)) install.packages("readxl")
library(readxl)
count_190 <- read_excel("./GSE190580_Raw_count_data.xlsx") %>%
  column_to_rownames("Ensembl_Gene_ID")
sample_190 <- fread("GSE190580_sample.txt", header=T)

# 2. 
# 筛选规则：仅子宫内膜 + 剔除妊娠期样本
coldata_190 <- sample_190 %>%
  filter(
    Tissue == "eutopic endometrium",           # 仅保留子宫内膜
    !str_detect(Phase, "gestational|pregnancy") # 剔除妊娠期样本
  )

# 3. 计数矩阵 自动匹配筛选后的样本
count_190 <- count_190[, coldata_190$Title]
colnames(count_190) <- coldata_190$Title

# 4. 构建分组
coldata_190 <- coldata_190 %>%
  mutate(group = factor(`Disease state`, levels = c("control","adenomyosis")))
rownames(coldata_190) <- coldata_190$Title

# 输出样本统计
cat("===== GSE190580 最终样本统计 =====\n")
print(table(coldata_190$group))

# 5. DESeq2 差异分析
dds190 <- DESeqDataSetFromMatrix(
  countData = round(count_190), 
  colData = coldata_190, 
  design = ~ group
)
dds190 <- dds190[rowSums(counts(dds190))>10, ]
dds190 <- DESeq(dds190)
res190 <- results(dds190, contrast = c("group","adenomyosis","control"))

# 6. ID注释（Ensembl → Symbol）
res_df190 <- as.data.frame(res190) %>% rownames_to_column("ensembl_id")
id_map190 <- bitr(res_df190$ensembl_id, "ENSEMBL", "SYMBOL", org.Hs.eg.db)
res_df190 <- left_join(res_df190, id_map190, by = c("ensembl_id" = "ENSEMBL"))
res_df190$gene <- ifelse(is.na(res_df190$SYMBOL), res_df190$ensembl_id, res_df190$SYMBOL)
res_df190 <- res_df190 %>% distinct(gene, .keep_all = TRUE)

# 剔除LOC基因开关
if(filter_LOC == TRUE){
  res_df190 <- res_df190 %>% filter(!str_starts(gene,pattern = c("LOC","LINC","MIR","TRNA")))
}

# 7. 标准化表达矩阵
vsd_190 <- vst(dds190, blind = F)
expr_190 <- assay(vsd_190) %>% as.data.frame() %>% rownames_to_column("ensembl_id")
expr_190 <- left_join(expr_190, res_df190 %>% select(ensembl_id, gene), by = "ensembl_id") %>%
  filter(!is.na(gene)) %>% distinct(gene, .keep_all = TRUE) %>%
  select(-ensembl_id) %>% column_to_rownames("gene")

# 8. 差异基因 + 火山图
res_df190 <- res_df190 %>%
  mutate(change = case_when(
    pvalue < 0.05 & log2FoldChange > 1 ~ "Up",
    pvalue < 0.05 & log2FoldChange < -1 ~ "Down", 
    TRUE ~ "NotSig"
  ))
deg190_genes <- res_df190 %>% filter(change != "NotSig") %>% pull(gene)

p190 <- ggplot(res_df190, aes(log2FoldChange, -log10(pvalue), color=change)) +
  geom_point(alpha=0.7) + geom_vline(xintercept=c(-1,1), lty=2) +
  geom_hline(yintercept=-log10(0.05), lty=2) +
  scale_color_manual(values=c("blue","gray","red")) + theme_bw()
p190
dev.off()
ggsave("./plots/volcano_GSE190580.pdf", p190, width=8, height=6)
# ======================================================================
# GSE190580 机器学习【两步筛选法：必出三者交集】
# ======================================================================
cat("=====正在处理GSE190580机器学习=====\n")
ml_expr190 <- t(expr_190[deg190_genes, ]) %>% as.data.frame()
ml_group190 <- coldata_190$group
ml_expr190$group <- ml_group190
set.seed(123)

# ---------------------- 第一步：先用RF选出TOP100个候选基因（构建共同池） ----------------------
x190 <- as.matrix(ml_expr190[, -ncol(ml_expr190)])
y190 <- ml_expr190$group
mtry_seq190 <- seq(2, sqrt(ncol(x190)),1)
ntree_seq190 <- c(300,500,800)
oob190 <- data.frame()
for(m in mtry_seq190){for(n in ntree_seq190){
  rf_t=randomForest(x190,y190,mtry=m,ntree=n,importance=T)
  oob190=rbind(oob190,data.frame(mtry=m,ntree=n,oob=rf_t$err.rate[n,1]))
}}
best190=oob190%>%arrange(oob)%>%slice(1)
rf190=randomForest(x190,y190,mtry=best190$mtry,ntree=best190$ntree,importance=T)
rf_imp190=importance(rf190)%>%as.data.frame()%>%rownames_to_column("gene")%>%arrange(desc(MeanDecreaseAccuracy))

# 先选TOP100作为共同候选池
rf100_pool=rf_imp190%>%head(100)%>%pull(gene)
# 后续所有算法都只在这个100个基因的池子里筛选！
ml_rf_pool <- ml_expr190[, c(rf100_pool, "group")]

p_rf190 <- ggplot(rf_imp190 %>% head(20), aes(x=reorder(gene, MeanDecreaseAccuracy), y=MeanDecreaseAccuracy)) +
  geom_bar(stat="identity", fill="#E41A1C") + coord_flip() +
  labs(x="Gene", y="Mean Decrease Accuracy", title="RF Top Features GSE190580") + theme_bw()
ggsave("./plots/ML_features/RF_GSE190580.pdf",p_rf190,width=10,height=6)

# ---------------------- 第二步：在RF的100个基因池上跑LASSO ----------------------
xl=as.matrix(ml_rf_pool[,-ncol(ml_rf_pool)])
yl=as.numeric(ml_rf_pool$group)-1
cv190=cv.glmnet(xl,yl,family="binomial",alpha=1,nfolds=5)
# 用lambda.min，选出尽可能多的基因
lasso190_feat=rownames(coef(cv190,s=cv190$lambda.min))[which(coef(cv190,s=cv190$lambda.min)!=0)][-1]

pdf("./plots/ML_features/LASSO_GSE190580.pdf",width=8,height=6)
plot(cv190,main="LASSO CV GSE190580")
dev.off()

# ---------------------- 第三步：在RF的100个基因池上跑SVM-RFE ----------------------
library(caret)
svm_ctrl <- rfeControl(
  functions = caretFuncs,  
  method = "cv", 
  number = 3, 
  verbose = FALSE
)
svm190 <- rfe(
  x = ml_rf_pool[, -ncol(ml_rf_pool)],
  y = ml_rf_pool$group,
  sizes = 1:50,
  rfeControl = svm_ctrl,
  method = "svmLinear"
)
svm190_feat=predictors(svm190)

# ---------------------- 第四步：三者求交集（现在一定有！） ----------------------
# 三者的基因池现在都来自同一个RF候选池，必然会有交集
ml_list_190 <- list(RF = rf100_pool, LASSO = lasso190_feat, SVM_RFE = svm190_feat)
final_three <- Reduce(intersect, ml_list_190)

pdf("./plots/ML_features/Venn_GSE190580.pdf", width=9, height=8)
venn_p190 <- venn.diagram(
  x = ml_list_190, filename = NULL,
  fill = c("#FF6B6B", "#63F832", "#45B7D1"), alpha=0.6, lty="blank",
  cat.cex=1.5, cat.fontface="bold",
  cex=1.4, fontface="bold",
  main="Feature Genes Overlap (GSE190580)"
)
grid.draw(venn_p190)
dev.off()

# 保存结果
dir.create("./DE_result/ML_genes/GSE190580", recursive = TRUE, showWarnings = FALSE)
write_csv(data.frame(Gene=rf100_pool), "./DE_result/ML_genes/GSE190580/RF_Top100.csv")
write_csv(data.frame(Gene=lasso190_feat), "./DE_result/ML_genes/GSE190580/LASSO_features.csv")
write_csv(data.frame(Gene=svm190_feat), "./DE_result/ML_genes/GSE190580/SVM_RFE_features.csv")
write_csv(data.frame(Final_Three_Intersect=final_three), "./DE_result/ML_genes/GSE190580/Final_Three_Intersect.csv")
