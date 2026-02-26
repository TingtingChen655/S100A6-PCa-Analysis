library(dplyr)
library(limma)
library(ggplot2)
library(ggrepel)
library(stringr)
setwd("")
expr<-read.table("gene.FPKM.matrix.annotation.txt",header = T,stringsAsFactors = F,sep = "\t")
expr<-expr[-which(expr$Gene.name==""),]
expr<-expr[-unique(which(is.na(expr),arr.ind = T)[,1]),]
expr <- expr %>%
  mutate(across(-Gene.name, as.numeric))
expr_mean <- aggregate(. ~ Gene.name, data = expr, FUN = mean)
rownames(expr_mean)<-expr_mean$Gene.name
expr_mean<-expr_mean[,-1]
expr_mean<-expr_mean[,c(1:6,7:12)]
# expr_mean<-expr_mean[,-c(2,11)]
expr_mean <- log2(expr_mean + 1)
expr_mean<-expr_mean[-which(apply(expr_mean, 1, function(x){length(which(x==0))})==12),]
group <- factor(c(rep("Stimulated", 6), rep("Control", 6)))
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
contrast.matrix <- makeContrasts(Stimulated - Control, levels = design)
fit <- lmFit(expr_mean, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
res <- topTable(fit2, coef = 1, n = Inf, sort.by = "P")
DEGs<-rownames(res)[which(res$P.Value < 0.05 & abs(res$logFC) > 1)]
#绘制差异表达热图
DEGs_expr<-expr_mean[DEGs,]
DEGs_expr<-log2(DEGs_expr+1)
pheatmap::pheatmap(as.matrix(DEGs_expr),scale = "row")

res$Gene <- rownames(res)
logFC_cutoff <- 1
P_cutoff <- 0.05
res <- res %>%
  mutate(change = case_when(
    P.Value < P_cutoff & logFC > logFC_cutoff ~ "UP",
    P.Value < P_cutoff & logFC < -logFC_cutoff ~ "DOWN",
    TRUE ~ "NOT"
  ))
top_genes <- res %>%
  filter(change != "NOT") %>%
  arrange(P.Value) %>%
  head(30) # 标记前20个，避免画面太满
res$label <- ifelse(res$Gene %in% top_genes$Gene, res$Gene, "")

ggplot(res, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(data = filter(res, change == "NOT"), 
             color = "#DCDCDC", alpha = 0.4, size = 0.5) +
  geom_point(data = filter(res, change == "DOWN"), 
             color = "#105186", alpha = 0.8, size = 3.0) +
  geom_point(data = filter(res, change == "UP"), 
             color = "#FF4500", alpha = 0.8, size = 3.0) +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), lty = "dashed", color = "#666666", lwd = 0.3) +
  geom_hline(yintercept = -log10(P_cutoff), lty = "dashed", color = "#666666", lwd = 0.3) +
  geom_text_repel(aes(label = label),
                  size = 3.8,
                  fontface = "bold.italic",
                  box.padding = unit(0.6, "lines"),
                  point.padding = unit(0.5, "lines"),
                  segment.color = "black",
                  segment.size = 0.3,
                  max.overlaps = Inf) +
  scale_x_continuous(limits = c(-max(abs(res$logFC)), max(abs(res$logFC)))) + # 让 0 对齐中心
  labs(title = "Differential Expression Analysis",
       x = expression(log[2](Fold~Change)),
       y = expression(-log[10](P-value))) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11, color = "black")
  )
#############
#功能富集分析
###########
#GO
library(clusterProfiler)
res_GO<-enrichGO(gene = DEGs,OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont = "ALL",pAdjustMethod = "none",qvalueCutoff = 1)
a<-as.data.frame(res_GO)
plot_data <- as.data.frame(res_GO) %>%
  group_by(ONTOLOGY) %>%
  slice_max(order_by = -pvalue, n = 10) %>% 
  ungroup() %>%
  mutate(logP = -log10(pvalue)) %>%
  # 自动处理过长的 GO 名称（超过 50 个字符自动换行，防止挤出图片）
  # mutate(Description = str_wrap(Description, width = 50)) %>%
  arrange(ONTOLOGY, logP) %>%
  mutate(Description = factor(Description, levels = unique(Description)))

ggplot(plot_data, aes(x = logP, y = Description)) +
  geom_segment(aes(x = 0, xend = max(logP) * 1.1, y = Description, yend = Description),
               # color = "#F1F8E9", size = 12, lineend = "round") +
               color = "#F0F7FF", size = 10, lineend = "round") +
  geom_segment(aes(x = 0, xend = logP, y = Description, yend = Description, color = logP),
               size = 10, lineend = "round",alpha = 0.6) +
  geom_text(aes(x = 0.1, y = Description, label = Description),
            hjust = 0, size = 3.2, color = "black", fontface = "bold.italic") +
  # scale_color_gradientn(colors = c("#F7FCB9", "#ADDD8E", "#78C679", "#238B45", "#00441B")) +
  scale_color_gradientn(colors = c("#E7F5FF", "#D0EBFF", "#A5D8FF", "#74C0FC", "#4DABF7")) +
  facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "-log10 (P-value)", y = NULL) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),          
    axis.text.y = element_blank(),          
    axis.ticks.y = element_blank(),
    axis.line.x = element_line(color = "black"),
    # 分面标签（BP/CC/MF）的样式微调
    strip.background = element_rect(fill = "#0288D1", color = NA),
    # strip.background = element_rect(fill = "#238B45", color = NA),
    strip.text = element_text(color = "white", face = "bold", size = 10),
    panel.spacing = unit(1, "lines"),      
    legend.position = "none") +coord_cartesian(clip = "off")
############
#KEGG
gene_convert <- bitr(DEGs, 
                     fromType = "SYMBOL",   # 输入是基因名
                     toType = "ENTREZID",   # 输出是 Entrez ID
                     OrgDb = "org.Hs.eg.db")
res_kegg <- enrichKEGG(gene = gene_convert$ENTREZID,
                       organism = 'hsa',
                       pvalueCutoff = 0.05,pAdjustMethod = "none",qvalueCutoff = 1)
a<-as.data.frame(res_kegg)
##############
#REACTOME
my_gmt <- read.gmt("c2.cp.reactome.v2026.1.Hs.symbols.gmt") 
res_custom <- enricher(gene = DEGs,
                       TERM2GENE = my_gmt,
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 1,minGSSize = 1,pAdjustMethod = "none")
a<-as.data.frame(res_custom)
plot_data <- as.data.frame(res_custom) %>%
  top_n(22, wt = -pvalue) %>% 
  mutate(logP = -log10(pvalue)) %>%
  arrange(logP) %>%
  mutate(Description = factor(Description, levels = Description))
win.graph(width = 60,height = 60)
#棒棒糖图
ggplot(plot_data, aes(x = logP, y = reorder(Description, logP))) +
  geom_segment(aes(x = 0, xend = logP, y = Description, yend = Description), 
               color = "grey90", size = 1) +
  geom_point(aes(color = logP), size = 10, alpha = 0.3) +
  geom_point(aes(color = logP), size = 6) +
  scale_color_gradientn(colors = c("#BFDEFF", "#85C1E9", "#3498DB", "#2874A6", "#154360"))+
  # scale_color_gradientn(colors = c("#F7FCB9", "#ADDD8E", "#78C679", "#238B45", "#00441B")) +
  theme_test() +
  labs(x = "-log10(P-value)", y = NULL, color = "Significance") +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 11, color = "black", face = "bold"),
    panel.grid = element_blank())

#柱状图
ggplot(plot_data, aes(x = logP, y = reorder(Description, logP), fill = logP)) +
  geom_col(width = 0.7, color = "white", size = 0.3, alpha = 0.9) +
  scale_fill_gradientn(colors = c("#F7FCB9", "#ADDD8E", "#78C679", "#238B45", "#00441B")) +
  scale_color_gradientn(colors = c("#F7FCB9", "#ADDD8E", "#78C679", "#238B45", "#00441B")) +
  labs(x = "-log10 (P-value)", y = NULL, fill = "-log10(P)") +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", linetype = "dashed"), 
    axis.text.y = element_text(size = 11, color = "black", face = "italic"),
    axis.line.x = element_line(color = "black"),
    legend.position = "right",
    plot.margin = margin(20, 20, 20, 20)) +
  coord_cartesian(clip = "off")
#文字内嵌胶囊图
ggplot(plot_data, aes(x = logP, y = Description)) +
  geom_segment(aes(x = 0, xend = max(logP) * 1.1, y = Description, yend = Description),
               color = "#F5F5F5", size = 12, lineend = "round") +
  geom_segment(aes(x = 0, xend = logP, y = Description, yend = Description, color = logP),
               size = 12, lineend = "round",alpha = 0.6) +
  geom_text(aes(x = 0.1, y = Description, label = Description),
            hjust = 0, size = 3.8, color = "black", fontface = "bold.italic") +
  scale_color_gradientn(colors = c("#F7FCB9", "#ADDD8E", "#78C679", "#238B45", "#00441B")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "-log10 (P-value)", y = NULL) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),         
    axis.text.y = element_blank(),         
    axis.ticks.y = element_blank(),
    axis.line.x = element_line(color = "black"),
    legend.position = "none",
    plot.margin = margin(20, 20, 20, 20))
##############
#PCA主成分分析
pca_input <- t(log2(expr_mean+1))
pca_input <- t(expr_mean)
sample_info <- data.frame(
  sample = rownames(pca_input),
  group = factor(c(rep("PC3_S", 6), rep("PC3_con", 6))) 
)
pca_res <- prcomp(pca_input, scale. = TRUE)

# 提取主成分得分 (Coordinates)
pca_data <- as.data.frame(pca_res$x)
pca_data$sample <- rownames(pca_data)
pca_data$group <- sample_info$group

# 计算每个主成分的方差贡献率 (%)
var_explained <- pca_res$sdev^2 / sum(pca_res$sdev^2) * 100
pc1_var <- round(var_explained[1], 2)
pc2_var <- round(var_explained[2], 2)
my_colors <- c("PC3_S" = "#E67E22", "PC3_con" = "#5DADE2")
my_shapes <- c("PC3_S" = 21, "PC3_con" = 24) # 21是实心圆，24是实心三角

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, fill = group, shape = group)) +
  # 绘制中心虚线
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  # 绘制点：使用 fill 填充颜色，color 设置边框
  geom_point(size = 4, color = "black", stroke = 0.5) +
  # 添加样本标签
  geom_text_repel(aes(label = sample), 
                  size = 3.5, 
                  color = "black",
                  box.padding = 0.3,
                  show.legend = FALSE) +
  # 设置颜色和形状
  scale_fill_manual(values = my_colors) +
  scale_shape_manual(values = my_shapes) +
  # 设置坐标轴标签（包含贡献率）
  labs(x = paste0("PC1(", pc1_var, "%)"),
       y = paste0("PC2(", pc2_var, "%)"),
       title = "PCA analysis") +
  # 主题美化：白底、黑框、无背景网格
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    legend.title = element_blank(),
    legend.position = c(0.9, 0.9), # 图例放在右上角
    legend.background = element_blank()
  ) +
  # 调整图例中点的外观
  guides(fill = guide_legend(override.aes = list(shape = c(21, 24), size = 4)))

# 显示图片
print(p_pca)
#############
#表达热图
library(ggtree)
library(ggtreeExtra)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize) # 用于颜色映射
mat_scaled <- t(scale(t(DEGs_expr))) 
col_fun = colorRamp2(c(-2, 0, 2), c("#67B3DA", "white", "#105186"))
ann_col = data.frame(Group = c(rep("Stimulated", 6), rep("Control", 6)))
col_ann = list(Group = c("Stimulated" = "#E67E22", "Control" = "#5DADE2"))
column_ha = HeatmapAnnotation(df = ann_col, col = col_ann, show_annotation_name = FALSE)
res_sorted <- res[DEGs,] %>% arrange(desc(logFC)) # 按 logFC 从大到小排序
mat_final <- mat_scaled[res_sorted$Gene, ]
library(ComplexHeatmap)
library(circlize)
col_red_blue = colorRamp2(c(-2, 0, 2),c("#2166AC", "#FFFFFF", "#B2182B"))
Heatmap(mat_final, 
        name = "Z-score",
        cluster_rows = FALSE,     
        cluster_columns = FALSE,  
        col = col_red_blue,
        rect_gp = gpar(col = "white", lwd = 1), 
        right_annotation = rowAnnotation(
          log2FC = anno_barplot(
            res_sorted$logFC, 
            gp = gpar(fill = ifelse(res_sorted$logFC > 0, "#B2182B", "#2166AC"), col = "white"), 
            width = unit(2, "cm"))),
        column_title = "Genes ordered by log2FC",
        row_names_gp = gpar(fontsize = 10, fontface = "bold.italic"),
        column_names_rot = 45,
        heatmap_legend_param = list(
          title = "Z-score",
          at = c(-2, 0, 2),
          labels = c("Low", "0", "High"),
          grid_height = unit(6, "mm")))