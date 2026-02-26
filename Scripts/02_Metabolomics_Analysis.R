#代谢分析
meta<-read.table("all_metab_data1.txt",header = T,stringsAsFactors = F,sep = "\t", quote = "",
                 fill = T,check.names = F,row.names = 1)
meta <- log2(meta + 1)
meta<-meta[-which(apply(meta, 1, function(x){length(which(x==0))})==12),]
group <- factor(c(rep("Control", 6), rep("Stimulated", 6)))
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
contrast.matrix <- makeContrasts(Stimulated - Control, levels = design)
fit <- lmFit(meta, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
res <- topTable(fit2, coef = 1, n = Inf, sort.by = "P")
DEMs<-rownames(res)[which(res$P.Value < 0.05 & abs(res$logFC) > 1)]
DEMs<-DEMs[-which(grepl(pattern = "metab_",x = DEMs,ignore.case = T))]
#绘制差异表达热图
DEMs_expr<-meta[DEMs,]
DEMs_expr<-log2(DEMs_expr+1)
mat_scaled <- t(scale(t(DEMs_expr))) 
col_fun = colorRamp2(c(-2, 0, 2), c("#67B3DA", "white", "#105186"))
ann_col = data.frame(Group = c( rep("Control", 6),rep("Stimulated", 6)))
col_ann = list(Group = c("Stimulated" = "#E67E22", "Control" = "#5DADE2"))
column_ha = HeatmapAnnotation(df = ann_col, col = col_ann, show_annotation_name = FALSE)
res_sorted <- res[DEMs,] %>% arrange(desc(logFC)) # 按 logFC 从大到小排序
res_sorted$Gene<-rownames(res_sorted)
mat_final <- mat_scaled[res_sorted$Gene, ]
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

res$Gene <- rownames(res)
logFC_cutoff <- 1
P_cutoff <- 0.05
res <- res %>%
  mutate(change = case_when(
    P.Value < P_cutoff & logFC > logFC_cutoff ~ "UP",
    P.Value < P_cutoff & logFC < -logFC_cutoff ~ "DOWN",
    TRUE ~ "NOT"
  ))
top_genes <-DEMs
res$label <- ifelse(res$Gene %in% top_genes, res$Gene, "")

ggplot(res, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(data = filter(res, change == "NOT"), 
             color = "#DCDCDC", alpha = 0.4, size = 0.5) +
  geom_point(data = filter(res, change == "DOWN"), 
             color = "#105186", alpha = 0.8, size = 3.0) +
  geom_point(data = filter(res, change == "UP"), 
             color = "#FF4500", alpha = 0.8, size = 3.0) +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), lty = "dashed", color = "#666666", lwd = 0.3) +
  geom_hline(yintercept = -log10(P_cutoff), lty = "dashed", color = "#666666", lwd = 0.3) +
  # geom_text_repel(aes(label = label),
  #                 size = 3.8,
  #                 fontface = "bold.italic",
  #                 box.padding = unit(0.6, "lines"),
  #                 point.padding = unit(0.5, "lines"),
  #                 segment.color = "black",
  #                 segment.size = 0.3,
  #                 max.overlaps = Inf) +
  scale_x_continuous(limits = c(-max(abs(res$logFC)), max(abs(res$logFC)))) + # 让 0 对齐中心
  labs(title = "Differential Expression Analysis of metabolites",
       x = expression(log[2](Fold~Change)),
       y = expression(-log[10](P-value))) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11, color = "black")
  )
#代谢物功能富集分析
library(enrichmet)
ann<-read.table("all_metab_data.txt",header = T,stringsAsFactors = F,sep = "\t", quote = "",
                fill = T,check.names = F,row.names = 1)
ann_name<-ann[-which(grepl(pattern = "metab_",x = ann$Metabolite,ignore.case = T)),]
ann_kegg<-ann[-which(ann$`KEGG Compound ID`=="-"),]
meta_data<-ann_kegg
rownames(meta_data)<-meta_data$Metabolite
DEM_da<-meta_data[DEMs,]
inter_DEM<-DEM_da$`KEGG Compound ID`[-which(is.na(DEM_da$`KEGG Compound ID`))]
inter_DEM <- str_extract_all(inter_DEM, "C\\d+")
inter_DEM<- unique(unlist(inter_DEM))
backgroundList <- str_extract_all(ann_kegg$`KEGG Compound ID`, "C\\d+")
backgroundList <- unique(unlist(backgroundList))



library(KEGGREST)
library(clusterProfiler)
library(dplyr)
link_data <- keggLink("pathway", "compound")
kegg_map <- data.frame(
  Compound = names(link_data),
  Pathway = link_data
)
kegg_map$Compound <- gsub("cpd:", "", kegg_map$Compound)
kegg_map$Pathway <- gsub("path:", "", kegg_map$Pathway)
term2gene <- kegg_map[, c("Pathway", "Compound")]
pathway_names <- keggList("pathway")
term2name <- data.frame(
  Pathway = names(pathway_names),
  Name = pathway_names
)
term2name$Pathway <- gsub("path:", "", term2name$Pathway)
# 运行富集分析
enrich_res <- enricher(
  gene = inter_DEM,              
  universe = NULL,               
  TERM2GENE = term2gene,        
  TERM2NAME = term2name,         
  pvalueCutoff = 0.05,          
  qvalueCutoff = 1   ,
  pAdjustMethod = "none"
)

plot_data <- as.data.frame(enrich_res)
top_n <- min(15, nrow(plot_data))
plot_data <- plot_data %>% 
  arrange(pvalue) %>% 
  head(top_n)
plot_data$logP <- -log10(plot_data$pvalue)
plot_data <- plot_data %>% arrange(logP)
plot_data$Description <- factor(plot_data$Description, levels = plot_data$Description)
p<-ggplot(plot_data, aes(x = logP, y = reorder(Description, logP))) +
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
win.graph(width = 90,height = 60)
p

ggplot(plot_data, aes(x = logP, y = Description)) +
  geom_segment(aes(x = 0, xend = max(logP) * 1.1, y = Description, yend = Description),
               color = "#F5F5F5", size = 12, lineend = "round") +
  geom_segment(aes(x = 0, xend = logP, y = Description, yend = Description, color = logP),
               size = 12, lineend = "round",alpha = 0.6) +
  geom_text(aes(x = 0.1, y = Description, label = Description),
            hjust = 0, size = 3.8, color = "black", fontface = "bold.italic") +
  scale_color_gradientn(colors = c("#BFDEFF", "#85C1E9", "#3498DB", "#2874A6", "#154360"))+
  # scale_color_gradientn(colors = c("#F7FCB9", "#ADDD8E", "#78C679", "#238B45", "#00441B")) +
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
