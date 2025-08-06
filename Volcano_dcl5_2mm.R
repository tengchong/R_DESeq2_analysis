# 安装和加载所需包（如果未安装）
install.packages("ggplot2")
install.packages("ggrepel")

library(ggplot2)
library(ggrepel)

# 读取数据（修改为你的文件名）
deg <- read.csv("/Users/tc-macpro/Downloads/DEG_dcl5_2_fertile_vs_ms_2_0mm.csv", row.names = 1)

# 替换 NA 的 padj 为 1（避免绘图报错）
deg$padj[is.na(deg$padj)] <- 1

# 添加显著性标签
deg$significance <- "Not Sig"
deg$significance[deg$padj < 0.05 & deg$log2FoldChange > 1] <- "Up"
deg$significance[deg$padj < 0.05 & deg$log2FoldChange < -1] <- "Down"

# 火山图绘制
ggplot(deg, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Volcano Plot",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted p-value",
    color = "Significance"
  )

# 选 padj 最小的前10个显著基因进行标注
top10 <- deg[deg$padj < 0.05 & abs(deg$log2FoldChange) > 1, ]
top10 <- top10[order(top10$padj), ][1:min(10, nrow(top10)), ]

ggplot(deg, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_text_repel(data = top10, aes(label = rownames(top10)), size = 4) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Volcano Plot (Top Genes Labeled)",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted p-value",
    color = "Significance"
  )

# 替换 NA，避免报错
res_2_0mm$padj[is.na(res_2_0mm$padj)] <- 1

# 筛选显著差异：padj < 0.05 且 |log2FC| > 1
deg <- subset(res_2_0mm, padj < 0.05 & abs(log2FoldChange) > 1)

# 标记方向
deg$direction <- ifelse(deg$log2FoldChange > 1, "Up",
                        ifelse(deg$log2FoldChange < -1, "Down", "NS"))

# 统计上下调数量
table(deg$direction)

# （参考你之前的变量名：deg）
deg_filtered <- deg[, c("log2FoldChange", "padj", "direction")]
deg_filtered$gene_id <- rownames(deg_filtered)  # 把行名变为列

# 保存为过滤过的 DEG 表
write.csv(deg_filtered, file = "/Users/tc-macpro/Downloads/filtered_DEG_dcl5_2_fertile_vs_ms_2_0mm.csv", row.names = FALSE)
4

# 读取 GFF3（过滤 gene 层级）
gff_lines <- readLines("/Users/tc-macpro/Downloads/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3")
gff_lines <- gff_lines[!grepl("^#", gff_lines)]
gff_split <- strsplit(gff_lines, "\t")
gff_df <- do.call(rbind, gff_split)
colnames(gff_df) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")

gff_genes <- gff_df[gff_df[, "type"] == "gene", ]

extract_attr <- function(attr_vec, key) {
  sapply(attr_vec, function(x) {
    match <- regmatches(x, regexpr(paste0(key, "=([^;]+)"), x))
    if (length(match) == 0 || match == "") {
      return(NA)
    } else {
      return(sub(paste0(key, "="), "", match))
    }
  })
}
gene_annot <- data.frame(
  gene_id      = extract_attr(gff_genes[, "attributes"], "ID"),
  alias        = extract_attr(gff_genes[, "attributes"], "Alias"),
  uniprot_id   = extract_attr(gff_genes[, "attributes"], "UNIPROT_ID"),
  uniprot_name = extract_attr(gff_genes[, "attributes"], "UNIPROT_NAME"),
  go_id        = extract_attr(gff_genes[, "attributes"], "GO_ID"),
  go_name      = extract_attr(gff_genes[, "attributes"], "GO_NAME"),
  stringsAsFactors = FALSE
)

deg_filtered_df <- as.data.frame(deg_filtered)

deg_annot <- merge(deg_filtered_df, gene_annot, by = "gene_id", all.x = TRUE)
write.csv(deg_annot, file = "/Users/tc-macpro/Downloads/filtered_DEG_annot_dcl5_2_fertile_vs_ms_2_0mm.csv", row.names = FALSE)
