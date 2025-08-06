# install DESeq2
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

library(DESeq2)
library(readr)

# 读取 counts 和 metadata
counts <- read.csv("/Users/tc-macpro/Downloads/gene_count_matrix_dcl5_2.csv", row.names = 1)
colnames(counts)

metadata <- read.csv("/Users/tc-macpro/Downloads/metadata_dcl5_2.csv", header = TRUE, row.names = 1)

# 确保列名和 metadata row 名匹配
all(colnames(counts) == rownames(metadata))  # 应该是 TRUE

# 创建新的 condition 列，用于 PCA 和 DEG
metadata$condition <- paste(metadata$fertility, metadata$stage, sep = "_")

# 检查
head(metadata)

# 构建 DESeq2 对象
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ condition)

# 设置参考组为 fertile_0.5mm
dds$condition <- relevel(dds$condition, ref = "fertile_2_0mm")

# 差异分析
dds <- DESeq(dds)

# 比较 1: fertile_2_0mm vs ms_2_0mm
res_2_0mm <- results(dds, contrast = c("condition", "ms_2_0mm", "fertile_2_0mm"))
write.csv(as.data.frame(res_2_0mm), file = "/Users/tc-macpro/Downloads/DEG_dcl5_2_fertile_vs_ms_2_0mm.csv")

# 比较 2: fertile_1_5mm vs ms_1_5mm
res_1_5mm <- results(dds, contrast = c("condition", "ms_1_5mm", "fertile_1_5mm"))
write.csv(as.data.frame(res_1_5mm), file = "/Users/tc-macpro/Downloads/DEG_dcl5_2_fertile_vs_ms_1_5mm.csv")

# 比较 3: fertile_2_5mm vs ms_2_5mm
res_2_5mm <- results(dds, contrast = c("condition", "ms_2_5mm", "fertile_2_5mm"))
write.csv(as.data.frame(res_2_5mm), file = "/Users/tc-macpro/Downloads/DEG_dcl5_2_fertile_vs_ms_2_5mm.csv")

# 比较 1: fertile_3_0mm vs ms_3_0mm
res_3_0mm <- results(dds, contrast = c("condition", "ms_3_0mm", "fertile_3_0mm"))
write.csv(as.data.frame(res_3_0mm), file = "/Users/tc-macpro/Downloads/DEG_dcl5_2_fertile_vs_ms_3_0mm.csv")





