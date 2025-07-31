# install DESeq2
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

library(DESeq2)
library(readr)

# 读取 counts 和 metadata
counts <- read.csv("/Users/tc-macpro/Downloads/gene_count_matrix.csv", row.names = 1)
colnames(counts)

metadata <- read.csv("/Users/tc-macpro/Downloads/metadata_mtm.csv", header = TRUE, row.names = 1)

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
dds$condition <- relevel(dds$condition, ref = "fertile_0_5mm")

# 差异分析
dds <- DESeq(dds)

# 比较 1: fertile_0.5mm vs ms_0.5mm
res_0.5mm <- results(dds, contrast = c("condition", "ms_0_5mm", "fertile_0_5mm"))
write.csv(as.data.frame(res_0.5mm), file = "DEG_fertile_vs_ms_0_5mm.csv")

# 比较 2: fertile_1.0mm vs ms_1.0mm
res_1.0mm <- results(dds, contrast = c("condition", "ms_1_0mm", "fertile_1_0mm"))
write.csv(as.data.frame(res_1.0mm), file = "DEG_fertile_vs_ms_1.0mm.csv")

# 比较 3: fertile_1.5mm vs ms_1.5mm
res_1.5mm <- results(dds, contrast = c("condition", "ms_1_5mm", "fertile_1_5mm"))
write.csv(as.data.frame(res_1.5mm), file = "DEG_fertile_vs_ms_1.5mm.csv")