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

dds <- dds[rowSums(counts(dds)) > 10, ]  # 删除低表达基因

vsd <- vst(dds, blind = TRUE)  # 可替换为 rlog()

# 3. 提取 PCA 数据
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

library(ggplot2)
library(RColorBrewer)

# 分组定义
pcaData$fertility <- ifelse(grepl("ms", pcaData$group), "ms", "fertile")
pcaData$stage <- gsub(".*_(.*)", "\\1", pcaData$group)

# 自定义颜色映射（可根据具体需要手动调）
library(RColorBrewer)

# 拿出更深的蓝色和紫色，跳过最浅的两档
fertile_colors <- brewer.pal(9, "Reds")[c(4, 6, 8, 9)]   # 对应 0.5mm, 1.0mm, 1.5mm
ms_colors      <- brewer.pal(9, "Purples")[c(4, 6, 8, 9)]

# 手动设置颜色对应 group（确保顺序对应 group levels）
group_levels <- sort(unique(pcaData$group))  # 例如 fertile_0.5mm, ..., ms_1.5mm
group_colors <- c(fertile_colors, ms_colors)
names(group_colors) <- group_levels
names(group_colors) <- group_levels
ggplot(pcaData, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  scale_color_manual(values = group_colors) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )