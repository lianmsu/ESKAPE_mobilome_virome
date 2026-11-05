# ===========================================================================
# ESKAPE_mobilome_virome -- analyses / 03_gene_annotation
# ARG prevalence, drug classes and ARG/MRG/VFG/DS relationships
# Manuscript      : Fig. 1C-D
# Original path   : 2ESKAPE_MGE/3genome_gene_annotation/1_ARG.R
#
# Inputs are read from CS_DATA (default <repo>/data) and figures are
# written to CS_WORK (default <repo>/results):
#     Sys.setenv(CS_DATA = "/path/to/data", CS_WORK = "/path/to/results")
# ===========================================================================

# ---------------------------------------------------------------------------
CS_HOME <- normalizePath(Sys.getenv("CS_HOME", unset = "."), mustWork = FALSE)
CS_DATA <- normalizePath(Sys.getenv("CS_DATA", unset = file.path(CS_HOME, "data")), mustWork = FALSE)
CS_WORK <- normalizePath(Sys.getenv("CS_WORK", unset = file.path(CS_HOME, "results")), mustWork = FALSE)
CS_EXT  <- normalizePath(Sys.getenv("CS_EXT", unset = file.path(CS_HOME, "data", "external")), mustWork = FALSE)
dir.create(CS_WORK, recursive = TRUE, showWarnings = FALSE)
# ---------------------------------------------------------------------------

setwd(file.path(CS_DATA, "3genome_gene_annotation"))

#Data prepare##########################################################################
ESKAPE_genome_info <- read.delim2("ESKAPE_genome_info.tsv")

ARG_num <- read.csv("rgiargs_4_type.csv",sep = "\t")
VFG_num <- read.csv("diamond_VFG_out3_type.csv",sep = "\t")
MRG_num <- read.csv("diamond_MRG_out3_type.csv",sep = "\t")

DF_info <- read.delim2("all_4538_ESKAPE_genomes_renamed_defense_finder_systems.tsv")

library(dplyr)
library(tidyr)
DF_info <- DS_info %>%
  mutate(
    genome = gsub("-.*","",DS_info$sys_id)
  )
DS_info <- DF_info[DF_info$activity == "Defense",]

DS_num <- table(DS_info$genome) %>% as.data.frame()
colnames(DS_num) <- c("genome","DS")

ESKAPE_genome_info$genome <- gsub("\\.","_",ESKAPE_genome_info$genome)
ESKAPE_genome_info <- ESKAPE_genome_info %>%
  left_join(ARG_num, by = "genome") %>%
  left_join(VFG_num, by = "genome") %>%
  left_join(MRG_num, by = "genome") %>%
  left_join(DS_num, by = "genome")
ESKAPE_genome_info <- ESKAPE_genome_info %>%
  mutate(across(19:22, ~ replace_na(., 0)))

ESKAPE_genome_info$ARG_density <- ESKAPE_genome_info$ARG / ESKAPE_genome_info$Size *1000000
ESKAPE_genome_info$MRG_density <- ESKAPE_genome_info$MRG / ESKAPE_genome_info$Size *1000000
ESKAPE_genome_info$VFG_density <- ESKAPE_genome_info$VFG / ESKAPE_genome_info$Size *1000000
ESKAPE_genome_info$DS_density <- ESKAPE_genome_info$DS / ESKAPE_genome_info$Size *1000000

write.table(ESKAPE_genome_info,
            "ESKAPE_genome_info3.tsv",
            sep = "\t",           # TSV使用制表符分隔
            row.names = FALSE,    # 不写入行名
            col.names = TRUE,     # 写入列名
            quote = FALSE        # 不添加引号
          )            # NA值的表示方式

#ARG density#########################################################################
library(ggplot2)
mytheme <-   theme(
  panel.background = element_blank(),  # 设置背景为透明
  panel.grid.major = element_blank(),  # 去掉主要网格线
  panel.grid.minor = element_blank(),  # 去掉次要网格线
  panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
  axis.title = element_text( size = 12, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
  axis.text = element_text(size = 10, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
  axis.text.x = element_text(angle = 0, vjust = 0.8, hjust = 0.5),  # 横轴标签旋转
  axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
  legend.position = 'none',
  #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
  legend.title = element_text( size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
  legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
  legend.background = element_blank(),  # 设置图例背景为透明
  legend.box.background = element_rect(color = "grey80", linewidth = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
  legend.key = element_blank()  # 去掉图例键的背景
)

p_ARG_MRG <- 
ggplot(ESKAPE_genome_info, aes(x = ARG, y = MRG, color = Tax)) +
  geom_point(alpha = 0.6, size = 2) +
  #geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "dashed") +
  scale_color_manual(values = c("Ef"="#9f86c0","Sa"="#FFBE7D",
                               "Kp"="#76c893","Ab"="#F0E442","Pa"="#ff99c8",
                               "En"="#56B4E9")
  ) + mytheme
p_ARG_MRG

p_ARG_VFG <- 
  ggplot(ESKAPE_genome_info, aes(x = ARG, y = VFG, color = Tax)) +
  geom_point(alpha = 0.6, size = 2) +
  #geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "dashed") +
  scale_color_manual(values = c("Ef"="#9f86c0","Sa"="#FFBE7D",
                                "Kp"="#76c893","Ab"="#F0E442","Pa"="#ff99c8",
                                "En"="#56B4E9")
  ) + mytheme
p_ARG_VFG

p_ARG_DS <- 
  ggplot(ESKAPE_genome_info, aes(x = ARG, y = DS, color = Tax)) +
  geom_point(alpha = 0.6, size = 2) +
  #geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "dashed") +
  scale_color_manual(values = c("Ef"="#9f86c0","Sa"="#FFBE7D",
                                "Kp"="#76c893","Ab"="#F0E442","Pa"="#ff99c8",
                                "En"="#56B4E9")
  ) + mytheme
p_ARG_DS



# 计算Spearman相关系数
cor_spearman <- cor(ESKAPE_genome_info$ARG, 
                    ESKAPE_genome_info$VFG, 
                    use = "complete.obs", method = "spearman")

# 显著性检验
cor_test_spearman <- cor.test(ESKAPE_genome_info$ARG, 
                              ESKAPE_genome_info$VFG, 
                              method = "spearman")

cat("Spearman相关系数:", round(cor_spearman, 3), "\n")
cat("P值:", format(cor_test_spearman$p.value, scientific = TRUE, digits = 3), "\n")

library(cowplot)
g1 <- plot_grid(p_ARG_VFG,p_ARG_MRG,labels = LETTERS[2], nrow = 1)
g1

setwd(file.path(CS_DATA, "3genome_gene_annotation"))
ggsave(plot = g1,'1_ARG2VFGMRG.pdf',width = 5,height = 2.2,units = 'in')

p_ARG_VFG_density <- 
  ggplot(ESKAPE_genome_info[ESKAPE_genome_info$ARG_density<15,], 
         aes(x = ARG_density, y = VFG_density, color = Tax)) +
  geom_point(alpha = 0.6, size = 2) +
  #geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "dashed") +
  scale_color_manual(values = c("Ef"="#9f86c0","Sa"="#FFBE7D",
                                "Kp"="#76c893","Ab"="#F0E442","Pa"="#ff99c8",
                                "En"="#56B4E9")
  ) + mytheme
p_ARG_VFG_density

p_ARG_MRG_density <- 
  ggplot(ESKAPE_genome_info[ESKAPE_genome_info$ARG_density<15,], 
         aes(x = ARG_density, y = MRG_density, color = Tax)) +
  geom_point(alpha = 0.6, size = 2) +
  #geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "dashed") +
  scale_color_manual(values = c("Ef"="#9f86c0","Sa"="#FFBE7D",
                                "Kp"="#76c893","Ab"="#F0E442","Pa"="#ff99c8",
                                "En"="#56B4E9")
  ) + mytheme
p_ARG_MRG_density

p_ARG_DS_density <- 
  ggplot(ESKAPE_genome_info[ESKAPE_genome_info$ARG_density<15,], 
         aes(x = ARG_density, y = DS_density, color = Tax)) +
  geom_point(alpha = 0.6, size = 2) +
  #geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "dashed") +
  scale_color_manual(values = c("Ef"="#9f86c0","Sa"="#FFBE7D",
                                "Kp"="#76c893","Ab"="#F0E442","Pa"="#ff99c8",
                                "En"="#56B4E9")
  ) + mytheme
p_ARG_DS_density

library(cowplot)
g2 <- plot_grid(p_ARG_VFG_density,p_ARG_MRG_density,p_ARG_DS_density,
                labels = LETTERS[2], nrow = 1)
g2

setwd(file.path(CS_DATA, "3genome_gene_annotation"))
ggsave(plot = g2,'1_ARG2VFGMRG_density.pdf',width = 7.87,height = 2.5,units = 'in')

#DS frequency###########################################

# 创建计数矩阵
DS_count_matrix <- DS_info %>%
  count(genome, type) %>%  # 计算每个genome-type组合的出现次数
  complete(genome, type, fill = list(n = 0)) %>%  # 填充缺失组合为0
  pivot_wider(
    names_from = type,
    values_from = n,
    values_fill = 0
  ) %>%
  as.data.frame()

ESKAPE_DS_distribution <- ESKAPE_genome_info %>%
  select(genome, Tax)

ESKAPE_DS_distribution <- ESKAPE_DS_distribution %>%
  left_join(.,DS_count_matrix, suffix = c("genome","genome"))

ESKAPE_DS_distribution[is.na(ESKAPE_DS_distribution)] <- 0
ESKAPE_DS_distribution <- ESKAPE_DS_distribution[,-1]

#统计col2列每个元素的数量
element_counts <- ESKAPE_DS_distribution %>%
  count(Tax)

#分组统计
result <- function(df) {
  output <- df %>%
    group_by(Tax) %>%
    summarise(across(where(is.numeric), ~mean(. != 0, na.rm = TRUE) *100 ))
  return(output)
}
result_data <- result(ESKAPE_DS_distribution) %>% as.data.frame()

#过滤数据
library(dplyr)

rownames(result_data) <- result_data$Tax
result_data <- result_data[,-1]
#过滤列#
#删除全为0的列
result_data <- result_data[, !colSums(result_data) == 0]
#筛选主要的DSs

ds_filter <- colSums(result_data) %>% as.data.frame() %>% setNames("Propotion_sum")
ds_filter$DS <- rownames(ds_filter)
ds_filter <- ds_filter[ds_filter$Propotion_sum>=50,]

result_data <- result_data[,colnames(result_data) %in% ds_filter$DS]

#绘制热图
library(pheatmap)

# 设置颜色
#colors <- colorRampPalette(c('white', colorRampPalette('#F9B14F', '#8F2731')))(1000)
colors <- c('white', colorRampPalette(c('#F1EDE4', '#3C4D61'))(9999))
p <- pheatmap(result_data,
              #color = colorRampPalette(colors)(1000),
              #color = colorRampPalette(c('white','lightpink','red'))(1000),
              color = colors,
              cluster_cols = F,
              cluster_rows = F,
              angle_col = 90
)
p

library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)

# 转换为长格式
result_long <- as.data.frame(result_data) %>%
  rownames_to_column("Species") %>%
  pivot_longer(cols = -Species, names_to = "Defense_System", values_to = "Percentage")

result_long$Species <- factor(result_long$Species, levels = c("Ef","Sa","Kp","Ab","Pa","En"))

# 创建热图
p_DS_frequency <- 
ggplot(result_long, aes(x = Species, y = Defense_System, 
                        fill = Species, alpha = Percentage)) +
  geom_tile(color = "white", size = 0) +
  # 添加数值标签
  geom_text(aes(label = sprintf("%.1f", Percentage)), 
            color = "black", size = 2.5, fontface = "bold") +
  scale_fill_manual(values = c("Ef"="#9f86c0","Sa"="#FFBE7D",
                               "Kp"="#76c893","Ab"="#F0E442","Pa"="#ff99c8",
                               "En"="#56B4E9")
  )+
  scale_alpha_continuous(range = c(0.2, 1),  # 透明度范围
                         guide = "none") +   # 隐藏alpha图例
  scale_x_discrete(position = "top") +
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_blank(),
    #axis.title = element_text( size = 12, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 9, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5),  # 横轴标签旋转
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_blank(), #设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_blank(), #设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", linewidth = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
p_DS_frequency


#ARG frequency try1############################################################################

ARG_detail_info <- read.delim2("rgiargs_detail_info.txt")
ARG_detail_info <- ARG_detail_info %>%
  mutate(
    genome = gsub("\\-.*","",ORF_ID)
  )

antibiotics_type <- table(ARG_detail_info$Antibiotic) %>% as.data.frame()
drug_class <- table(ARG_detail_info$Drug.Class) %>% as.data.frame()
gene_family <- table(ARG_detail_info$AMR.Gene.Family) %>% as.data.frame()
gene_name <- table(ARG_detail_info$Best_Hit_ARO) %>% as.data.frame()

# 创建计数矩阵
ARG_count_matrix <- ARG_detail_info %>%
  count(genome, Best_Hit_ARO) %>%  # 计算每个genome-type组合的出现次数
  complete(genome, Best_Hit_ARO, fill = list(n = 0)) %>%  # 填充缺失组合为0
  pivot_wider(
    names_from = Best_Hit_ARO,
    values_from = n,
    values_fill = 0
  ) %>%
  as.data.frame()

ESKAPE_ARG_distribution <- ESKAPE_genome_info %>%
  select(genome, Tax)

ESKAPE_ARG_distribution <- ESKAPE_ARG_distribution %>%
  left_join(.,ARG_count_matrix, suffix = c("genome","genome"))

ESKAPE_ARG_distribution[is.na(ESKAPE_ARG_distribution)] <- 0
ESKAPE_ARG_distribution <- ESKAPE_ARG_distribution[,-1]

#统计col2列每个元素的数量
element_counts <- ESKAPE_ARG_distribution %>%
  count(Tax)

#分组统计
result <- function(df) {
  output <- df %>%
    group_by(Tax) %>%
    summarise(across(where(is.numeric), ~mean(. != 0, na.rm = TRUE) * 100))
  return(output)
}
ARG_result_data <- result(ESKAPE_ARG_distribution) %>% as.data.frame()

#过滤数据
library(dplyr)

rownames(ARG_result_data) <- ARG_result_data$Tax
ARG_result_data <- ARG_result_data[,-1]
#过滤列#
#删除全为0的列
ARG_result_data <- ARG_result_data[, !colSums(ARG_result_data) == 0]
#筛选主要的DSs

ARG_filter <- colSums(ARG_result_data) %>% 
  as.data.frame() %>% setNames("Propotion_sum")
ARG_filter$ARG <- rownames(ARG_filter)
ARG_filter <- ARG_filter[ARG_filter$Propotion_sum>=100,]

ARG_result_data <- ARG_result_data[,colnames(ARG_result_data) %in% 
                                     ARG_filter$ARG]


#绘制热图
library(pheatmap)

# 设置颜色
#colors <- colorRampPalette(c('white', colorRampPalette('#F9B14F', '#8F2731')))(1000)
colors <- c('white', colorRampPalette(c('#F1EDE4', '#3C4D61'))(9999))
p <- pheatmap(ARG_result_data,
              #color = colorRampPalette(colors)(1000),
              #color = colorRampPalette(c('white','lightpink','red'))(1000),
              color = colors,
              cluster_cols = F,
              cluster_rows = F,
              angle_col = 90
)
p

library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)

# 转换为长格式
ARG_result_long <- as.data.frame(ARG_result_data) %>%
  rownames_to_column("Species") %>%
  pivot_longer(cols = -Species, names_to = "ARG", values_to = "Percentage")

ARG_result_long$Species <- factor(ARG_result_long$Species, levels = c("Ef","Sa","Kp","Ab","Pa","En"))

# 创建热图
p_ARG_frequency <- 
  ggplot(ARG_result_long, aes(x = ARG, y = Species, 
                          fill = Species, alpha = Percentage)) +
  geom_tile(color = "white", size = 0) +
  # 添加数值标签
  geom_text(aes(label = sprintf("%.1f", Percentage)), 
            color = "black", size = 3, fontface = "bold") +
  scale_fill_manual(values = c("Ef"="#9f86c0","Sa"="#FFBE7D",
                               "Kp"="#76c893","Ab"="#F0E442","Pa"="#ff99c8",
                               "En"="#56B4E9")
  )+
  scale_alpha_continuous(range = c(0.1, 1),  # 透明度范围
                         guide = "none") +   # 隐藏alpha图例
  scale_x_discrete(position = "top") +
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_blank(),
    #axis.title = element_text( size = 12, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 30, vjust = 0, hjust = 0),  # 横轴标签旋转
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_blank(), #设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_blank(), #设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", linewidth = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
p_ARG_frequency

#ARG frequency try2############################################################################

ARG_detail_info <- read.delim2("rgiargs_detail_info.txt")
ARG_detail_info <- ARG_detail_info %>%
  mutate(
    genome = gsub("\\-.*","",ORF_ID)
  )

antibiotics_type <- table(ARG_detail_info$Antibiotic) %>% as.data.frame()
drug_class <- table(ARG_detail_info$Drug.Class) %>% as.data.frame()
gene_family <- table(ARG_detail_info$AMR.Gene.Family) %>% as.data.frame()
gene_name <- table(ARG_detail_info$Best_Hit_ARO) %>% as.data.frame()

# 拆分Drug.Class列
ARG_detail_expanded <- ARG_detail_info[ARG_detail_info$Drug.Class!="",]
ARG_detail_expanded <- ARG_detail_expanded %>%
  separate_rows(Drug.Class, sep = ";")

# 去除空白字符
ARG_detail_expanded <- ARG_detail_expanded %>%
  mutate(Drug.Class = trimws(Drug.Class))  # 去除首尾空格

# 创建计数矩阵
ARG_count_matrix <- ARG_detail_expanded %>%
  count(genome, Drug.Class) %>%  # 计算每个genome-type组合的出现次数
  complete(genome, Drug.Class, fill = list(n = 0)) %>%  # 填充缺失组合为0
  pivot_wider(
    names_from = Drug.Class,
    values_from = n,
    values_fill = 0
  ) %>%
  as.data.frame()

ESKAPE_ARG_distribution <- ESKAPE_genome_info %>%
  select(genome, Tax)

ESKAPE_ARG_distribution <- ESKAPE_ARG_distribution %>%
  left_join(.,ARG_count_matrix, suffix = c("genome","genome"))

ESKAPE_ARG_distribution[is.na(ESKAPE_ARG_distribution)] <- 0
ESKAPE_ARG_distribution <- ESKAPE_ARG_distribution[,-1]

#统计col2列每个元素的数量
element_counts <- ESKAPE_ARG_distribution %>%
  count(Tax)

#分组统计
result <- function(df) {
  output <- df %>%
    group_by(Tax) %>%
    summarise(across(where(is.numeric), ~mean(. !=0,na.rm = TRUE) *100 ))
  return(output)
}
ARG_result_data <- result(ESKAPE_ARG_distribution) %>% as.data.frame()

#过滤数据
library(dplyr)

rownames(ARG_result_data) <- ARG_result_data$Tax
ARG_result_data <- ARG_result_data[,-1]
#过滤列#
#删除全为0的列
ARG_result_data <- ARG_result_data[, !colSums(ARG_result_data) == 0]
#筛选主要的DSs

ARG_filter <- colSums(ARG_result_data) %>% 
  as.data.frame() %>% setNames("Propotion_sum")
ARG_filter$ARG <- rownames(ARG_filter)
ARG_filter <- ARG_filter[ARG_filter$Propotion_sum>=0,]

ARG_result_data <- ARG_result_data[,colnames(ARG_result_data) %in% 
                                     ARG_filter$ARG]


#绘制热图
library(pheatmap)

# 设置颜色
#colors <- colorRampPalette(c('white', colorRampPalette('#F9B14F', '#8F2731')))(1000)
colors <- c('white', colorRampPalette(c('#F1EDE4', '#3C4D61'))(9999))
p <- pheatmap(ARG_result_data,
              #color = colorRampPalette(colors)(1000),
              #color = colorRampPalette(c('white','lightpink','red'))(1000),
              color = colors,
              cluster_cols = F,
              cluster_rows = F,
              angle_col = 90
)
p

library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)

# 转换为长格式
ARG_result_long <- as.data.frame(ARG_result_data) %>%
  rownames_to_column("Species") %>%
  pivot_longer(cols = -Species, names_to = "ARG", values_to = "Percentage")

ARG_result_long$Species <- factor(ARG_result_long$Species, levels = c("Ef","Sa","Kp","Ab","Pa","En"))

ARG_result_long$ARG <- gsub(" antibiotic","", ARG_result_long$ARG)
ARG_result_long$ARG <- gsub("disinfecting agents and antiseptics",
                            "DAA", ARG_result_long$ARG)
ARG_result_long$ARG <- gsub("penicillin",
                            "P", ARG_result_long$ARG)

# 创建热图
p_ARG_frequency <- 
  ggplot(ARG_result_long, aes(x = Species, y = ARG, 
                              fill = Species, alpha = Percentage)) +
  geom_tile(color = "white", size = 0) +
  # 添加数值标签
  geom_text(aes(label = sprintf("%.1f", Percentage)), 
            color = "black", size = 2.5, fontface = "bold") +
  scale_fill_manual(values = c("Ef"="#9f86c0","Sa"="#FFBE7D",
                               "Kp"="#76c893","Ab"="#F0E442","Pa"="#ff99c8",
                               "En"="#56B4E9")
  )+
  scale_alpha_continuous(range = c(0.2, 1),  # 透明度范围
                         guide = "none") +   # 隐藏alpha图例
  scale_x_discrete(position = "top") +
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_blank(),
    #axis.title = element_text( size = 12, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 9, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),  # 横轴标签旋转
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_blank(), #设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_blank(), #设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", linewidth = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
p_ARG_frequency

#ARG frequency try3############################################################################

ARG_detail_info <- read.delim2("rgiargs_detail_info.txt")
ARG_detail_info <- ARG_detail_info %>%
  mutate(
    genome = gsub("\\-.*","",ORF_ID)
  )

antibiotics_type <- table(ARG_detail_info$Antibiotic) %>% as.data.frame()
drug_class <- table(ARG_detail_info$Drug.Class) %>% as.data.frame()
gene_family <- table(ARG_detail_info$AMR.Gene.Family) %>% as.data.frame()
gene_name <- table(ARG_detail_info$Best_Hit_ARO) %>% as.data.frame()

# 拆分antibiotics_type列
ARG_detail_expanded <- ARG_detail_info[ARG_detail_info$Antibiotic!="",]
ARG_detail_expanded <- ARG_detail_expanded %>%
  separate_rows(Antibiotic, sep = ";")

# 去除空白字符
ARG_detail_expanded <- ARG_detail_expanded %>%
  mutate(Antibiotic = trimws(Antibiotic))  # 去除首尾空格

# 创建计数矩阵
ARG_count_matrix <- ARG_detail_expanded %>%
  count(genome, Antibiotic) %>%  # 计算每个genome-type组合的出现次数
  complete(genome, Antibiotic, fill = list(n = 0)) %>%  # 填充缺失组合为0
  pivot_wider(
    names_from = Antibiotic,
    values_from = n,
    values_fill = 0
  ) %>%
  as.data.frame()

ESKAPE_ARG_distribution <- ESKAPE_genome_info %>%
  select(genome, Tax)

ESKAPE_ARG_distribution <- ESKAPE_ARG_distribution %>%
  left_join(.,ARG_count_matrix, suffix = c("genome","genome"))

ESKAPE_ARG_distribution[is.na(ESKAPE_ARG_distribution)] <- 0
ESKAPE_ARG_distribution <- ESKAPE_ARG_distribution[,-1]

#统计col2列每个元素的数量
element_counts <- ESKAPE_ARG_distribution %>%
  count(Tax)

#分组统计
result <- function(df) {
  output <- df %>%
    group_by(Tax) %>%
    summarise(across(where(is.numeric), ~mean(., na.rm = TRUE) ))
  return(output)
}
ARG_result_data <- result(ESKAPE_ARG_distribution) %>% as.data.frame()

#过滤数据
library(dplyr)

rownames(ARG_result_data) <- ARG_result_data$Tax
ARG_result_data <- ARG_result_data[,-1]
#过滤列#
#删除全为0的列
ARG_result_data <- ARG_result_data[, !colSums(ARG_result_data) == 0]
#筛选主要的DSs

ARG_filter <- colSums(ARG_result_data) %>% 
  as.data.frame() %>% setNames("Propotion_sum")
ARG_filter$ARG <- rownames(ARG_filter)
ARG_filter <- ARG_filter[ARG_filter$Propotion_sum>=10,]

ARG_result_data <- ARG_result_data[,colnames(ARG_result_data) %in% 
                                     ARG_filter$ARG]


#绘制热图
library(pheatmap)

# 设置颜色
#colors <- colorRampPalette(c('white', colorRampPalette('#F9B14F', '#8F2731')))(1000)
colors <- c('white', colorRampPalette(c('#F1EDE4', '#3C4D61'))(9999))
p <- pheatmap(ARG_result_data,
              #color = colorRampPalette(colors)(1000),
              #color = colorRampPalette(c('white','lightpink','red'))(1000),
              color = colors,
              cluster_cols = F,
              cluster_rows = F,
              angle_col = 90
)
p

library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)

# 转换为长格式
ARG_result_long <- as.data.frame(ARG_result_data) %>%
  rownames_to_column("Species") %>%
  pivot_longer(cols = -Species, names_to = "ARG", values_to = "Percentage")

ARG_result_long$Species <- factor(ARG_result_long$Species, levels = c("Ef","Sa","Kp","Ab","Pa","En"))

# 创建热图
p_ARG_frequency <- 
  ggplot(ARG_result_long, aes(x = Species, y = ARG, 
                              fill = Species, alpha = Percentage)) +
  geom_tile(color = "white", size = 0) +
  # 添加数值标签
  geom_text(aes(label = sprintf("%.1f", Percentage)), 
            color = "black", size = 2.5, fontface = "bold") +
  scale_fill_manual(values = c("Ef"="#9f86c0","Sa"="#FFBE7D",
                               "Kp"="#76c893","Ab"="#F0E442","Pa"="#ff99c8",
                               "En"="#56B4E9")
  )+
  scale_alpha_continuous(range = c(0.2, 1),  # 透明度范围
                         guide = "none") +   # 隐藏alpha图例
  scale_x_discrete(position = "top") +
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_blank(),
    #axis.title = element_text( size = 12, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),  # 横轴标签旋转
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_blank(), #设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_blank(), #设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", linewidth = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
p_ARG_frequency

library(cowplot)
g3 <- plot_grid(p_ARG_frequency,p_DS_frequency,
                #labels = LETTERS[2], 
                rel_widths = c(2, 1.7),
                nrow = 1)
g3

setwd(file.path(CS_DATA, "3genome_gene_annotation"))
ggsave(plot = g3,'1_ARG_DS_heatmap3.pdf',width = 5.1,height = 3.75,units = 'in')

save.image("./1_ARG.RData")

load("./1_ARG.RData")


