# ===========================================================================
# ESKAPE_mobilome_virome -- analyses / 02_phage_catalogue
# Genome features of complete phages (AAU, pI, length)
# Manuscript      : Fig. 3
# Original path   : 2ESKAPE_MGE/2virus_genome_info/2_phage_genome_features.R
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


#比较不同ESKAPE病毒的GC length AAU pI

setwd(file.path(CS_DATA, "2virus_genome_info"))

library(tidyr)
library(dplyr)
library(stringr)

load("1_virus_genome_info.RData")

#AAU####
complete_phage <- 
  ESKAPE_all_phage_info_tax[ESKAPE_all_phage_info_tax$checkv_quality %in% 
                              c("Complete"),]

table(complete_phage$Host_tax)


AAU_matrix_phage <- read.csv("AAU_matrix.csv")
complete_phage_host <- select(complete_phage,ID,Host_tax)
complete_phage_host$ID <- gsub("\\|", "__", complete_phage_host$ID)
complete_phage_host_AAU <- 
  left_join(complete_phage_host, AAU_matrix_phage,
                                     by = c("ID"="Genome"))
complete_phage_host_AAU <- select(complete_phage_host_AAU, -ID)

# 按 Host_tax 分组，并对所有其他数值列求平均值
AAU_result <- complete_phage_host_AAU %>%
  group_by(Host_tax) %>%
  # 修改点：使用 ~ mean(.x, na.rm = TRUE)
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = 'drop') %>%
  as.data.frame()


library(pheatmap)
plot_data <- AAU_result
rownames(plot_data) <- plot_data$Host_tax  # 设置行名
plot_data$Host_tax <- NULL                 # 删除原始的 Host_tax 列，只留数值
plot_data <- plot_data*100

# --- 3. 绘制美化热图 ---
p_AAU_heatmap <- 
pheatmap(plot_data,
         #color = color_palette,          # 使用自定义配色
         cluster_rows = TRUE,            # 行聚类
         cluster_cols = FALSE,           # 列不聚类 (保持氨基酸顺序 A-Y)
         cutree_rows = 4,                # 切割为 4 类
         
         # --- 边框与网格 ---
         border_color = "grey60",        # 边框颜色稍深一点，更清晰
         cellwidth = 15,                 # 增加单元格宽度
         cellheight = 12,                # 增加单元格高度
         
         # --- 字体与标签 ---
         fontsize = 10,                  # 基础字体大小
         fontsize_row = 9,               # 行名字体 (宿主)
         fontsize_col = 8,               # 列名字体 (氨基酸)
         angle_col = 0,                 # 列名旋转 45 度，防止重叠
         display_numbers = FALSE,        # 是否显示数值 (格子小则不显示，设为 TRUE 可显示)
         number_format = "%.2f",         # 数值格式 (如果显示的话)
         number_color = "grey30",        # 数值颜色
         
         # --- 聚类树美化 ---
         treeheight_row = 50,            # 行聚类树高度
         treeheight_col = 0,             # 列聚类树高度 (因为没聚类，设为0)
         clustering_distance_rows = "euclidean", # 距离算法: euclidean, correlation 等
         clustering_method = "complete",       # 聚类方法: complete, average, ward.D2
         
         # --- 标题与图例 ---
         main = "Amino Acid Usage Frequency by complete ESKAPE phages", # 主标题
         legend = TRUE,
         legend_breaks = c(0, 5, 10), # 图例刻度 (可选)
         legend_labels = c("0", "5%", "10%"), # 图例标签 (可选)
         
         # --- 其他 ---
         na_col = "grey90",              # 如果有 NA，显示的颜色
         silent = FALSE                  # 是否静默运行
)


#pI####

pI_distribution_phage <- read.csv("pI_distribution_matrix.csv")

complete_phage_host_pI <- 
  left_join(complete_phage_host, pI_distribution_phage,
            by = c("ID"="Genome"))
complete_phage_host_pI <- select(complete_phage_host_pI, -ID)

# 按 Host_tax 分组，并对所有其他数值列求平均值
pI_result <- complete_phage_host_pI %>%
  group_by(Host_tax) %>%
  # 修改点：使用 ~ mean(.x, na.rm = TRUE)
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = 'drop') %>%
  as.data.frame()

library(ggplot2)
library(tidyr)
library(dplyr)
library(tibble)

# 1. 数据预处理
# 假设你的数据框名为 pI_result
pI_result2 <- pI_result

# 步骤 1: 保存 Host_tax 值
host_tax_values <- pI_result2$Host_tax

# 步骤 2: 只保留数值列进行转置
numeric_pI_result2 <- pI_result2[, !colnames(pI_result2) %in% "Host_tax"] 

# 步骤 3: 转置，并将原始 Host_tax 值作为列名
transposed_pI_result2 <- t(numeric_pI_result2)
colnames(transposed_pI_result2) <- host_tax_values # 这里是关键！

# 步骤 4: 转换为数据框
transposed_pI_result2 <- as.data.frame(transposed_pI_result2)

# 步骤 5: 将行名（pI 值）转换为一列
transposed_pI_result2$pI_Value <- 
  as.numeric(gsub("pI_", "", rownames(transposed_pI_result2)))

# 步骤 6: 长格式转换
pI_plot_data <- transposed_pI_result2 %>%
  select(pI_Value, everything()) %>% # 把 pI_Value 放在第一列
  pivot_longer(
    cols = -pI_Value,             # 除了 pI_Value 外的所有列
    names_to = "Host_tax",        # 新列名：宿主分类（来自原来的列名，即 Host_tax 值）
    values_to = "Frequency"       # 新列名：频率值
  ) %>%
  mutate(
    Host_tax = factor(Host_tax, levels = unique(host_tax_values)), # 保持原有顺序
    Frequency_Percent = Frequency * 100                            # 转换为百分比
  )

# 2. 绘制折线图
p_pI_distribution <- ggplot(pI_plot_data, 
            aes(x = pI_Value, 
                y = Frequency_Percent, 
                color = Host_tax, group = Host_tax)) +
  geom_line(linewidth = 0.5, alpha = 1) +          # 绘制折线，线宽 1，透明度 0.8
  geom_point(size = 0.5, alpha = 0.6) +            # 添加数据点 (可选，类似 Python 默认标记)
  
  # 标签与标题
  labs(
  #  title = "Relative Frequency Distribution of pI Values",
  #  subtitle = "Grouped by Host Taxonomy",
    x = "pI",              # X 轴标签
    y = "Relative frequency (%)",                  # Y 轴标签
  #  color = "Host Taxonomy"                        # 图例标题
  ) +
  scale_color_manual(values = c("Ef"="#9f86c0","Sa"="#FFBE7D",
                               "Kp"="#76c893","Ab"="#F0E442","Pa"="#ff99c8",
                               "En"="#56B4E9", "En,Kp"="black")
  ) +
  # 主题美化 (类似 matplotlib 默认风格但更现代)
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, linewidth = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_text( size = 12, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0, vjust = 0.8, hjust = 0.5),  # 横轴标签旋转
    axis.ticks = element_line(color = "grey20", linewidth = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_text( size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", linewidth = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )+
  # 坐标轴调整 (根据数据范围自动调整，也可手动设定)
  scale_x_continuous(breaks = seq(0, 14, by = 1)) + # X 轴刻度每 1 单位一个
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) # Y 轴顶部留白 5%

# 3. 显示图形
print(p_pI_distribution)

#GC含量与宿主比较####

phage_characteristics <- read.delim2("Molecular_Final_Full_Stats.txt")
phage_characteristics$Sample <- gsub(".fa","",phage_characteristics$Sample)

phage_characteristics <- phage_characteristics %>%
  left_join(complete_phage_host, by = c("Sample"="ID"))

# 自动检测并转换所有合适的列为数值类型
phage_characteristics <- type.convert(phage_characteristics, as.is = TRUE)

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

# 1. 整理数据：将需要比较的指标转为长格式
# 假设你的数据框名为 merged_phage_compare
plot_data_phage_charac <- phage_characteristics %>%
  select(Host_tax,Size_bp, GC, AvgMW, CARSC, NARSC, SARSC) %>%
  pivot_longer(cols = -Host_tax, 
               names_to = "Metric", 
               values_to = "Value")

plot_data_phage_charac$Host_tax <- factor(plot_data_phage_charac$Host_tax, 
                                          levels = c("Ef","Sa","Kp","Ab","Pa","En","En,Kp"))
# 2. 绘图
# 使用 facet_wrap 分面显示不同指标，scales = "free_y" 让每个图有自己的 Y 轴范围
p_phage_characteristics_comparison <- 
  ggplot(plot_data_phage_charac, 
         aes(x = Host_tax, y = Value, color = Host_tax)) +
  #  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
  geom_boxplot(alpha = 0.7, 
               outlier.shape = 1,          # 实心圆点
               outlier.size = 1,            # 离群点大小
               outlier.alpha = 0.7,         # 离群点透明度
               width = 0.6,
               size = 0.4                   # 边框线粗细
  ) +
  stat_summary(
    fun = mean,  # 计算均值
    aes(fill = Host_tax),
    geom = "point",  # 绘制点
    shape = 21,  # 菱形形状
    size = 1.5,    # 点的大小
    color = "black"    # 边框颜色
  )+
  #  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  facet_wrap(~Metric, scales = "free_y", ncol = 3) +
  theme_bw() +
  # 自动添加统计检验（例如 Kruskal-Wallis）
  #  stat_compare_means(label = "p.signif", label.x = 1.5) +
  # labs(title = "Molecular Characteristics Comparison across Ecosystems",
  #      x = "Ecosystem Type",
  #      y = "Measured Value") +
  #scale_color_brewer(palette = "Set1") +
  #scale_fill_brewer(palette = "Set1") +
  scale_color_manual(values = c("Ef"="#9f86c0","Sa"="#FFBE7D",
                                "Kp"="#76c893","Ab"="#F0E442","Pa"="#ff99c8",
                                "En"="#56B4E9", "En,Kp"="black")
  ) +
  scale_fill_manual(values = c("Ef"="#9f86c0","Sa"="#FFBE7D",
                                "Kp"="#76c893","Ab"="#F0E442","Pa"="#ff99c8",
                                "En"="#56B4E9", "En,Kp"="black")
  ) +
  # theme(
  #   strip.background = element_rect(fill = "gray90"),
  #   strip.text = element_text(face = "bold"),
  #   axis.text.x = element_text(angle = 45, hjust = 1),
  #   legend.position = "none"
  # )+
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "plain"),  # 标题居中，加粗，字体大小为16
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),  # 设置图形边距
    plot.background = element_blank(),  # 去除图形背景
    panel.background = element_blank(),  # 去除面板背景
    panel.border = element_rect(linewidth = 0.6),  # 设置面板边框
    panel.grid = element_blank(),  # 去除网格线
    #axis.title.x = element_text(size = 10, face = "plain",color = 'black'),  # 设置 x 轴标题
    #axis.title.y = element_text(size = 10, face = "plain",color = 'black'),  # 设置 y 轴标题
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    #axis.ticks.x = element_line(linewidth = 0.5),  # 设置 x 轴刻度线
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(linewidth = 0.5),  # 设置 y 轴刻度线
    axis.ticks.length.x = unit(0.1, "cm"),  # 设置 x 轴刻度长度
    axis.ticks.length.y = unit(0.1, "cm"),  # 设置 y 轴刻度长度
    #axis.text.x = element_text(angle = 30, size = 10, hjust = 1, vjust = 1,color = 'black'),  # 设置 x 轴文本
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 10, color = 'black'),  # 设置 y 轴文本
    legend.position = "none"  # 去除图例
  )
p_phage_characteristics_comparison



ggsave(plot = p_phage_characteristics_comparison,
       '2_p_phage_characteristics_comparison.pdf',
       width = 9.4,height = 4,units = 'in')
ggsave(plot = p_AAU_heatmap,
       '2_p_AAU_heatmap.pdf',
       width = 5.9,height = 2,units = 'in')
ggsave(plot = p_pI_distribution,
       '2_p_pI_distribution.pdf',
       width = 3.5,height = 2,units = 'in')



#图片组合
library(Cairo)
library(cowplot)
# g_complete_phage_features <- plot_grid(p_AAU_heatmap,
#                                       p_pI_distribution,
#                                       p_phage_characteristics_comparison,
#                                       ncol = 2, nrow = 2
#                                       # labels = LETTERS[1:2],
#                                       # rel_widths = c(1, 2)
#                                       # label_x = 0.09,            # 标签x坐标（1为最右侧）
#                                       #  label_y = 0.9             # 标签y坐标（1为最上方）
#                           )
# g_complete_phage_features

ggsave(plot = g_MGE_spacer_competition,'1_g_MGE_spacer_competition2.pdf',
       width = 2.7,height = 2.7,units = 'in')

save.image("./2_phage_genome_features.RData")






