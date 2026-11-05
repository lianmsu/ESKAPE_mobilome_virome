# ===========================================================================
# ESKAPE_mobilome_virome -- analyses / 06_phage_defense
# Defense and anti-defense genes of the phages
# Manuscript      : Fig. 3F
# Original path   : 2ESKAPE_MGE/6_0phage_DF/1_phage_DF.R
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

setwd(file.path(CS_DATA, "6_0phage_DF"))

library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

ESKAPE_all_phage_info <- read.delim2("ESKAPE_all_phage_info_filter.tsv")

ESKAPE_all_phage_DF <- read.delim2("ESKAPE_all_phage_info_filter_defense_finder_systems.tsv")
ESKAPE_all_phage_DF <- distinct(ESKAPE_all_phage_DF)
table(ESKAPE_all_phage_DF$activity)

#不同MGE上的Anti-defense的比例####
load(file.path(CS_DATA, "4genome_MGE/1_ESKAPE_MGE.RData"))

MGE_ds_num <- all_mge_ds %>%
  count(MGE_type, name = "count")
MGE_ds_num$type <- "DS"

MGE_antids_num <- all_mge_antids %>%
  count(MGE_type, name = "count")
MGE_antids_num$type <- "Anti-DS"
#MGE_antids_num <- MGE_antids_num[-2,]

noprovirus_phage_DF <- 
  ESKAPE_all_phage_DF[!grepl("^GCA", ESKAPE_all_phage_DF$sys_id), ]
noprovirus_phage_DF <- distinct(noprovirus_phage_DF)
table(noprovirus_phage_DF$activity)


noprovirus_phage_DF_num <- 
  data.frame(MGE_type = "Phage", count = c(782,1428), type = c("DS","Anti-DS"))

all_MGE_DF_num <- rbind(noprovirus_phage_DF_num, MGE_ds_num, MGE_antids_num)


# 准备数据：计算每个MGE_type中DS和Anti-DS的比例
all_MGE_DF_prop <- all_MGE_DF_num %>%
  group_by(MGE_type) %>%
  mutate(
    total = sum(count),  # 每个MGE_type的总数
    prop = count / total,  # 计算比例
    percent = round(prop * 100, 1)  # 转换为百分比
  ) %>%
  ungroup()

# 查看处理后的数据
print(all_MGE_DF_prop)

all_MGE_DF_prop$MGE_type <- 
  factor(all_MGE_DF_prop$MGE_type, levels = c("Phage","Prophage","Plasmid","ICEs","Integron","Composite_transposon"))

# 绘制横向比例柱状图
p_MGE_DS_propotion <-
ggplot(all_MGE_DF_prop, 
       aes(x = prop, y = reorder(MGE_type, percent), fill = type)) +
  geom_col(position = "fill", width = 0.7, alpha = 0.8) +
  geom_text(
    aes(label = paste0( percent, "%")),
    position = position_fill(vjust = 0.5),
    size = 3.5,
    color = "white",
    fontface = "bold"
  ) +
  scale_x_continuous(
    labels = scales::percent_format(),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_fill_manual(
    values = c("DS" = "#1f77b4", "Anti-DS" = "#ff7f0e"),
    name = "Type"
  ) +
  theme(
    panel.spacing = unit(0.1, "lines"),  # 将分面间距设为0
    panel.background = element_blank(),  # 设置背景为透明
    # panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, linewidth = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_blank(),
    #axis.title = element_text( size = 12, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),  # 横轴标签旋转
    axis.ticks = element_line(color = "grey20", linewidth = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_blank(), #设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_blank(), #设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", linewidth = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
p_MGE_DS_propotion

#不同ESKAPE的phage上携带的AntiDS####
noprovirus_phage_DF <- noprovirus_phage_DF %>%
  mutate(
    contig = str_replace(sys_beg, "_\\d+$", "")
  )
noprovirus_phage_antiDS <- 
  noprovirus_phage_DF[noprovirus_phage_DF$activity == "Antidefense",]
noprovirus_phage_antiDS <- 
  noprovirus_phage_antiDS %>%
  left_join(select(ESKAPE_all_phage_info,seq_name,Host_tax), 
            by = c("contig"="seq_name"))
# 使用separate_rows拆分逗号分隔的值
noprovirus_phage_antiDS_expanded <- noprovirus_phage_antiDS %>%
  separate_rows(Host_tax, sep = ",") %>%
  mutate(Host_tax = trimws(Host_tax))  # 去除可能存在的空格
noprovirus_phage_antiDS_expanded <- distinct(noprovirus_phage_antiDS_expanded)
noprovirus_phage_antiDS_expanded_num <- 
  noprovirus_phage_antiDS_expanded %>%
  count(type, Host_tax, name = "count")
noprovirus_phage_antiDS_expanded_num$Host_tax <- 
  factor(noprovirus_phage_antiDS_expanded_num$Host_tax, levels = c("Ef","Sa","Kp","Ab","Pa","En"))

noprovirus_phage_info <- 
  ESKAPE_all_phage_info[!grepl("^GCA", ESKAPE_all_phage_info$seq_name), ]
table(noprovirus_phage_info$Host_tax)
noprovirus_phage_info_expanded <- noprovirus_phage_info %>%
  separate_rows(Host_tax, sep = ",") %>%
  mutate(Host_tax = trimws(Host_tax))  # 去除可能存在的空格
noprovirus_phage_info_expanded_num <- 
  table(noprovirus_phage_info_expanded$Host_tax) %>%
  as.data.frame()
colnames(noprovirus_phage_info_expanded_num) <- c("Host_tax", "Phage_num")

# 添加归一化列
noprovirus_phage_antiDS_normalized <- noprovirus_phage_antiDS_expanded_num %>%
  left_join(noprovirus_phage_info_expanded_num, by = "Host_tax") %>%
  mutate(
    count_normalized = count / Phage_num,
    percent = (count / Phage_num) * 100  # 可选：转换为百分比
  ) %>%
  arrange(type, Host_tax)

library(scales)
p_dot_heatmap_noprovirus_phage_antiDS <-
  ggplot(noprovirus_phage_antiDS_normalized, aes(x = Host_tax, y = type)) +
  geom_point(aes(size = count_normalized, color = Host_tax), alpha = 0.7) +
  scale_size_continuous(
    name = "count_normalized",
    range = c(1, 10),  # 控制点的大小范围
    breaks = pretty_breaks(n = 5)
  ) +
  scale_color_manual(values = c("Ef"="#9f86c0","Sa"="#FFBE7D",
                                "Kp"="#76c893","Ab"="#F0E442","Pa"="#ff99c8",
                                "En"="#56B4E9")
  )+
#  facet_grid(~ MGE_type, scales = "free", space = "free") +
  #facet_wrap(~ MGE_type, scales = "free", ncol = 5) +
  theme(
    panel.spacing = unit(0.1, "lines"),  # 将分面间距设为0
    panel.background = element_blank(),  # 设置背景为透明
    # panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
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
p_dot_heatmap_noprovirus_phage_antiDS

ggsave(plot = p_dot_heatmap_noprovirus_phage_antiDS,
       '1_heatmap_noprovirus_phage_antiDS_legend.pdf',
       width = 5,height = 4.66,units = 'in')

library(Cairo)
library(cowplot)
g_noprovirus_phage_antiDS <- 
  plot_grid(
    p_MGE_DS_propotion,
    p_dot_heatmap_noprovirus_phage_antiDS,
    labels = LETTERS[3:4], 
    ncol = 1,
    align = "v",  # 垂直对齐
    axis = "lr",  # 左和右对齐
    rel_widths = c(1, 1),
    rel_heights = c(1, 1.2)  # 根据实际高度调整
  )

g_noprovirus_phage_antiDS
ggsave(plot = g1,'1_MGEnum_antiDS.pdf',width = 8,height = 2.5,units = 'in')

save.image("./1_phage_DF.RData")

#不同ESKAPE的phage上携带的DS####
noprovirus_phage_DF <- noprovirus_phage_DF %>%
  mutate(
    contig = str_replace(sys_beg, "_\\d+$", "")
  )
noprovirus_phage_DS <- 
  noprovirus_phage_DF[noprovirus_phage_DF$activity == "Defense",]
noprovirus_phage_DS <- 
  noprovirus_phage_DS %>%
  left_join(select(ESKAPE_all_phage_info,seq_name,Host_tax), 
            by = c("contig"="seq_name"))
# 使用separate_rows拆分逗号分隔的值
noprovirus_phage_DS_expanded <- noprovirus_phage_DS %>%
  separate_rows(Host_tax, sep = ",") %>%
  mutate(Host_tax = trimws(Host_tax))  # 去除可能存在的空格
noprovirus_phage_DS_expanded <- distinct(noprovirus_phage_DS_expanded)
noprovirus_phage_DS_expanded_num <- 
  noprovirus_phage_DS_expanded %>%
  count(type, Host_tax, name = "count")
noprovirus_phage_DS_expanded_num$Host_tax <- 
  factor(noprovirus_phage_DS_expanded_num$Host_tax, levels = c("Ef","Sa","Kp","Ab","Pa","En"))

noprovirus_phage_info <- 
  ESKAPE_all_phage_info[!grepl("^GCA", ESKAPE_all_phage_info$seq_name), ]
table(noprovirus_phage_info$Host_tax)
noprovirus_phage_info_expanded <- noprovirus_phage_info %>%
  separate_rows(Host_tax, sep = ",") %>%
  mutate(Host_tax = trimws(Host_tax))  # 去除可能存在的空格
noprovirus_phage_info_expanded_num <- 
  table(noprovirus_phage_info_expanded$Host_tax) %>%
  as.data.frame()
colnames(noprovirus_phage_info_expanded_num) <- c("Host_tax", "Phage_num")

# 添加归一化列
noprovirus_phage_DS_normalized <- noprovirus_phage_DS_expanded_num %>%
  left_join(noprovirus_phage_info_expanded_num, by = "Host_tax") %>%
  mutate(
    count_normalized = count / Phage_num,
    percent = (count / Phage_num) * 100  # 可选：转换为百分比
  ) %>%
  arrange(type, Host_tax)

library(scales)
p_dot_heatmap_noprovirus_phage_DS <-
  ggplot(noprovirus_phage_DS_normalized, aes(x = Host_tax, y = type)) +
  geom_point(aes(size = count_normalized, color = Host_tax), alpha = 0.7) +
  scale_size_continuous(
    name = "count_normalized",
    range = c(1, 10),  # 控制点的大小范围
    breaks = pretty_breaks(n = 5)
  ) +
  scale_color_manual(values = c("Ef"="#9f86c0","Sa"="#FFBE7D",
                                "Kp"="#76c893","Ab"="#F0E442","Pa"="#ff99c8",
                                "En"="#56B4E9")
  )+
  coord_flip()+
  #  facet_grid(~ MGE_type, scales = "free", space = "free") +
  #facet_wrap(~ MGE_type, scales = "free", ncol = 5) +
  theme(
    panel.spacing = unit(0.1, "lines"),  # 将分面间距设为0
    panel.background = element_blank(),  # 设置背景为透明
    # panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_blank(),
    #axis.title = element_text( size = 12, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  # 横轴标签旋转
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    #legend.position = 'none',
    legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_blank(), #设置图例标题字体为粗体，大小为12，颜色为深灰色
    #legend.text = element_blank(), #设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", linewidth = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
p_dot_heatmap_noprovirus_phage_DS

ggsave(plot = p_dot_heatmap_noprovirus_phage_DS,
       '1_p_dot_heatmap_noprovirus_phage_DS_legend.pdf',
       width = 11,height = 3,units = 'in')

#不同ESKAPE的phage上携带的ARG####
load("./1_phage_DF.RData")

All_phage_ARG <- read.delim2("rgiargs_detail_info.txt")
All_phage_ARG <- All_phage_ARG %>%
  mutate(
    ORF = gsub(" \\#.*","",ORF_ID),
    contig = gsub("_\\d+$","",ORF)
  )

All_phage_ARG <- All_phage_ARG %>%
  mutate(
    # 按"#"分割字符串
    parts = str_split_fixed(ORF_ID, " # ", n = 5),
    ORF_start = as.numeric(parts[, 2]),  # 第二部分是起始位置
    ORF_end = as.numeric(parts[, 3])     # 第三部分是结束位置
  ) %>%
  select(-parts)  # 移除中间列

noprovirus_phage_ARG_info <- 
  All_phage_ARG[!grepl("^GCA", All_phage_ARG$contig), ]
table(noprovirus_phage_ARG_info$Drug.Class)

# 拆分Drug.Class列
noprovirus_phage_ARG_expanded <- 
  noprovirus_phage_ARG_info[noprovirus_phage_ARG_info$Drug.Class!="",] %>%
  separate_rows(Drug.Class, sep = ";") %>%
  mutate(Drug.Class = trimws(Drug.Class))  # 去除首尾空格

noprovirus_phage_ARG_expanded <-
  noprovirus_phage_ARG_expanded %>% 
  left_join(select(ESKAPE_all_phage_info,seq_name,Host_tax), 
            by = c("contig"="seq_name"))

noprovirus_phage_ARG_expanded_num <- 
  noprovirus_phage_ARG_expanded %>%
  count(Drug.Class, Host_tax, name = "count")
noprovirus_phage_ARG_expanded_num$Host_tax <- 
  factor(noprovirus_phage_ARG_expanded_num$Host_tax, 
         levels = c("Ef","Sa","Kp","Ab","Pa","En"))

# 添加归一化列
noprovirus_phage_ARG_expanded_num_normalized <- 
  noprovirus_phage_ARG_expanded_num %>%
  left_join(noprovirus_phage_info_expanded_num, by = "Host_tax") %>%
  mutate(
    count_normalized = count / Phage_num,
    percent = (count / Phage_num) * 100  # 可选：转换为百分比
  ) 
noprovirus_phage_ARG_expanded_num_normalized$Host_tax <- 
  factor(noprovirus_phage_ARG_expanded_num_normalized$Host_tax, 
         levels = c("Ef","Sa","Kp","Ab","Pa","En"))


library(scales)
p_dot_heatmap_noprovirus_phage_ARG <-
  ggplot(noprovirus_phage_ARG_expanded_num_normalized, 
         aes(x = Host_tax, y = Drug.Class)) +
  geom_point(aes(size = count_normalized, color = Host_tax), alpha = 0.7) +
  scale_size_continuous(
    name = "count_normalized",
    range = c(1, 10),  # 控制点的大小范围
    breaks = pretty_breaks(n = 5)
  ) +
  scale_color_manual(values = c("Ef"="#9f86c0","Sa"="#FFBE7D",
                                "Kp"="#76c893","Ab"="#F0E442","Pa"="#ff99c8",
                                "En"="#56B4E9")
  )+
  #  facet_grid(~ MGE_type, scales = "free", space = "free") +
  #facet_wrap(~ MGE_type, scales = "free", ncol = 5) +
  theme(
    panel.spacing = unit(0.1, "lines"),  # 将分面间距设为0
    panel.background = element_blank(),  # 设置背景为透明
    # panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, linewidth = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_blank(),
    #axis.title = element_text( size = 12, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),  # 横轴标签旋转
    axis.ticks = element_line(color = "grey20", linewidth = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_blank(), #设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_blank(), #设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", linewidth = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
p_dot_heatmap_noprovirus_phage_ARG

ggsave(plot = p_dot_heatmap_noprovirus_phage_ARG,
       '1_p_dot_heatmap_noprovirus_phage_ARG_legend.pdf',
       width = 5,height = 4.66,units = 'in')

library(Cairo)
library(cowplot)
g_noprovirus_phage_ARG <- plot_grid(g_noprovirus_phage_antiDS,
                                    p_dot_heatmap_noprovirus_phage_ARG,
                #labels = LETTERS[2], 
                rel_widths = c(1, 1),
                nrow = 1)
g_noprovirus_phage_ARG

ggsave(plot = g_noprovirus_phage_ARG,
       '1_g_noprovirus_phage_ARG.pdf',
       width = 6.98,height = 4.16,units = 'in')

save.image("./1_phage_DF.RData")


#ESKAPE病毒的lifestyle比较####

##所有病毒整体情况####

bacphip <- read.delim2("ESKAPE_all_phage_info_filter.fa.bacphlip")
colnames(bacphip)[1] <- "seq_name"
bacphip <- bacphip %>%
  mutate(Lifestyle = case_when(
    Virulent > 0.9 ~ "Virulent",
    Temperate > 0.9 ~ "Temperate",
    TRUE ~ "Unknown"
  ))
table(bacphip$Lifestyle)

ESKAPE_all_phage_info2 <- ESKAPE_all_phage_info %>%
  left_join(bacphip, by = c("contig_id"="seq_name"))

ESKAPE_all_phage_info2$Lifestyle <- ifelse(
  grepl("^GCA", 
        ESKAPE_all_phage_info2$contig_id, ignore.case = TRUE),
  "Temperate",
  ESKAPE_all_phage_info2$Lifestyle
)
table(ESKAPE_all_phage_info2$Lifestyle)

##病毒生活方式统计###
all_phage_lifestyle <- table(ESKAPE_all_phage_info2$Lifestyle)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

# 准备数据
df_lifestyle <- as.data.frame(all_phage_lifestyle)
colnames(df_lifestyle) <- c("Lifestyle", "Count")

# 计算百分比和标签位置
df_lifestyle <- df_lifestyle %>%
  mutate(
    Percentage = Count / sum(Count) * 100,
    Label = paste0(Count, " (", round(Percentage, 1), "%)"),
    ymax = cumsum(Percentage),
    ymin = c(0, head(ymax, n = -1))
  )

# 创建环形图
p_lifestyle_pie <- ggplot(df_lifestyle, aes(ymax = ymax, ymin = ymin, 
                                        xmax = 4, xmin = 3, 
                                        fill = Lifestyle)) +
  geom_rect() +
  coord_polar(theta = "y") +
  xlim(c(2, 4)) +  # 控制环形的宽度
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right"
  ) +
  scale_fill_brewer(palette = "Set2") +  # 使用漂亮的颜色
  geom_text(aes(x = 3.5, y = (ymin + ymax) / 2, 
                label = paste0(Lifestyle, "\n", Label)),
            size = 3.5, color = "black")

print(p_lifestyle_pie)

ggsave(file.path(CS_WORK, "1_p_lifestyle_pie.pdf"),p_lifestyle_pie,
       width = 3.52,height = 2.54,units = 'in')

##不同宿主病毒的裂解比例####
table(ESKAPE_all_phage_info2$Host_tax)
table(ESKAPE_all_phage_info2$Lifestyle)

# 1. 创建交叉表（计数）
lifestyle_by_host <- table(
  ESKAPE_all_phage_info2$Host_tax,
  ESKAPE_all_phage_info2$Lifestyle
)

library(ggplot2)
library(dplyr)
library(scales)

# 1. 条形图（堆叠比例）
p_lifestyle_by_host <- ESKAPE_all_phage_info2 %>%
  filter(!Host_tax %in% c("Ab,Kp", "Ab,Pa")) %>%
  mutate(Host_tax = factor(Host_tax, 
                           levels = c("Ef", "Sa", "Kp", "Ab", "Pa", "En", "En,Kp"))) %>%
  count(Host_tax, Lifestyle) %>%
  group_by(Host_tax) %>%
  mutate(prop = n / sum(n)) %>%
  ggplot(aes(x = Host_tax, y = prop, fill = Lifestyle)) +
  geom_col(position = "fill") +
  geom_text(aes(label = paste0(round(prop * 100, 1), "%")),
            position = position_fill(vjust = 0.5),
            size = 3.5, color = "black",fontface = "bold") +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values = c("Temperate" = "#8da0cb", 
                               "Unknown" = "#CCCCCC", 
                               "Virulent" = "#fc8d62")) +
  theme(
    panel.spacing = unit(0.1, "lines"),  # 将分面间距设为0
    panel.background = element_blank(),  # 设置背景为透明
    # panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, linewidth = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_blank(),
    #axis.title = element_text( size = 12, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),  # 横轴标签旋转
    axis.ticks = element_line(color = "grey20", linewidth = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_blank(), #设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_blank(), #设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", linewidth = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )

print(p_lifestyle_by_host)

ggsave(file.path(CS_WORK, "1_p_lifestyle_by_host.pdf"),p_lifestyle_by_host,
       width = 2.5,height = 2.54,units = 'in')

# 2. 条形图（计数）
# p_lifestyle_by_host <- ESKAPE_all_phage_info2 %>%
#   count(Host_tax, Lifestyle) %>%
#   ggplot(aes(x = Host_tax, y = n, fill = Lifestyle)) +
#   geom_col(position = "dodge") +
#   geom_text(aes(label = n), 
#             position = position_dodge(width = 0.9),
#             vjust = -0.5, size = 3) +
#   scale_fill_manual(values = c("Temperate" = "#2E8B57", 
#                                "Unknown" = "#808080", 
#                                "Virulent" = "#CD5C5C")) +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     plot.title = element_text(hjust = 0.5, face = "bold")
#   )
# 
# print(p_lifestyle_by_host)
#ESKAPE病毒的分子特征比较####

save.image("./1_phage_DF.RData")

load("./1_phage_DF.RData")

#比较不同lifestyle的病毒中的DS/anti-defense比例####
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

ESKAPE_all_phage_DF2 <- ESKAPE_all_phage_DF
ESKAPE_all_phage_DF2 <- ESKAPE_all_phage_DF2 %>%
  mutate(
    contig = str_replace(sys_beg, "_\\d+$", "")
  )

ESKAPE_all_phage_info2

ESKAPE_all_phage_DF2 <- ESKAPE_all_phage_DF2 %>%
  left_join(
    ESKAPE_all_phage_info2 %>% 
      select(ID, Lifestyle),
    by = c("contig" = "ID")
  )

ESKAPE_all_phage_DF2 <- distinct(ESKAPE_all_phage_DF2)
table(ESKAPE_all_phage_DF2$activity)

DS_temperate_lytic <- 
  ESKAPE_all_phage_DF2[ESKAPE_all_phage_DF2$Lifestyle %in% c("Temperate", "Virulent"),]
length(unique(DS_temperate_lytic[DS_temperate_lytic$Lifestyle=="Temperate",]$contig))
length(unique(DS_temperate_lytic[DS_temperate_lytic$Lifestyle=="Virulent",]$contig))


Lifestyle_activity_count <- ESKAPE_all_phage_DF2 %>%
  count(Lifestyle, activity)

plot_df <- Lifestyle_activity_count %>%
  filter(Lifestyle %in% c("Temperate", "Virulent")) %>%
  group_by(Lifestyle) %>%
  mutate(Proportion = n / sum(n))

ggplot(plot_df, aes(x = Lifestyle,
                    y = Proportion,
                    fill = activity)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c(
    "Defense" = "#4DBBD5",
    "Antidefense" = "#E64B35"
  ))+
  labs(x = NULL,
       y = "Proportion",
       fill = "Activity") +
  theme_classic(base_size = 14)

chisq.test(
  table(
    ESKAPE_all_phage_DF2$Lifestyle,
    ESKAPE_all_phage_DF2$activity
  )[c("Temperate","Virulent"), ]
)


library(ggplot2)

p_DS_antiDS_proportion_by_lifestyle <-
  ggplot(plot_df,
            aes(x = Lifestyle,
                y = Proportion,
                fill = activity)) +
  geom_bar(stat = "identity",
           width = 0.65,
           color = "black"
           ) +
  # geom_text(aes(label = Percent_label),
  #           position = position_stack(vjust = 0.5),
  #           size = 5,
  #           family = "Arial") +
  scale_fill_manual(values = c(
    "Defense" = "#1f77b4",
    "Antidefense" = "#ff7f0e"
  )) +
  scale_y_continuous(
    limits = c(0,1),
    expand = c(0,0),
    labels = scales::percent_format()
  ) +
  labs(x = NULL,
       y = "Proportion",
       fill = NULL) +
  theme(
    panel.spacing = unit(0.1, "lines"),  # 将分面间距设为0
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    #panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, linewidth = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_blank(),
    #axis.title = element_text( size = 12, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 12, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),  # 横轴标签旋转
    axis.ticks = element_line(color = "grey20", linewidth = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_blank(), #设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_blank(), #设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", linewidth = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
p_DS_antiDS_proportion_by_lifestyle

ggsave(file.path(CS_WORK, "1_p_DS_antiDS_proportion_by_lifestyle.pdf"),p_DS_antiDS_proportion_by_lifestyle,
       width = 2.5,height = 3,units = 'in')





