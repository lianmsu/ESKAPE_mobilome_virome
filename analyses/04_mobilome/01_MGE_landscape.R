# ===========================================================================
# ESKAPE_mobilome_virome -- analyses / 04_mobilome
# MGE counts per lineage, ARG and anti-defense cargo, nested MGE cases
# Manuscript      : Fig. 2
# Original path   : 2ESKAPE_MGE/4genome_MGE/1_ESKAPE_MGE.R
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

library(dplyr)
library(stringr)
library(tidyr)
library(fuzzyjoin)
library(ggplot2)
library(scales)
setwd(file.path(CS_DATA, "4genome_MGE"))

#prophage###

genomad_prophage <- read.delim2("all_4538_ESKAPE_genomes_renamed_virus_summary.tsv")
# 过滤不包含plasmid的行
genomad_prophage <- genomad_prophage %>%
  filter(!grepl("plasmid", seq_name, ignore.case = TRUE))
# 提取 | 前面的部分
genomad_prophage$contig <- sub("\\|.*$", "", genomad_prophage$seq_name)
genomad_prophage$genome <- gsub("-.*","",genomad_prophage$contig)
genomad_prophage <- genomad_prophage %>%
  mutate(
    position = str_extract(seq_name, "(?<=provirus_)\\d+_\\d+"),       # 提取"数字_数字"格式的位置信息
    begin = as.numeric(str_extract(position, "^\\d+")),    # 提取_前的数字作为起始位置
    end = as.numeric(str_extract(position, "\\d+$"))       # 提取_后的数字作为结束位置
  ) %>% 
  select(-position)

genomad_prophage <- genomad_prophage %>%
  mutate(
    begin = ifelse(topology != "Provirus", 1, begin),
    end = ifelse(topology != "Provirus", length, end)
  )

genomad_prophage$MGE_type <- "Prophage"

#Plasmid##
plasmid_info <- read.delim2("plasmid_info.tsv")
plasmid_info <- plasmid_info %>%
  mutate(
    genome = gsub("-.*","",plasmid_info$contig)
  )

plasmid_info$MGE_type <- "Plasmid"

#IS##
IS <- read.delim2("4538genomes_filtered_completeIS.tsv")
IS <- IS %>% 
  mutate(
    genome = gsub("-.*","",IS$seqID),
    contig = seqID
  )

IS$MGE_type <- "IS"

#integron##
integrons <- read.delim2("all_4538_ESKAPE_genomes_renamed.integrons")
complete_integrons <- integrons[integrons$type == "complete",]
# 处理数据
complete_integrons_clean <- complete_integrons %>%
  group_by(ID_replicon) %>%          # 按ID_replicon分组
  summarise(
    min_pos_beg = min(pos_beg),      # 最小起始位置
    max_pos_end = max(pos_end),      # 最大结束位置
    length = max_pos_end - min_pos_beg +1,
    .groups = "drop"                 # 取消分组
  ) %>%
  arrange(ID_replicon)               # 按ID排序

complete_integrons_clean <- complete_integrons_clean %>%
  mutate(
    genome = gsub("-.*","", ID_replicon),
    contig = ID_replicon
  )

complete_integrons_clean$MGE_type <- "Integron"

#ICE###
ice <- read.delim2("processed_ice_sequences_ID.txt",header = F)
colnames(ice) <- c("seq_name","range","length")
ice$seq_name <- gsub("\\.","_",ice$seq_name)
ice <- separate(ice, range, into = c("num1", "num2"), sep = "-", convert = TRUE) %>%
  mutate(
    start = pmin(num1, num2),  # 取较小值作为 start
    end = pmax(num1, num2)     # 取较大值作为 end
  ) %>%
  select(-num1, -num2)  # 移除中间列

ice$length <- ice$end - ice$start + 1
ice$genome <- gsub("-.*","",ice$seq_name)

ice$MGE_type <- "ICEs"

#MGE数量展示####

plasmid_info$start <- 1
plasmid_info$end <- plasmid_info$length
plasmid_n <- select(plasmid_info, genome, contig, start, end, length, MGE_type)

genomad_prophage$start <- genomad_prophage$begin
Prophage_n <- select(genomad_prophage, genome, contig, start, end, length, MGE_type)

IS$start <- IS$isBegin
IS$end <- IS$isEnd
IS$length <- IS$isLen
IS_n <- select(IS, genome, contig, start, end, length, MGE_type)

library(dplyr)

# 按基因组和contig分组，排序
IS_sorted <- IS_n %>%
  arrange(genome, contig, start)

# 找出邻近的IS对
find_composite_transposons <- function(df, max_distance = 100000) {
  results <- data.frame()
  
  # 按contig分组
  for(ctg in unique(df$contig)) {
    ctg_data <- df %>% filter(contig == ctg) %>% arrange(start)
    
    if(nrow(ctg_data) >= 2) {
      for(i in 1:(nrow(ctg_data)-1)) {
        for(j in (i+1):nrow(ctg_data)) {
          # 计算两个IS之间的距离
          distance <- ctg_data$start[j] - ctg_data$end[i]
          
          # 如果距离在阈值内，可能是复合转座子
          if(distance > 0 & distance <= max_distance) {
            # 确定方向
            direction <- ifelse(
              ctg_data$start[i] < ctg_data$start[j] & 
                ctg_data$end[i] < ctg_data$end[j],
              "same", "opposite"
            )
            
            results <- rbind(results, data.frame(
              genome = ctg_data$genome[i],
              contig = ctg,
              IS1_start = ctg_data$start[i],
              IS1_end = ctg_data$end[i],
              IS2_start = ctg_data$start[j],
              IS2_end = ctg_data$end[j],
              distance = distance,
              direction = direction,
              total_length = (ctg_data$end[j] - ctg_data$start[i]),
              IS1_length = ctg_data$length[i],
              IS2_length = ctg_data$length[j]
            ))
          }
        }
      }
    }
  }
  return(results)
}

# 应用函数
composite_candidates <- find_composite_transposons(IS_n, max_distance = 100000)
composite_transposons_n <- 
  select(composite_candidates, genome, contig, IS1_start, IS2_end, total_length)
composite_transposons_n$MGE_type <- "Composite_transposon"
colnames(composite_transposons_n) <- c("genome","contig","start","end","length","MGE_type")



complete_integrons_clean$start <- complete_integrons_clean$min_pos_beg
complete_integrons_clean$end <- complete_integrons_clean$max_pos_end
integron_n <- select(complete_integrons_clean, genome, contig, start, end, length, MGE_type)

write.table(integron_n,
            "ESKAPE_integron_n.tsv",
            sep = "\t",           # TSV使用制表符分隔
            row.names = FALSE,    # 不写入行名
            col.names = TRUE,     # 写入列名
            quote = FALSE        # 不添加引号
            )            # NA值的表示方式

ice$contig <- ice$seq_name
ice_n <- select(ice, genome, contig, start, end, length, MGE_type)

all_MGE <- rbind(plasmid_n, Prophage_n, IS_n, integron_n, ice_n)

ESKAPE_genome_info <- read.delim2("ESKAPE_genome_info.tsv")
ESKAPE_genome_info$genome <- gsub("\\.","_",ESKAPE_genome_info$genome)

# 假设两个数据框的匹配列都叫 "genome"
all_MGE <- left_join(all_MGE, 
                     select(ESKAPE_genome_info, genome, Tax),
                     by = "genome")

mge_wide <- all_MGE %>%
  count(genome, MGE_type, name = "count") %>%
  pivot_wider(
    id_cols = genome,
    names_from = MGE_type,
    values_from = count,
    values_fill = 0  # 关键：将NA填充为0
  )
ESKAPE_genome_info <- left_join(ESKAPE_genome_info, 
                        mge_wide,
                        by = "genome")
ESKAPE_genome_info <- ESKAPE_genome_info %>%
  mutate(across(19:23, ~ replace(., is.na(.), 0)))

write.table(ESKAPE_genome_info,
            "ESKAPE_genome_info2.tsv",
            sep = "\t",           # TSV使用制表符分隔
            row.names = FALSE,    # 不写入行名
            col.names = TRUE,     # 写入列名
            quote = FALSE        # 不添加引号
          )            # NA值的表示方式

mge_matrix <- select(ESKAPE_genome_info, genome, Prophage, 
                     Plasmid, ICEs, Integron, IS, Tax)
mge_long <- mge_matrix %>%
  pivot_longer(
    cols = c(Prophage, Plasmid, ICEs, Integron, IS),  # 指定要转换的列
    names_to = "MGE_type",       # 新列名：存储原来的列名
    values_to = "count"          # 新列名：存储原来的值
  )

##绘图####
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ggplot2)
library(dplyr)

# 1. 按Tax和MGE_type汇总
tax_summary <- mge_long %>%
  group_by(Tax, MGE_type) %>%
  summarise(
    mean_count = mean(count),
    median_count = median(count),
    total_count = sum(count),
    n_genomes = n(),
    .groups = "drop"
  )
tax_summary$Tax <- factor(tax_summary$Tax, levels = c("Ef","Sa","Kp","Ab","Pa","En"))
# 2. 绘制聚合热图
p_mge_num_heatmap <- 
ggplot(tax_summary, aes(x = MGE_type, y = Tax)) +
  geom_tile(aes(fill = mean_count), color = "white", size = 1) +
  geom_text(aes(label = round(mean_count, 1)), 
            color = "white", size = 4, fontface = "plain") +
  scale_fill_gradient2(
    low = "white",
    mid = "#1982c4",
    high = "#ff9f1c",
    midpoint = 2,
    name = "Mean Count"
  ) +
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
p_mge_num_heatmap

##绘图####
library(rstatix)
md_p_val2 <- mge_long %>%
  # 先过滤掉全为0的分组组合
  group_by(Tax, MGE_type) %>%
  filter(!all(count == 0)) %>%
  ungroup() %>%
  # 再进行检验
  group_by(MGE_type) %>%
  wilcox_test(count ~ Tax) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj") %>%
  add_xy_position(x = "Tax", dodge = 0.5, group = "MGE_type")

library(ggplot2)
library(ggsignif)#geom_signif
library(ggpubr)#stat_pvalue_manual stat_compare_means
library(gghalves)
library(ggsci)
library(ggbeeswarm)
p_mge_tax_box <- ggplot(mge_long, 
                        aes(x = Tax, y = count)) +
  geom_quasirandom(
    aes(
      x = Tax,
      y = count,
      color = Tax
    ),
    dodge.width = 0,
    width = 0.25,   # 控制横向展开
    varwidth = FALSE,
    size = 0.5
  )+
  geom_boxplot(
    aes(fill = Tax),
    color = "black",
    linewidth = 0.4,
    alpha = 0.5,          # ← 整体透明度
    outlier.shape = NA,
    width = 0.8
  )+

#  stat_pvalue_manual(md_p_val2,label = "p.adj.signif",label.size=5,hide.ns = T)+
  facet_wrap(.~MGE_type,nrow=5, scales = "free_y")+ #分⾯操作
  # 定义颜⾊
  # scale_fill_npg()+
  # scale_color_npg()+
  scale_fill_manual(values = c("Ef"="#9f86c0","Sa"="#FFBE7D",
                                "Kp"="#76c893","Ab"="#F0E442","Pa"="#ff99c8",
                                "En"="#56B4E9")
  )+
  scale_color_manual(values = c("Ef"="#9f86c0","Sa"="#FFBE7D",
                                "Kp"="#76c893","Ab"="#F0E442","Pa"="#ff99c8",
                                "En"="#56B4E9")
  )+
  labs(x="River",y="log10 (DS Density +1) × 1e-3")+
  theme(
    panel.spacing = unit(0.1, "lines"),  # 将分面间距设为0
    strip.background = element_blank(),      # 移除分面背景
    strip.text = element_blank(),            # 移除分面标题文字
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
p_mge_tax_box

ggsave(plot = p_mge_tax_box,'1_p_mge_tax_box.pdf',width = 2.28,height = 4.37,units = 'in')


#MGE与ARG关系####
ARG_detail_info <- read.delim2("rgiargs_detail_info.txt")
ARG_detail_info <- ARG_detail_info %>%
  mutate(
    genome = gsub("\\-.*","",ORF_ID),
    # 方法A：使用str_split
    contig = sapply(strsplit(ORF_ID, "_", fixed = TRUE), 
                    function(x) paste(x[1:3], collapse = "_"))
  )

ARG_detail_info <- ARG_detail_info %>%
  mutate(
    # 按"#"分割字符串
    parts = str_split_fixed(ORF_ID, " # ", n = 5),
    ORF_start = as.numeric(parts[, 2]),  # 第二部分是起始位置
    ORF_end = as.numeric(parts[, 3])     # 第三部分是结束位置
  ) %>%
  select(-parts)  # 移除中间列

##Prophage上的ARG####

#对ARG_detail_info进行过滤，只保留contig列与genomad_prophage$contig一致，
#且ORF_start大于等于genomad_prophage$begin，且ORF_end小于等于genomad_prophage$end的行

# 进行连接和过滤
prophage_ARG <- ARG_detail_info %>%
  inner_join(genomad_prophage, 
             by = c("contig" = "contig")) %>%
  filter(
    ORF_start >= begin,
    ORF_end <= end
  ) 
# 拆分Drug.Class列
prophage_ARG <- prophage_ARG[prophage_ARG$Drug.Class!="",]
prophage_ARG <- prophage_ARG %>%
  separate_rows(Drug.Class, sep = ";")

# 去除空白字符
prophage_ARG <- prophage_ARG %>%
  mutate(Drug.Class = trimws(Drug.Class))  # 去除首尾空格

prophage_ARG <- select(prophage_ARG, genome.x, Drug.Class, MGE_type)

colnames(prophage_ARG)[1] <- "genome"
prophage_ARG <- left_join(prophage_ARG, 
                        select(ESKAPE_genome_info, genome, Tax),
                        by = "genome")
##plasmid上的ARG####
# 进行连接和过滤
plasmid_ARG <- ARG_detail_info %>%
  inner_join(plasmid_n, 
             by = c("contig" = "contig")) %>%
  filter(
    ORF_start >= start,
    ORF_end <= end
  )
# 拆分Drug.Class列
plasmid_ARG <- plasmid_ARG[plasmid_ARG$Drug.Class!="",]
plasmid_ARG <- plasmid_ARG %>%
  separate_rows(Drug.Class, sep = ";")

# 去除空白字符
plasmid_ARG <- plasmid_ARG %>%
  mutate(Drug.Class = trimws(Drug.Class))  # 去除首尾空格

plasmid_ARG <- select(plasmid_ARG, genome.x, Drug.Class, MGE_type)

colnames(plasmid_ARG)[1] <- "genome"
plasmid_ARG <- left_join(plasmid_ARG, 
                          select(ESKAPE_genome_info, genome, Tax),
                          by = "genome")

##IS上的ARG####
# 进行连接和过滤
IS_ARG <- ARG_detail_info %>%
  inner_join(IS_n, 
             by = c("contig" = "contig")) %>%
  filter(
    ORF_start >= start,
    ORF_end <= end
  )
# 拆分Drug.Class列
IS_ARG <- IS_ARG[IS_ARG$Drug.Class!="",]
IS_ARG <- IS_ARG %>%
  separate_rows(Drug.Class, sep = ";")

# 去除空白字符
IS_ARG <- IS_ARG %>%
  mutate(Drug.Class = trimws(Drug.Class))  # 去除首尾空格

IS_ARG <- select(IS_ARG, genome.x, Drug.Class, MGE_type)

colnames(IS_ARG)[1] <- "genome"
IS_ARG <- left_join(IS_ARG, 
                         select(ESKAPE_genome_info, genome, Tax),
                         by = "genome")

##IS上的ARG####
# 进行连接和过滤
composite_transposons_ARG <- ARG_detail_info %>%
  inner_join(composite_transposons_n, 
             by = c("contig" = "contig")) %>%
  filter(
    ORF_start >= start,
    ORF_end <= end
  )
# 拆分Drug.Class列
composite_transposons_ARG <- 
  composite_transposons_ARG[composite_transposons_ARG$Drug.Class!="",]
composite_transposons_ARG <- composite_transposons_ARG %>%
  separate_rows(Drug.Class, sep = ";")

# 去除空白字符
composite_transposons_ARG <- composite_transposons_ARG %>%
  mutate(Drug.Class = trimws(Drug.Class))  # 去除首尾空格

composite_transposons_ARG <- 
  select(composite_transposons_ARG, genome.x, Drug.Class, MGE_type)

colnames(composite_transposons_ARG)[1] <- "genome"
composite_transposons_ARG <- left_join(composite_transposons_ARG, 
                    select(ESKAPE_genome_info, genome, Tax),
                    by = "genome")

##Integron上的ARG####
# 进行连接和过滤
integron_ARG <- ARG_detail_info %>%
  inner_join(integron_n, 
             by = c("contig" = "contig")) %>%
  filter(
    ORF_start >= start,
    ORF_end <= end
  )
# 拆分Drug.Class列
integron_ARG <- integron_ARG[integron_ARG$Drug.Class!="",]
integron_ARG <- integron_ARG %>%
  separate_rows(Drug.Class, sep = ";")

# 去除空白字符
integron_ARG <- integron_ARG %>%
  mutate(Drug.Class = trimws(Drug.Class))  # 去除首尾空格

integron_ARG <- select(integron_ARG, genome.x, Drug.Class, MGE_type)

colnames(integron_ARG)[1] <- "genome"
integron_ARG <- left_join(integron_ARG, 
                    select(ESKAPE_genome_info, genome, Tax),
                    by = "genome")

##ICE上的ARG####
# 进行连接和过滤
# ice_ARG <- ARG_detail_info %>%
#   inner_join(ice_n, 
#              by = c("contig" = "contig")) %>%
#   filter(
#     ORF_start >= start,
#     ORF_end <= end
#   )
ice_ARG <- ARG_detail_info %>%
  fuzzy_inner_join(
    ice_n,
    by = c("contig" = "contig"),
    match_fun = function(x, y) str_detect(x, fixed(y))
  ) %>%
  filter(
    ORF_start >= start,
    ORF_end <= end
  )

# 拆分Drug.Class列
ice_ARG <- ice_ARG[ice_ARG$Drug.Class!="",]
ice_ARG <- ice_ARG %>%
  separate_rows(Drug.Class, sep = ";")

# 去除空白字符
ice_ARG <- ice_ARG %>%
  mutate(Drug.Class = trimws(Drug.Class))  # 去除首尾空格

ice_ARG <- select(ice_ARG, genome.x, Drug.Class, MGE_type)

colnames(ice_ARG)[1] <- "genome"
ice_ARG <- left_join(ice_ARG, 
                          select(ESKAPE_genome_info, genome, Tax),
                          by = "genome")

##绘图####
all_mge_arg <- 
  rbind(plasmid_ARG, prophage_ARG, composite_transposons_ARG, 
        ice_ARG, integron_ARG)
all_mge_arg <- all_mge_arg[,-1]

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
# 1. 计算计数
dot_heatmap_data <- all_mge_arg %>%
  count(Drug.Class, MGE_type, Tax, name = "count")

dot_heatmap_data$Tax <- factor(dot_heatmap_data$Tax, levels = c("Ef","Sa","Kp","Ab","Pa","En"))

dot_heatmap_data$Drug.Class <- gsub(" antibiotic","", dot_heatmap_data$Drug.Class)
dot_heatmap_data$Drug.Class <- gsub("disinfecting agents and antiseptics",
                            "DAA", dot_heatmap_data$Drug.Class)
dot_heatmap_data$Drug.Class <- gsub("penicillin",
                            "P", dot_heatmap_data$Drug.Class)
dot_heatmap_data$MGE_type[dot_heatmap_data$MGE_type == "Composite_transposon"] <- "IS"


# 2. 检查数据
print(head(dot_heatmap_data))
print(paste("Unique Drug.Class:", n_distinct(dot_heatmap_data$Drug.Class)))
print(paste("Unique MGE_type:", n_distinct(dot_heatmap_data$MGE_type)))
print(paste("Unique Tax:", n_distinct(dot_heatmap_data$Tax)))

# 添加归一化列
dot_heatmap_data_normalized <- dot_heatmap_data %>%
  left_join(tax_summary, by = c("MGE_type", "Tax")) %>%
  mutate(
    count_normalized = count / total_count,
    percent = (count / total_count) * 100  # 可选：转换为百分比
  ) %>%
  arrange(Drug.Class, MGE_type)

# 3. 创建分面热图
p_dot_heatmap_arg_normalized <-
ggplot(dot_heatmap_data_normalized, aes(x = Tax, y = Drug.Class)) +
  geom_point(aes(size = count_normalized, color = Tax), alpha = 0.7) +
  scale_size_continuous(
    name = "Count",
    range = c(1, 10),  # 控制点的大小范围
    breaks = pretty_breaks(n = 5)
  ) +
  scale_color_manual(values = c("Ef"="#9f86c0","Sa"="#FFBE7D",
                               "Kp"="#76c893","Ab"="#F0E442","Pa"="#ff99c8",
                               "En"="#56B4E9")
  )+
  facet_grid(~ MGE_type, scales = "free", space = "free") +
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
p_dot_heatmap_arg_normalized

p_dot_heatmap_arg_normalized_legend <-
  ggplot(dot_heatmap_data_normalized, aes(x = Tax, y = Drug.Class)) +
  geom_point(aes(size = count_normalized, color = Tax), alpha = 0.7) +
  scale_size_continuous(
    name = "Count",
    range = c(1, 10),  # 控制点的大小范围
    breaks = pretty_breaks(n = 5)
  ) +
  scale_color_manual(values = c("Ef"="#9f86c0","Sa"="#FFBE7D",
                                "Kp"="#76c893","Ab"="#F0E442","Pa"="#ff99c8",
                                "En"="#56B4E9")
  )+
  facet_grid(~ MGE_type, scales = "free", space = "free") +
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
    #legend.position = 'none',
    legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    #legend.title = element_blank(), #设置图例标题字体为粗体，大小为12，颜色为深灰色
    #legend.text = element_blank(), #设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", linewidth = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
p_dot_heatmap_arg_normalized_legend

ggsave(plot = p_dot_heatmap_arg_normalized_legend,
       '1_p_dot_heatmap_arg_normalized_legend.pdf',
       width = 5.7,height = 6.2,units = 'in')

#MGE与DS关系####

DF_system_info <- read.delim2("all_4538_ESKAPE_genomes_renamed_defense_finder_systems.tsv")
DF_gene_info <- read.delim2("all_4538_ESKAPE_genomes_renamed_defense_finder_genes.tsv")

ORFgff <- read.delim2("DF_gene_gff.gff", header = F, sep = '\t')
ORFgff <- ORFgff[c(1:107553),]
library(stringr)
ORFgff <- ORFgff %>%
  mutate(
    # 提取ID=后面的内容，然后取最后一个_后的数字
    gene_num = str_extract(V9, "ID=[^;]+") %>% 
      str_remove("ID=") %>% 
      str_extract("[^_]+$"),
    hit_id = paste0(V1, "_", gene_num)
  ) %>%
  select(-gene_num)  # 移除临时列

DF_gene_info <- left_join(DF_gene_info, ORFgff, by = "hit_id")
DF_gene_info <- select(DF_gene_info, replicon, hit_id, V4, V5, gene_name, type, subtype, activity)

#用DF_defense的sys_beg去匹配DG数据框的hit_id列，匹配到的行的V4列的值作为DF_defense对应行的begin
#用DF_defense的sys_end去匹配DG数据框的hit_id列，匹配到的行的V5列的值作为DF_defense对应行的end

# 首先为DF_defense创建begin和end列（如果不存在）
if(!"begin" %in% names(DF_system_info)) DF_system_info$begin <- NA
if(!"DSend" %in% names(DF_system_info)) DF_system_info$DSend <- NA
# 创建从hit_id到V4/V5的映射
v4_mapping <- setNames(DF_gene_info$V4, DF_gene_info$hit_id)
v5_mapping <- setNames(DF_gene_info$V5, DF_gene_info$hit_id)
# 用sys_beg匹配并填充begin列
DF_system_info$begin <- ifelse(
  DF_system_info$sys_beg %in% names(v4_mapping),
  v4_mapping[as.character(DF_system_info$sys_beg)],
  DF_system_info$begin
)
# 用sys_end匹配并填充end列
DF_system_info$DSend <- ifelse(
  DF_system_info$sys_end %in% names(v5_mapping),
  v5_mapping[as.character(DF_system_info$sys_end)],
  DF_system_info$DSend
)



##统计质粒上的DS####
head(plasmid_n)
head(DF_defense)

DF_defense <- DF_system_info[DF_system_info$activity == "Defense",]

DF_defense <- DF_defense %>%
  mutate(
    contig = str_replace(sys_beg, "_\\d+$", "")  
  )

plasmid_DS <- DF_defense %>%
  inner_join(plasmid_n, 
             by = c("contig" = "contig")) %>%
  filter(
    begin >= start,
    DSend <= end
  )
plasmid_DS <- select(plasmid_DS, genome, type, subtype, MGE_type)

plasmid_DS <- left_join(plasmid_DS, 
                     select(ESKAPE_genome_info, genome, Tax),
                     by = "genome")

##prophage上的DS####

prophage_DS <- DF_defense %>%
  inner_join(Prophage_n, 
             by = c("contig" = "contig")) %>%
  filter(
    begin >= start,
    DSend <= end
  )

prophage_DS <- select(prophage_DS, genome, type, subtype, MGE_type)

prophage_DS <- left_join(prophage_DS, 
                        select(ESKAPE_genome_info, genome, Tax),
                        by = "genome")
##is上的DS####

IS_DS <- DF_defense %>%
  inner_join(IS_n, 
             by = c("contig" = "contig")) %>%
  filter(
    begin >= start,
    DSend <= end
  )

IS_DS <- select(IS_DS, genome, type, subtype, MGE_type)

IS_DS <- left_join(IS_DS, 
                         select(ESKAPE_genome_info, genome, Tax),
                         by = "genome")
##composite_transposons上的DS####

composite_transposons_DS <- DF_defense %>%
  inner_join(composite_transposons_n, 
             by = c("contig" = "contig")) %>%
  filter(
    begin >= start,
    DSend <= end
  )

composite_transposons_DS <- 
  select(composite_transposons_DS, genome, type, subtype, MGE_type)

composite_transposons_DS <- left_join(composite_transposons_DS, 
                   select(ESKAPE_genome_info, genome, Tax),
                   by = "genome")



##Integron上的DS####

integron_DS <- DF_defense %>%
  inner_join(integron_n, 
             by = c("contig" = "contig")) %>%
  filter(
    begin >= start,
    DSend <= end
  )

integron_DS <- select(integron_DS, genome, type, subtype, MGE_type)

integron_DS <- left_join(integron_DS, 
                         select(ESKAPE_genome_info, genome, Tax),
                         by = "genome")
##ICE上的DS####

ice_DS <- DF_defense %>%
  fuzzy_inner_join(
    ice_n,
    by = c("contig" = "contig"),
    match_fun = function(x, y) str_detect(x, fixed(y))
  ) %>%
  filter(
    begin >= start,
    DSend <= end
  )

ice_DS <- select(ice_DS, genome, type, subtype, MGE_type)

ice_DS <- left_join(ice_DS, 
                         select(ESKAPE_genome_info, genome, Tax),
                         by = "genome")

##绘图####
all_mge_ds <- rbind(plasmid_DS, prophage_DS, composite_transposons_DS, ice_DS, integron_DS)
all_mge_ds <- all_mge_ds[,-1]


# 1. 计算计数
dot_heatmap_ds <- all_mge_ds %>%
  count(type, MGE_type, Tax, name = "count")

dot_heatmap_ds$Tax <- factor(dot_heatmap_ds$Tax, levels = c("Ef","Sa","Kp","Ab","Pa","En"))

type_totals <- dot_heatmap_ds %>%
  group_by(type) %>%
  summarise(total_count = sum(count, na.rm = TRUE))
type_totals <- type_totals[type_totals$total_count >110,]

dot_heatmap_ds <- dot_heatmap_ds[dot_heatmap_ds$type %in% type_totals$type,]

dot_heatmap_ds$MGE_type[dot_heatmap_ds$MGE_type == "Composite_transposon"] <- "IS"

# 添加归一化列
dot_heatmap_ds_normalized <- dot_heatmap_ds %>%
  left_join(tax_summary, by = c("MGE_type", "Tax")) %>%
  mutate(
    count_normalized = count / total_count,
    percent = (count / total_count) * 100  # 可选：转换为百分比
  ) %>%
  arrange(type, MGE_type)

# 3. 创建分面热图
p_dot_heatmap_ds_normalized <-
  ggplot(dot_heatmap_ds_normalized, aes(x = Tax, y = type)) +
  geom_point(aes(size = count_normalized, color = Tax), alpha = 0.7) +
  scale_size_continuous(
    name = "Count",
    range = c(1, 10),  # 控制点的大小范围
    breaks = pretty_breaks(n = 5)
  ) +
  scale_color_manual(values = c("Ef"="#9f86c0","Sa"="#FFBE7D",
                                "Kp"="#76c893","Ab"="#F0E442","Pa"="#ff99c8",
                                "En"="#56B4E9")
  )+
  facet_grid(~ MGE_type, scales = "free", space = "free") +
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
p_dot_heatmap_ds_normalized

p_dot_heatmap_ds_normalized_legend <-
  ggplot(dot_heatmap_ds_normalized, aes(x = Tax, y = type)) +
  geom_point(aes(size = count_normalized, color = Tax), alpha = 0.7) +
  scale_size_continuous(
    name = "Count",
    range = c(1, 10),  # 控制点的大小范围
    breaks = pretty_breaks(n = 5)
  ) +
  scale_color_manual(values = c("Ef"="#9f86c0","Sa"="#FFBE7D",
                                "Kp"="#76c893","Ab"="#F0E442","Pa"="#ff99c8",
                                "En"="#56B4E9")
  )+
  facet_grid(~ MGE_type, scales = "free", space = "free") +
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
    #legend.position = 'none',
    legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    #legend.title = element_blank(), #设置图例标题字体为粗体，大小为12，颜色为深灰色
    #legend.text = element_blank(), #设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", linewidth = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
p_dot_heatmap_ds_normalized_legend

ggsave(plot = p_dot_heatmap_ds_normalized_legend,
       '1_p_dot_heatmap_ds_normalized_legend.pdf',
       width = 6.9,height = 5.5,units = 'in')



#MGE与anti-DS关系####
##统计质粒上的antiDS####

DF_antidefense <- DF_system_info[DF_system_info$activity == "Antidefense",]

DF_antidefense <- DF_antidefense %>%
  mutate(
    contig = str_replace(sys_beg, "_\\d+$", "")  
  )

plasmid_antiDS <- DF_antidefense %>%
  inner_join(plasmid_n, 
             by = c("contig" = "contig")) %>%
  filter(
    begin >= start,
    DSend <= end
  )
plasmid_antiDS <- select(plasmid_antiDS, genome, type, subtype, MGE_type)

plasmid_antiDS <- left_join(plasmid_antiDS, 
                        select(ESKAPE_genome_info, genome, Tax),
                        by = "genome")

##prophage上的DS####

prophage_antiDS <- DF_antidefense %>%
  inner_join(Prophage_n, 
             by = c("contig" = "contig")) %>%
  filter(
    begin >= start,
    DSend <= end
  )

prophage_antiDS <- select(prophage_antiDS, genome, type, subtype, MGE_type)

prophage_antiDS <- left_join(prophage_antiDS, 
                         select(ESKAPE_genome_info, genome, Tax),
                         by = "genome")
##is上的DS####

IS_antiDS <- DF_antidefense %>%
  inner_join(IS_n, 
             by = c("contig" = "contig")) %>%
  filter(
    begin >= start,
    DSend <= end
  )

IS_antiDS <- select(IS_antiDS, genome, type, subtype, MGE_type)

IS_antiDS <- left_join(IS_antiDS, 
                   select(ESKAPE_genome_info, genome, Tax),
                   by = "genome")
##composite_transposons上的DS####

composite_transposons_antiDS <- DF_antidefense %>%
  inner_join(composite_transposons_n, 
             by = c("contig" = "contig")) %>%
  filter(
    begin >= start,
    DSend <= end
  )

composite_transposons_antiDS <- 
  select(composite_transposons_antiDS, genome, type, subtype, MGE_type)

composite_transposons_antiDS <- left_join(composite_transposons_antiDS, 
                       select(ESKAPE_genome_info, genome, Tax),
                       by = "genome")

##Integron上的DS####

integron_antiDS <- DF_antidefense %>%
  inner_join(integron_n, 
             by = c("contig" = "contig")) %>%
  filter(
    begin >= start,
    DSend <= end
  )

integron_antiDS <- select(integron_antiDS, genome, type, subtype, MGE_type)

integron_antiDS <- left_join(integron_antiDS, 
                         select(ESKAPE_genome_info, genome, Tax),
                         by = "genome")
##ICE上的DS####

ice_antiDS <- DF_antidefense %>%
  fuzzy_inner_join(
    ice_n,
    by = c("contig" = "contig"),
    match_fun = function(x, y) str_detect(x, fixed(y))
  ) %>%
  filter(
    begin >= start,
    DSend <= end
  )

ice_antiDS <- select(ice_antiDS, genome, type, subtype, MGE_type)

ice_antiDS <- left_join(ice_antiDS, 
                    select(ESKAPE_genome_info, genome, Tax),
                    by = "genome")

##绘图####
all_mge_antids <- rbind(plasmid_antiDS, prophage_antiDS, composite_transposons_antiDS, 
                        ice_antiDS, integron_antiDS)
all_mge_antids <- all_mge_antids[,-1]


# 1. 计算计数
dot_heatmap_antids <- all_mge_antids %>%
  count(type, MGE_type, Tax, name = "count")

dot_heatmap_antids$Tax <- factor(dot_heatmap_antids$Tax, levels = c("Ef","Sa","Kp","Ab","Pa","En"))

antitype_totals <- dot_heatmap_antids %>%
  group_by(type) %>%
  summarise(total_count = sum(count, na.rm = TRUE))

dot_heatmap_antids$MGE_type[dot_heatmap_antids$MGE_type == "Composite_transposon"] <- "IS"

# 添加归一化列
dot_heatmap_antids_normalized <- dot_heatmap_antids %>%
  left_join(tax_summary, by = c("MGE_type", "Tax")) %>%
  mutate(
    count_normalized = count / total_count,
    percent = (count / total_count) * 100  # 可选：转换为百分比
  ) %>%
  arrange(type, MGE_type)

# 3. 创建分面热图
p_dot_heatmap_antids_normalized <-
  ggplot(dot_heatmap_antids_normalized, aes(x = Tax, y = type)) +
  geom_point(aes(size = count_normalized, color = Tax), alpha = 0.7) +
  scale_size_continuous(
    name = "Count",
    range = c(1, 10),  # 控制点的大小范围
    breaks = pretty_breaks(n = 5)
  ) +
  scale_color_manual(values = c("Ef"="#9f86c0","Sa"="#FFBE7D",
                                "Kp"="#76c893","Ab"="#F0E442","Pa"="#ff99c8",
                                "En"="#56B4E9")
  )+
  facet_grid(~ MGE_type, scales = "free", space = "free") +
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
p_dot_heatmap_antids_normalized

p_dot_heatmap_antids_normalized_legend <-
  ggplot(dot_heatmap_antids_normalized, aes(x = Tax, y = type)) +
  geom_point(aes(size = count_normalized, color = Tax), alpha = 0.7) +
  scale_size_continuous(
    name = "Count",
    range = c(1, 10),  # 控制点的大小范围
    breaks = pretty_breaks(n = 5)
  ) +
  scale_color_manual(values = c("Ef"="#9f86c0","Sa"="#FFBE7D",
                                "Kp"="#76c893","Ab"="#F0E442","Pa"="#ff99c8",
                                "En"="#56B4E9")
  )+
  facet_grid(~ MGE_type, scales = "free", space = "free") +
  #facet_wrap(~ MGE_type, scales = "free", ncol = 5) +
  theme(
    panel.spacing = unit(0.1, "lines"),  # 将分面间距设为0
    panel.background = element_blank(),  # 设置背景为透明
    #panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_blank(),
    #axis.title = element_text( size = 12, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),  # 横轴标签旋转
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    #legend.position = 'none',
    legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    #legend.title = element_blank(), #设置图例标题字体为粗体，大小为12，颜色为深灰色
    #legend.text = element_blank(), #设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", linewidth = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
p_dot_heatmap_antids_normalized_legend

ggsave(plot = p_dot_heatmap_antids_normalized_legend,
       '1_p_dot_heatmap_antids_normalized_legend.pdf',
       width = 5.5,height = 6.2,units = 'in')


library(Cairo)
library(cowplot)

g1 <- plot_grid(p_mge_num_heatmap,p_dot_heatmap_antids_normalized,
                labels = LETTERS[1:2], ncol = 2,
                rel_widths = c(1, 2)
                # label_x = 0.09,            # 标签x坐标（1为最右侧）
                #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g1

ggsave(plot = g1,'1_MGEnum_antiDS2.pdf',width = 8,height = 2.5,units = 'in')


g2 <- plot_grid(p_dot_heatmap_arg_normalized,p_dot_heatmap_ds_normalized,
                labels = LETTERS[3:4], ncol = 1
                #rel_widths = c(3, 2)
                # label_x = 0.09,            # 标签x坐标（1为最右侧）
                #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g2
ggsave(plot = g2,'1_ARG_DS_inMGE2.pdf',width = 5.7,height = 9,units = 'in')

save.image("1_ESKAPE_MGE.RData")
#load("1_ESKAPE_MGE.RData")


#MGE之间的重叠，携带ARG####
#prophage与plasmid不考虑
#MGE interaction，使得ARG从胞内转移变得有胞外转移的可能/传递到别的细菌的可能
setwd(file.path(CS_DATA, "4genome_MGE"))
load("1_ESKAPE_MGE.RData")

##plasmid-ice####
ice_n2 <- plasmid_n %>%
  fuzzy_inner_join(
    ice_n,
    by = c("contig" = "contig"),
    match_fun = function(x, y) str_detect(x, fixed(y))
  )
ice_n2 <- select(ice_n2, contig.x, contig.y)

ice_n3 <- ice_n
# 找到匹配的索引
match_idx <- match(ice_n3$contig, ice_n2$contig.y)
# 更新contig列
ice_n3 <- ice_n3 %>%
  mutate(
    contig = ifelse(!is.na(match_idx), 
                    ice_n2$contig.x[match_idx], 
                    contig)
  )

ice_in_plasmid <- ice_n3[ice_n3$contig,]
ice_in_plasmid <- ice_n3[grepl("plasmid", ice_n3$contig, ignore.case = TRUE), ]

ice_in_plasmid_ARG <- ARG_detail_info %>%
  inner_join(ice_in_plasmid, 
             by = c("contig" = "contig")) %>%
  filter(
    ORF_start >= start,
    ORF_end <= end
  )

##plasmid-IS####
IS_in_plasmid <- IS_n[grepl("plasmid", IS_n$contig, ignore.case = TRUE), ]
IS_in_plasmid_ARG <- ARG_detail_info %>%
  inner_join(IS_in_plasmid, 
             by = c("contig" = "contig")) %>%
  filter(
    ORF_start >= start,
    ORF_end <= end
  )

##plasmid-integron####
integron_in_plasmid <- integron_n[grepl("plasmid", integron_n$contig, ignore.case = TRUE), ]
integron_in_plasmid_ARG <- ARG_detail_info %>%
  inner_join(integron_in_plasmid, 
             by = c("contig" = "contig")) %>%
  filter(
    ORF_start >= start,
    ORF_end <= end
  )

##integron-in-plasmid与ice-in-plasmid####
integron_in_plasmid
ice_in_plasmid

# 假设两个数据框都有: contig, start, end 列
plasmid_integron_ice_overlap_rows <- fuzzy_inner_join(
  integron_in_plasmid,
  ice_in_plasmid,
  by = c("contig", "start", "end"),
  match_fun = list(`==`, function(x, y) x <= y, function(x, y) x >= y)
) %>%
  filter(start.x <= end.y & end.x >= start.y)  # 确保范围确实重叠

plasmid_integron_ice_overlap_ARG <- ARG_detail_info %>%
  inner_join(plasmid_integron_ice_overlap_rows, 
             by = c("contig" = "contig.x")) %>%
  filter(
    ORF_start >= start.x,
    ORF_end <= end.x
  )

plasmid_integron_ice_overlap_ARG_plasmid <- plasmid_n[plasmid_n$contig %in% 
                                                        plasmid_integron_ice_overlap_ARG$contig,]
plasmid_integron_ice_overlap_ARG_plasmid <- left_join(plasmid_integron_ice_overlap_ARG_plasmid, 
                                                      select(ESKAPE_genome_info, genome, Tax),
                                                      by = "genome")
#GCA_015287605_1-2-plasmid GCA_051091715_1-2-plasmid, 质粒上携带integron和ice

##plasmid selection####

genomad_plasmid <- read.delim2("all_4538_ESKAPE_genomes_renamed_plasmid_summary.tsv")
genomad_plasmid_complete <- genomad_plasmid[genomad_plasmid$topology %in% c("DTR","ITR"),]

selected_plasmid_withice <- ice_in_plasmid_ARG[ice_in_plasmid_ARG$contig %in% genomad_plasmid_complete$seq_name,]
selected_plasmid_withintegron <- integron_in_plasmid_ARG[integron_in_plasmid_ARG$contig %in% genomad_plasmid_complete$seq_name,]

#GCA_052906805_1-4-plasmid
#GCA_009387815_1-2-plasmid



##prophage-integron####

# 假设两个数据框都有: contig, start, end 列
library(data.table)

setDT(integron_n)
setDT(Prophage_n)

# 为 foverlaps 准备区间列（必须是 start <= end）
setkey(integron_n, contig, start, end)
setkey(Prophage_n, contig, start, end)

integron_in_prophage <- foverlaps(
  integron_n,
  Prophage_n,
  by.x = c("contig", "start", "end"),
  by.y = c("contig", "start", "end"),
  type = "within",      # 任意区间重叠
  nomatch = 0L       # 等价于 inner join
)


integron_in_prophage_ARG <- ARG_detail_info %>%
  inner_join(integron_in_prophage, 
             by = c("contig" = "contig")) %>%
  filter(
    ORF_start >= i.start,
    ORF_end <= i.end
  )

##prophage-ice####

# 假设两个数据框都有: contig, start, end 列
library(data.table)

setDT(ice_n)
setDT(Prophage_n)

# 为 foverlaps 准备区间列（必须是 start <= end）
setkey(ice_n, contig, start, end)
setkey(Prophage_n, contig, start, end)

ice_in_prophage <- foverlaps(
  ice_n,
  Prophage_n,
  by.x = c("contig", "start", "end"),
  by.y = c("contig", "start", "end"),
  type = "within",      # 任意区间重叠
  nomatch = 0L       # 等价于 inner join
)


ice_in_prophage_ARG <- ARG_detail_info %>%
  inner_join(ice_in_prophage, 
             by = c("contig" = "contig")) %>%
  filter(
    ORF_start >= i.start,
    ORF_end <=  i.end
  )

##phage selection####

#GCA_026073395_1-1
#GCA_006304565_1-1

# 提取包含 "GCA_026073395_1-1" 的行
genomad_prophage_clean %>%
  filter(if_any(everything(), ~ str_detect(., fixed("GCA_026073395_1-1"))))

#GCA_026073395_1-1|provirus_2569287_2630967

genomad_prophage_clean %>%
  filter(if_any(everything(), ~ str_detect(., fixed("GCA_006304565_1-1"))))

#GCA_006304565_1-1|provirus_3098378_3120159

save.image("1_ESKAPE_MGE.RData")


#MGE trade-off ?####
head(mge_matrix)

mytheme <-   theme(
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
)

p_plasmid_prophage <- 
  ggplot(mge_matrix, aes(x = Prophage, y = Plasmid, color = Tax)) +
  geom_point(alpha = 0.6, size = 2) +
  #geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "dashed") +
  scale_color_manual(values = c("Ef"="#9f86c0","Sa"="#FFBE7D",
                                "Kp"="#76c893","Ab"="#F0E442","Pa"="#ff99c8",
                                "En"="#56B4E9")
  ) + mytheme
p_plasmid_prophage

#统计各种基因之间的关系
library(psych)

corr_test <- corr.test(mge_matrix[,c(-1,-7)], method = "spearman")
corr_test$r  # 相关系数矩阵
corr_test$p  # p值矩阵
library(corrplot)
corrplot(corr_test$r, 
         method = "color",
         type = "upper",
         tl.col = "black",
         addCoef.col = "black",
         number.cex = 0.7)


df.count <- mge_matrix %>%
  mutate(Prophage = factor(Prophage, levels = 0:15)) %>%
  mutate(group = case_when(Prophage == 1 ~ "group1", TRUE ~ "group2"))
df.msd <- mge_matrix %>%
  group_by(Prophage) %>%
  summarise_at(vars(Plasmid), list(max = max, mean = mean, sd = sd, median = median, skew = moments::skewness)) %>% 
  as.data.frame() %>%
  mutate(Prophage = as.numeric(Prophage))
df.msd
# Visualize
library(scales)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(cowplot)
set.seed(12345)
p_plasmid_prophage <- ggplot() +
  geom_violin(data = df.count, aes(x = Prophage, y = Plasmid, fill = Prophage),
              position = position_dodge(0.5), linewidth = 0.5,drop = F) +
  geom_jitter(data = df.count, aes(x = Prophage, y = Plasmid), size = 0.01, color = "grey") +
  # scale_x_discrete(breaks = 0:25, labels = c(0:25, ">=25")) +
  scale_y_continuous(trans = log2_trans(), breaks = c(1,2,4,8,16,32,64), labels = c(1,2,4,8,16,32,"≥64")) +
  scale_fill_viridis_d(option = 'C')+
  # scale_fill_manual(values = c("1" = "#FFEFCC", "2" = "#FECCCC", "3" = "#99CCCC", "4" = "#CCCCFF", "5" = "#3399CC"))+
  labs(title = "", x = "Prophage count", y = "Plasmid count", fill = "") +
  theme_bw() +
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0), "cm"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = 0.5),
        panel.grid = element_blank(),
        #  axis.line.x = element_line(linewidth = 0.5),
        #  axis.line.y = element_line(linewidth = 0.5),
        axis.title.x = element_text(size = 12, face = "plain"),
        axis.title.y = element_text(size = 12, face = "plain"),
        axis.ticks.x = element_line(linewidth = 0.5),
        axis.ticks.y = element_line(linewidth = 0.5),
        axis.ticks.length.x = unit(0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        axis.text.x = element_text(angle = 90,size = 10,hjust = 0.5,vjust = 0.5),
        axis.text.y = element_text(size = 10),
        legend.position = "none")
p_plasmid_prophage


p_plasmid_prophage_2 <- ggscatter(df.msd, x = "Prophage", y = "mean", add = "reg.line", add.params = list(color = "blue", fill = NA, linetype = "dashed")) +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
           label.x = 2, label.y = 8) +
  stat_regline_equation(label.x = 2, label.y = 6) +
  labs(title = "", x = "Prophage count", y = "mean Plasmid count", fill = "") +
  theme_bw() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = 0.5),
        panel.grid = element_blank(),
        # axis.line.x = element_line(linewidth = 0.5),
        # axis.line.y = element_line(linewidth = 0.5),
        axis.title.x = element_text(size = 12, face = "plain"),
        axis.title.y = element_text(size = 12, face = "plain"),
        axis.ticks.x = element_line(linewidth = 0.5),
        axis.ticks.y = element_line(linewidth = 0.5),
        axis.ticks.length.x = unit(0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "none")
p_plasmid_prophage_2

save.image("1_ESKAPE_MGE.RData")

load("1_ESKAPE_MGE.RData")

