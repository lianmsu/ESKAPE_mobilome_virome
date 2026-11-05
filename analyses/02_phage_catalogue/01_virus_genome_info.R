# ===========================================================================
# ESKAPE_mobilome_virome -- analyses / 02_phage_catalogue
# Quality, taxonomy, host and lifestyle of the retained viruses
# Manuscript      : Fig. 3A-E
# Original path   : 2ESKAPE_MGE/2virus_genome_info/1_virus_genome_info_o.R
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


setwd(file.path(CS_DATA, "2virus_genome_info"))
library(tidyr)
library(dplyr)
library(stringr)

#病毒宿主注释####

ESKAPE_all_phage_info <- read.delim2("phage_detailed_info_trim.tsv")

VHDB_info <- read.delim2("VHDB_ESKAPE_virus.tsv")

VHDB_info[c(25:30),]
head(VHDB_info, n = 30)


# 拆分 refseq.id 列
VHDB_expanded <- VHDB_info %>%
  separate_rows(refseq.id, sep = ",\\s*") %>%
  mutate(refseq.id = trimws(refseq.id))  # 去除空白字符

# # 最基础的方法
# write.table(VHDB_expanded, 
#             "VHDB_expanded.tsv",
#             sep = "\t",           # TSV使用制表符分隔
#             row.names = FALSE,    # 不写入行名
#             col.names = TRUE,     # 写入列名
#             quote = FALSE        # 不添加引号
#             )            # NA值的表示方式

VHDB_expanded <- read.delim2("VHDB_expanded.tsv")



ESKAPE_all_phage_info <- ESKAPE_all_phage_info %>%
  mutate(
    Host_genome = str_remove(seq_name, "\\|prophage.*")
  )

# 方法1B：只添加Host_tax列，保持其他列不变
ESKAPE_all_phage_info <- ESKAPE_all_phage_info %>%
  mutate(Host_tax = VHDB_expanded$Host[match(Host_genome, VHDB_expanded$refseq.id)])

ESKAPE_all_phage_info <- ESKAPE_all_phage_info %>%
  select(-Host_genome)
#
library(stringr)
ESKAPE_all_phage_info <- ESKAPE_all_phage_info %>%
  mutate(
    Host_genome = ifelse(
      str_starts(seq_name, "GCA"),
      str_remove(seq_name, "-.*"),
      seq_name
    )
  )

ESKAPE_genome_info <- read.delim2("ESKAPE_genome_info.tsv")
ESKAPE_genome_info$genome <- gsub("\\.", "_", ESKAPE_genome_info$genome)

# 找到ESKAPE_all_phage_info$Host_genome在ESKAPE_genome_info$genome中的位置
matches <- match(ESKAPE_all_phage_info$Host_genome, ESKAPE_genome_info$genome)
# 使用match结果更新Host_tax
ESKAPE_all_phage_info <- ESKAPE_all_phage_info %>%
  mutate(
    Host_tax = ifelse(
      !is.na(matches),  # 如果匹配成功
      ESKAPE_genome_info$Tax[matches],  # 使用匹配的Tax值
      Host_tax  # 否则保持原样
    )
  )


ESKAPE_all_phage_info <- ESKAPE_all_phage_info %>%
  select(-Host_genome)

ESKAPE_all_phage_info <- ESKAPE_all_phage_info %>%
  mutate(
    seq_name_tmp = 
      str_remove(seq_name, "\\|prophage.*")
    )

#IMGVR处理
IMGVR_spacer_info <- read.delim2("spacers_IMGVR.blastn.1mis100cov.out2",
                                 header = F)
colnames(IMGVR_spacer_info) <- c("spacer","seq_id")
IMGVR_spacer_info <- IMGVR_spacer_info %>%
  mutate(
    Host_genome = str_remove(spacer, "-.*")
  )
# 
matches <- match(IMGVR_spacer_info$Host_genome, ESKAPE_genome_info$genome)
# 
IMGVR_spacer_info <- IMGVR_spacer_info %>%
  mutate(
    Host_tax = ifelse(
      !is.na(matches),  # 如果匹配成功
      ESKAPE_genome_info$Tax[matches],  # 使用匹配的Tax值
      Host_tax  # 否则保持原样
    )
  )
IMGVR_spacer_cleaned <- IMGVR_spacer_info %>%
  group_by(seq_id) %>%
  summarise(
    seq_name = first(seq_id),  # 去冗余，每个seq_id只保留一个
    Host_tax = paste(unique(Host_tax), collapse = ","),  # 合并不同的Host_tax
    .groups = 'drop'
  ) %>%
  select(seq_name, Host_tax)  # 选择需要的列

IMGVR_spacer_cleaned <- IMGVR_spacer_cleaned %>%
  mutate(
    Host_tax = ifelse(
      str_detect(Host_tax, ","),
      sapply(str_split(Host_tax, ","), function(x) paste(sort(x), collapse = ",")),
      Host_tax
    )
  )

#VIRE处理
VIRE_spacer_info <- read.delim2("spacers_VIRE.blastn.1mis100cov.out2",
                                header = F)
colnames(VIRE_spacer_info) <- c("spacer","seq_id")
VIRE_spacer_info <- VIRE_spacer_info %>%
  mutate(
    Host_genome = str_remove(spacer, "-.*")
  )
# 
matches <- match(VIRE_spacer_info$Host_genome, ESKAPE_genome_info$genome)
# 
VIRE_spacer_info <- VIRE_spacer_info %>%
  mutate(
    Host_tax = ifelse(
      !is.na(matches),  # 如果匹配成功
      ESKAPE_genome_info$Tax[matches],  # 使用匹配的Tax值
      Host_tax  # 否则保持原样
    )
  )
VIRE_spacer_cleaned <- VIRE_spacer_info %>%
  group_by(seq_id) %>%
  summarise(
    seq_name = first(seq_id),  # 去冗余，每个seq_id只保留一个
    Host_tax = paste(unique(Host_tax), collapse = ","),  # 合并不同的Host_tax
    .groups = 'drop'
  ) %>%
  select(seq_name, Host_tax)  # 选择需要的列

VIRE_spacer_cleaned <- VIRE_spacer_cleaned %>%
  mutate(
    Host_tax = ifelse(
      str_detect(Host_tax, ","),
      sapply(str_split(Host_tax, ","), function(x) paste(sort(x), collapse = ",")),
      Host_tax
    )
  )


spacer2virus_tax <- rbind(IMGVR_spacer_cleaned,VIRE_spacer_cleaned)

#合并到病毒metadata中

# 
matches <- match(ESKAPE_all_phage_info$seq_name_tmp, spacer2virus_tax$seq_name)
# 使用match结果更新Host_tax
ESKAPE_all_phage_info <- ESKAPE_all_phage_info %>%
  mutate(
    Host_tax = ifelse(
      !is.na(matches),  # 如果匹配成功
      spacer2virus_tax$Host_tax[matches],  # 使用匹配的Tax值
      Host_tax  # 否则保持原样
    )
  )


##iPHoP过滤结果#####
iphop_res_phage_id <- read.delim2("ESKAPE_phage.id",header = F)
colnames(iphop_res_phage_id) <- "Phage_id"

iphop_res_detail_info <- read.delim2("iphop_phage2hostgenome.tsv")


ESKAPE_all_phage_info_other <- 
  ESKAPE_all_phage_info[!ESKAPE_all_phage_info$ID %in% iphop_res_phage_id$Phage_id,]


ESKAPE_all_phage_info <- ESKAPE_all_phage_info %>%
  mutate(contain_method = case_when(
    # 条件1: ID列以GCA开头
    grepl("^GCA", ID) ~ "Prophage",
    
    # 条件2: ID在iphop_res_detail_info$Virus中
    ID %in% iphop_res_detail_info$Virus ~ 
      # 查找对应的Main.method
      iphop_res_detail_info$Main.method[match(ID, iphop_res_detail_info$Virus)],
    
    # 条件3: topology列为DTR或ITR
    topology %in% c("DTR", "ITR") ~ topology,
    
    # 默认情况
    TRUE ~ "only_my_crispr"
  ))

ESKAPE_all_phage_info_filter <- 
  ESKAPE_all_phage_info[ESKAPE_all_phage_info$contain_method != "only_my_crispr",]

write.table(ESKAPE_all_phage_info_filter,
            "ESKAPE_all_phage_info_filter.tsv",
            sep = "\t",           # TSV使用制表符分隔
            row.names = FALSE,    # 不写入行名
            col.names = TRUE,     # 写入列名
            quote = FALSE        # 不添加引号
          )            # NA值的表示方式

ESKAPE_all_phage_info_filter <- read.delim2("ESKAPE_all_phage_info_filter.tsv")


##1. 物种桑基图 ###########################
table(ESKAPE_all_phage_info_filter$Host_tax)

# 拆分 classification 列
ESKAPE_all_phage_info_tax <- ESKAPE_all_phage_info_filter %>%
  separate(
    taxonomy,
    into = c("Virus", "Realm", "Kingdom", "Phylum", "Class", "Order", "Family"),
    sep = ";",
    remove = FALSE  # 如果保留原始列，设为 FALSE；否则设为 TRUE
  )

library(tidyverse)
library(ggsankeyfier) 
library(MetBrewer)

df <- ESKAPE_all_phage_info_tax

# 统计 k 列中每个k的计数
realm_counts <- df %>%
  count(Realm)

# 统计 phylum 列中每个门的计数
p_counts <- df %>%
  count(Phylum)
# 统计 class 列中每个门的计数
c_counts <- df %>%
  count(Class)

# 统计 order 列中每个门的计数
o_counts <- df %>%
  count(Order)

# 统计 family 列中每个门的计数
f_counts <- df %>%
  count(Family)

#df <- df[,c(12:18)]
df <- df[df$Virus != "Unclassified",]
df <- df[,c(27:33)]

library(tidyverse)
library(ggsankeyfier) 
library(MetBrewer)
df1 <- df %>% select(1,2) %>% group_by(Virus,Realm) %>% count() %>%
  pivot_stages_longer(.,stages_from = c("Virus", "Realm"),
                      values_from = "n")
df2 <- df %>% select(2,4) %>% group_by(Realm,Phylum) %>% count() %>%
  pivot_stages_longer(.,stages_from = c("Realm", "Phylum"),
                      values_from = "n")
df3 <- df %>% select(4,5) %>% group_by(Phylum,Class) %>% count() %>%
  pivot_stages_longer(.,stages_from = c("Phylum", "Class"),
                      values_from = "n")

df4 <- df %>% select(5,6) %>% group_by(Class,Order) %>% count() %>%
  pivot_stages_longer(.,stages_from = c("Class", "Order"),
                      values_from = "n")
df5 <- df %>% select(6,7) %>% group_by(Order,Family) %>% count() %>%
  pivot_stages_longer(.,stages_from = c("Order", "Family"),
                      values_from = "n")

p_phage_sankey_4level <- 
  ggplot(data=df1,aes(x = stage,y =n,group = node,
                      edge_id = edge_id,connector = connector))+
  # 绘制第 1，2 层级
  geom_sankeyedge(aes(fill = node),
                  position = position_sankey(order ="ascending",v_space="auto",
                                             width = 0.01))+
  geom_sankeynode(aes(fill=node,color=node),
                  position = position_sankey(order = "ascending",v_space ="auto",
                                             width = 0.05))+
  geom_text(data=df1 %>% filter(connector=="from"),
            aes(label = node),stat = "sankeynode",
            position = position_sankey(v_space ="auto",order="ascending",nudge_x=0.05),
            hjust=0,size=3.5,fontface="plain",color="black")+
  # 绘制第 2，3 层级
  geom_sankeyedge(data=df2,aes(fill = node),
                  position = position_sankey(order = "ascending",v_space ="auto",
                                             width = 0.01))+
  geom_sankeynode(data=df2,aes(fill=node,color=node),
                  position = position_sankey(order = "ascending",v_space ="auto",
                                             width = 0.05))+
  geom_text(data=df1 %>% filter(connector=="to"),
            aes(label = node),stat = "sankeynode",angle=0,
            position = position_sankey(v_space ="auto",order="ascending", nudge_x=0.05),
            hjust=0,size=3.5,vjust=0.5,color="black",fontface="plain")+
  geom_text(data=df2 %>% filter(connector=="to"),
            aes(label = node),stat = "sankeynode",angle=0,
            position = position_sankey(v_space ="auto",order="ascending",nudge_x=0.05),
            hjust=0,size=3.5,color="black",fontface="plain")+
  # 绘制第 3，4 层级
  geom_sankeyedge(data=df3,aes(fill = node),
                  position = position_sankey(order = "ascending",v_space ="auto",
                                             width = 0.01))+
  geom_sankeynode(data=df3,aes(fill=node,color=node),
                  position = position_sankey(order = "ascending",v_space ="auto",
                                             width = 0.05))+
  geom_text(data=df3 %>% filter(connector=="to"),
            aes(label = node),stat = "sankeynode",angle=0,
            position = position_sankey(v_space ="auto",order="ascending", nudge_x=0.05),
            hjust=0,size=3.5,vjust=0.5,color="black",fontface="plain")+
  # 绘制第 4，5 层级
  # geom_sankeyedge(data=df4,aes(fill = node),
  #                 position = position_sankey(order = "ascending",v_space ="auto",
  #                                            width = 0.01))+
  # geom_sankeynode(data=df4,aes(fill=node,color=node),
  #                 position = position_sankey(order = "ascending",v_space ="auto",
  #                                            width = 0.05))+
  # geom_text(data=df4 %>% filter(connector=="to"),
  #           aes(label = node),stat = "sankeynode",angle=0,
  #           position = position_sankey(v_space ="auto",order="ascending", nudge_x=0.05),
  #           hjust=0,size=3,vjust=0.5,color="black",fontface="plain")+
# 
coord_cartesian(clip="off")+
  scale_x_discrete(position = "bottom")+
  # scale_fill_manual(values = met.brewer("Nizami")) +
  # scale_color_manual(values = met.brewer("Nizami")) +
  theme_void()+
  theme(plot.margin = margin(0,0,0,0,unit = "in"),
        axis.text.x=element_text(color="black",face="plain",size=10,
                                 margin = margin(b=0.1,unit = "cm")),
        legend.position="none")

p_phage_sankey_4level



##2. 不同宿主的分布##################################################
library(ggplot2)
library(dplyr)

# 准备数据：统计每个Host_tax中不同Class的数量
ESKAPE_all_phage_info_tax <- ESKAPE_all_phage_info_tax %>%
  mutate(
    Class = ifelse(Class %in% c("Caudoviricetes", "Faserviricetes"), Class, "Others")
  ) 

plot_data <- ESKAPE_all_phage_info_tax %>%
  group_by(Host_tax, Class) %>%
  summarise(count = n(), .groups = 'drop') %>%
  ungroup()

plot_data[is.na(plot_data)] <- "Unclassified"

plot_data1 <- plot_data[!plot_data$Host_tax %in% c("Ab,Pa","Ab,Kp","Ab,En,Kp","Ab,En"),]

plot_data1$Host_tax <- factor(plot_data1$Host_tax, 
                             levels = c("Ef","Sa","Kp","Ab","Pa","En","En,Kp"))

library(ggplot2)
# 绘制堆叠柱状图
p_major_phage <- 
ggplot(plot_data1, aes(x = Host_tax, y = count, fill = Class)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  # 添加数值标签
  geom_text(aes(label = count), 
            position = position_stack(vjust = 0.5),
            size = 3, color = "black", fontface = "bold") +
  # 配色方案
  #scale_fill_brewer(palette = "Set3") +
  scale_fill_manual(values = c("Ef"="#9f86c0","Faserviricetes"="#A58AFF",
                               "Kp"="#76c893","Ab"="#F0E442","Pa"="#ff99c8",
                               "Caudoviricetes"="#FB61D7")
  ) +
  # 美化主题
  theme(
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
p_major_phage

plot_data2 <- plot_data[plot_data$Host_tax %in% c("Ab,Pa","Ab,Kp","Ab,En,Kp","Ab,En"),]

p_other_phage <- 
  ggplot(plot_data2, aes(x = Host_tax, y = count, fill = Class)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  # 添加数值标签
  geom_text(aes(label = count), 
            position = position_stack(vjust = 0.5),
            size = 3, color = "black", fontface = "bold") +
  # 配色方案
  #scale_fill_brewer(palette = "Set3") +
  coord_flip() +  # 关键：将坐标轴翻转
  scale_fill_manual(values = c("Ef"="#9f86c0","Faserviricetes"="#A58AFF",
                               "Kp"="#76c893","Ab"="#F0E442","Pa"="#ff99c8",
                               "Caudoviricetes"="#FB61D7")
  ) +
  # 美化主题
  theme(
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
p_other_phage

##3. 病毒checkv序列的质量####

library(ggplot2)
library(ggbeeswarm)  # 用于更好的散点排列
# 计算每个类别的数量
quality_counts <- ESKAPE_all_phage_info_filter %>%
  count(checkv_quality, name = "count") %>%
  mutate(label = paste0(checkv_quality, " (n=", count, ")"))

quality_counts <- ESKAPE_all_phage_info_filter %>%
  count(checkv_quality, name = "count") %>%
  mutate(label = paste0(checkv_quality))


# 合并计数回原始数据
ESKAPE_all_phage_info_filter <- ESKAPE_all_phage_info_filter %>%
  left_join(quality_counts, by = "checkv_quality")

# 确保有序因子（包含计数标签）
ESKAPE_all_phage_info_filter$label <- factor(
  ESKAPE_all_phage_info_filter$label,
  levels = c("Complete (n=893)","High-quality (n=8252)",
             "Medium-quality (n=3001)","Low-quality (n=2350)"),
  ordered = TRUE
)

# 确保有序因子（包含计数标签）
ESKAPE_all_phage_info_filter$label <- factor(
  ESKAPE_all_phage_info_filter$label,
  levels = c("Complete","High-quality",
             "Medium-quality","Low-quality"),
  ordered = TRUE
)

# # 确保checkv_quality是有序因子
# checkv_phage_quality$checkv_quality <- factor(
#   checkv_phage_quality$checkv_quality,
#   levels = c("Not-determined", "Low-quality", "Medium-quality", "High-quality", "Complete"),
#   ordered = TRUE
# )
# 创建颜色映射

# 绘制横向箱线图
p_phage_quality <- 
  ggplot(ESKAPE_all_phage_info_filter, 
         aes(x = contig_length, 
             y = label, 
             fill = checkv_quality,
             color = checkv_quality)) +
  scale_y_discrete(limits = rev)+
  # 添加箱线图
  geom_boxplot(alpha = 0.7, 
               outlier.shape = 1,          # 实心圆点
               outlier.size = 1,            # 离群点大小
               outlier.alpha = 0.7,         # 离群点透明度
               width = 0.6,
               size = 0.4                   # 边框线粗细
  ) +
  # 设置颜色和填充（使用相同颜色方案）
  # scale_fill_manual(values = quality_colors) +
  # scale_color_manual(values = quality_colors) +
  # 添加散点表示所有数据点
  # geom_quasirandom(
  #   aes(color = checkv_quality),  # 使用颜色区分
  #   alpha = 0.5,                 # 半透明
  #   size = 1.5,                  # 点大小
  #   groupOnX = FALSE             # 在y轴方向排列
  # ) +
  
# 使用viridis颜色方案
scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.8) +
  scale_color_viridis_d(option = "D", begin = 0.2, end = 0.8) +
  
  # 对数变换x轴
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  # 
  # 添加标签和主题
  labs(
    title = "Contig length distribution by CheckV quality",
    x = "Contig length (bp, log10 scale)",
    y = "CheckV quality level",
    #caption = "Each point represents one contig"
  ) +
  
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "plain"),  # 标题居中，加粗，字体大小为16
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),  # 设置图形边距
    plot.background = element_blank(),  # 去除图形背景
    panel.background = element_blank(),  # 去除面板背景
    panel.border = element_rect(linewidth = 0.6),  # 设置面板边框
    panel.grid = element_blank(),  # 去除网格线
    axis.title.x = element_text(size = 10, face = "plain",color = 'black'),  # 设置 x 轴标题
    axis.title.y = element_text(size = 10, face = "plain",color = 'black'),  # 设置 y 轴标题
    axis.ticks.x = element_line(linewidth = 0.5),  # 设置 x 轴刻度线
    axis.ticks.y = element_line(linewidth = 0.5),  # 设置 y 轴刻度线
    axis.ticks.length.x = unit(0.1, "cm"),  # 设置 x 轴刻度长度
    axis.ticks.length.y = unit(0.1, "cm"),  # 设置 y 轴刻度长度
    axis.text.x = element_text(angle = 0, size = 10, hjust = 0.5, vjust = 0.5,color = 'black'),  # 设置 x 轴文本
    axis.text.y = element_text(size = 10, color = 'black'),  # 设置 y 轴文本
    legend.position = "none"  # 去除图例
  )
p_phage_quality

ggsave(file.path(CS_WORK, "1_p_phage_quality2.pdf"),
       p_phage_quality ,width = 3,height = 3,units = 'in')


ggsave(file.path(CS_WORK, "1_p_phage_sankey_4level.pdf"),
       p_phage_sankey_4level ,width = 3,height = 3,units = 'in')

ggsave(file.path(CS_WORK, "1_p_major_phage.pdf"),
       p_major_phage ,width = 3,height = 3,units = 'in')

ggsave(file.path(CS_WORK, "1_p_other_phage.pdf"),
       p_other_phage ,width = 1.5,height = 1.2,units = 'in')

library(Cairo)
library(cowplot)
#cairo_pdf(filename = file.path(CS_WORK, "1_genome_metadata.pdf") , width = 8,height = 4,unit = 'in')
g1 <- plot_grid(p_phage_quality, p_phage_sankey_4level,p_major_phage,
                labels = LETTERS[1:3], ncol = 3,
                rel_widths = c(1, 1,1)
                # label_x = 0.09,            # 标签x坐标（1为最右侧）
                #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g1


#ggsave(file.path(CS_WORK, "1_virus_genome_info.pdf"),g1,width = 8,height = 3.2,units = 'in')

p_other_phage
#ggsave(file.path(CS_WORK, "1_p_other_phage.pdf"),p_other_phage,width = 2.24,height = 2.74,units = 'in')

##长度分布密度图####

library(ggplot2)
library(dplyr)

# 计算统计信息
stats <- ESKAPE_all_phage_info_tax %>%
  summarise(
    mean_len = mean(length, na.rm = TRUE),
    median_len = median(length, na.rm = TRUE),
    min_len = min(length, na.rm = TRUE),
    max_len = max(length, na.rm = TRUE),
    n = n()
  )

p_phage_length_density <-
ggplot(ESKAPE_all_phage_info_tax, aes(x = length)) +
  geom_density(fill = "#4DBBD5", alpha = 0.6, color = "black", linewidth = 0.6) +
  # 添加均值和分位数线
  geom_vline(aes(xintercept = mean_len), data = stats, 
             color = "red", linewidth = 0.6, linetype = "dashed") +
  geom_vline(aes(xintercept = median_len), data = stats,
             color = "blue", linewidth = 0.6, linetype = "dotted") +
  # 添加统计信息文本
  annotate("text", x = Inf, y = Inf, 
           label = paste0("N = ", stats$n, "\n",
                          "Mean = ", round(stats$mean_len, 0), "\n",
                          "Median = ", stats$median_len),
           hjust = 1.1, vjust = 1.1, size = 4) +
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
  )
p_phage_length_density

ggsave(file.path(CS_WORK, "1_p_phage_length_density.pdf"),p_phage_length_density,
       width = 4.85,height = 3.8,units = 'in')


#complete phage############################

complete_phage <- ESKAPE_all_phage_info_tax[ESKAPE_all_phage_info_tax$topology %in% c("DTR","ITR"),]

table(complete_phage$Host_tax)


save.image("1_virus_genome_info.RData")

load("1_virus_genome_info.RData")










