# ===========================================================================
# ESKAPE_mobilome_virome -- analyses / 03_gene_annotation
# CRISPR-Cas subtype distribution across the ESKAPE lineages
# Manuscript      : Fig. 1C
# Original path   : 2ESKAPE_MGE/3genome_gene_annotation/1CasTyper/1_CRISPRCas_tax.R
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


setwd(file.path(CS_DATA, "3genome_gene_annotation/1CasTyper"))

CRISPR_Cas <- read.delim2("CRISPR_Cas.tab")

library(stringr)
CRISPR_Cas$genome <- str_extract(CRISPR_Cas$Contig, "^[^-]*")

all_4538genome_metadata <- read.delim2("ESKAPE_genome_info.tsv")
all_4538genome_metadata$genome <- gsub("\\.","_",all_4538genome_metadata$genome)
CRISPR_Cas_genome <- all_4538genome_metadata[all_4538genome_metadata$genome %in% CRISPR_Cas$genome,]

# 按照CRISPR_Cas$genome的顺序重新排序
CRISPR_Cas_genome2 <- CRISPR_Cas_genome[
  match(CRISPR_Cas$genome, CRISPR_Cas_genome$genome),
]

CRISPR_Cas_genome2 <- merge(CRISPR_Cas_genome, CRISPR_Cas, by.x = "genome",by.y = "genome")

#I-F挑选
table(CRISPR_Cas$Prediction)

I_A <- CRISPR_Cas_genome2[CRISPR_Cas_genome2$Prediction=="I-A",]
I_B <- CRISPR_Cas_genome2[CRISPR_Cas_genome2$Prediction=="I-B",]
I_C <- CRISPR_Cas_genome2[CRISPR_Cas_genome2$Prediction=="I-C",]
I_D <- CRISPR_Cas_genome2[CRISPR_Cas_genome2$Prediction=="I-D",]
I_E <- CRISPR_Cas_genome2[CRISPR_Cas_genome2$Prediction=="I-E",]
I_F <- CRISPR_Cas_genome2[CRISPR_Cas_genome2$Prediction=="I-F",]
I_G <- CRISPR_Cas_genome2[CRISPR_Cas_genome2$Prediction=="I-G",]

II_A <- CRISPR_Cas_genome2[CRISPR_Cas_genome2$Prediction=="II-A",]
II_B <- CRISPR_Cas_genome2[CRISPR_Cas_genome2$Prediction=="II-B",]
II_C <- CRISPR_Cas_genome2[CRISPR_Cas_genome2$Prediction=="II-C",]

III_A <- CRISPR_Cas_genome2[CRISPR_Cas_genome2$Prediction=="III-A",]
III_B <- CRISPR_Cas_genome2[CRISPR_Cas_genome2$Prediction=="III-B",]
III_C <- CRISPR_Cas_genome2[CRISPR_Cas_genome2$Prediction=="III-C",]
III_D <- CRISPR_Cas_genome2[CRISPR_Cas_genome2$Prediction=="III-D",]

V_A <- CRISPR_Cas_genome2[CRISPR_Cas_genome2$Prediction=="V-A",]
V_B1 <- CRISPR_Cas_genome2[CRISPR_Cas_genome2$Prediction=="V-B1",]
V_M <- CRISPR_Cas_genome2[CRISPR_Cas_genome2$Prediction=="V-M",]

IV_A3 <- CRISPR_Cas_genome2[CRISPR_Cas_genome2$Prediction=="IV-A3",]

Hybrid <- CRISPR_Cas_genome2[grepl("Hybrid", CRISPR_Cas_genome2$Prediction), ]

library(dplyr)
IV_A3 <- IV_A3 %>%
  mutate(genomic_context = case_when(
    grepl("plasmid", Contig, ignore.case = TRUE) ~ "plasmid",
    TRUE ~ "chromosome"  # 默认值
  ))

IV_A3_propotion <- table(IV_A3$genomic_context) %>% as.data.frame()

# 先计算百分比
IV_A3_propotion <- IV_A3_propotion %>%
  mutate(
    percentage = Freq / sum(Freq) * 100,
    label = paste0(Freq, "\n(", round(percentage, 1), "%)")
  )

library(ggplot2)
# 绘制环形图
p_IV_A3 <- 
ggplot(IV_A3_propotion, aes(x = 2, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", width = 0.8, color = "white") +
  coord_polar(theta = "y", start = 0) +
  xlim(0.5, 2.5) +  # 控制环的厚度（0.5-2.5之间）
  scale_fill_manual(
    values = c("chromosome" = "#2E86AB", "plasmid" = "#A23B72"),
    name = ""
  ) +
  # 在扇形中心添加数值标签
  # 在扇形中心添加数值标签
  geom_text( aes(label = label), 
             position = position_stack(vjust = 0.5),
             color = "white", # 文本颜色
             size = 4) + # 文本大小
  # 在中心添加总计
  annotate("text", x = 0.5, y = 0, 
           label = paste("Totol\n", sum(IV_A3_propotion$Freq)),
           size = 5, fontface = "bold", color = "gray30") +
  theme_void() +  # 移除背景元素
  theme(legend.position = "none",  # 移除默认图例
        axis.text.x = element_blank(),
        axis.ticks = element_blank())  # 移除坐标轴标签和刻度
p_IV_A3

####################
I_E <- I_E %>%
  mutate(genomic_context = case_when(
    grepl("plasmid", Contig, ignore.case = TRUE) ~ "plasmid",
    TRUE ~ "chromosome"  # 默认值
  ))

I_E_propotion <- table(I_E$genomic_context) %>% as.data.frame()

# 先计算百分比
I_E_propotion <- I_E_propotion %>%
  mutate(
    percentage = Freq / sum(Freq) * 100,
    label = paste0(Freq, "\n(", round(percentage, 1), "%)")
  )

library(ggplot2)
# 绘制环形图
p_I_E <- 
  ggplot(I_E_propotion, aes(x = 2, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", width = 0.8, color = "white") +
  coord_polar(theta = "y", start = 0) +
  xlim(0.5, 2.5) +  # 控制环的厚度（0.5-2.5之间）
  scale_fill_manual(
    values = c("chromosome" = "#2E86AB", "plasmid" = "#A23B72"),
    name = ""
  ) +
  # 在扇形中心添加数值标签
  # 在扇形中心添加数值标签
  geom_text( aes(label = label), 
             position = position_stack(vjust = 0.5),
             color = "white", # 文本颜色
             size = 4) + # 文本大小
  # 在中心添加总计
  annotate("text", x = 0.5, y = 0, 
           label = paste("Totol\n", sum(I_E_propotion$Freq)),
           size = 5, fontface = "bold", color = "gray30") +
  theme_void() +  # 移除背景元素
  theme(legend.position = "none",  # 移除默认图例
        axis.text.x = element_blank(),
        axis.ticks = element_blank())  # 移除坐标轴标签和刻度
p_I_E

##################
I_F <- I_F %>%
  mutate(genomic_context = case_when(
    grepl("plasmid", Contig, ignore.case = TRUE) ~ "plasmid",
    TRUE ~ "chromosome"  # 默认值
  ))

I_F_propotion <- table(I_F$genomic_context) %>% as.data.frame()

# 先计算百分比
I_F_propotion <- I_F_propotion %>%
  mutate(
    percentage = Freq / sum(Freq) * 100,
    label = paste0(Freq, "\n(", round(percentage, 1), "%)")
  )

library(ggplot2)
# 绘制环形图
p_I_F <- 
  ggplot(I_F_propotion, aes(x = 2, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", width = 0.8, color = "white") +
  coord_polar(theta = "y", start = 0) +
  xlim(0.5, 2.5) +  # 控制环的厚度（0.5-2.5之间）
  scale_fill_manual(
    values = c("chromosome" = "#2E86AB", "plasmid" = "#A23B72"),
    name = ""
  ) +
  # 在扇形中心添加数值标签
  geom_text( aes(label = label), 
             position = position_stack(vjust = 0.5),
             color = "white", # 文本颜色
             size = 4) + # 文本大小
  # 在中心添加总计
  annotate("text", x = 0.5, y = 0, 
           label = paste("Totol\n", sum(I_F_propotion$Freq)),
           size = 5, fontface = "bold", color = "gray30") +
  theme_void() +  # 移除背景元素
  theme(legend.position = "none",  # 移除默认图例
        axis.text.x = element_blank(),
        axis.ticks = element_blank())  # 移除坐标轴标签和刻度
  
p_I_F

ggsave(plot = p_IV_A3,'p_IV_A3.pdf',width = 3,height = 3,units = 'in')
ggsave(plot = p_I_E,'p_I_E.pdf',width = 3,height = 3,units = 'in')
ggsave(plot = p_I_F,'p_I_F.pdf',width = 3,height = 3,units = 'in')




