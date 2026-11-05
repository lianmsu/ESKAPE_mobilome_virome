# ===========================================================================
# ESKAPE_mobilome_virome -- analyses / 04_mobilome
# CRISPR spacer targeting of the different MGE types
# Manuscript      : Fig. 2B
# Original path   : 2ESKAPE_MGE/4_1MGE_spacer_competition/1_MGE_competition.R
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

setwd(file.path(CS_DATA, "4_1MGE_spacer_competition"))

ESKAPE_MGE_num <- read.delim2("ESKAPE_genome_info2.tsv")

ESKAPE_DS_num <- read.delim2("ESKAPE_genome_info3.tsv")

library(dplyr)
library(ggplot2)
ESKAPE_MGE_num <- ESKAPE_MGE_num %>%
  left_join(ESKAPE_DS_num[,c(1, 19:26)], by = c("genome"))

MGE_spacer_competition <- readxl::read_excel("MGE_spacer_competition_res.xlsx")


# 先计算百分比
plasmid_propotion <- 
  MGE_spacer_competition[MGE_spacer_competition$Source == "Plasmid",]
plasmid_propotion <- plasmid_propotion %>%
  mutate(
    percentage = number / sum(number) * 100,
    label = paste0(number, "\n(", round(percentage, 1), "%)")
  )
# 绘制环形图
p_plasmid_propotion <- 
  ggplot(plasmid_propotion, aes(x = 2, y = number, fill = target)) +
  geom_bar(stat = "identity", width = 0.8, color = "white") +
  coord_polar(theta = "y", start = 0) +
  xlim(0.5, 2.5) +  # 控制环的厚度（0.5-2.5之间）
  scale_fill_manual(
    values = c("Plasmid" = "#9A01FD", "Prophage" = "#FA482E",
               "Integron" = "#FCBE6E", "ICE" = "#93B4E9"),
    name = ""
  ) +
  # 在扇形中心添加数值标签
  # 在扇形中心添加数值标签
  geom_text( aes(label = label), 
             position = position_stack(vjust = 0.5),
             color = "white", # 文本颜色
             size = 4) + # 文本大小
  # 在中心添加总计
  # annotate("text", x = 0.5, y = 0, 
  #          label = paste("Totol\n", sum(plasmid_propotion$number)),
  #          size = 5, fontface = "bold", color = "gray30") +
  theme_void() +  # 移除背景元素
  theme(legend.position = "none",  # 移除默认图例
        axis.text.x = element_blank(),
        axis.ticks = element_blank())  # 移除坐标轴标签和刻度
p_plasmid_propotion

# 先计算百分比
prophage_propotion <- 
  MGE_spacer_competition[MGE_spacer_competition$Source == "Prophage",]
prophage_propotion <- prophage_propotion %>%
  mutate(
    percentage = number / sum(number) * 100,
    label = paste0(number, "\n(", round(percentage, 1), "%)")
  )
# 绘制环形图
p_prophage_propotion <- 
  ggplot(prophage_propotion, aes(x = 2, y = number, fill = target)) +
  geom_bar(stat = "identity", width = 0.8, color = "white") +
  coord_polar(theta = "y", start = 0) +
  xlim(0.5, 2.5) +  # 控制环的厚度（0.5-2.5之间）
  scale_fill_manual(
    values = c("Plasmid" = "#9A01FD", "Prophage" = "#FA482E",
               "Integron" = "#FCBE6E", "ICE" = "#93B4E9"),
    name = ""
  ) +
  # 在扇形中心添加数值标签
  # 在扇形中心添加数值标签
  geom_text( aes(label = label), 
             position = position_stack(vjust = 0.5),
             color = "white", # 文本颜色
             size = 4) + # 文本大小
  # 在中心添加总计
  # annotate("text", x = 0.5, y = 0, 
  #          label = paste("Totol\n", sum(plasmid_propotion$number)),
  #          size = 5, fontface = "bold", color = "gray30") +
  theme_void() +  # 移除背景元素
  theme(legend.position = "none",  # 移除默认图例
        axis.text.x = element_blank(),
        axis.ticks = element_blank())  # 移除坐标轴标签和刻度
p_prophage_propotion

# 先计算百分比
integron_propotion <- 
  MGE_spacer_competition[MGE_spacer_competition$Source == "Integron",]
integron_propotion <- integron_propotion %>%
  mutate(
    percentage = number / sum(number) * 100,
    label = paste0(number, "\n(", round(percentage, 1), "%)")
  )
# 绘制环形图
p_integron_propotion <- 
  ggplot(integron_propotion, aes(x = 2, y = number, fill = target)) +
  geom_bar(stat = "identity", width = 0.8, color = "white") +
  coord_polar(theta = "y", start = 0) +
  xlim(0.5, 2.5) +  # 控制环的厚度（0.5-2.5之间）
  scale_fill_manual(
    values = c("Plasmid" = "#9A01FD", "Prophage" = "#FA482E",
               "Integron" = "#FCBE6E", "ICE" = "#93B4E9"),
    name = ""
  ) +
  # 在扇形中心添加数值标签
  # 在扇形中心添加数值标签
  geom_text( aes(label = label), 
             position = position_stack(vjust = 0.5),
             color = "white", # 文本颜色
             size = 4) + # 文本大小
  # 在中心添加总计
  # annotate("text", x = 0.5, y = 0, 
  #          label = paste("Totol\n", sum(plasmid_propotion$number)),
  #          size = 5, fontface = "bold", color = "gray30") +
  theme_void() +  # 移除背景元素
  theme(legend.position = "none",  # 移除默认图例
        axis.text.x = element_blank(),
        axis.ticks = element_blank())  # 移除坐标轴标签和刻度
p_integron_propotion

# 先计算百分比
ICE_propotion <- 
  MGE_spacer_competition[MGE_spacer_competition$Source == "ICE",]
ICE_propotion <- ICE_propotion %>%
  mutate(
    percentage = number / sum(number) * 100,
    label = paste0(number, "\n(", round(percentage, 1), "%)")
  )
# 绘制环形图
p_ICE_propotion <- 
  ggplot(ICE_propotion, aes(x = 2, y = number, fill = target)) +
  geom_bar(stat = "identity", width = 0.8, color = "white") +
  coord_polar(theta = "y", start = 0) +
  xlim(0.5, 2.5) +  # 控制环的厚度（0.5-2.5之间）
  scale_fill_manual(
    values = c("Plasmid" = "#9A01FD", "Prophage" = "#FA482E",
               "Integron" = "#FCBE6E", "ICE" = "#93B4E9"),
    name = ""
  ) +
  # 在扇形中心添加数值标签
  # 在扇形中心添加数值标签
  geom_text( aes(label = label), 
             position = position_stack(vjust = 0.5),
             color = "white", # 文本颜色
             size = 4) + # 文本大小
  # 在中心添加总计
  # annotate("text", x = 0.5, y = 0, 
  #          label = paste("Totol\n", sum(plasmid_propotion$number)),
  #          size = 5, fontface = "bold", color = "gray30") +
  theme_void() +  # 移除背景元素
  theme(legend.position = "none",  # 移除默认图例
        axis.text.x = element_blank(),
        axis.ticks = element_blank())  # 移除坐标轴标签和刻度
p_ICE_propotion

library(Cairo)
library(cowplot)

g_MGE_spacer_competition <- plot_grid(p_prophage_propotion,
                                      p_plasmid_propotion,
                                      p_ICE_propotion,
                                      p_integron_propotion,
                 ncol = 2, nrow = 2
               # labels = LETTERS[1:2],
               # rel_widths = c(1, 2)
                # label_x = 0.09,            # 标签x坐标（1为最右侧）
                #  label_y = 0.9             # 标签y坐标（1为最上方）
                )
g_MGE_spacer_competition

ggsave(plot = g_MGE_spacer_competition,'1_g_MGE_spacer_competition2.pdf',
       width = 2.7,height = 2.7,units = 'in')

save.image("./1_MGE_competition.RData")

