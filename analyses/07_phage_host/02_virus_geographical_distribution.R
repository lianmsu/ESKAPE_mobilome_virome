# ===========================================================================
# ESKAPE_mobilome_virome -- analyses / 07_phage_host
# Ecosystem and geographical distribution of the IMG/VR and VIRE viruses
# Manuscript      : Fig. 4A-B
# Original path   : 2ESKAPE_MGE/7_1partial_virus_source/1_phage_geographical_distribution.R
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


setwd(file.path(CS_DATA, "7_1partial_virus_source"))
library(dplyr)


ESKAPE_IMGVR_phage <- read.delim2("ESKAPE_IMGVR_phage2.tsv")

IMG_genomeCart_metadata <- readxl::read_excel("IMG_genomecart_metadata.xlsx")

ESKAPE_IMGVR_phage <- ESKAPE_IMGVR_phage %>%
  left_join(IMG_genomeCart_metadata %>% select(`IMG Genome ID`, 
                                               Ecosystem,`Ecosystem Category`, Latitude,
                                               Longitude),
            by = c("Taxon_oid" = "IMG Genome ID"))

table(ESKAPE_IMGVR_phage$Ecosystem)
table(ESKAPE_IMGVR_phage$`Ecosystem Category`)

ESKAPE_IMGVR_phage$`Ecosystem Category` <- ifelse(
  is.na(ESKAPE_IMGVR_phage$`Ecosystem Category`), 
  "Unknown", 
  ESKAPE_IMGVR_phage$`Ecosystem Category`
)


category_summary <- ESKAPE_IMGVR_phage %>%
  count(`Ecosystem Category`) %>%
  arrange(desc(n))
print(category_summary)

library(ggplot2)
library(dplyr)

nature_colors <- c("#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#79706E",
                   "#B07AA1", "#FFBE7D", "#8CD17D", "#86BCB6", "#499894", "#D37295")

# 水平条形图版本
p_IMGVR_phage_ecotype <- 
ggplot(category_summary, aes(x = n, y = reorder(`Ecosystem Category`, n), 
                             fill = `Ecosystem Category`)) +
  geom_bar(stat = "identity", width = 0.7, color = "black",linewidth = 0.1) +
  geom_text(aes(label = n), hjust = -0.2, size = 3.5, color = "black") +
  labs(
    y = "",
    x = "Number of ESKAPE phages identified from IMG/VR"
  ) +
  scale_fill_manual(values = c("Built environment"="#9b5de5","Terrestrial"="#FFBE7D",
                               "Plants"="#59A14F","Aquatic"="#00b4d8","Unknown"="#79706E",
                               "Host_associated"="#E15759","Mammals: Human"="#D37295"
                               )) +
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, linewidth = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_text( size = 12, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 12, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
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
p_IMGVR_phage_ecotype
##地理信息##########


self_selected_phage_4 <- 
  ESKAPE_IMGVR_phage[!ESKAPE_IMGVR_phage$Latitude %in% c("none","NA","n/a"),]
self_selected_phage_4 <- 
  self_selected_phage_4[!self_selected_phage_4$Longitude %in% c("none","NA","n/a"),]
self_selected_phage_4 <- self_selected_phage_4 %>%
  mutate(
    Latitude  = as.numeric(Latitude),
    Longitude = as.numeric(Longitude)
  )

redundancy_count <- self_selected_phage_4 %>%
  select(`Ecosystem Category`, Latitude, Longitude) %>%   # 只保留最后三列
  group_by(`Ecosystem Category`, Latitude, Longitude) %>% # 按三列分组
  summarise(virus_number = n(), .groups = "drop")             # 计数

#参考https://cloud.tencent.com/developer/article/2206435
library(ggplot2)
world<-map_data("world")
world.map <- 
ggplot() +
  geom_polygon(data=world,aes(x=long,y=lat,group=group),
               fill="#e9ecef")+
  geom_point(data=redundancy_count, 
             aes(x=Longitude,y=Latitude, 
                 fill = `Ecosystem Category`, size = virus_number),
             shape = 21,      # 带边框的填充圆
             colour = "black",# 边框颜色
             stroke = 0.05     # 边框宽度（mm）
             )+
  scale_size(range = c(2,6))+
  scale_fill_manual(values = c("Built environment"="#9b5de5","Terrestrial"="#FFBE7D",
                               "Plants"="#59A14F","Aquatic"="#00b4d8","Unknown"="#79706E",
                               "Host_associated"="#E15759","Mammals: Human"="#D37295"
  )) +
  labs( x = "Longitude (°)", y = "Latitude (°)") +
  theme_bw()+
  scale_y_continuous(expand = expansion(mult=c(0,0)))+
  scale_x_continuous(expand = expansion(add=c(0,0)))+
  theme(legend.position = 'none')
world.map


world.map2 <- 
  ggplot() +
  geom_polygon(data=world,aes(x=long,y=lat,group=group),
               fill="#e9ecef")+
  geom_point(data=redundancy_count, 
             aes(x=Longitude,y=Latitude, 
                 color = `Ecosystem Category`, size = virus_number),
  )+
  scale_color_manual(values = c("Built environment"="#9b5de5","Terrestrial"="#FFBE7D",
                               "Plants"="#59A14F","Aquatic"="#00b4d8","Unknown"="#79706E",
                               "Host_associated"="#E15759","Mammals: Human"="#D37295"
  )) +
  theme_bw()+
  scale_y_continuous(expand = expansion(mult=c(0,0)))+
  scale_x_continuous(expand = expansion(add=c(0,0)))+
  theme(legend.position = 'none')
world.map2

library(Cairo)
library(cowplot)
#cairo_pdf(filename = file.path(CS_WORK, "1_genome_metadata.pdf") , width = 8,height = 4,unit = 'in')
g1 <- plot_grid(world.map,p_virus_unclassified_MGE,
                labels = letters[2], ncol = 2,
                rel_widths = c(2, 1)
                # label_x = 0.09,            # 标签x坐标（1为最右侧）
                #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g1

ggsave(file.path(CS_WORK, "1_phage_geographical_distribution.pdf"),g1,width = 8.46,height = 3.66,units = 'in')


ggsave(file.path(CS_WORK, "1_p_IMGVR_phage_ecotype.pdf"),
       p_IMGVR_phage_ecotype,
       width = 8,height = 5,units = 'in')

ggsave(file.path(CS_WORK, "1_p_world.map.pdf"),
       world.map,
       width = 5.6,height = 3.4,units = 'in')


#vire病毒####

ESKAPE_vire_phage <- read.delim2("ESKAPE_vire_phage.tsv")

vire_database_info <- read.delim2("vire_matched_rows.tsv")

# 加载必要的包
library(dplyr)
library(tidyr)
library(stringr)
# 创建主要类别分类
classify_ecosystem_main <- function(microntology) {
  case_when(
    # 人类相关
    str_detect(microntology, "human host") & 
      str_detect(microntology, "digestive tract|intestine") ~ "Human_Gut",
    
    str_detect(microntology, "human host") & 
      str_detect(microntology, "skin") ~ "Human_Skin",
    
    str_detect(microntology, "human host") & 
      str_detect(microntology, "airways|respiratory") ~ "Human_Respiratory",
    
    str_detect(microntology, "human host") & 
      str_detect(microntology, "urogenital") ~ "Human_Urogenital",
    
    str_detect(microntology, "human host") & 
      str_detect(microntology, "mouth|oral") ~ "Human_Oral",
    
    str_detect(microntology, "human host") & 
      !str_detect(microntology, "digestive tract|intestine|skin|airways|urogenital|mouth") ~ "Human_Other",
    
    # 其他动物相关
    str_detect(microntology, "mammalian host|animal host") & 
      !str_detect(microntology, "human") ~ "Other_Animals",
    
    # 植物相关
    str_detect(microntology, "plant host") ~ "Plants",
    
    # 水生环境
    str_detect(microntology, "aquatic:marine") ~ "Marine",
    str_detect(microntology, "aquatic:fresh water") ~ "Freshwater",
    str_detect(microntology, "aquatic:brackish") ~ "Brackish",
    str_detect(microntology, "aquatic") ~ "Aquatic_Other",
    
    # 陆地环境
    str_detect(microntology, "terrestrial:soil") ~ "Soil",
    str_detect(microntology, "terrestrial:forest") ~ "Forest",
    str_detect(microntology, "terrestrial:grassland") ~ "Grassland",
    str_detect(microntology, "terrestrial") ~ "Terrestrial_Other",
    
    # 人工环境
    str_detect(microntology, "built environment") ~ "Built_Environment",
    str_detect(microntology, "wastewater") ~ "Wastewater",
    str_detect(microntology, "agriculture") ~ "Agricultural",
    str_detect(microntology, "food") ~ "Food_Related",
    str_detect(microntology, "anthropogenic") ~ "Anthropogenic_Other",
    
    # 其他
    str_detect(microntology, "hydrothermal") ~ "Hydrothermal",
    str_detect(microntology, "contaminated") ~ "Contaminated",
    str_detect(microntology, "ancient") ~ "Ancient",
    
    TRUE ~ "Unclassified"
  )
}

# 应用主要类别分类
vire_database_info <- vire_database_info %>%
  mutate(ecosystem_main = classify_ecosystem_main(microntology))

# 查看主要类别统计
cat("\n=== 主要生态系统类别统计 ===\n")
main_summary <- table(vire_database_info$ecosystem_main, useNA = "ifany") %>%
  sort(decreasing = TRUE)
print(main_summary)

category_summary_vire <- vire_database_info %>%
  count(ecosystem_main) %>%
  arrange(desc(n))
print(category_summary_vire)

nature_colors <- c("#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#79706E",
                   "#B07AA1", "#FFBE7D", "#8CD17D", "#86BCB6", "#499894", "#D37295")

# 水平条形图版本
p_vire_phage_ecotype <- 
  ggplot(category_summary_vire, aes(x = n, y = reorder(ecosystem_main, n), 
                               fill = ecosystem_main)) +
  geom_bar(stat = "identity", width = 0.7, color = "black",linewidth = 0.1) +
  geom_text(aes(label = n), hjust = -0.2, size = 3.5, color = "black") +
  labs(
    y = "",
    x = "Number of ESKAPE phages identified from VIRE"
  ) +
  scale_fill_manual(values = c("Built environment"="#9b5de5","Terrestrial_Other"="#FFBE7D",
                               "Plants"="#59A14F","Freshwater"="#00b4d8",
                               "Aquatic_Other"="#00b4d8",
                               "Unknown"="#79706E",
                               "Other_Animals"="#E15759","Human_Gut"="#D37295"
  )) +
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, linewidth = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_text( size = 12, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 12, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
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
p_vire_phage_ecotype

ggsave(file.path(CS_WORK, "1_p_vire_phage_ecotype.pdf"),
       p_vire_phage_ecotype,
       width = 3.24,height = 5,units = 'in')

save.image("./1_phage_geographical_distribution.RData")
