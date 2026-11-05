# ===========================================================================
# ESKAPE_mobilome_virome -- analyses / 07_phage_host
# Functional types of the viral genes targeted by host spacers
# Manuscript      : Fig. 4C
# Original path   : 2ESKAPE_MGE/7_2protospacer_gene/1_gene_annotation.R
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


setwd(file.path(CS_DATA, "7_2protospacer_gene"))

eggnog_res <- readxl::read_excel("out.emapper.annotations.xlsx")

# 加载必要的包
library(dplyr)
library(stringr)
library(tidyr)

# 创建gene_type分类函数
classify_gene_type <- function(description, pfams = NA) {
  # 合并Description和PFAMs列进行搜索
  desc <- tolower(as.character(description))
  pfam <- tolower(as.character(pfams))
  combined <- paste(desc, pfam, sep = " ")
  
  case_when(
    # ============= 1. 尾丝蛋白相关 (tail fiber) =============
    str_detect(combined, "tail fiber|tail fibre|fibronectin|tailfiber|tail fibre protein|fiber protein") ~ "Tail fiber",
    
    # ============= 2. 尾管蛋白 (tail tube) =============
    str_detect(combined, "tail tube|tail tubular|tube protein|tail_tube") ~ "Tail tube",
    
    # ============= 3. 尾鞘蛋白 (tail sheath) =============
    str_detect(combined, "tail sheath|sheath protein|tail_sheath") ~ "Tail sheath",
    
    # ============= 4. 尾丝蛋白组件 (tail component) - 您特别指定的 =============
    str_detect(combined, "tail component") ~ "Tail component",
    
    # ============= 5. 尾丝蛋白测量蛋白 (tape measure) =============
    str_detect(combined, "tape measure|tape measure protein|tail tape measure") ~ "Tape measure",
    
    # ============= 6. 其他尾丝相关蛋白 =============
    str_detect(combined, "minor tail|tail|putative tail|tail collar|tail-collar|tail collar domain|collar domain") ~ "Tail associated",
    
    # ============= 7. 基板蛋白 (baseplate) =============
    str_detect(combined, "baseplate|phage-baseplate|base plate|base-plate|baseplate j|baseplate assembly") ~ "Baseplate",
    
    # ============= 8. 噬菌体裂解相关 =============
    str_detect(combined, "lysozyme|lytic transglycosylase|transglycosylase|lytic|endolysin|lysin|holin|spanin|lysis|transglycosylase|peptidoglycan|amidase|lytic|muramidase|n-acetylmuramoyl|lytic enzyme|cell wall hydrolase|nlpc/p60|glycosyl hydrolase") ~ "Lysis",
    
    # ============= 9. 噬菌体结构蛋白 =============
    str_detect(combined, "capsid|head protein|head morphogenesis|scaffolding|portal protein|phage portal|head-tail connector|head completion|coat protein|major capsid") ~ "Structural capsid",
    
    # ============= 10. 噬菌体组装/形态发生 =============
    str_detect(combined, "morphogenesis|assembly|chaperone|virion assembly|head-tail joining") ~ "Assembly",
    
    # ============= 11. 噬菌体DNA包装 =============
    str_detect(combined, "terminase|terminase-like|terminase large|terminase small|dna packaging|packaging|gpA|gpP|nu1|terminase_") ~ "DNA packaging",
    
    # ============= 12. 噬菌体DNA代谢 =============
    str_detect(combined, "helicase|helicases|primase|polymerase|primase-helicase|toprim-like|toprim_2|dna polymerase|dna primase|replication|replicase|dna replication|dna binding|dna metabolic|helicase activity|single-strand binding") ~ "DNA replication",
    
    # ============= 13. 重组/修复 =============
    str_detect(combined, "recombinase|exodeoxyribonuclease|recombination|rec|reca|recb|recd|rect|rece|recombinational|ruv|holiday|rusa|yqaj|exonuclease|endonuclease|nuclease|integrase|integration|recombinase") ~ "Recombination/repair",
    
    # ============= 14. 甲基化相关 =============
    str_detect(combined, "methylase|Methyltransferase|methyltransferase|methyltransf|dna methylase|dam|dcm|methylase|methyltransferase|methylase activity|methyltransferase activity|methylates") ~ "Methylation",
    
    # ============= 15. 转录调控 =============
    str_detect(combined, "transcription factor|transcription|regulation|transcriptional regulator|repressor|antirepressor|activator|cro|cI|regulatory|helix-turn-helix|hth domain|transcription regulation") ~ "Transcriptional regulation",
    
    # ============= 16. 抗转录终止 =============
    str_detect(combined, "antitermination|antiterminator|antitermination protein|q protein") ~ "Antitermination",
    
    # ============= 17. 毒素-抗毒素系统 =============
    str_detect(combined, "toxin|antitoxin|addiction|kila|toxin-antitoxin|ydat|toxin") ~ "Toxin-antitoxin",
    
    # ============= 18. 细胞壁降解酶 =============
    str_detect(combined, "cellulase|chitinase|chitin catabolic|chitinase activity|glycosyl hydrolase family|cellulase activity") ~ "Cell wall degradation",
    
    # ============= 19. 激酶/磷酸酶 =============
    str_detect(combined, "kinase|phosphatase|phosphorylase|phosphoesterase|atpase|gtpase") ~ "Kinase/phosphatase",
    
    # ============= 20. 未知功能结构域 =============
    str_detect(combined, "duf[0-9]+|domain of unknown function|uncharacterized protein|conserved protein|unknown function|hypothetical protein") ~ "unknown_function",
    
    # ============= 21. 噬菌体整合酶 =============
    str_detect(combined, "integrase|int family|phage integrase|lambda integrase|arm dna-binding") ~ "Integrase",
    
    # ============= 22. SOS反应相关 =============
    str_detect(combined, "sos response|lexa|reca|sos|damage-inducible") ~ "SOS response",
    
    # ============= 23. 其他酶类 =============
    str_detect(combined, "dehydrogenase|oxidoreductase|synthase|synthetase|transferase|hydrolase|ligase|isomerase|lyase|protease|peptidase") ~ "Enzyme other",
    
    # ============= 24. 膜蛋白/转运 =============
    str_detect(combined, "membrane protein|transmembrane|integral membrane|transporter|atp-binding protein|abc transporter|secretion") ~ "Membrane transport",
    
    # ============= 25. 其他噬菌体相关 =============
    str_detect(combined, "phage related|Siphovirus|phage protein|prophage|bacteriophage") ~ "Phage other",
    
    # ============= 26. 染色体分离 =============
    str_detect(combined, "chromosome segregation|nuclear chromosome segregation") ~ "Chromosome segregation",
    
    # ============= 26. 染色体分离 =============
    str_detect(combined, "-") ~ "unannotated",
    
    # ============= 27. 默认分类 =============
    TRUE ~ "Other"
  )
}

#eggnog_res <- eggnog_res[eggnog_res$Description!="-",]

# 应用分类函数
eggnog_res <- eggnog_res %>%
  mutate(gene_type = classify_gene_type(Description, PFAMs))

# 查看分类结果
cat("=== 基因类型分类统计 ===\n")
type_summary <- table(eggnog_res$gene_type, useNA = "ifany") %>%
  sort(decreasing = TRUE)
print(type_summary)

# 计算百分比
type_percent <- prop.table(type_summary) * 100
cat("\n=== 百分比分布 ===\n")
print(round(type_percent, 2))


library(ggplot2)
library(forcats)

# 准备数据用于绘图
plot_data <- as.data.frame(type_summary) %>%
  rename(GeneType = Var1, Count = Freq) %>%
  # 过滤掉低频类别以便于可视化
  filter(Count > 10) %>%
  mutate(Percentage = Count / sum(Count) * 100)

plot_data <- plot_data[!plot_data$GeneType %in% c("unannotated","unknown_function"),]


# 条形图
p_phage_gene_targeted <-
ggplot(plot_data, aes(x = fct_reorder(GeneType, Count), y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(
       x = "Gene type",
       y = "Number of phage gene targeted by CRISPR spacer") +
  theme(
    panel.spacing = unit(0.1, "lines"),  # 将分面间距设为0
    panel.background = element_blank(),  # 设置背景为透明
    # panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, linewidth = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    #axis.title = element_blank(),
    axis.title = element_text( size = 12, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
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
p_phage_gene_targeted

ggsave(plot = p_phage_gene_targeted,'1_p_phage_gene_targeted.pdf',
       width = 6.22,height = 4.88,units = 'in')

save.image("./1_gene_annotation.RData")

load("./1_gene_annotation.RData")


ESKAPE_phage_methylase <- eggnog_res[eggnog_res$gene_type == "Methylation",]

write.table(ESKAPE_phage_methylase, 
            "ESKAPE_phage_methylase.tsv", 
            sep = "\t",           # TSV使用制表符分隔
            row.names = FALSE,    # 不写入行名
            col.names = TRUE,     # 写入列名
            quote = FALSE        # 不添加引号
            )





