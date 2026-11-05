# ===========================================================================
# ESKAPE_mobilome_virome -- analyses / 06_phage_defense
# Functional annotation of the phage genes (hhsuite2 -> KEGG/COG)
# Manuscript      : Fig. 3F
# Original path   : 2ESKAPE_MGE/6_1phage_gene_annotation/1_phage_gene_annotation.R
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


setwd(file.path(CS_DATA, "6_1phage_gene_annotation"))

all_phage_gff <- read.delim2("selected_phage_treated.gff",
                                header = F, sep = '\t')

library(dplyr)
library(stringr)
all_phage_gff <- all_phage_gff %>%
  mutate(
    # 提取ID=后面的内容，然后取最后一个_后的数字
    gene_num = str_extract(V9, "ID=[^;]+") %>% 
      str_remove("ID=") %>% 
      str_extract("[^_]+$"),
    hit_id = paste0(V1, "_", gene_num)
  ) %>%
  select(-gene_num)  # 移除临时列
all_phage_gff <- all_phage_gff[,c(1,3,4,5,7,10)]

colnames(all_phage_gff) <- c("contig","type","start","stop","strand",
                           "name")

#type和name可以改为需要的内容（type决定颜色，name决定标签）

write.table(all_phage_gff,"all_phage_gff.tsv",sep = "\t",quote = F)


primary_anno <- read.delim2("phage_annotation_summary.tsv")

all_phage_gff$seq <- sprintf("seq%d", 1:nrow(all_phage_gff))


all_phage_gff <- all_phage_gff %>%
  left_join(primary_anno, by = "seq")

write.table(all_phage_gff,"all_phage_gff_anno.tsv",sep = "\t",quote = F)
save.image("./1_phage_gene_annotation.RData")
################################################################################

all_phage_gff_manually <- read.delim2("all_phage_gff_anno.tsv")

library(dplyr)
library(stringr)

classified_gff <- all_phage_gff_manually %>%
  mutate(
    # 保留原始type列信息
    original_type = type,
    
    # 开始根据name列的内容创建新的type列
    type = case_when(
      # ========== 1. Lysis-related proteins ==========
      # 裂解相关：溶菌酶、孔蛋白、内溶素、跨膜蛋白等
      str_detect(name, regex("lysozyme|endolysin|holin|lysis|spanin|transglycosylase|amidase|L-alanyl-D-glutamate peptidase|peptidase|murein", ignore_case = TRUE)) ~ "Lysis",
      
      # ========== 2. Virion structure and morphogenesis ==========
      # 病毒颗粒结构和形态发生：衣壳、尾部、纤维、门户、组装蛋白等
      str_detect(name, regex("capsid|tail|fiber|portal|terminase|virion|morphogenesis|head|sheath|baseplate|ATPase|tape|neck|assembly|chaperone|adaptor|completion|decoration|structural protein|spike|wedge|hub|tube|major tail|minor tail|major head|minor head|tail terminator|head-tail adaptor|head maturation protease|head closure|baseplate|tail length tape|tail assembly|head decoration|scaffolding", ignore_case = TRUE)) ~ "Virion structure and morphogenesis",
      
      # ========== 3. DNA replication, recombination, repair and packaging ==========
      # DNA复制、重组、修复和包装
      str_detect(name, regex("replication|replicase|helicase|primase|polymerase|ligase|recombination|recombinase|Rec|integrase|excisionase|exonuclease|nuclease|endonuclease|terminase|packaging|Dna|parition|Par|partition|segregation|condensation|Scp|SMC|RusA|Ruv|Holliday|resolution|resolvase|ssb|single.*strand.*DNA|DNA.*binding|DNA.*methyl|methyltransferase|DNA.*repair|Rad|Sbc|nucleotidyltransferase|initiation|recombination|antirecombination", ignore_case = TRUE)) ~ "DNA replication, recombination, repair and packaging",
      
      # ========== 4. Host interaction and anti-defense systems ==========
      # 宿主互作和抗防御系统
      str_detect(name, regex("anti.*restriction|antirestriction|anti.*defense|anti.*CRISPR|DarB|Acr|immunity|superinfection|toxin|antitoxin|addiction|antibiotic|resistance|virulence|vir|Vfr|regulator|transcriptional|repressor|activator|ArsR|LysR|TetR|LacI|MerR|AraC|sigma|anti-termination|Q|CII|CIII|regulatory|FleQ|ComR|chromatin|remodeler|dps|oxidative|stress|HigB|toxin|antidote", ignore_case = TRUE)) ~ "Host interaction and anti-defense systems",
      
      # ========== 5. Metabolism and biosynthesis ==========
      # 代谢和生物合成
      str_detect(name, regex("phosphofructokinase|pyruvate.*dehydrogenase|kinase|dehydrogenase|transferase|reductase|synthase|biosynthesis|metabolism|metabolic|glycosyltransferase|methyltransferase|acetyltransferase|phosphatase|phosphorylase|peptidase|dipeptidase|dipeptidyl.*peptidase|protease|phosphoadenosine", ignore_case = TRUE)) ~ "Metabolism and biosynthesis",
      
      # ========== 6. Membrane and transport proteins ==========
      # 膜蛋白和转运蛋白
      str_detect(name, regex("membrane.*protein|periplasmic|lipoprotein|transport|transporter|secretion|Tat|YifL|channel|pore|porin|efflux", ignore_case = TRUE)) ~ "Membrane and transport proteins",
      
      # ========== 7. DNA modification and restriction-modification ==========
      # DNA修饰和限制-修饰系统
      str_detect(name, regex("restriction.*endonuclease|modification.*methylase|HhaI|RsrI|DNA.*methylase|methylase|restriction|modification", ignore_case = TRUE)) ~ "DNA modification and restriction-modification",
      
      # ========== 8. Unannotated/hypothetical proteins ==========
      # 未注释/假设蛋白
      is.na(name) | name == "" | str_detect(name, regex("^[0-9A-Z_]+$|^gp[0-9]+$|^orf[0-9]+$|hypothetical|unknown|unannotated|predicted|putative|conserved.*protein|domain.*protein|family.*protein|^protein.*$|^[A-Za-z]+[0-9]+$", ignore_case = TRUE)) ~ "Unannotated/hypothetical proteins",
      
      # ========== 9. Other and miscellaneous ==========
      # 其他和杂项
      TRUE ~ "Other function"
    ),
    
    # 添加更详细的子分类（可选）
    subtype = case_when(
      # Lysis 子分类
      str_detect(name, regex("lysozyme|endolysin|transglycosylase|amidase|peptidase", ignore_case = TRUE)) & type == "Lysis" ~ "Endolysin",
      str_detect(name, regex("holin|spanin", ignore_case = TRUE)) & type == "Lysis" ~ "Membrane disruption",
      
      # Virion 子分类
      str_detect(name, regex("capsid|head|major capsid", ignore_case = TRUE)) & type == "Virion structure and morphogenesis" ~ "Capsid/Head",
      str_detect(name, regex("tail|sheath|tube|fiber", ignore_case = TRUE)) & type == "Virion structure and morphogenesis" ~ "Tail",
      str_detect(name, regex("portal|terminase", ignore_case = TRUE)) & type == "Virion structure and morphogenesis" ~ "DNA packaging",
      str_detect(name, regex("assembly|chaperone|adaptor|completion", ignore_case = TRUE)) & type == "Virion structure and morphogenesis" ~ "Assembly/Chaperone",
      
      # DNA 子分类
      str_detect(name, regex("replication|primase|helicase|polymerase", ignore_case = TRUE)) & type == "DNA replication, recombination, repair and packaging" ~ "Replication",
      str_detect(name, regex("recombination|recombinase|integrase|excisionase", ignore_case = TRUE)) & type == "DNA replication, recombination, repair and packaging" ~ "Recombination/integration",
      str_detect(name, regex("repair|nuclease|exonuclease|endonuclease|ligase", ignore_case = TRUE)) & type == "DNA replication, recombination, repair and packaging" ~ "Repair",
      str_detect(name, regex("partition|Par|segregation", ignore_case = TRUE)) & type == "DNA replication, recombination, repair and packaging" ~ "Partition/segregation",
      
      # 其他情况的默认值
      TRUE ~ NA_character_
    )
  )

# 查看分类结果统计
category_summary <- classified_gff %>%
  count(type, sort = TRUE)

print("基因功能分类统计:")
print(category_summary)

# 查看每个类别中的示例基因
examples_by_category <- classified_gff %>%
  group_by(type) %>%
  slice_head(n = 5) %>%
  select(name, type, subtype) %>%
  arrange(type)

print("各分类示例基因:")
print(examples_by_category, n = Inf)

# 验证一些关键基因的分类
test_genes <- c(
  "Lysozyme; peptidoglycan dydrolase, cationic peptides, monomer",
  "Major Capsid",
  "Terminase large subunit",
  "DNA methyltransferase",
  "holin; lysis",
  "Integrase",
  "anti-restriction protein",
  "transcriptional regulator",
  "phosphofructokinase",
  "membrane protein"
)

test_results <- classified_gff %>%
  filter(name %in% test_genes) %>%
  select(name, type, subtype) %>%
  distinct()

print("测试基因分类结果:")
print(test_results)


write.table(classified_gff,"selected_phage_anno_withtype.tsv",sep = "\t",quote = F)

save.image("./1_phage_gene_annotation.RData")

#######

setwd(CS_EXT)

my_17_phage <- read.delim2("Candidate_17_new_phages.tsv")


phage_genomad <- read.delim2("nr_phage_rep_seq_virus_summary.tsv")

library(dplyr)
my_17_phage <- my_17_phage %>% 
  left_join(phage_genomad,by=c("seq_name"="seq_name"))


phage_protospacer <- read.delim2("filtered_candidate17.out",header = F)

phage_protospacer <- phage_protospacer[,c(1,2,15,16)]
colnames(phage_protospacer) <- c("host","contig","start","stop")
phage_protospacer$type <- "Protospacer"

library(stringr)
phage_protospacer <- phage_protospacer %>% 
  mutate(
    host      = str_remove(host, "-[^-]*$"),      # 删掉最后一个 “-” 及之后
    contig = str_remove(contig, "\\|.*")   # 删掉第一个 “|” 及之后
  )

phage_protospacer <- unique(phage_protospacer)
phage_protospacer$strand <- "+"

#host info

all_4649_genome_metadata <- read.delim2("all_4649_genome_metadata.tsv")

phage_protospacer <- phage_protospacer %>%
  left_join(all_4649_genome_metadata, by = c("host"="Genome"))

phage_protospacer <- phage_protospacer[,c(2:21,1)]

write.table(phage_protospacer,"phage_protospacer_anno.tsv",quote = F,
            sep = '\t')






