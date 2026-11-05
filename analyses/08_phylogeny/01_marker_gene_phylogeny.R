# ===========================================================================
# ESKAPE_mobilome_virome -- analyses / 08_phylogeny
# Marker-gene phylogenies (MCP, portal, TerL) of the phages
# Manuscript      : Fig. 6
# Original path   : 2ESKAPE_MGE/7_3marker_gene_phylogeny/1_phage_gene_annotation.R
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


setwd(file.path(CS_DATA, "7_3marker_gene_phylogeny"))

all_32phage_gff <- read.delim2("32_phage_treated.gff",
                                header = F, sep = '\t')

library(dplyr)
library(stringr)
all_32phage_gff <- all_32phage_gff %>%
  mutate(
    # 提取ID=后面的内容，然后取最后一个_后的数字
    gene_num = str_extract(V9, "ID=[^;]+") %>% 
      str_remove("ID=") %>% 
      str_extract("[^_]+$"),
    hit_id = paste0(V1, "_", gene_num)
  ) %>%
  select(-gene_num)  # 移除临时列
all_32phage_gff <- all_32phage_gff[,c(1,3,4,5,7,10)]

colnames(all_32phage_gff) <- c("contig","type","start","stop","strand",
                           "name")

#type和name可以改为需要的内容（type决定颜色，name决定标签）

#write.table(all_phage_gff,"all_phage_gff.tsv",sep = "\t",quote = F)


primary_anno <- read.delim2("32phage_annotation_summary.tsv")

all_32phage_gff$seq <- sprintf("seq%d", 1:nrow(all_32phage_gff))


all_32phage_gff <- all_32phage_gff %>%
  left_join(primary_anno, by = "seq")

write.table(all_32phage_gff,"all_32phage_gff_anno.tsv",sep = "\t",quote = F)
save.image("./1_phage_gene_annotation.RData")
################################################################################

#major capsid protein

# 忽略大小写进行匹配（推荐，因为注释中可能出现"Major capsid"等变体）
filtered_32phage_MCP <- 
  subset(all_32phage_gff, 
         grepl("major capsid", detail_anno, ignore.case = TRUE) &
         !grepl("head decoration", detail_anno, ignore.case = TRUE) &
         !grepl("head maturation protease", detail_anno, ignore.case = TRUE) &
         !grepl("5WK1_M", detail_anno, ignore.case = TRUE)
           )

filtered_32phage_TerL <- 
  subset(all_32phage_gff, 
         grepl("terminase large|large terminase|large subunit", 
               detail_anno, ignore.case = TRUE))

filtered_32phage_portal <- 
  subset(all_32phage_gff, 
         grepl("portal protein|portal", detail_anno, ignore.case = TRUE) &
           !grepl("Tail terminator|Tail tubular|Gene 7|Major tail protein|Minor capsid protein", detail_anno, ignore.case = TRUE)
         )

save.image("./1_phage_gene_annotation.RData")


write.table(filtered_32phage_MCP, 
            "filtered_32phage_MCP.tsv", 
            sep = "\t",           # TSV使用制表符分隔
            row.names = FALSE,    # 不写入行名
            col.names = TRUE,     # 写入列名
            quote = FALSE        # 不添加引号
            )
write.table(filtered_32phage_portal, 
            "filtered_32phage_portal.tsv", 
            sep = "\t",           # TSV使用制表符分隔
            row.names = FALSE,    # 不写入行名
            col.names = TRUE,     # 写入列名
            quote = FALSE        # 不添加引号
            )
write.table(filtered_32phage_TerL, 
            "filtered_32phage_TerL.tsv", 
            sep = "\t",           # TSV使用制表符分隔
            row.names = FALSE,    # 不写入行名
            col.names = TRUE,     # 写入列名
            quote = FALSE        # 不添加引号
            )

selected_7_phage <- read.delim2("selected_phage_anno_withtype.tsv")
selected_7_phage_MCP <- 
  subset(selected_7_phage, 
         grepl("major capsid", detail_anno, ignore.case = TRUE) &
           !grepl("head decoration", detail_anno, ignore.case = TRUE) &
           !grepl("head maturation protease", detail_anno, ignore.case = TRUE) &
           !grepl("5WK1_M", detail_anno, ignore.case = TRUE)
        )
write.table(selected_7_phage_MCP, 
            "selected_7_phage_MCP.tsv", 
            sep = "\t",           # TSV使用制表符分隔
            row.names = FALSE,    # 不写入行名
            col.names = TRUE,     # 写入列名
            quote = FALSE        # 不添加引号
        )
selected_7_phage_TerL <- 
  subset(selected_7_phage, 
         grepl("terminase large|large terminase|large subunit", 
               detail_anno, ignore.case = TRUE))
write.table(selected_7_phage_TerL, 
            "selected_7_phage_TerL.tsv", 
            sep = "\t",           # TSV使用制表符分隔
            row.names = FALSE,    # 不写入行名
            col.names = TRUE,     # 写入列名
            quote = FALSE        # 不添加引号
        )
selected_7_phage_portal <- 
  subset(selected_7_phage, 
         grepl("portal protein|portal", detail_anno, ignore.case = TRUE) &
           !grepl("Tail terminator|Tail tubular|Gene 7|Major tail protein|Minor capsid protein", detail_anno, ignore.case = TRUE)
  )

write.table(selected_7_phage_portal, 
            "selected_7_phage_portal.tsv", 
            sep = "\t",           # TSV使用制表符分隔
            row.names = FALSE,    # 不写入行名
            col.names = TRUE,     # 写入列名
            quote = FALSE        # 不添加引号
        )

save.image("./1_phage_gene_annotation.RData")


################################################################################
all_phage_gff <- read.delim2("ESKAPE_complete_phage_treated.gff",
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

#write.table(all_phage_gff,"all_phage_gff.tsv",sep = "\t",quote = F)


primary_anno <- read.delim2("phage_annotation_summary.tsv")

all_phage_gff$seq <- sprintf("seq%d", 1:nrow(all_phage_gff))


all_phage_gff <- all_phage_gff %>%
  left_join(primary_anno, by = "seq")

write.table(all_phage_gff,"all_phage_gff_anno.tsv",sep = "\t",quote = F)
save.image("./1_phage_gene_annotation.RData")
#######################

#major capsid protein

# 忽略大小写进行匹配（推荐，因为注释中可能出现"Major capsid"等变体）
filtered_phage_MCP <- 
  subset(all_phage_gff, 
         grepl("major capsid", detail_anno, ignore.case = TRUE) &
           !grepl("head decoration", detail_anno, ignore.case = TRUE) &
           !grepl("head maturation protease", detail_anno, ignore.case = TRUE) &
           !grepl("5WK1_M", detail_anno, ignore.case = TRUE)
  )

filtered_phage_TerL <- 
  subset(all_phage_gff, 
         grepl("terminase large|large terminase|large subunit", 
               detail_anno, ignore.case = TRUE))

filtered_phage_portal <- 
  subset(all_phage_gff, 
         grepl("portal protein|portal", detail_anno, ignore.case = TRUE) &
           !grepl("Tail terminator|Tail tubular|Gene 7|Major tail protein|Minor capsid protein", detail_anno, ignore.case = TRUE)
  )





save.image("./1_phage_gene_annotation.RData")

ESKAPE_phage_metadata <- read.delim2("ESKAPE_all_phage_info_filter.tsv")

filtered_phage_MCP <- filtered_phage_MCP %>%
  left_join(ESKAPE_phage_metadata, by = c("contig"="ID"))

filtered_phage_portal <- filtered_phage_portal %>%
  left_join(ESKAPE_phage_metadata, by = c("contig"="ID"))

filtered_phage_TerL <- filtered_phage_TerL %>%
  left_join(ESKAPE_phage_metadata, by = c("contig"="ID"))

write.table(filtered_phage_MCP, 
            "filtered_phage_MCP.tsv", 
            sep = "\t",           # TSV使用制表符分隔
            row.names = FALSE,    # 不写入行名
            col.names = TRUE,     # 写入列名
            quote = FALSE        # 不添加引号
)
write.table(filtered_phage_portal, 
            "filtered_phage_portal.tsv", 
            sep = "\t",           # TSV使用制表符分隔
            row.names = FALSE,    # 不写入行名
            col.names = TRUE,     # 写入列名
            quote = FALSE        # 不添加引号
)
write.table(filtered_phage_TerL, 
            "filtered_phage_TerL.tsv", 
            sep = "\t",           # TSV使用制表符分隔
            row.names = FALSE,    # 不写入行名
            col.names = TRUE,     # 写入列名
            quote = FALSE        # 不添加引号
)

