# ===========================================================================
# ESKAPE_mobilome_virome -- analyses / 05_phage_clustering
# vContact3 genome-sharing network and taxonomic assignment
# Manuscript      : Fig. 5A
# Original path   : 2ESKAPE_MGE/5phage_vcontact3/1_vcontact3_analysis.R
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

setwd(CS_EXT)

myxoco_all_phage_info <- read.delim2("phage_detailed_info_trim.tsv")

myxoco_complete_phage_info <- 
  myxoco_all_phage_info[myxoco_all_phage_info$checkv_quality %in% 
                                                      c("Complete"),]

vcontact3_result <- read.csv("part1.cyjs default node.csv")
length(unique(vcontact3_result$id))

potential_novel_complete_phage <- 
  myxoco_complete_phage_info[! myxoco_complete_phage_info$seq_name %in% vcontact3_result$id,]
head(potential_novel_complete_phage)


write.table(potential_novel_complete_phage, 
            "potential_novel_complete_phage_ID.tsv", 
            sep = "\t",           # TSV使用制表符分隔
            row.names = FALSE,    # 不写入行名
            col.names = TRUE,     # 写入列名
            quote = FALSE        # 不添加引号
            )






save.image("1_vcontact3_analysis.RData")

load("1_vcontact3_analysis.RData")




