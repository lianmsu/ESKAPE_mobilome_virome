# ===========================================================================
# ESKAPE_mobilome_virome -- analyses / 07_phage_host
# Multi-host phages and their CRISPR-based host links
# Manuscript      : Fig. 4
# Original path   : 2ESKAPE_MGE/7_0multi_host_phage/1_multi_host_phage.R
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


setwd(file.path(CS_DATA, "7_0multi_host_phage"))

ESKAPE_all_phage_info_filter <- read.delim2("ESKAPE_all_phage_info_filter.tsv")

EnKp_phage <- 
  ESKAPE_all_phage_info_filter[ESKAPE_all_phage_info_filter$Host_tax == "En,Kp",]
EnKp_phage_complete <- 
  EnKp_phage[EnKp_phage$checkv_quality=="Complete",]

En_phage <- 
  ESKAPE_all_phage_info_filter[ESKAPE_all_phage_info_filter$Host_tax == "En",]
En_phage_complete <-
  En_phage[En_phage$checkv_quality=="Complete",]

Kp_phage <- 
  ESKAPE_all_phage_info_filter[ESKAPE_all_phage_info_filter$Host_tax == "Kp",]
Kp_phage_complete <-
  Kp_phage[Kp_phage$checkv_quality=="Complete",]

# 102 47 443

#提取对应的complete的病毒信息
En_Kp_EnKp_phage_complete <- 
  rbind(En_phage_complete, Kp_phage_complete, EnKp_phage_complete)

write.table(En_Kp_EnKp_phage_complete,
            "En_Kp_EnKp_phage_complete.tsv",
            sep = "\t",           # TSV使用制表符分隔
            row.names = FALSE,    # 不写入行名
            col.names = TRUE,     # 写入列名
            quote = FALSE        # 不添加引号
            )            # NA值的表示方式

save.image("./1_multi_host_phage.RData")

#根据vcontact2或vcontact3的结果进行过滤，只保留目标病毒
vcontact2_res <- read.table("vcontact2_res.ntw",header = T)
vcontact2_res_filter <- 
  vcontact2_res[(vcontact2_res$Source %in% En_Kp_EnKp_phage_complete$ID) &
  (vcontact2_res$Target %in% En_Kp_EnKp_phage_complete$ID),
   ]
write.table(vcontact2_res_filter,
            "vcontact2_res_filter.tsv",
            sep = "\t",           # TSV使用制表符分隔
            row.names = FALSE,    # 不写入行名
            col.names = TRUE,     # 写入列名
            quote = FALSE        # 不添加引号
            )            # NA值的表示方式

load("./1_multi_host_phage.RData")


vcontact3_res <- read.csv("part1_cosmo_data.csv",sep = ";")
vcontact3_res_filter <- 
  vcontact3_res[(vcontact3_res$source %in% En_Kp_EnKp_phage_complete$ID) &
  (vcontact3_res$target %in% En_Kp_EnKp_phage_complete$ID),
  ]
write.table(vcontact3_res_filter,
            "vcontact3_res_filter.tsv",
            sep = "\t",           # TSV使用制表符分隔
            row.names = FALSE,    # 不写入行名
            col.names = TRUE,     # 写入列名
            quote = FALSE        # 不添加引号
          )            # NA值的表示方式


vcontact3_res_metadata <- read.csv("part1_cosmo_metadata.csv",sep = ";")

save.image("./1_multi_host_phage.RData")





