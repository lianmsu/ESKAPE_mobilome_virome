# ===========================================================================
# ESKAPE_mobilome_virome -- analyses / 01_genome_landscape
# Genome metadata, quality and the world map of the 4 538 genomes
# Manuscript      : Fig. 1A-B
# Original path   : 2ESKAPE_MGE/1genome/1_genome_info.R
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


setwd(file.path(CS_DATA, "1genome"))

Ef_NCBI_info <- read.delim2("data_summary_Ef.tsv")
Sa_NCBI_info <- read.delim2("data_summary_Sa.tsv")
Kp_NCBI_info <- read.delim2("data_summary_Kp.tsv")
Ab_NCBI_info <- read.delim2("data_summary_Ab.tsv")
Pa_NCBI_info <- read.delim2("data_summary_Pa.tsv")
En_NCBI_info <- read.delim2("data_summary_En.tsv")

Ef_NCBI_info$Tax <- "Ef"
Sa_NCBI_info$Tax <- "Sa"
Kp_NCBI_info$Tax <- "Kp"
Ab_NCBI_info$Tax <- "Ab"
Pa_NCBI_info$Tax <- "Pa"
En_NCBI_info$Tax <- "En"

NCBI_info <- rbind(Ef_NCBI_info, Sa_NCBI_info,
                   Kp_NCBI_info, Ab_NCBI_info,
                   Pa_NCBI_info, En_NCBI_info)

drep_res <- read.csv("Widb.csv")
drep_res$genome <- gsub(".fa","",drep_res$genome)
drep_res <- drep_res[,c(1,3,4)]

library(dplyr)
ESKAPE_info <- drep_res %>%
  left_join(NCBI_info,
            by = c("genome" = "Assembly.Accession"))

write.table(ESKAPE_info,
            "ESKAPE_genome_info.tsv",
            sep = "\t",           # TSV使用制表符分隔
            row.names = FALSE,    # 不写入行名
            col.names = TRUE,     # 写入列名
            quote = FALSE        # 不添加引号
            )            # NA值的表示方式

ESKAPE_info <- read.delim2("ESKAPE_genome_info.tsv")

######################################################################

library(stringr)
library(dplyr)
library(tidyr)

parse_biosample_tidy <- function(file_path) {
  # 读取文件
  content <- readLines(file_path)
  
  # 找到样本开始的位置
  sample_starts <- which(grepl("^\\d+: ", content))
  sample_ends <- c(sample_starts[-1] - 1, length(content))
  
  # 处理每个样本
  samples <- lapply(seq_along(sample_starts), function(i) {
    start <- sample_starts[i]
    end <- sample_ends[i]
    sample_text <- content[start:end]
    
    # 提取信息
    bio_sample <- str_extract(sample_text[str_detect(sample_text, "BioSample:")], 
                              "BioSample:\\s*([^;]+)") %>%
      str_remove("BioSample:\\s*") %>%
      trimws() %>%
      ifelse(length(.) == 0, NA, .)
    
    # 如果没有找到BioSample，尝试从Accession行提取
    if (is.na(bio_sample) || length(bio_sample) == 0) {
      bio_sample <- str_extract(sample_text[str_detect(sample_text, "Accession:")], 
                                "Accession:\\s*([^\\t]+)") %>%
        str_remove("Accession:\\s*") %>%
        trimws() %>%
        ifelse(length(.) == 0, NA, .)
    }
    
    # 提取latitude and longitude
    lat_lon_line <- sample_text[str_detect(sample_text, "/latitude and longitude")]
    if (length(lat_lon_line) > 0) {
      lat_lon <- str_extract(lat_lon_line, '="([^"]+)"') %>%
        str_remove_all('^="|"$')
    } else {
      lat_lon <- "not determined"
    }
    
    # 分离latitude和longitude（如果可能）
    if (lat_lon != "not determined") {
      # 尝试解析坐标格式
      if (str_detect(lat_lon, "[NSEW]")) {
        # 包含方向信息
        latitude <- ifelse(str_detect(lat_lon, "[NS]"), lat_lon, "not determined")
        longitude <- ifelse(str_detect(lat_lon, "[EW]"), lat_lon, "not determined")
      } else if (str_detect(lat_lon, ",")) {
        # 逗号分隔
        parts <- str_split(lat_lon, ",", simplify = TRUE)
        latitude <- trimws(parts[1])
        longitude <- ifelse(ncol(parts) > 1, trimws(parts[2]), "not determined")
      } else {
        latitude <- lat_lon
        longitude <- lat_lon
      }
    } else {
      latitude <- "not determined"
      longitude <- "not determined"
    }
    
    # 提取geographic location
    geo_line <- sample_text[str_detect(sample_text, "/geographic location")]
    if (length(geo_line) > 0) {
      geographic_location <- str_extract(geo_line, '="([^"]+)"') %>%
        str_remove_all('^="|"$')
    } else {
      geographic_location <- "not determined"
    }
    
    # 返回数据框
    data.frame(
      BioSample = ifelse(is.na(bio_sample), "not determined", bio_sample),
      latitude = latitude,
      longitude = longitude,
      geographic_location = geographic_location,
      stringsAsFactors = FALSE
    )
  })
  
  # 合并所有样本
  result_df <- bind_rows(samples)
  
  # 移除完全重复的行
  result_df <- distinct(result_df)
  
  return(result_df)
}

# 使用函数
result_df <- parse_biosample_tidy("biosample_result.txt")
print(result_df)


library(stringr)

format_coordinates_robust <- function(df) {
  
  # 预处理：清理数据
  df$coord_clean <- df$latitude
  
  # 标记无法解析的格式
  df$parsed <- FALSE
  df$lat_dec <- NA
  df$lon_dec <- NA
  df$coord_type <- NA
  
  for (i in 1:nrow(df)) {
    coord <- df$coord_clean[i]
    
    if (is.na(coord) || 
        coord %in% c("not determined", "missing", "-", "", "NA")) {
      df$coord_type[i] <- "missing"
      next
    }
    
    # 类型1: 纬度/经度格式 (-37.8462/144.10096)
    if (grepl("^\\-?\\d+\\.\\d+/\\-?\\d+\\.\\d+$", coord)) {
      parts <- str_split(coord, "/", simplify = TRUE)
      df$lat_dec[i] <- as.numeric(parts[1])
      df$lon_dec[i] <- as.numeric(parts[2])
      df$coord_type[i] <- "lat_lon_slash"
      df$parsed[i] <- TRUE
      next
    }
    
    # 类型2: 带方向的十进制格式 (0.0515 N 37.6456 E)
    # 匹配模式: 数字 方向 数字 方向
    dir_pattern <- "([\\-\\d\\.]+)\\s*([NS])\\s+([\\-\\d\\.]+)\\s*([EW])"
    dir_match <- str_match(coord, dir_pattern)
    
    if (!any(is.na(dir_match))) {
      lat_val <- as.numeric(dir_match[1,2])
      lat_dir <- dir_match[1,3]
      lon_val <- as.numeric(dir_match[1,4])
      lon_dir <- dir_match[1,5]
      
      # 根据方向调整
      if (lat_dir == "S") lat_val <- -lat_val
      if (lon_dir == "W") lon_val <- -lon_val
      
      df$lat_dec[i] <- lat_val
      df$lon_dec[i] <- lon_val
      df$coord_type[i] <- "decimal_with_direction"
      df$parsed[i] <- TRUE
      next
    }
    
    # 类型3: 可能缺少空格 (0.18 N 32.34 E)
    # 尝试更灵活的模式
    relaxed_pattern <- "([\\-\\d\\.]+)[\\s°]*([NS])?[\\s,]+([\\-\\d\\.]+)[\\s°]*([EW])?"
    relaxed_match <- str_match(coord, relaxed_pattern)
    
    if (!any(is.na(relaxed_match))) {
      lat_val <- as.numeric(relaxed_match[1,2])
      lat_dir <- relaxed_match[1,3]
      lon_val <- as.numeric(relaxed_match[1,4])
      lon_dir <- relaxed_match[1,5]
      
      # 如果有方向信息，调整
      if (!is.na(lat_dir) && lat_dir == "S") lat_val <- -lat_val
      if (!is.na(lon_dir) && lon_dir == "W") lon_val <- -lon_val
      
      # 验证合理性
      if (abs(lat_val) <= 90 && abs(lon_val) <= 180) {
        df$lat_dec[i] <- lat_val
        df$lon_dec[i] <- lon_val
        df$coord_type[i] <- "relaxed_decimal"
        df$parsed[i] <- TRUE
        next
      }
    }
    
    # 类型4: 尝试提取所有数字
    numbers <- as.numeric(str_extract_all(coord, "\\-?\\d+\\.?\\d*")[[1]])
    
    if (length(numbers) >= 2) {
      # 尝试不同的组合
      possible_pairs <- list(
        c(numbers[1], numbers[2]),
        c(numbers[2], numbers[1])
      )
      
      for (pair in possible_pairs) {
        lat_candidate <- pair[1]
        lon_candidate <- pair[2]
        
        # 检查是否在合理范围内
        if (abs(lat_candidate) <= 90 && abs(lon_candidate) <= 180) {
          df$lat_dec[i] <- lat_candidate
          df$lon_dec[i] <- lon_candidate
          df$coord_type[i] <- "extracted_numbers"
          df$parsed[i] <- TRUE
          break
        }
      }
      
      if (df$parsed[i]) next
    }
    
    # 无法解析
    df$coord_type[i] <- "unparsed"
  }
  
  # 统计信息
  type_counts <- table(df$coord_type)
  cat("坐标解析统计:\n")
  for (type in names(type_counts)) {
    cat(sprintf("  %-25s: %d (%.1f%%)\n", 
                type, type_counts[type], type_counts[type]/nrow(df)*100))
  }
  
  cat(sprintf("\n总成功解析: %d/%d (%.1f%%)\n",
              sum(df$parsed), nrow(df), sum(df$parsed)/nrow(df)*100))
  
  return(df)
}

# 使用函数
result_formatted <- format_coordinates_robust(result_df)

result_formatted <- result_formatted[,c(1,7,8,4)]

# 方法1B：根据某些条件选择重复项（例如，非NA值优先）
result_formatted_unique <- result_formatted %>%
  group_by(BioSample) %>%
  arrange(desc(!is.na(lat_dec)), desc(!is.na(lon_dec))) %>%
  slice(1) %>%
  ungroup()
# 
ESKAPE_info <- ESKAPE_info %>%
  left_join(result_formatted_unique, by = c("BioSample" = "BioSample"))

redundancy_count <- ESKAPE_info %>%
  select(Tax, lon_dec, lat_dec) %>%   # 只保留最后三列
  group_by(Tax, lon_dec, lat_dec) %>% # 按三列分组
  summarise(genome_number = n(), .groups = "drop")             # 计数

#参考https://cloud.tencent.com/developer/article/2206435
library(ggplot2)
world<-map_data("world")
world.map <- 
  ggplot() +
  geom_polygon(data=world,aes(x=long,y=lat,group=group),
               fill="#e9ecef")+
  geom_point(data=redundancy_count, 
             aes(x=lon_dec,y=lat_dec, 
                 fill = Tax, 
                 size = genome_number,
                 shape = Tax),
#             shape = 21,      # 带边框的填充圆
             colour = "black",# 边框颜色
             stroke = 0.1     # 边框宽度（mm）
  )+
  scale_size(range = c(2,6))+
  scale_fill_manual(values = c("Ef"="#9f86c0","Sa"="#FFBE7D",
                               "Kp"="#76c893","Ab"="#F0E442","Pa"="#ff99c8",
                               "En"="#56B4E9")
                    ) +
  scale_shape_manual(values = c(
    "Ef" = 21,  # 填充圆
    "Sa" = 22,  # 填充方形
    "Kp" = 23,  # 填充菱形
    "Ab" = 24,  # 填充三角形
    "Pa" = 25,  # 填充倒三角形
    "Host_associated" = 8,  # 星形
    "En" = 22   # 方形
  )) +
  theme_bw()+
  scale_y_continuous(expand = expansion(mult=c(0,0)))+
  scale_x_continuous(expand = expansion(add=c(0,0)))+
  theme(legend.position = 'none')
world.map
#2252个基因组的坐标是NA

# 绘制横条形图
MAG_number <- table(ESKAPE_info$Tax) %>% as.data.frame()

MAG_number$Var1 <- factor(MAG_number$Var1, levels = c("Ef","Sa","Kp","Ab","Pa","En"))

mytheme <-   theme(
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

p_MAGnumber <- ggplot(MAG_number, 
                              aes(x = Var1 ,y = Freq)) +
  geom_bar(stat = "identity", aes(fill = Var1),width = 0.5) +  # 使用stat = "identity"因为我们的值已经是计数或比例
  #coord_flip() +  # 翻转坐标轴，使条形图水平
  #scale_x_discrete(limits = rev) +  # X轴类别逆序
  # scale_fill_gradient(low = "#FFF1EB",
  #                     high = "#ACE0F9")+
  # scale_fill_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
  #                              "LCJ"="#FFE2CE","s"= "#C4D8E9",
  #                              "CJ"= "#BEBCDF","HH"="#FB8C62","Ⅶ"="#FFF3CA"
  # )) +  # 手动设置颜色
  scale_fill_manual(values = c("Ef"="#9f86c0","Sa"="#FFBE7D",
                               "Kp"="#76c893","Ab"="#F0E442","Pa"="#ff99c8",
                               "En"="#56B4E9")
                    )+
  # scale_fill_viridis_c(option = "inferno") + ggtitle("'inferno'")+
  #  geom_vline(xintercept = 4.5, linetype = "dashed", color = "black")+
  #  geom_vline(xintercept = 5.5, linetype = "dashed", color = "black")+
  geom_text(aes(label = round(Freq,2)),vjust = 1,hjust = 0.5,size = 4)+
  labs(x = "", y = "Number of DS type")+  # 设置标签和标题
  mytheme
p_MAGnumber

#######################################################################

##2.2 基因组大小及gc含量####
# genome_size 的复合图
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(ggsignif)#geom_signif
library(gghalves)
library(ggsci)

genome_checkm2 <- read.delim2("quality_report.tsv")
genome_checkm2$Name <- gsub("_renamed","",genome_checkm2$Name)
genome_checkm2$Name <- gsub("\\.","_",genome_checkm2$Name)

genome_checkm2$genome_size <- genome_checkm2$Genome_Size / 1000000
ESKAPE_info$genome <- gsub("\\.","_",ESKAPE_info$genome)

matches <- match(genome_checkm2$Name, ESKAPE_info$genome)
genome_checkm2 <- genome_checkm2 %>%
  mutate(
    Tax = ifelse(
      !is.na(matches),  # 如果匹配成功
      ESKAPE_info$Tax[matches],  # 使用匹配的Tax值
      Host_tax  # 否则保持原样
    )
  )

genome_checkm2$Tax <- factor(genome_checkm2$Tax, levels = c("Ef","Sa","Kp","Ab","Pa","En"))

# p1 <- ggplot(all_4649_genome, aes(x = phylum, y = genome_size, fill = phylum)) +
#   geom_violin(alpha = 0.6, trim = FALSE) +
#   geom_boxplot(width = 0.6, alpha = 0.8, outlier.shape = NA) +
p1 <-
  ggplot(genome_checkm2, 
         aes(x = Tax, y = genome_size)) +
  geom_half_boxplot(aes(fill = Tax), color ="black",
                    side = "l",
                    errorbar.draw = T,
                    outlier.shape = NA,
                    width=0.8
  )+
  geom_half_point(aes(color = Tax),
                  side = "r",size = 0.5,
                  transformation_params = list(height = 0,width = 0.001,seed = 2))+
  scale_y_continuous(labels = scales::comma) +  # 使用逗号分隔大数字
  labs(
    x = "Taxonomy",
    y = "Genome Size (Mbp)") +
  scale_fill_manual(values = c("Ef"="#9f86c0","Sa"="#FFBE7D",
                               "Kp"="#76c893","Ab"="#F0E442","Pa"="#ff99c8",
                               "En"="#56B4E9")
  )+
  scale_color_manual(values = c("Ef"="#9f86c0","Sa"="#FFBE7D",
                               "Kp"="#76c893","Ab"="#F0E442","Pa"="#ff99c8",
                               "En"="#56B4E9")
  )+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +  # 移除图例，因为x轴已有标签
  stat_compare_means(method = "kruskal.test", label = "p.format", 
                     label.x = 1, 
                     label.y = max(genome_checkm2$genome_size, na.rm = TRUE) * 0.9)+
  theme_minimal() +  # 使用简洁主题
  theme(legend.position = "none")+  # 图例放在底部
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    #panel.border = element_blank(),  # 去掉边界线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    #axis.line = element_line(color = "grey20", size = 0.5),  # 设置坐标轴线条为深灰色，更加细致
    axis.title = element_text(face = "plain", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.title.x = element_blank(),
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0,vjust = 1,hjust = 0.5), #横轴旋转
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = "none",  # 设置图例位置在右侧
    legend.title = element_text(face = "plain", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )

print(p1)

genome_checkm2$GC_Content <- as.numeric(genome_checkm2$GC_Content)
genome_checkm2$Coding_Density <- as.numeric(genome_checkm2$Coding_Density)
genome_checkm2$Completeness <- as.numeric(genome_checkm2$Completeness)
# gc_percentage 的复合图
p2 <- ggplot(genome_checkm2, 
             aes(x = Tax, y = GC_Content)) +
  geom_half_boxplot(aes(fill = Tax), color ="black",
                    side = "l",
                    errorbar.draw = T,
                    outlier.shape = NA,
                    width=0.8
  )+
  geom_half_point(aes(color = Tax),
                  side = "r",size = 0.5,
                  transformation_params = list(height = 0,width = 0.001,seed = 2))+
  labs(
    x = "Taxonomy",
    y = "GC content") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  stat_compare_means(method = "kruskal.test", label = "p.format",
                     label.x = 1, 
                     label.y = max(genome_checkm2$GC_Content, na.rm = TRUE) * 0.9)+
  scale_fill_manual(values = c("Ef"="#9f86c0","Sa"="#FFBE7D",
                               "Kp"="#76c893","Ab"="#F0E442","Pa"="#ff99c8",
                               "En"="#56B4E9")
  )+
  scale_color_manual(values = c("Ef"="#9f86c0","Sa"="#FFBE7D",
                                "Kp"="#76c893","Ab"="#F0E442","Pa"="#ff99c8",
                                "En"="#56B4E9")
  )+
  theme_minimal() +  # 使用简洁主题
  theme(legend.position = "none")+  # 图例放在底部
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    #panel.border = element_blank(),  # 去掉边界线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    #axis.line = element_line(color = "grey20", size = 0.5),  # 设置坐标轴线条为深灰色，更加细致
    axis.title = element_text(face = "plain", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.title.x = element_blank(),
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0,vjust = 1,hjust = 0.5), #横轴旋转
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = "none",  # 设置图例位置在右侧
    legend.title = element_text(face = "plain", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )

print(p2)

# gc_percentage 的复合图
p3 <- ggplot(genome_checkm2, 
             aes(x = Tax, y = Coding_Density)) +
  geom_half_boxplot(aes(fill = Tax), color ="black",
                    side = "l",
                    errorbar.draw = T,
                    outlier.shape = NA,
                    width=0.8
  )+
  geom_half_point(aes(color = Tax),
                  side = "r",size = 0.5,
                  transformation_params = list(height = 0,width = 0.001,seed = 2))+
  labs(
    x = "",
    y = "Coding density") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  stat_compare_means(method = "kruskal.test", label = "p.format",
                     label.x = 1, 
                     label.y = max(genome_checkm2$Coding_Density, na.rm = TRUE) * 0.9)+
  scale_fill_manual(values = c("Ef"="#9f86c0","Sa"="#FFBE7D",
                               "Kp"="#76c893","Ab"="#F0E442","Pa"="#ff99c8",
                               "En"="#56B4E9")
  )+
  scale_color_manual(values = c("Ef"="#9f86c0","Sa"="#FFBE7D",
                                "Kp"="#76c893","Ab"="#F0E442","Pa"="#ff99c8",
                                "En"="#56B4E9")
  )+
  theme_minimal() +  # 使用简洁主题
  theme(legend.position = "none")+  # 图例放在底部
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    #panel.border = element_blank(),  # 去掉边界线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    #axis.line = element_line(color = "grey20", size = 0.5),  # 设置坐标轴线条为深灰色，更加细致
    axis.title = element_text(face = "plain", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.title.x = element_blank(),
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0,vjust = 1,hjust = 0.5), #横轴旋转
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = "none",  # 设置图例位置在右侧
    legend.title = element_text(face = "plain", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
print(p3)


mytheme <-   theme(
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


library(Cairo)
library(cowplot)
#cairo_pdf(filename = file.path(CS_WORK, "1_genome_metadata.pdf") , width = 8,height = 4,unit = 'in')
g1 <- plot_grid(p1,p2,p3,p_MAGnumber,
                labels = LETTERS[3:6], ncol = 1
                #rel_heights = c(1, 1, 1.5),
                # label_x = 0.09,            # 标签x坐标（1为最右侧）
                #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g1


ggsave(file.path(CS_WORK, "1_genome_info.pdf"),g1,width = 3,height = 7,units = 'in')

ggsave(file.path(CS_WORK, "1_genome_globalmap.pdf"),world.map,width = 5,height = 3.5,units = 'in')











