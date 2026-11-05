# ESKAPE_mobilome_virome

本仓库是论文 **《Functional profiling of ESKAPE viromes uncovers resistance-limiting
phage-host dynamics》**（Li P., Liu Q., Deng C., Ni J.，bioRxiv 2026，
[doi:10.64898/2026.08.16.745060](https://doi.org/10.64898/2026.08.16.745060)）的配套代码。

仓库包含从公开基因组与宏基因组数据到论文图件的全部脚本：ESKAPE 基因组集的构建与
抗性组/防御组注释、可移动遗传元件（MGE）调查、ESKAPE 病毒组的构建与注释、基于 CRISPR
间隔序列的互作分析，以及产出统计结果和图件的 R 分析脚本。

> English version: **[README.md](README.md)**

---

## 研究概要

* **基因组集**：从 NCBI RefSeq 下载 11 947 个 ESKAPE 完整基因组，以 99.9% ANI（覆盖度 30%）
  去冗余，得到 **4 538 个非冗余基因组**，覆盖 *A. baumannii*、*E. faecium*、*Enterobacter*
  spp.、*K. pneumoniae*、*P. aeruginosa* 和 *S. aureus*。
* **可移动元件**：共鉴定 35 312 个 MGE——10 199 个质粒、18 411 个前噬菌体、4 799 个
  ICE/IME、1 594 个整合子、309 个插入序列和 84 个复合转座子。嵌套结构（转座子位于 ICE、
  ICE 位于质粒、整合子位于前噬菌体等）使低移动性的 ARG 获得了水平转移能力。
* **病毒组**：将基因组中的前噬菌体与 IMG/VR、VIRE 整合，得到 **14 496 个 ESKAPE 相关病毒**
  的非冗余目录，包含 4 个此前未报道的病毒类群。ESKAPE 噬菌体几乎不携带 ARG，并与质粒
  相互拮抗，即"限制抗性"的一面。
* **军备竞赛**：烈性噬菌体富集抗防御基因；新病毒携带的甲基转移酶本身又是宿主 CRISPR-Cas
  系统的靶标。

## 仓库结构

```
ESKAPE_mobilome_virome/
├── pipelines/          # 上游生物信息流程（SLURM 作业脚本与 python 辅助脚本）
│   ├── 01_genome_catalogue/            下载、重命名、dRep、CheckM2、GTDB-Tk
│   ├── 02_crispr_and_defense/          minced、CRISPRCasFinder、CRISPRCasTyper、DefenseFinder
│   ├── 03_resistome_and_bgc/           RGI/CARD、diamond-SARG、Macrel、antiSMASH
│   ├── 04_mobile_genetic_elements/     IntegronFinder、ICEfinder2、ISEScan、特殊 MGE
│   ├── 05_virus_discovery/             geNomad、CheckV、聚类、最终病毒目录
│   ├── 06_virus_annotation/            iPHoP/BACPHLIP、噬菌体特征、DefenseFinder、RGI
│   ├── 07_spacer_competition/          MGE 间隔序列竞争、spacer 比对 IMG/VR 与 VIRE
│   ├── 08_phage_clustering_and_host/   vContact3、vContact2、新病毒、多宿主噬菌体
│   ├── 09_phylogeny/                   标志基因（MCP、portal、TerL）与特殊基因系统发育
│   └── 10_experiments/                 湿实验噬菌体序列与宿主菌株
├── analyses/           # 下游统计与绘图（R，另含一个 Perl 脚本）
│   ├── 01_genome_landscape/            Fig. 1A-B
│   ├── 02_phage_catalogue/             Fig. 3A-E
│   ├── 03_gene_annotation/             Fig. 1C-D
│   ├── 04_mobilome/                    Fig. 2
│   ├── 05_phage_clustering/            Fig. 5A-B
│   ├── 06_phage_defense/               Fig. 3F
│   ├── 07_phage_host/                  Fig. 4
│   ├── 08_phylogeny/                   Fig. 6
│   └── 09_experiments/                 分离噬菌体的密码子使用
├── data/               # 小型基础表格（见 data/README.md）
└── docs/               # 流程、目录、文件对应、输入输出、软件环境
```

每个脚本都保留了它在服务器工作目录中的子目录名（例如
`pipelines/05_virus_discovery/04_final_catalogue/46final_virus/1virus_filter.sh`），因此每一步都能
回溯到原始工作目录。每个文件头部都写明了用途、对应的论文图件以及它在服务器上的原始位置。

## 文档

| 文档 | 内容 |
| --- | --- |
| [docs/01_pipeline_overview.md](docs/01_pipeline_overview.md) | 端到端流程，逐阶段列出输入与输出 |
| [docs/02_directory_layout.md](docs/02_directory_layout.md) | 目录说明、命名约定、环境变量 |
| [docs/03_file_mapping.md](docs/03_file_mapping.md) | 每个脚本：服务器原路径 → 新路径 |
| [docs/04_script_inputs_outputs.md](docs/04_script_inputs_outputs.md) | 各模块的输入、输出与对应图件 |
| [docs/05_software_and_environment.md](docs/05_software_and_environment.md) | 软件版本、conda 环境、R 包、数据库 |

## 复现

**1. 上游流程。** `pipelines/` 下的脚本就是实际使用的作业脚本（保留了 `#SBATCH` 头），
除路径外与生产运行一致：所有站点相关的绝对路径都改成了变量，每个脚本开头都有说明。

```bash
export PROJECT_ROOT=/path/to/08ESKAPE      # 项目目录
export HOME_ROOT=/path/to/home             # 同一账号下的其他项目
export DB_ROOT=/path/to/databases          # 参考数据库
```

按 [docs/01](docs/01_pipeline_overview.md) 给出的顺序提交；conda 环境见
[docs/05](docs/05_software_and_environment.md)。

**2. R 分析。** `analyses/` 下的脚本需要 R >= 4.2 以及
[docs/05](docs/05_software_and_environment.md) 中列出的 R 包。输入表通过 `CS_DATA` 定位，
图件写入 `CS_WORK`，默认指向仓库内目录：

```r
Sys.setenv(CS_DATA = "/path/to/data", CS_WORK = "/path/to/results")
source("analyses/01_genome_landscape/01_genome_metadata_and_map.R")
```

`CS_EXT` 指向从姊妹项目（捕食性细菌病毒组）导入的表格（`data/external`），部分章节会复用。

## 数据可用性

基因组目录、MGE 表格和 ESKAPE 病毒组以补充表格形式随预印本提供。前几个分析所需的小型
基础表已随仓库提供，见 [`data/`](data/README.md)；较大的中间结果（geNomad/CheckV 输出、
vContact3 结果、spacer 比对表）可用 `pipelines/` 中的脚本重新生成。原始数据均为公开数据：
基因组来自 NCBI RefSeq，病毒参考序列来自 IMG/VR、VIRE 和 Virus-Host DB。

## 引用

```bibtex
@article{Li2026eskapevirome,
  title   = {Functional profiling of ESKAPE viromes uncovers resistance-limiting
             phage-host dynamics},
  author  = {Li, Pengwei and Liu, Quan and Deng, Chunfang and Ni, Jinren},
  journal = {bioRxiv},
  year    = {2026},
  doi     = {10.64898/2026.08.16.745060}
}
```

使用本仓库脚本时请引用上述预印本。
