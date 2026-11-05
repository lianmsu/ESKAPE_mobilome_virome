# Software and environment

The production runs used the conda environments named in each script (`source activate ...`).
The table below lists the software, the environment it ran in, and where it is used.

## Bioinformatics software

| Tool | Environment (as activated in the scripts) | Used in |
| --- | --- | --- |
| dRep | `drep` | `01_genome_catalogue` |
| CheckM2 | `checkm2` | `01_genome_catalogue` |
| GTDB-Tk | `gtdbtk-2.5.2`, `gtdbtk2.3` (GTDB R226 / R214 databases) | `01_genome_catalogue`, `10_experiments` |
| minced | `minced` | `02_crispr_and_defense` |
| CRISPRCasFinder | `crisprcasfinder` | `02_crispr_and_defense` |
| CRISPRCasTyper | `cctyper` | `02_crispr_and_defense` |
| DefenseFinder | `defensefinder` | `02_crispr_and_defense`, `06_virus_annotation` |
| RGI / CARD | `RGI_CARD` | `03_resistome_and_bgc`, `06_virus_annotation` |
| diamond | `eggnog` | `03_resistome_and_bgc` |
| Prodigal | `prokka`, `instrain` | `03_resistome_and_bgc`, `06_virus_annotation` |
| Macrel | `macrel` | `03_resistome_and_bgc` |
| antiSMASH | `antiSMASH` | `03_resistome_and_bgc` |
| IntegronFinder | `integronfinder` | `04_mobile_genetic_elements` |
| ICEfinder2 | `icefinder2` | `04_mobile_genetic_elements` |
| ISEScan | `isescan` | `04_mobile_genetic_elements` |
| geNomad | `genomad` | `05_virus_discovery` |
| CheckV | `checkv` | `05_virus_discovery` |
| mmseqs2 | `mmseqs2` | `05_virus_discovery` |
| iPHoP | `iphop` | `06_virus_annotation` |
| BACPHLIP | `bacphlip` | `06_virus_annotation` |
| hhsuite2 | `hhsuite` | `06_virus_annotation`, `09_phylogeny` |
| eggNOG-mapper | `eggnog` | `07_spacer_competition` |
| BLAST+ | `blast` | `07_spacer_competition`, `09_phylogeny` |
| vContact3 / vContact2 | `vcontact3`, `vcontact2` | `08_phage_clustering_and_host` |
| comparem | `comparem` | `08_phage_clustering_and_host` |
| clinker | `clinker` | `08_phage_clustering_and_host` |
| ViPTreeGen | `viptree` | `08_phage_clustering_and_host` |
| IQ-TREE | `iqtree` | `09_phylogeny` |
| seqkit | `seqkit` | throughout |

Where a script does not pin a version, the version is the one installed in that conda
environment; the environment names above are exactly those used in the production runs.

## R

The analyses were run with R 4.3 and require:

```
ape            complexheatmap  dplyr          ggalluvial     ggbreak
ggcorrplot     ggforce         ggpubr         ggraph         ggplot2
ggrepel        ggsankey       ggsignif       ggvenn         grid
gridExtra      igraph         patchwork      pheatmap       plyr
RColorBrewer   readxl         reshape2       rstatix        scales
stringr        tibble         tidyr          tidyverse      vegan
viridis        wesanderson
```

Install them with `install.packages()` (Bioconductor for `ComplexHeatmap` if used) and check the
package list at the top of each script.

## Databases

| Database | Used for |
| --- | --- |
| NCBI RefSeq | the 11 947 complete ESKAPE genomes |
| GTDB (R226 metadata, R214 mash for GTDB-Tk 2.3) | taxonomy and trees |
| CARD | RGI annotation of ARGs |
| SARG | diamond annotation of ARGs |
| BacMet | metal resistance genes (MRG) |
| VFDB | virulence factors (VFG) |
| geNomad database | virus/plasmid prediction |
| CheckV database | virus quality |
| IMG/VR | spacer targets and virus source |
| VIRE | spacer targets and virus source |
| Virus-Host DB | host information of reference viruses |
| antiSMASH / Macrel models | BGC and AMP prediction |

## Notes

* The scripts were written for a SLURM cluster: the `#SBATCH` headers of the production runs are
  kept. Queue names (`cpu`, `fat`, `cu24`, `cu08`, ...) and walltimes are those of the original
  cluster and must be adapted.
* `2ESKAPE_MGE` (the R working copy) and the two cluster copies (`lfat1t01`, `pkusz_CPU`) overlap
  in a few places; the module headers record which copy a script came from.
* Files in directories whose name ended in `outdated` or `wrong` were excluded on purpose.
