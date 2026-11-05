# Pipeline overview

The analysis has four layers: a genome catalogue, its functional annotation (resistome,
defensome, MGEs), the construction and annotation of the virome, and the phage-host interaction
analyses. The R scripts turn the resulting tables into the figures.

## Stage by stage

### 1. Genome catalogue (`pipelines/01_genome_catalogue`)

| Step | Script | Input | Output |
| --- | --- | --- | --- |
| Download and rename | `01_dereplication/1genome/1rename_genome.sh` | NCBI RefSeq assemblies | `all_4538_ESKAPE_genomes_renamed.fa` |
| Dereplication | `01_dereplication/2drep/drep_MAGs.sh`, `2drep/drep_MAGs.sh` | 11 947 genomes | 4 538 genomes (99.9 % ANI, 50 % coverage) |
| Quality | `02_quality_and_taxonomy/3checkm2/checkm2.sh` | dereplicated genomes | `quality_report.tsv` |
| Taxonomy | `02_quality_and_taxonomy/4GTDB/GTDB_infer.sh`, `3gtdbtk/01GTDB_classify.sh` | dereplicated genomes | GTDB-Tk summary and tree |

### 2. CRISPR arrays and defense systems (`pipelines/02_crispr_and_defense`)

| Step | Script | Purpose |
| --- | --- | --- |
| Spacer extraction | `01_minced/22minced/1minced.sh`, `25level4_minced/1minced.sh` | minced CRISPR spacers on the genomes and on the level-4 arrays |
| CRISPR-Cas detection | `02_crisprcasfinder/23CRISPRCasFinder/1CRISPRCasFinder.sh` | CRISPRCasFinder `-levelMin 4` |
| Subtyping | `03_crisprcaster/24CRISPRCasTyper/1CRISPRTyper.sh` | CRISPRCasTyper on the Cas-containing contigs |
| Defense systems | `04_defensefinder/30Defensefinder/01DefenseFinder.sh` | DefenseFinder on the genomes |

### 3. Resistome and secondary metabolites (`pipelines/03_resistome_and_bgc`)

| Step | Script | Purpose |
| --- | --- | --- |
| ARG annotation | `01_resistome_annotation/31ARG/1RGI_CARD.sh`, `2mag_orf_diamond_li.sh` | RGI/CARD and diamond against SARG |
| AMP prediction | `02_amp_macrel/32AMP/*.sh` | Prodigal + Macrel |
| BGC prediction | `03_antismash/33antismash/*.sh` | antiSMASH on the merged genomes |

### 4. Mobile genetic elements (`pipelines/04_mobile_genetic_elements`)

| Step | Script | Purpose |
| --- | --- | --- |
| Integrons | `01_ice_integron_isescan/4integronfinder/01integron_finder.sh` | IntegronFinder |
| ICEs | `01_ice_integron_isescan/5ICEfinder/icefinder_*.sh` | ICEfinder2 (four batches + full set) |
| ISs | `01_ice_integron_isescan/6ISEscan/ISEScan_4538MAGs.sh` | ISEScan |
| Special MGEs | `02_special_mge/34special_MGE/*`, `5_ICEfinder_special_MGE/*` | the nested MGEs used for the case studies |

### 5. Virus discovery (`pipelines/05_virus_discovery`)

| Step | Script | Purpose |
| --- | --- | --- |
| Virus prediction | `01_genomad/21genomad/01genomad_lpw.sh`, `45virus_genomad/…` | geNomad on the genomes |
| Quality | `02_checkv/42genomad_checkv/01checkv_lpw.sh`, `45virus_checkv/…` | CheckV |
| Clustering | `03_clustering/44virus_checkv_cluster/0checkv_cluster.sh`, `44virus_mmseqs/2mmseqs.sh` | ANI clustering and mmseqs2 |
| Catalogue | `04_final_catalogue/46final_virus/1virus_filter.sh` | length/score/hallmark rules that define the 14 496 viruses |

### 6. Virus annotation (`pipelines/06_virus_annotation`)

| Step | Script | Purpose |
| --- | --- | --- |
| Host and lifestyle | `01_host_prediction/46iPHoP/1_iPHoP_lpw.sh`, `3_bacphlip_lpw.sh` | iPHoP, BACPHLIP |
| Features | `02_phage_features/46phage_features/1_molecular_characteristics.sh` (+ python helpers) | GC, NARSC, length |
| Location | `03_location_and_source/46phage_location_source/…` | ecosystem and geography |
| Defense | `04_virus_defense/47virus_defensefinder/01DefenseFinder.sh` | DefenseFinder on the viruses |
| Resistome | `05_virus_resistome/48virus_RGI/1RGI_CARD.sh` | RGI/CARD on the viruses |
| References | `06_virus_and_host_databases/*` | Virus-Host DB and VIRE subsets |
| Selected viruses | `07_selected_virus_annotation/49selected_virus_annotation/*` | hhsuite2 annotation of the novel viruses |

### 7. CRISPR spacer competition (`pipelines/07_spacer_competition`)

| Step | Script | Purpose |
| --- | --- | --- |
| MGE competition | `01_mge_crispr_spacers/35MGE_spacer_competition/1CRISPRTyper*.sh` | spacers per MGE type |
| IMG/VR | `02_spacer_to_imgvr/40spacer2IMGVR/*` | makeblastdb + blastn + filtering |
| VIRE | `03_spacer_to_vire/41spacer2VIRE/*` | same for VIRE |
| Target genes | `04_spacer_target_genes/54spacer_target_gene/*` | gene extraction, eggNOG-mapper, methylase tree |

### 8. Phage clustering and host range (`pipelines/08_phage_clustering_and_host`)

| Step | Script | Purpose |
| --- | --- | --- |
| vContact3 | `01_vcontact3/50vcontact3/1vcontact3*.sh` | genome-sharing network |
| vContact2 | `02_vcontact2/51vcontact2/*` | companion clustering |
| Novel viruses | `03_novel_virus/52potential_novel_virus/*` | comparem AAI, clinker maps, ViPTree |
| Multi-host phages | `04_multi_host_phage/53multi_host_phage/*` | host range and links |

### 9. Phylogeny (`pipelines/09_phylogeny`)

| Step | Script | Purpose |
| --- | --- | --- |
| Marker genes | `01_marker_genes/60marker_gene_phylogeny/**` | hhsuite2 annotation, IQ-TREE for MCP, portal and TerL |
| Special genes | `02_special_genes/61special_gene_phylogeny/*` | blast + IQ-TREE for e.g. VanX |
| Methylases | `03_methylase/54spacer_target_gene/methylase_phylogeny/*` | methyltransferase tree |

### 10. Experiments (`pipelines/10_experiments`)

Codon usage and sequence analysis of the isolated phages (`apc.pl`/`apc.py`), GTDB-Tk trees of
the host strains, and the defense/prophage content of the experimental hosts.

## From tables to figures

| Figure | Content | Analyses |
| --- | --- | --- |
| Fig. 1 | Genome map, genome traits, ARG/DS prevalence, trait relationships | `analyses/01_genome_landscape`, `analyses/03_gene_annotation` |
| Fig. 2 | MGE counts, spacer competition, ARG/anti-defense cargo, nested MGEs | `analyses/04_mobilome` |
| Fig. 3 | Virome quality, taxonomy, host, lifestyle, defense vs anti-defense | `analyses/02_phage_catalogue`, `analyses/06_phage_defense` |
| Fig. 4 | IMG/VR and VIRE distribution, spacer-targeted viral genes | `analyses/07_phage_host` |
| Fig. 5 | Whole-proteome phylogeny and AAI of the novel groups | `analyses/05_phage_clustering` |
| Fig. 6 | Genome comparison, gene phylogenies, structural models | `analyses/08_phylogeny` |
