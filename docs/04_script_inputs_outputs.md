# Script inputs, outputs and section map

Cluster paths below are relative to `${PROJECT_ROOT}` (the 08ESKAPE directory); the R paths are
relative to `CS_DATA` (`data/` by default).

## pipelines

### `01_genome_catalogue`

* **Input** — NCBI RefSeq assemblies (`ncbi_*` folders, `all_genomes`), GTDB metadata.
* **Output** — `all_4538_ESKAPE_genomes_renamed.fa` (4 538 genomes), `res_ESKAPE`
  (dRep), `quality_report.tsv` (CheckM2), `gtdbtk.bac120.summary.tsv` and the GTDB tree.

### `02_crispr_and_defense`

* **Input** — the renamed genomes, the level-4 CRISPR arrays.
* **Output** — `minced_res/*.gff`, `Crisprs_REPORT.tsv`, `all_CRISPR_level4.fa`,
  `CCTyper_res*` (subtypes), `res1/*defense_finder_genes.tsv`.

### `03_resistome_and_bgc`

* **Input** — Prodigal ORFs of the genomes (`4538genomes_prodigal.faa`), CARD and SARG databases.
* **Output** — `0ARG_RGI/*.txt` (RGI), `diamond_ARG_out*` and `*_gene_out*_table.csv`
  (diamond + `turn_to_table.py`), Macrel peptide output, antiSMASH `res1/`.

### `04_mobile_genetic_elements`

* **Input** — the renamed genomes and the plasmid/prophage subsets.
* **Output** — IntegronFinder `res/`, ICEfinder2 `ICEres/`, ISEScan `ISEscan_res/`, and the
  annotation of the nested MGEs used in Fig. 2E-H.

### `05_virus_discovery`

* **Input** — genomes (`genomad`), the IMG/VR and VIRE virus sets.
* **Output** — `*_virus_summary.tsv` (geNomad), CheckV quality tables, ANI/mmseqs clusters, and
  the final non-redundant catalogue (`46final_virus`).

### `06_virus_annotation`

* **Input** — the final virus catalogue.
* **Output** — iPHoP host predictions, BACPHLIP lifestyle, per-phage features
  (GC/NARSC/length), DefenseFinder and RGI results on viruses, hhsuite2 annotation of the
  selected viruses, virus subsets of Virus-Host DB and VIRE.

### `07_spacer_competition`

* **Input** — CRISPR spacers (minced), MGE sequences, IMG/VR and VIRE nucleotide databases.
* **Output** — per-MGE spacer tables, BLAST hits against IMG/VR and VIRE, the filtered
  spacer-protospacer pairs, and the extracted/annotated target genes
  (`1_extract_gene.py`, eggNOG-mapper).

### `08_phage_clustering_and_host`

* **Input** — the virus catalogue proteins.
* **Output** — vContact3 and vContact2 networks, `aai_summary*.tsv` (comparem), clinker maps,
  ViPTree trees, multi-host phage tables.

### `09_phylogeny`

* **Input** — marker gene sequences (MCP, portal, TerL) and the special genes.
* **Output** — hhsuite2 annotation, alignments and IQ-TREE trees.

### `10_experiments`

* **Input** — the isolated phage genomes and the host strains.
* **Output** — codon usage (`apc.pl`/`apc.py`), GTDB-Tk trees and defense/prophage content of the
  hosts.

## analyses

### `01_genome_landscape/01_genome_metadata_and_map.R`

* **Input** — `1genome/{ESKAPE_genome_info.tsv, data_summary_*.tsv, quality_report.tsv, Widb.csv}`
* **Output** — `1_genome_info.pdf`, `1_genome_globalmap.pdf` (CS_WORK), updated
  `ESKAPE_genome_info.tsv` (CS_DATA)
* **Figure** — Fig. 1A-B

### `02_phage_catalogue`

| Script | Input | Figure |
| --- | --- | --- |
| `01_virus_genome_info.R` | `2virus_genome_info/{ESKAPE_phage.id, ESKAPE_all_phage_info*.tsv, iphop_phage2hostgenome.tsv}` | Fig. 3A-E |
| `02_phage_genome_features.R` | `2virus_genome_info/{AAU_matrix.csv, pI_distribution_matrix.csv, per_phage_stats.tsv}` | Fig. 3 |

### `03_gene_annotation`

| Script | Input | Figure |
| --- | --- | --- |
| `01_CRISPR_Cas_taxonomy.R` | `3genome_gene_annotation/1CasTyper/ESKAPE_CasType.tsv` | Fig. 1C |
| `02_ARG_resistome.R` | `3genome_gene_annotation/{rgiargs_4_type.csv, diamond_MRG_out3_type.csv, diamond_VFG_out3_type.csv}` | Fig. 1C-D |

### `04_mobilome`

| Script | Input | Figure |
| --- | --- | --- |
| `01_MGE_landscape.R` | `4genome_MGE/{ESKAPE_genome_info2.tsv, plasmid_info.tsv, ICE_contig_info.tsv, ESKAPE_integron_n.tsv, *virus_summary.tsv}` | Fig. 2A, 2C-D |
| `02_spacer_competition.R` | `4_1MGE_spacer_competition/ESKAPE_genome_info3.tsv` + the spacers | Fig. 2B |
| `03_special_MGE_annotation.R` | `4genome_MGE/special_MGE/*` | Fig. 2E-H |

### `05_phage_clustering`

| Script | Input | Figure |
| --- | --- | --- |
| `01_vcontact3_analysis.R` | vContact3 output, `selected_phage_ID.tsv` | Fig. 5A |
| `02_AAI_heatmap.R` | `6_3comparem/aai_summary*.tsv` | Fig. 5B |

### `06_phage_defense`

| Script | Input | Figure |
| --- | --- | --- |
| `01_phage_defense_systems.R` | `6_0phage_DF/*defense_finder_systems.tsv`, `rgiargs_detail_info.txt` | Fig. 3F |
| `02_phage_gene_annotation.R` | `6_1phage_gene_annotation/{all_phage_gff_anno.tsv, phage_annotation_summary.tsv}` | Fig. 3F |

### `07_phage_host`

| Script | Input | Figure |
| --- | --- | --- |
| `01_multi_host_phage.R` | `7_0multi_host_phage/{En_Kp_EnKp_phage_complete.tsv, vcontact2_res_filter.tsv}` | Fig. 4 |
| `02_virus_geographical_distribution.R` | `7_1partial_virus_source/{ESKAPE_IMGVR_phage2.tsv, ESKAPE_vire_phage.tsv, vire_matched_rows.tsv}` | Fig. 4A-B |
| `03_spacer_target_genes.R` | `7_2protospacer_gene/ESKAPE_phage_methylase.tsv` + annotate tables | Fig. 4C |

### `08_phylogeny/01_marker_gene_phylogeny.R`

* **Input** — `7_3marker_gene_phylogeny/{filtered_*_{MCP,portal,TerL}.tsv, all_32phage_gff_anno.tsv}`
* **Figure** — Fig. 6

### `09_experiments/01_phage_codon_usage.pl`

* **Input** — the isolated phage genomes (`all_4_novel_phage.fa`, GFF files).
* **Output** — codon usage tables.
