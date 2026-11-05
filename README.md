# ESKAPE_mobilome_virome

Code accompanying **"Functional profiling of ESKAPE viromes uncovers resistance-limiting
phage-host dynamics"** — Li P., Liu Q., Deng C., Ni J.
bioRxiv (2026). [doi:10.64898/2026.08.16.745060](https://doi.org/10.64898/2026.08.16.745060)

The repository contains the scripts that go from public genome and metagenome data to the
figures of the preprint: the genome catalogue and its resistome/defensome annotation, the
mobile genetic element (MGE) survey, the construction and annotation of the ESKAPE virome, the
CRISPR spacer-based interaction analyses, and the R analyses that produce the statistics and
the figures.

> 中文版说明见 **[README_zh.md](README_zh.md)**。

---

## Study in brief

* **Genome catalogue.** 11,947 complete ESKAPE genomes were downloaded from NCBI RefSeq and
  dereplicated at 99.9 % ANI / 30 % coverage, giving **4,538 non-redundant genomes** covering
  *A. baumannii*, *E. faecium*, *Enterobacter* spp., *K. pneumoniae*, *P. aeruginosa* and
  *S. aureus*.
* **Mobilome.** 35,312 MGEs were identified: 10,199 plasmids, 18,411 prophages, 4,799 ICEs/IMEs,
  1,594 integrons, 309 insertion sequences and 84 composite transposons. Nested MGE
  architectures (transposon inside ICE inside plasmid, integron inside prophage, ...) give
  low-mobility ARGs a route to horizontal transfer.
* **Virome.** Integration of the genome-encoded prophages with IMG/VR and VIRE produced a
  non-redundant catalogue of **14,496 ESKAPE-associated viruses**, including four previously
  unreported viral clades. ESKAPE phages carry almost no ARGs and antagonise plasmids, which is
  the "resistance-limiting" side of the phage-host relationship.
* **Arms race.** Virulent phages are enriched in anti-defense genes, and the novel viruses
  encode methyltransferases that are themselves targeted by host CRISPR-Cas systems.

## What is in here

```
ESKAPE_mobilome_virome/
├── pipelines/          # upstream bioinformatics (SLURM job scripts, python helpers)
│   ├── 01_genome_catalogue/            download, rename, dRep, CheckM2, GTDB-Tk
│   ├── 02_crispr_and_defense/          minced, CRISPRCasFinder, CRISPRCasTyper, DefenseFinder
│   ├── 03_resistome_and_bgc/           RGI/CARD, diamond-SARG, Macrel, antiSMASH
│   ├── 04_mobile_genetic_elements/     IntegronFinder, ICEfinder2, ISEScan, special MGEs
│   ├── 05_virus_discovery/             geNomad, CheckV, clustering, final catalogue
│   ├── 06_virus_annotation/            iPHoP/BACPHLIP, phage features, DefenseFinder, RGI
│   ├── 07_spacer_competition/          MGE spacer competition, spacer -> IMG/VR and VIRE
│   ├── 08_phage_clustering_and_host/   vContact3, vContact2, novel viruses, multi-host phages
│   ├── 09_phylogeny/                   marker genes (MCP, portal, TerL), special genes
│   └── 10_experiments/                 wet-lab phage sequences and host strains
├── analyses/           # downstream statistics and figures (R, one Perl helper)
│   ├── 01_genome_landscape/            Fig. 1A-B
│   ├── 02_phage_catalogue/             Fig. 3A-E
│   ├── 03_gene_annotation/             Fig. 1C-D
│   ├── 04_mobilome/                    Fig. 2
│   ├── 05_phage_clustering/            Fig. 5A-B
│   ├── 06_phage_defense/               Fig. 3F
│   ├── 07_phage_host/                  Fig. 4
│   ├── 08_phylogeny/                   Fig. 6
│   └── 09_experiments/                 codon usage of the isolated phages
├── data/               # small base tables (see data/README.md)
└── docs/               # workflow, layout, file mapping, inputs/outputs, software
```

Each script keeps the original sub-directory name in its path (for example
`pipelines/05_virus_discovery/04_final_catalogue/46final_virus/1virus_filter.sh`), so a step can
always be traced back to the working copy it came from. Every file carries a header with its
purpose, the manuscript figure, and its original location on the cluster.

## Documentation

| Document | Content |
| --- | --- |
| [docs/01_pipeline_overview.md](docs/01_pipeline_overview.md) | End-to-end workflow, stage by stage, with inputs and outputs |
| [docs/02_directory_layout.md](docs/02_directory_layout.md) | What lives where, naming conventions, environment variables |
| [docs/03_file_mapping.md](docs/03_file_mapping.md) | Every script: original path → new path |
| [docs/04_script_inputs_outputs.md](docs/04_script_inputs_outputs.md) | Inputs, outputs and the figure produced by each module |
| [docs/05_software_and_environment.md](docs/05_software_and_environment.md) | Software versions, conda environments, R packages, databases |

## Reproducing the analyses

**1. Upstream pipelines.** The scripts under `pipelines/` are the SLURM job scripts used for the
production runs (`#SBATCH` headers included). They are reproduced as executed, except that the
site-specific absolute paths were replaced by variables. Each script lists them in its header:

```bash
export PROJECT_ROOT=/path/to/08ESKAPE      # the project directory
export HOME_ROOT=/path/to/home             # home of the companion projects
export DB_ROOT=/path/to/databases          # reference databases
```

Submit the stages in the order given in [docs/01](docs/01_pipeline_overview.md); the conda
environments are listed in [docs/05](docs/05_software_and_environment.md).

**2. R analyses.** The scripts under `analyses/` require R >= 4.2 and the packages in
[docs/05](docs/05_software_and_environment.md). Inputs are read from `CS_DATA` and figures are
written to `CS_WORK`, both defaulting to the repository layout:

```r
Sys.setenv(CS_DATA = "/path/to/data", CS_WORK = "/path/to/results")
source("analyses/01_genome_landscape/01_genome_metadata_and_map.R")
```

`CS_EXT` points at the tables imported from the companion predatory-bacteria project
(`data/external`), which a few sections re-use.

## Data availability

The genome catalogue, the MGE tables and the ESKAPE virome are provided as supplementary tables
with the preprint. The small base tables used by the first analyses are included under
[`data/`](data/README.md); the larger intermediates (geNomad/CheckV outputs, vContact3 results,
spacer BLAST tables) are regenerated by `pipelines/`. All primary data are public: the genomes
come from NCBI RefSeq, the virus references from IMG/VR, VIRE and Virus-Host DB.

## Citation

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

If you use the scripts, please cite the preprint above.
