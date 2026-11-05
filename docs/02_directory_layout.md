# Directory layout

```
ESKAPE_mobilome_virome/
├── pipelines/       upstream work, numbered in execution order
│   ├── 01_genome_catalogue/
│   ├── 02_crispr_and_defense/
│   ├── 03_resistome_and_bgc/
│   ├── 04_mobile_genetic_elements/
│   ├── 05_virus_discovery/
│   ├── 06_virus_annotation/
│   ├── 07_spacer_competition/
│   ├── 08_phage_clustering_and_host/
│   ├── 09_phylogeny/
│   └── 10_experiments/
├── analyses/        statistics and figures (R, one Perl helper)
│   ├── 01_genome_landscape/
│   ├── 02_phage_catalogue/
│   ├── 03_gene_annotation/
│   ├── 04_mobilome/
│   ├── 05_phage_clustering/
│   ├── 06_phage_defense/
│   ├── 07_phage_host/
│   ├── 08_phylogeny/
│   └── 09_experiments/
├── data/            small base tables, see data/README.md
└── docs/
```

## Conventions

* **Two layers.** `pipelines/` holds everything that turns raw data into a table (SLURM job
  scripts, python and perl helpers); `analyses/` holds everything that turns a table into a
  statistic or a figure (R).
* **Provenance in the path.** Inside a module the original working-copy directory is kept as a
  sub-directory, e.g. `pipelines/05_virus_discovery/04_final_catalogue/46final_virus/`. This
  makes it possible to compare a script with the cluster copy at any time. `00_notes/` collects
  the small `readme.sh` command logs of a stage.
* **Numbering.** `pipelines/` is numbered by stage of the workflow (01 → 10). `analyses/` is
  numbered by topic, and each module lists the figure it produces.
* **Headers.** Every file starts with a header giving the module, a one-line purpose, the
  manuscript figure, and the original path on the cluster or in the working copy (`lfat1t01:…`,
  `pkusz_CPU:…`, `2ESKAPE_MGE/…`).
* **No personal information.** All paths containing a user name, account, mail address or
  personal directory were removed. The remaining site-specific parts are variables.

## Environment variables

| Variable | Used by | Default | Meaning |
| --- | --- | --- | --- |
| `PROJECT_ROOT` | shell scripts | `/path/to/08ESKAPE` | the 08ESKAPE project directory |
| `HOME_ROOT` | shell scripts | `/path/to/home` | home of the companion projects (`06mobilome_predators`, `04Qinghai_Tibet_Plateau`, …) |
| `SHARED_HOME` | shell scripts | `/path/to/shared_home` | shared home on the second cluster |
| `DB_ROOT` | shell scripts | `/path/to/databases` | reference databases (geNomad, CheckM2, GTDB, CARD, …) |
| `CS_HOME` | R scripts | `.` | repository root |
| `CS_DATA` | R scripts | `$CS_HOME/data` | root of the input tables |
| `CS_WORK` | R scripts | `$CS_HOME/results` | where figures and tables are written |
| `CS_EXT` | R scripts | `$CS_HOME/data/external` | tables imported from the companion project |

## Where the outputs go

* `pipelines/` writes into `PROJECT_ROOT/<original directory>`, exactly as the production runs
  did; nothing is written inside the repository.
* R figures go to `CS_WORK` (git-ignored, created on the fly); data tables that a later section
  of the same script re-reads stay in `CS_DATA`, as in the original notebooks.
* Cluster logs (`job.*.out`, `nohup_*.log`) are not part of the repository.
