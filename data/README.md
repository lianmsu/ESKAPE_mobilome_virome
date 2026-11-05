# `data/` — base tables for the R analyses

```
data/
├── 1genome/
│   ├── ESKAPE_genome_info.tsv        # the 4 538 non-redundant genomes and their metadata
│   ├── data_summary_{Ab,Ef,En,Kp,Pa,Sa}.tsv   # per-species NCBI assembly summaries
│   ├── quality_report.tsv            # CheckM2 completeness/contamination
│   └── Widb.csv                      # isolation metadata used for the world map
├── 2virus_genome_info/
│   ├── ESKAPE_phage.id               # identifiers of the retained viruses
│   ├── per_phage_stats.tsv           # length, GC, coding density per virus
│   └── Molecular_Final_Full_Stats.txt
└── external/                         # (not distributed) tables from the companion project
```

The R scripts resolve this directory through `CS_DATA` and, for the few tables imported from the
companion predatory-bacteria project, through `CS_EXT` (default `data/external`).

Everything else that the analyses read is produced by the pipelines and is therefore not
distributed here:

* `3genome_gene_annotation/` — RGI/CARD and diamond tables (`rgiargs_4_type.csv`,
  `diamond_MRG_out3_type.csv`, `diamond_VFG_out3_type.csv`) from `pipelines/03`.
* `4genome_MGE/`, `4_1MGE_spacer_competition/` — MGE tables from `pipelines/04` and `07`.
* `6_0phage_DF/`, `6_1phage_gene_annotation/`, `6_3comparem/` — virus annotation from
  `pipelines/06` and `08`.
* `7_*` — host range, distribution and spacer-target tables from `pipelines/07` and `08`.
