# File mapping (original working copies -> this repository)

Three sources were merged: the local R working copy (`2ESKAPE_MGE`), the cluster copy on
`lfat1t01` (`/share/home/<user>/08ESKAPE`) and the cluster copy on `pkusz_CPU`
(`/home/<user>/08ESKAPE`). Directories whose name ended in `outdated` or `wrong`
are not part of the repository.

| New path | Original path | Purpose |
| --- | --- | --- |
| `analyses/01_genome_landscape/01_genome_metadata_and_map.R` | `2ESKAPE_MGE/1genome/1_genome_info.R` | Genome metadata, quality and the world map of the 4 538 genomes |
| `analyses/02_phage_catalogue/01_virus_genome_info.R` | `2ESKAPE_MGE/2virus_genome_info/1_virus_genome_info_o.R` | Quality, taxonomy, host and lifestyle of the retained viruses |
| `analyses/02_phage_catalogue/02_phage_genome_features.R` | `2ESKAPE_MGE/2virus_genome_info/2_phage_genome_features.R` | Genome features of complete phages (AAU, pI, length) |
| `analyses/03_gene_annotation/01_CRISPR_Cas_taxonomy.R` | `2ESKAPE_MGE/3genome_gene_annotation/1CasTyper/1_CRISPRCas_tax.R` | CRISPR-Cas subtype distribution across the ESKAPE lineages |
| `analyses/03_gene_annotation/02_ARG_resistome.R` | `2ESKAPE_MGE/3genome_gene_annotation/1_ARG.R` | ARG prevalence, drug classes and ARG/MRG/VFG/DS relationships |
| `analyses/04_mobilome/01_MGE_landscape.R` | `2ESKAPE_MGE/4genome_MGE/1_ESKAPE_MGE.R` | MGE counts per lineage, ARG and anti-defense cargo, nested MGE cases |
| `analyses/04_mobilome/02_spacer_competition.R` | `2ESKAPE_MGE/4_1MGE_spacer_competition/1_MGE_competition.R` | CRISPR spacer targeting of the different MGE types |
| `analyses/04_mobilome/03_special_MGE_annotation.R` | `2ESKAPE_MGE/4genome_MGE/special_MGE/1_MGE_genome_annotation.R` | Annotation of the representative nested MGEs |
| `analyses/05_phage_clustering/01_vcontact3_analysis.R` | `2ESKAPE_MGE/5phage_vcontact3/1_vcontact3_analysis.R` | vContact3 genome-sharing network and taxonomic assignment |
| `analyses/05_phage_clustering/02_AAI_heatmap.R` | `2ESKAPE_MGE/6_3comparem/1_aai_heatmap.R` | Pairwise AAI heat map of the new viral groups |
| `analyses/06_phage_defense/01_phage_defense_systems.R` | `2ESKAPE_MGE/6_0phage_DF/1_phage_DF.R` | Defense and anti-defense genes of the phages |
| `analyses/06_phage_defense/02_phage_gene_annotation.R` | `2ESKAPE_MGE/6_1phage_gene_annotation/1_phage_gene_annotation.R` | Functional annotation of the phage genes (hhsuite2 -> KEGG/COG) |
| `analyses/07_phage_host/01_multi_host_phage.R` | `2ESKAPE_MGE/7_0multi_host_phage/1_multi_host_phage.R` | Multi-host phages and their CRISPR-based host links |
| `analyses/07_phage_host/02_virus_geographical_distribution.R` | `2ESKAPE_MGE/7_1partial_virus_source/1_phage_geographical_distribution.R` | Ecosystem and geographical distribution of the IMG/VR and VIRE viruses |
| `analyses/07_phage_host/03_spacer_target_genes.R` | `2ESKAPE_MGE/7_2protospacer_gene/1_gene_annotation.R` | Functional types of the viral genes targeted by host spacers |
| `analyses/08_phylogeny/01_marker_gene_phylogeny.R` | `2ESKAPE_MGE/7_3marker_gene_phylogeny/1_phage_gene_annotation.R` | Marker-gene phylogenies (MCP, portal, TerL) of the phages |
| `analyses/09_experiments/01_phage_codon_usage.pl` | `2ESKAPE_MGE/8_experiments/apc.pl` | Codon usage analysis of the isolated phages (Perl helper) |
| `pipelines/01_genome_catalogue/00_notes/1genome/0readme.sh` | `lfat1t01:1genome/0readme.sh` | Genome renaming, GTDB metadata matching and plasmid extraction commands |
| `pipelines/01_genome_catalogue/01_dereplication/1genome/0readme.sh` | `pkusz_CPU:1genome/0readme.sh` | Download and rename the complete ESKAPE genomes from NCBI RefSeq |
| `pipelines/01_genome_catalogue/01_dereplication/1genome/1rename_genome.sh` | `pkusz_CPU:1genome/1rename_genome.sh` | Download and rename the complete ESKAPE genomes from NCBI RefSeq |
| `pipelines/01_genome_catalogue/01_dereplication/2drep/0readme.sh` | `pkusz_CPU:2drep/0readme.sh` | dRep dereplication of the 11 947 ESKAPE genomes |
| `pipelines/01_genome_catalogue/01_dereplication/2drep/drep_MAGs.sh` | `pkusz_CPU:2drep/drep_MAGs.sh` | dRep dereplication of the 11 947 ESKAPE genomes |
| `pipelines/01_genome_catalogue/02_quality_and_taxonomy/3checkm2/checkm2.sh` | `lfat1t01:1genome/3checkm2/checkm2.sh` | Genome quality with CheckM2 |
| `pipelines/01_genome_catalogue/02_quality_and_taxonomy/3gtdbtk/01GTDB_classify.sh` | `pkusz_CPU:3gtdbtk/01GTDB_classify.sh` | GTDB-Tk classification of the dereplicated genomes |
| `pipelines/01_genome_catalogue/02_quality_and_taxonomy/4GTDB/GTDB_infer.sh` | `lfat1t01:1genome/4GTDB/GTDB_infer.sh` | GTDB-Tk taxonomy and reference tree of the dereplicated genomes |
| `pipelines/01_genome_catalogue/02_quality_and_taxonomy/align/GTDB_infer.sh` | `pkusz_CPU:3gtdbtk/align/GTDB_infer.sh` | GTDB-Tk classification of the dereplicated genomes |
| `pipelines/02_crispr_and_defense/00_notes/23CRISPRCasFinder/readme.sh` | `lfat1t01:23CRISPRCasFinder/readme.sh` | CRISPR-Cas detection with CRISPRCasFinder (level 4 arrays) |
| `pipelines/02_crispr_and_defense/00_notes/DG_gff/readme.sh` | `lfat1t01:30Defensefinder/res1/DG_gff/readme.sh` | Defense system annotation of the ESKAPE genomes with DefenseFinder |
| `pipelines/02_crispr_and_defense/00_notes/res1/readme.sh` | `lfat1t01:30Defensefinder/res1/readme.sh` | Defense system annotation of the ESKAPE genomes with DefenseFinder |
| `pipelines/02_crispr_and_defense/01_minced/22minced/1minced.sh` | `lfat1t01:22minced/1minced.sh` | CRISPR spacer extraction with minced on the ESKAPE genomes |
| `pipelines/02_crispr_and_defense/01_minced/25level4_minced/1minced.sh` | `lfat1t01:25level4_minced/1minced.sh` | minced spacer extraction on the level-4 CRISPR arrays |
| `pipelines/02_crispr_and_defense/01_minced/minced_res/0readme.sh` | `lfat1t01:22minced/minced_res/0readme.sh` | CRISPR spacer extraction with minced on the ESKAPE genomes |
| `pipelines/02_crispr_and_defense/02_crisprcasfinder/23CRISPRCasFinder/1CRISPRCasFinder.sh` | `lfat1t01:23CRISPRCasFinder/1CRISPRCasFinder.sh` | CRISPR-Cas detection with CRISPRCasFinder (level 4 arrays) |
| `pipelines/02_crispr_and_defense/03_crisprcaster/24CRISPRCasTyper/0readme.sh` | `lfat1t01:24CRISPRCasTyper/0readme.sh` | CRISPR-Cas subtyping with CRISPRCasTyper |
| `pipelines/02_crispr_and_defense/03_crisprcaster/24CRISPRCasTyper/1CRISPRTyper.sh` | `lfat1t01:24CRISPRCasTyper/1CRISPRTyper.sh` | CRISPR-Cas subtyping with CRISPRCasTyper |
| `pipelines/02_crispr_and_defense/03_crisprcaster/24CRISPRCasTyper/2readme.sh` | `lfat1t01:24CRISPRCasTyper/2readme.sh` | CRISPR-Cas subtyping with CRISPRCasTyper |
| `pipelines/02_crispr_and_defense/04_defensefinder/30Defensefinder/01DefenseFinder.sh` | `lfat1t01:30Defensefinder/01DefenseFinder.sh` | Defense system annotation of the ESKAPE genomes with DefenseFinder |
| `pipelines/03_resistome_and_bgc/01_resistome_annotation/0ARG_RGI/0readme.sh` | `lfat1t01:31ARG/0ARG_RGI/0readme.sh` | ARG annotation: RGI/CARD and diamond against SARG |
| `pipelines/03_resistome_and_bgc/01_resistome_annotation/31ARG/0readme.sh` | `lfat1t01:31ARG/0readme.sh` | ARG annotation: RGI/CARD and diamond against SARG |
| `pipelines/03_resistome_and_bgc/01_resistome_annotation/31ARG/1RGI_CARD.sh` | `lfat1t01:31ARG/1RGI_CARD.sh` | ARG annotation: RGI/CARD and diamond against SARG |
| `pipelines/03_resistome_and_bgc/01_resistome_annotation/31ARG/2mag_orf_diamond_li.sh` | `lfat1t01:31ARG/2mag_orf_diamond_li.sh` | ARG annotation: RGI/CARD and diamond against SARG |
| `pipelines/03_resistome_and_bgc/01_resistome_annotation/31ARG/3diamond_res_pretreatment.sh` | `lfat1t01:31ARG/3diamond_res_pretreatment.sh` | ARG annotation: RGI/CARD and diamond against SARG |
| `pipelines/03_resistome_and_bgc/01_resistome_annotation/31ARG/4turn_to_table.py` | `lfat1t01:31ARG/4turn_to_table.py` | ARG annotation: RGI/CARD and diamond against SARG |
| `pipelines/03_resistome_and_bgc/02_amp_macrel/32AMP/01prodigal.sh` | `lfat1t01:32AMP/01prodigal.sh` | Antimicrobial peptide prediction with Macrel |
| `pipelines/03_resistome_and_bgc/02_amp_macrel/32AMP/02macrel.sh` | `lfat1t01:32AMP/02macrel.sh` | Antimicrobial peptide prediction with Macrel |
| `pipelines/03_resistome_and_bgc/02_amp_macrel/32AMP/0readme.sh` | `lfat1t01:32AMP/0readme.sh` | Antimicrobial peptide prediction with Macrel |
| `pipelines/03_resistome_and_bgc/03_antismash/33antismash/01prodigal.sh` | `lfat1t01:33antismash/01prodigal.sh` | Biosynthetic gene cluster prediction with antiSMASH |
| `pipelines/03_resistome_and_bgc/03_antismash/33antismash/02antiSMASH.sh` | `lfat1t01:33antismash/02antiSMASH.sh` | Biosynthetic gene cluster prediction with antiSMASH |
| `pipelines/03_resistome_and_bgc/03_antismash/33antismash/0readme.sh` | `lfat1t01:33antismash/0readme.sh` | Biosynthetic gene cluster prediction with antiSMASH |
| `pipelines/04_mobile_genetic_elements/00_notes/integron_filter/readme.sh` | `pkusz_CPU:4integronfinder/integron_filter/readme.sh` | Integron detection with IntegronFinder |
| `pipelines/04_mobile_genetic_elements/01_ice_integron_isescan/4integronfinder/01integron_finder.sh` | `pkusz_CPU:4integronfinder/01integron_finder.sh` | Integron detection with IntegronFinder |
| `pipelines/04_mobile_genetic_elements/01_ice_integron_isescan/5ICEfinder/00readme.sh` | `pkusz_CPU:5ICEfinder/00readme.sh` | ICE detection with ICEfinder2 (batches of 1000 and the full set) |
| `pipelines/04_mobile_genetic_elements/01_ice_integron_isescan/5ICEfinder/icefinder_1000.sh` | `pkusz_CPU:5ICEfinder/icefinder_1000.sh` | ICE detection with ICEfinder2 (batches of 1000 and the full set) |
| `pipelines/04_mobile_genetic_elements/01_ice_integron_isescan/5ICEfinder/icefinder_2000.sh` | `pkusz_CPU:5ICEfinder/icefinder_2000.sh` | ICE detection with ICEfinder2 (batches of 1000 and the full set) |
| `pipelines/04_mobile_genetic_elements/01_ice_integron_isescan/5ICEfinder/icefinder_3000.sh` | `pkusz_CPU:5ICEfinder/icefinder_3000.sh` | ICE detection with ICEfinder2 (batches of 1000 and the full set) |
| `pipelines/04_mobile_genetic_elements/01_ice_integron_isescan/5ICEfinder/icefinder_4538.sh` | `pkusz_CPU:5ICEfinder/icefinder_4538.sh` | ICE detection with ICEfinder2 (batches of 1000 and the full set) |
| `pipelines/04_mobile_genetic_elements/01_ice_integron_isescan/6ISEscan/00readme.sh` | `pkusz_CPU:6ISEscan/00readme.sh` | Insertion sequence detection with ISEScan |
| `pipelines/04_mobile_genetic_elements/01_ice_integron_isescan/6ISEscan/ISEScan_4538MAGs.sh` | `pkusz_CPU:6ISEscan/ISEScan_4538MAGs.sh` | Insertion sequence detection with ISEScan |
| `pipelines/04_mobile_genetic_elements/01_ice_integron_isescan/ICEres/0readme.sh` | `pkusz_CPU:5ICEfinder/ICEres/0readme.sh` | ICE detection with ICEfinder2 (batches of 1000 and the full set) |
| `pipelines/04_mobile_genetic_elements/01_ice_integron_isescan/ICEres/process_iceid.py` | `pkusz_CPU:5ICEfinder/ICEres/process_iceid.py` | ICE detection with ICEfinder2 (batches of 1000 and the full set) |
| `pipelines/04_mobile_genetic_elements/01_ice_integron_isescan/ISEscan_res/0readme.sh` | `pkusz_CPU:6ISEscan/ISEscan_res/0readme.sh` | Insertion sequence detection with ISEScan |
| `pipelines/04_mobile_genetic_elements/02_special_mge/34special_MGE/0readme.sh` | `lfat1t01:34special_MGE/0readme.sh` | Selection of the special MGEs used in the case studies |
| `pipelines/04_mobile_genetic_elements/02_special_mge/5_ICEfinder_special_MGE/0readme.sh` | `pkusz_CPU:5_ICEfinder_special_MGE/0readme.sh` | ICEfinder2 on the selected plasmid/prophage set |
| `pipelines/05_virus_discovery/01_genomad/21genomad/01genomad_lpw.sh` | `lfat1t01:21genomad/01genomad_lpw.sh` | geNomad virus and plasmid detection on the ESKAPE genomes |
| `pipelines/05_virus_discovery/01_genomad/45virus_genomad/01genomad_lpw.sh` | `lfat1t01:45virus_genomad/01genomad_lpw.sh` | geNomad on the second virus batch |
| `pipelines/05_virus_discovery/02_checkv/42genomad_checkv/01checkv_lpw.sh` | `lfat1t01:42genomad_checkv/01checkv_lpw.sh` | CheckV quality control of the geNomad virus candidates |
| `pipelines/05_virus_discovery/02_checkv/42genomad_checkv/02readme.sh` | `lfat1t01:42genomad_checkv/02readme.sh` | CheckV quality control of the geNomad virus candidates |
| `pipelines/05_virus_discovery/02_checkv/45virus_checkv/01checkv_lpw.sh` | `lfat1t01:45virus_checkv/01checkv_lpw.sh` | CheckV on the second virus batch |
| `pipelines/05_virus_discovery/03_clustering/43combined_all_virus/0readme.sh` | `lfat1t01:43combined_all_virus/0readme.sh` | Combine the virus candidates from all sources |
| `pipelines/05_virus_discovery/03_clustering/44virus_checkv_cluster/0checkv_cluster.sh` | `lfat1t01:44virus_checkv_cluster/0checkv_cluster.sh` | CheckV + ANI clustering of the combined virus set |
| `pipelines/05_virus_discovery/03_clustering/44virus_checkv_cluster/anicalc.py` | `lfat1t01:44virus_checkv_cluster/anicalc.py` | CheckV + ANI clustering of the combined virus set |
| `pipelines/05_virus_discovery/03_clustering/44virus_checkv_cluster/aniclust.py` | `lfat1t01:44virus_checkv_cluster/aniclust.py` | CheckV + ANI clustering of the combined virus set |
| `pipelines/05_virus_discovery/03_clustering/44virus_mmseqs/2mmseqs.sh` | `lfat1t01:44virus_mmseqs/2mmseqs.sh` | mmseqs2 clustering of the virus set |
| `pipelines/05_virus_discovery/04_final_catalogue/46final_virus/0readme.sh` | `lfat1t01:46final_virus/0readme.sh` | Filtering rules that define the final non-redundant ESKAPE virome |
| `pipelines/05_virus_discovery/04_final_catalogue/46final_virus/1virus_filter.sh` | `lfat1t01:46final_virus/1virus_filter.sh` | Filtering rules that define the final non-redundant ESKAPE virome |
| `pipelines/06_virus_annotation/01_host_prediction/46iPHoP/1_iPHoP_lpw.sh` | `lfat1t01:46iPHoP/1_iPHoP_lpw.sh` | Host prediction with iPHoP and lifestyle with BACPHLIP |
| `pipelines/06_virus_annotation/01_host_prediction/46iPHoP/2_iPHoP_res_filter.sh` | `lfat1t01:46iPHoP/2_iPHoP_res_filter.sh` | Host prediction with iPHoP and lifestyle with BACPHLIP |
| `pipelines/06_virus_annotation/01_host_prediction/46iPHoP/3_bacphlip_lpw.sh` | `lfat1t01:46iPHoP/3_bacphlip_lpw.sh` | Host prediction with iPHoP and lifestyle with BACPHLIP |
| `pipelines/06_virus_annotation/01_host_prediction/iphop_output_refdb/res_filter.sh` | `lfat1t01:46iPHoP/iphop_output_refdb/res_filter.sh` | Host prediction with iPHoP and lifestyle with BACPHLIP |
| `pipelines/06_virus_annotation/02_phage_features/46phage_features/0readme.sh` | `lfat1t01:46phage_features/0readme.sh` | Molecular features of the phages: GC, NARSC, length, quality |
| `pipelines/06_virus_annotation/02_phage_features/46phage_features/1_molecular_characteristics.sh` | `lfat1t01:46phage_features/1_molecular_characteristics.sh` | Molecular features of the phages: GC, NARSC, length, quality |
| `pipelines/06_virus_annotation/02_phage_features/46phage_features/compute_features_v1.py` | `lfat1t01:46phage_features/compute_features_v1.py` | Molecular features of the phages: GC, NARSC, length, quality |
| `pipelines/06_virus_annotation/02_phage_features/46phage_features/get_gc_and_narsc.py` | `lfat1t01:46phage_features/get_gc_and_narsc.py` | Molecular features of the phages: GC, NARSC, length, quality |
| `pipelines/06_virus_annotation/02_phage_features/46phage_features/get_gc_and_narsc_python3.py` | `lfat1t01:46phage_features/get_gc_and_narsc_python3.py` | Molecular features of the phages: GC, NARSC, length, quality |
| `pipelines/06_virus_annotation/03_location_and_source/46phage_location_source/0readme.sh` | `lfat1t01:46phage_location_source/0readme.sh` | Location and source metadata of the retained viruses |
| `pipelines/06_virus_annotation/04_virus_defense/47virus_defensefinder/01DefenseFinder.sh` | `lfat1t01:47virus_defensefinder/01DefenseFinder.sh` | Defense system annotation of the viruses |
| `pipelines/06_virus_annotation/05_virus_resistome/0ARG_RGI/0readme.sh` | `lfat1t01:48virus_RGI/0ARG_RGI/0readme.sh` | ARG annotation of the viruses with RGI/CARD |
| `pipelines/06_virus_annotation/05_virus_resistome/0ARG_RGI2/0readme.sh` | `lfat1t01:48virus_RGI/0ARG_RGI2/0readme.sh` | ARG annotation of the viruses with RGI/CARD |
| `pipelines/06_virus_annotation/05_virus_resistome/48virus_RGI/1RGI_CARD.sh` | `lfat1t01:48virus_RGI/1RGI_CARD.sh` | ARG annotation of the viruses with RGI/CARD |
| `pipelines/06_virus_annotation/06_virus_and_host_databases/1virushostdb/0readme.sh` | `lfat1t01:2virus_database/1virushostdb/0readme.sh` | ESKAPE virus subset of Virus-Host DB |
| `pipelines/06_virus_annotation/06_virus_and_host_databases/1virushostdb/1prodigal.sh` | `lfat1t01:2virus_database/1virushostdb/1prodigal.sh` | ESKAPE virus subset of Virus-Host DB |
| `pipelines/06_virus_annotation/06_virus_and_host_databases/2VIRE/0download.sh` | `lfat1t01:2virus_database/2VIRE/0download.sh` | VIRE download and ORF prediction for the ESKAPE virus references |
| `pipelines/06_virus_annotation/06_virus_and_host_databases/2VIRE/0readme.sh` | `lfat1t01:2virus_database/2VIRE/0readme.sh` | VIRE download and ORF prediction for the ESKAPE virus references |
| `pipelines/06_virus_annotation/06_virus_and_host_databases/2VIRE/1prodigal.sh` | `lfat1t01:2virus_database/2VIRE/1prodigal.sh` | VIRE download and ORF prediction for the ESKAPE virus references |
| `pipelines/06_virus_annotation/07_selected_virus_annotation/49selected_virus_annotation/0readme.sh` | `lfat1t01:49selected_virus_annotation/0readme.sh` | hhsuite2 annotation of the selected (novel) viruses |
| `pipelines/06_virus_annotation/07_selected_virus_annotation/49selected_virus_annotation/1hhsuite1.sh` | `lfat1t01:49selected_virus_annotation/1hhsuite1.sh` | hhsuite2 annotation of the selected (novel) viruses |
| `pipelines/06_virus_annotation/07_selected_virus_annotation/49selected_virus_annotation/2_anno_treatment.py` | `lfat1t01:49selected_virus_annotation/2_anno_treatment.py` | hhsuite2 annotation of the selected (novel) viruses |
| `pipelines/06_virus_annotation/07_selected_virus_annotation/49selected_virus_annotation/3_anno_summary.py` | `lfat1t01:49selected_virus_annotation/3_anno_summary.py` | hhsuite2 annotation of the selected (novel) viruses |
| `pipelines/07_spacer_competition/01_mge_crispr_spacers/1MGE_spacer/0readme.sh` | `lfat1t01:35MGE_spacer_competition/1MGE_spacer/0readme.sh` | Spacer-based competition between MGE types (CRISPRCasTyper on MGE sets) |
| `pipelines/07_spacer_competition/01_mge_crispr_spacers/2blastn_db/0readme.sh` | `lfat1t01:35MGE_spacer_competition/2blastn_db/0readme.sh` | Spacer-based competition between MGE types (CRISPRCasTyper on MGE sets) |
| `pipelines/07_spacer_competition/01_mge_crispr_spacers/35MGE_spacer_competition/0readme.sh` | `lfat1t01:35MGE_spacer_competition/0readme.sh` | Spacer-based competition between MGE types (CRISPRCasTyper on MGE sets) |
| `pipelines/07_spacer_competition/01_mge_crispr_spacers/35MGE_spacer_competition/1CRISPRTyper.sh` | `lfat1t01:35MGE_spacer_competition/1CRISPRTyper.sh` | Spacer-based competition between MGE types (CRISPRCasTyper on MGE sets) |
| `pipelines/07_spacer_competition/01_mge_crispr_spacers/35MGE_spacer_competition/1CRISPRTyper_ICE.sh` | `lfat1t01:35MGE_spacer_competition/1CRISPRTyper_ICE.sh` | Spacer-based competition between MGE types (CRISPRCasTyper on MGE sets) |
| `pipelines/07_spacer_competition/01_mge_crispr_spacers/35MGE_spacer_competition/1CRISPRTyper_integron.sh` | `lfat1t01:35MGE_spacer_competition/1CRISPRTyper_integron.sh` | Spacer-based competition between MGE types (CRISPRCasTyper on MGE sets) |
| `pipelines/07_spacer_competition/01_mge_crispr_spacers/35MGE_spacer_competition/1CRISPRTyper_prophage.sh` | `lfat1t01:35MGE_spacer_competition/1CRISPRTyper_prophage.sh` | Spacer-based competition between MGE types (CRISPRCasTyper on MGE sets) |
| `pipelines/07_spacer_competition/01_mge_crispr_spacers/3spacer2db/0readme.sh` | `lfat1t01:35MGE_spacer_competition/3spacer2db/0readme.sh` | Spacer-based competition between MGE types (CRISPRCasTyper on MGE sets) |
| `pipelines/07_spacer_competition/02_spacer_to_imgvr/40spacer2IMGVR/0readme.sh` | `lfat1t01:40spacer2IMGVR/0readme.sh` | BLAST of the host spacers against the IMG/VR database |
| `pipelines/07_spacer_competition/02_spacer_to_imgvr/40spacer2IMGVR/2makeblastdb.sh` | `lfat1t01:40spacer2IMGVR/2makeblastdb.sh` | BLAST of the host spacers against the IMG/VR database |
| `pipelines/07_spacer_competition/02_spacer_to_imgvr/40spacer2IMGVR/3blastn_IMGVR.sh` | `lfat1t01:40spacer2IMGVR/3blastn_IMGVR.sh` | BLAST of the host spacers against the IMG/VR database |
| `pipelines/07_spacer_competition/02_spacer_to_imgvr/40spacer2IMGVR/4res_filter.sh` | `lfat1t01:40spacer2IMGVR/4res_filter.sh` | BLAST of the host spacers against the IMG/VR database |
| `pipelines/07_spacer_competition/03_spacer_to_vire/41spacer2VIRE/0readme.sh` | `lfat1t01:41spacer2VIRE/0readme.sh` | BLAST of the host spacers against the VIRE database |
| `pipelines/07_spacer_competition/03_spacer_to_vire/41spacer2VIRE/2makeblastdb.sh` | `lfat1t01:41spacer2VIRE/2makeblastdb.sh` | BLAST of the host spacers against the VIRE database |
| `pipelines/07_spacer_competition/03_spacer_to_vire/41spacer2VIRE/3blastn_VIRE.sh` | `lfat1t01:41spacer2VIRE/3blastn_VIRE.sh` | BLAST of the host spacers against the VIRE database |
| `pipelines/07_spacer_competition/03_spacer_to_vire/41spacer2VIRE/4res_filter.sh` | `lfat1t01:41spacer2VIRE/4res_filter.sh` | BLAST of the host spacers against the VIRE database |
| `pipelines/07_spacer_competition/04_spacer_target_genes/54spacer_target_gene/0_readme.sh` | `lfat1t01:54spacer_target_gene/0_readme.sh` | Extraction and annotation of the genes targeted by CRISPR spacers |
| `pipelines/07_spacer_competition/04_spacer_target_genes/54spacer_target_gene/1_extract_gene.py` | `lfat1t01:54spacer_target_gene/1_extract_gene.py` | Extraction and annotation of the genes targeted by CRISPR spacers |
| `pipelines/07_spacer_competition/04_spacer_target_genes/54spacer_target_gene/1_extract_gene_vire.py` | `lfat1t01:54spacer_target_gene/1_extract_gene_vire.py` | Extraction and annotation of the genes targeted by CRISPR spacers |
| `pipelines/07_spacer_competition/04_spacer_target_genes/54spacer_target_gene/2_eggnog_mapper.sh` | `lfat1t01:54spacer_target_gene/2_eggnog_mapper.sh` | Extraction and annotation of the genes targeted by CRISPR spacers |
| `pipelines/07_spacer_competition/04_spacer_target_genes/methylase_phylogeny/0_readme.sh` | `lfat1t01:54spacer_target_gene/methylase_phylogeny/0_readme.sh` | Extraction and annotation of the genes targeted by CRISPR spacers |
| `pipelines/07_spacer_competition/04_spacer_target_genes/methylase_phylogeny/1_iqtree.sh` | `lfat1t01:54spacer_target_gene/methylase_phylogeny/1_iqtree.sh` | Extraction and annotation of the genes targeted by CRISPR spacers |
| `pipelines/08_phage_clustering_and_host/01_vcontact3/50vcontact3/1vcontact3.sh` | `lfat1t01:50vcontact3/1vcontact3.sh` | vContact3 clustering of the ESKAPE viruses |
| `pipelines/08_phage_clustering_and_host/01_vcontact3/50vcontact3/1vcontact3_2.sh` | `lfat1t01:50vcontact3/1vcontact3_2.sh` | vContact3 clustering of the ESKAPE viruses |
| `pipelines/08_phage_clustering_and_host/01_vcontact3/50vcontact3/1vcontact3_3.sh` | `lfat1t01:50vcontact3/1vcontact3_3.sh` | vContact3 clustering of the ESKAPE viruses |
| `pipelines/08_phage_clustering_and_host/02_vcontact2/51vcontact2/0readme.sh` | `lfat1t01:51vcontact2/0readme.sh` | vContact2 clustering for the companion comparison |
| `pipelines/08_phage_clustering_and_host/03_novel_virus/52potential_novel_virus/0readme.sh` | `lfat1t01:52potential_novel_virus/0readme.sh` | Identification and comparison of the potential novel viruses |
| `pipelines/08_phage_clustering_and_host/03_novel_virus/52potential_novel_virus/1comparem.sh` | `lfat1t01:52potential_novel_virus/1comparem.sh` | Identification and comparison of the potential novel viruses |
| `pipelines/08_phage_clustering_and_host/03_novel_virus/clinker_virus/0readme.sh` | `lfat1t01:52potential_novel_virus/clinker_virus/0readme.sh` | Identification and comparison of the potential novel viruses |
| `pipelines/08_phage_clustering_and_host/03_novel_virus/clinker_virus/1vcontact3_3.sh` | `lfat1t01:52potential_novel_virus/clinker_virus/1vcontact3_3.sh` | Identification and comparison of the potential novel viruses |
| `pipelines/08_phage_clustering_and_host/04_multi_host_phage/53multi_host_phage/0_readme.sh` | `lfat1t01:53multi_host_phage/0_readme.sh` | Multi-host phages and their host ranges |
| `pipelines/08_phage_clustering_and_host/04_multi_host_phage/53multi_host_phage/1vcontact3_3.sh` | `lfat1t01:53multi_host_phage/1vcontact3_3.sh` | Multi-host phages and their host ranges |
| `pipelines/08_phage_clustering_and_host/04_multi_host_phage/clinker/0readme.sh` | `lfat1t01:53multi_host_phage/clinker/0readme.sh` | Multi-host phages and their host ranges |
| `pipelines/09_phylogeny/01_marker_genes/60marker_gene_phylogeny/0_readme.sh` | `lfat1t01:60marker_gene_phylogeny/0_readme.sh` | hhsuite2 annotation and IQ-TREE phylogenies of MCP / portal / TerL |
| `pipelines/09_phylogeny/01_marker_genes/60marker_gene_phylogeny/1hhsuite1.sh` | `lfat1t01:60marker_gene_phylogeny/1hhsuite1.sh` | hhsuite2 annotation and IQ-TREE phylogenies of MCP / portal / TerL |
| `pipelines/09_phylogeny/01_marker_genes/60marker_gene_phylogeny/1hhsuite_32phage.sh` | `lfat1t01:60marker_gene_phylogeny/1hhsuite_32phage.sh` | hhsuite2 annotation and IQ-TREE phylogenies of MCP / portal / TerL |
| `pipelines/09_phylogeny/01_marker_genes/60marker_gene_phylogeny/2_anno_treatment.py` | `lfat1t01:60marker_gene_phylogeny/2_anno_treatment.py` | hhsuite2 annotation and IQ-TREE phylogenies of MCP / portal / TerL |
| `pipelines/09_phylogeny/01_marker_genes/60marker_gene_phylogeny/2_anno_treatment_32phage.py` | `lfat1t01:60marker_gene_phylogeny/2_anno_treatment_32phage.py` | hhsuite2 annotation and IQ-TREE phylogenies of MCP / portal / TerL |
| `pipelines/09_phylogeny/01_marker_genes/60marker_gene_phylogeny/3_anno_summary.py` | `lfat1t01:60marker_gene_phylogeny/3_anno_summary.py` | hhsuite2 annotation and IQ-TREE phylogenies of MCP / portal / TerL |
| `pipelines/09_phylogeny/01_marker_genes/60marker_gene_phylogeny/3_anno_summary_32phage.py` | `lfat1t01:60marker_gene_phylogeny/3_anno_summary_32phage.py` | hhsuite2 annotation and IQ-TREE phylogenies of MCP / portal / TerL |
| `pipelines/09_phylogeny/01_marker_genes/MCP/0_readme.sh` | `lfat1t01:60marker_gene_phylogeny/MCP/0_readme.sh` | hhsuite2 annotation and IQ-TREE phylogenies of MCP / portal / TerL |
| `pipelines/09_phylogeny/01_marker_genes/MCP/1_iqtree_39phageMCP.sh` | `lfat1t01:60marker_gene_phylogeny/MCP/1_iqtree_39phageMCP.sh` | hhsuite2 annotation and IQ-TREE phylogenies of MCP / portal / TerL |
| `pipelines/09_phylogeny/01_marker_genes/MCP/1_iqtree_completephageMCP.sh` | `lfat1t01:60marker_gene_phylogeny/MCP/1_iqtree_completephageMCP.sh` | hhsuite2 annotation and IQ-TREE phylogenies of MCP / portal / TerL |
| `pipelines/09_phylogeny/01_marker_genes/TerL/0_readme.sh` | `lfat1t01:60marker_gene_phylogeny/TerL/0_readme.sh` | hhsuite2 annotation and IQ-TREE phylogenies of MCP / portal / TerL |
| `pipelines/09_phylogeny/01_marker_genes/TerL/1_iqtree_39phageTerL.sh` | `lfat1t01:60marker_gene_phylogeny/TerL/1_iqtree_39phageTerL.sh` | hhsuite2 annotation and IQ-TREE phylogenies of MCP / portal / TerL |
| `pipelines/09_phylogeny/01_marker_genes/TerL/1_iqtree_completephageTerL.sh` | `lfat1t01:60marker_gene_phylogeny/TerL/1_iqtree_completephageTerL.sh` | hhsuite2 annotation and IQ-TREE phylogenies of MCP / portal / TerL |
| `pipelines/09_phylogeny/01_marker_genes/portal/0_readme.sh` | `lfat1t01:60marker_gene_phylogeny/portal/0_readme.sh` | hhsuite2 annotation and IQ-TREE phylogenies of MCP / portal / TerL |
| `pipelines/09_phylogeny/01_marker_genes/portal/1_iqtree_39phagePortal.sh` | `lfat1t01:60marker_gene_phylogeny/portal/1_iqtree_39phagePortal.sh` | hhsuite2 annotation and IQ-TREE phylogenies of MCP / portal / TerL |
| `pipelines/09_phylogeny/01_marker_genes/portal/1_iqtree_completephagePortal.sh` | `lfat1t01:60marker_gene_phylogeny/portal/1_iqtree_completephagePortal.sh` | hhsuite2 annotation and IQ-TREE phylogenies of MCP / portal / TerL |
| `pipelines/09_phylogeny/02_special_genes/61special_gene_phylogeny/0_readme.sh` | `lfat1t01:61special_gene_phylogeny/0_readme.sh` | Phylogenies of the special genes (e.g. VanX) found in the new viruses |
| `pipelines/09_phylogeny/02_special_genes/61special_gene_phylogeny/1_blast_nr.sh` | `lfat1t01:61special_gene_phylogeny/1_blast_nr.sh` | Phylogenies of the special genes (e.g. VanX) found in the new viruses |
| `pipelines/09_phylogeny/02_special_genes/61special_gene_phylogeny/1_blast_nr_VanX.sh` | `lfat1t01:61special_gene_phylogeny/1_blast_nr_VanX.sh` | Phylogenies of the special genes (e.g. VanX) found in the new viruses |
| `pipelines/09_phylogeny/02_special_genes/61special_gene_phylogeny/2_iqtree.sh` | `lfat1t01:61special_gene_phylogeny/2_iqtree.sh` | Phylogenies of the special genes (e.g. VanX) found in the new viruses |
| `pipelines/09_phylogeny/02_special_genes/61special_gene_phylogeny/2_iqtree_VanX2.sh` | `lfat1t01:61special_gene_phylogeny/2_iqtree_VanX2.sh` | Phylogenies of the special genes (e.g. VanX) found in the new viruses |
| `pipelines/09_phylogeny/02_special_genes/61special_gene_phylogeny/2_iqtree_VanX3.sh` | `lfat1t01:61special_gene_phylogeny/2_iqtree_VanX3.sh` | Phylogenies of the special genes (e.g. VanX) found in the new viruses |
| `pipelines/10_experiments/01_phage_sequence/01phage_sequence/0_readme.sh` | `lfat1t01:70experiments/01phage_sequence/0_readme.sh` | Codon usage and sequence analysis of the isolated phages |
| `pipelines/10_experiments/01_phage_sequence/01phage_sequence/apc.pl` | `lfat1t01:70experiments/01phage_sequence/apc.pl` | Codon usage and sequence analysis of the isolated phages |
| `pipelines/10_experiments/01_phage_sequence/01phage_sequence/apc.py` | `lfat1t01:70experiments/01phage_sequence/apc.py` | Codon usage and sequence analysis of the isolated phages |
| `pipelines/10_experiments/02_host_phylogeny/02host_CRBC_phylogeny/01GTDB_classify.sh` | `lfat1t01:70experiments/02host_CRBC_phylogeny/01GTDB_classify.sh` | Phylogeny of the CRBC host strains (GTDB-Tk) |
| `pipelines/10_experiments/02_host_phylogeny/02host_CRBC_phylogeny/0_readme.sh` | `lfat1t01:70experiments/02host_CRBC_phylogeny/0_readme.sh` | Phylogeny of the CRBC host strains (GTDB-Tk) |
| `pipelines/10_experiments/02_host_phylogeny/02host_CRBC_phylogeny/GTDB_infer.sh` | `lfat1t01:70experiments/02host_CRBC_phylogeny/GTDB_infer.sh` | Phylogeny of the CRBC host strains (GTDB-Tk) |
| `pipelines/10_experiments/03_host_defense/03host_DF/0_readme.sh` | `lfat1t01:70experiments/03host_DF/0_readme.sh` | Defense system content of the experimental host strains |
| `pipelines/10_experiments/04_strain_prophage/04strain_prophage/0_readme.sh` | `lfat1t01:70experiments/04strain_prophage/0_readme.sh` | Prophage content of the experimental host strains |
| `pipelines/10_experiments/05_extra_validation/70_experiments/0_readme.sh` | `pkusz_CPU:70_experiments/0_readme.sh` | Additional GTDB-Tk validation of the experimental strains |
| `pipelines/10_experiments/05_extra_validation/70_experiments/GTDB_infer.sh` | `pkusz_CPU:70_experiments/GTDB_infer.sh` | Additional GTDB-Tk validation of the experimental strains |
| `pipelines/10_experiments/05_extra_validation/70_experiments/dl_GTDB_classify.sh` | `pkusz_CPU:70_experiments/dl_GTDB_classify.sh` | Additional GTDB-Tk validation of the experimental strains |
| `pipelines/10_experiments/05_extra_validation/align/GTDB_infer.sh` | `pkusz_CPU:70_experiments/res_EcoliMG1655/align/GTDB_infer.sh` | Additional GTDB-Tk validation of the experimental strains |

## What changed

* **Paths.** `/lustre/home/<user>`, `/share/home/<user>` and `/home/<user>` were replaced by
  `${PROJECT_ROOT}`, `${HOME_ROOT}`, `${SHARED_HOME}` and `${DB_ROOT}`; the `--mail-user`
  SLURM directive was removed. R scripts resolve inputs through `CS_DATA` (`CS_EXT` for the
  companion project) and write figures to `CS_WORK`.
* **Headers.** Each file gained a header with its module, purpose, manuscript figure and
  original location; nothing else in the shell/python/perl code was rewritten.
* **Excluded.** Every `*_outdated` and `*_wrong` directory, the vendored copies of
  CRISPRCasFinder and ViPTreeGen, the ICEfinder2 python bundle and the per-genome result
  files were left out; only the scripts are kept.
* **R scripts.** The two manuscripts' figure code is unchanged apart from the `setwd()`
  calls, `load()` paths and the `ggsave()`/`pdf()` destinations.
