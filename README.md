# twomarinedemospongesGenome
Genome assemblies and annotations of two marine demosponges
- `1_assembly_polishing`: retained Flye assembly, mixed-read and Illumina-only NextPolish runs, purge_dups, Kraken2, Tiara and BlobToolKit workflows.
- `2_contamination_curation`: multi-evidence curation, small decision lists and final NJ-lan review filtering code.
- `3_repeat_annotation`: precursor RepeatModeler records, final-assembly RepeatMasker and repeat-landscape processing.
- `4_structural_annotation`: final-input preparation, BRAKER4 configurations and annotation consistency audit.
- `5_functional_annotation`: eggNOG-mapper, InterProScan, Swiss-Prot and gene-level functional summaries.
- `6_quality_control`: BUSCO, QUAST, ONT coverage and Merqury workflows.
- `7_figures`: source-table preparation, manuscript plotting scripts.
