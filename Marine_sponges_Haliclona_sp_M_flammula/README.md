# MarineSpongesGenome

Workflow code used to assemble, curate, annotate and validate the genomes of
*Haliclona* sp. and *Mycale flammula*.

| Species | Sample ID | Historical aliases in scripts and files |
|---|---|---|
| *Haliclona* sp. | `LHT-04` | `lht`, `LHT`, `LHT04` |
| *Mycale flammula* | `NJ-lan` | `lan`, `LAN`, `Lanm`, `njlan` |

## Reproducing the pipeline

Scripts are organized by analysis stage. Before running them, adapt the cluster
paths, modules, conda environments, database locations and input files to the
target system. Raw reads, assemblies, annotation releases, third-party software
and reference databases are not included.

The retained workflow starts from ONT and Illumina reads and covers assembly,
polishing, contamination curation, annotation, final validation and figure
generation. LHT-04 used a valid length-ranked 11.8-Gb ONT subset. For NJ-lan,
the retained `sub118` directory and `Lanm.longest_11.8Gb.fq.gz` filename are
legacy labels: that file contains all 10,186,340 ONT reads (27,612,556,815 bp),
not an 11.8-Gb subset. The retained workflow used Flye, NextPolish and
purge_dups; Pilon was not used.

The initial NextPolish pass used paired Illumina reads plus the complete ONT
reads; the final pass used paired Illumina reads only. The recovered `run.cfg`
files refer to site-local `sgs.fofn` and `lgs.fofn` files, which should list
`LHT_R1.fq.gz`/`LHT_R2.fq.gz` and `LHT04.ont.fq.gz` for LHT-04, or
`3_R1.fq.gz`/`3_R2.fq.gz` and `Lanm.ont.fq.gz` for NJ-lan.

## Repository structure

- `1_assembly_polishing`: retained Flye assembly, mixed-read and Illumina-only NextPolish runs, purge_dups, Kraken2, Tiara and BlobToolKit workflows.
- `2_contamination_curation`: multi-evidence curation, small decision lists and final NJ-lan review/filtering code.
- `3_repeat_annotation`: precursor RepeatModeler records, final-assembly RepeatMasker and repeat-landscape processing.
- `4_structural_annotation`: final-input preparation, BRAKER4 configurations and annotation consistency audit.
- `5_functional_annotation`: eggNOG-mapper, InterProScan, Swiss-Prot and gene-level functional summaries.
- `6_quality_control`: BUSCO, QUAST, ONT coverage and Merqury workflows.
- `7_figures`: source-table preparation, audit code and manuscript plotting scripts.

## Important provenance and limitations

- The locked final assemblies contain 844 records and 115,734,075 bp for LHT-04, and 7,117 records and 291,889,606 bp for NJ-lan. Their SHA-256 values are `026d8b2e3e9888de0763d42ce226fcf7c63062a145e6ad8478e2e8c1725b2499` and `cea9b34517164f02103111a3a1dfa7cfd65ec95a61c224da1901036a4cca120a`, respectively.
- The LHT-04 RepeatModeler library was built from the 856-contig curation precursor; the NJ-lan library was built from the 7,366-contig precursor. RepeatMasker was then rerun on both locked final assemblies.
- The LHT-04 856-to-844 read-supported curation retained 844 primary records, placed 10 records in a non-host sidecar, placed four in quarantine and split two source contigs. The final builder and sidecar files were not preserved with the archived code; the preceding 952-to-856 evidence builder and decision lists are included.
- The archived records do not preserve the FCS-GX version, GX database build, expected-host division or action report. They also do not preserve the final BLAST+ version, NCBI nt snapshot, search date or exact BLASTN command; these details must not be inferred from another sample set.
- InterProScan completed with the included Python 3.11 retry script. Earlier failed launchers are intentionally excluded.
- Merqury v1.3 with Meryl v1.4.1 used canonical 19-mers. QV was 35.81 for LHT-04 and 23.77 for NJ-lan. Reported reliable-k-mer completeness used multiplicity `>7` and `>3`, respectively; automatically inferred thresholds 150 and 101 were not used for completeness reporting.

The scripts retain historical `lht` and `lan` paths so that the executed
analysis remains traceable. Review and adapt all site-specific paths before
execution.

## License

Project-authored code and documentation are released under the MIT License;
see `LICENSE`. Third-party software and databases retain their providers'
terms.
