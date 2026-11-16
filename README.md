# Common transposable-element insertions show purifying selection but no climate signal in wild teosinte and domesticated maize

Analysis code and processed transposable-element (TE) resources for a study testing whether
population-level TE-insertion polymorphisms track climate in wild teosinte (*Zea mays* ssp.
*mexicana*, 17 populations) and domesticated maize (*Z. mays* ssp. *mays*, 14 sympatric
populations), benchmarked against genome-wide SNPs in an identical, structure-aware framework.

## Data sources (not redistributed here)

- Whole-genome resequencing reads: NCBI BioProject PRJNA657016 (HILO project; Calfee et al. 2021).
- SNP allele frequencies and ancestry estimates: Figshare doi:10.6084/m9.figshare.16641799
  (per-population `.mafs.gz` files, `ind_metadata_and_genomewide_ancestry.csv`,
  `K3_alphas_by_symp_pop.txt`; B73 RefGen_V4 coordinates).
- Reference genome: B73 RefGen_V5 (Zm-B73-REFERENCE-NAM-5.0, MaizeGDB) with its structural TE
  annotation, for TE calling; B73 RefGen_V4 gene annotation (Ensembl Plants release 50) for the
  SNP-coordinate analyses; B73 RefGen_V5 gene annotation (Ensembl Plants release 57) for the
  TE-to-gene distances.
- TE library: Maize TE Consortium (MTEC) consensus library (maizeTE02052020) distributed with EDTA.
- Climate: WorldClim v2.1 bioclimatic variables (2.5 arc-minute); elevation: SRTM 30 m
  (Open Topo Data API, dataset srtm30m).
- Known inversion coordinates: HILO project repository (`knownInv_v4_coord.txt`).

## Repository layout

```
scripts/teosinte/   upstream processing of the 17 mexicana populations (HPC cluster, SLURM)
scripts/maize/      the same pipeline for the 14 sympatric maize populations
scripts/shared/     SNP matrices, ancestry covariates and site elevation from public data
analysis/           genotype-environment association, figures and tables (run from the repository root)
data/te_teosinte/   teosinte TE-insertion calls and frequency matrices (the reusable resource)
data/te_maize/      maize TE-insertion frequency matrices
data/metadata/      per-population climate, elevation, ancestry/introgression and sample sheets
data/tables/        Tables 1 and 2 and the supplementary tables
data/               site lists of the 80,000 analysed SNPs (teosinte and maize)
figures/            Figures 1 to 4 and Figures S1 to S6
```

### Upstream pipeline (scripts/teosinte, scripts/maize)

| Step | Script | Tool |
|---|---|---|
| environment | `00_setup_env.sh`, `environment.yml` | conda |
| sample discovery and sample sheet | `01_find_data.sh`, `02_build_sample_table.sh` (teosinte); `00_build_sample_sheet.py` (maize) | Entrez Direct, ENA |
| climate | `03_extract_climate.R` | terra, geodata (WorldClim v2.1) |
| reference | `04_get_reference.sh`, `05_build_te_reference.sh` | BEDTools maskfasta, bwa index |
| download, QC, alignment | `06_download_align.sh` (teosinte); `01_download.sh`, `02_align.sh` (maize) | SRA Toolkit, fastp, bwa mem -M, SAMtools |
| population pools | `07_merge_pools.sh`; `03_merge_pools.sh` | samtools merge |
| TE-insertion frequencies | `08_ppileup.sh`, `09_te_frequency.sh`; `04_ppileup.sh`, `05_te_frequency.sh` | PoPoolationTE2 v1.10.04 |
| frequency matrix | `10_build_te_matrix.R`; `06_build_te_matrix.R` | R data.table |
| initial GEA of the teosinte matrix | `11_gea.R` | vegan, LEA |

PoPoolationTE2 parameters: ppileup map-qual 15; subsamplePpileup target coverage 10
(teosinte) or 8 (maize); identifySignatures mode separate, min-count 2; filterSignatures
max-otherte-count 2, max-structvar-count 2; pairupSignatures min-distance -200, max-distance 300.
Calls of the same TE family within 100 bp are clustered into one locus; satellite classes
(knob, Cent, rDNA, subtelomere) are removed; undetected populations are scored 0; polymorphic loci
have non-zero variance and are detected in at least two populations.

The cluster scripts read the project directory from the `PROJ` environment variable (default: the current directory).

### Public-data inputs (scripts/shared)

- `build_snp_matrix.py` builds the SNP frequency matrices from the Figshare `.mafs.gz` files:
  the non-reference allele frequency at sites present in every population (374,335 teosinte;
  421,868 maize), subsampled to 80,000 sites (the maize matrix restricted to sites polymorphic
  across populations; the teosinte matrix built without that filter, `--keep-monomorphic`, so that
  1,742 of its sites are monomorphic and are dropped by the variance filters of the analyses).  The
  site lists `data/snp_sites_teosinte.tsv` and `data/snp_sites_maize.tsv` reproduce the analysed
  subsets exactly (`--sites`).
- `ancestry_covariates.py` builds the per-population maize-ancestry (teosinte) and
  mexicana-introgression (maize) covariates from the Figshare ancestry files.
- `site_elevation.py` extracts the SRTM 30 m elevation at each sampling coordinate.

### Analysis (analysis/)

All scripts are run from the repository root (`python3 analysis/<script>.py`) and share
`analysis/common.py` (data loaders, partial RDA with permutation test, per-locus partial
correlation, Benjamini-Hochberg FDR, genomic inflation factor).

| Script | Content |
|---|---|
| `01_teosinte_te_rda.py` | TE insertions versus climate PCs and temperature, full and fully observed sets, per-locus FDR |
| `02_teosinte_snp_rda.py` | SNP benchmark: positive control, SNP versus climate, per-SNP FDR |
| `03_snp_count_matched.py` | SNP RDA on random subsets matched to the TE locus count |
| `04_snp_robustness.py` | chr9 region excluded or alone; leave-one-population-out |
| `05_snp_peaks.py` | genomic distribution of significant SNPs; enrichment in chr9:120-150 Mb and inv4m |
| `06_snp_noisematch.py` | noise-matched SNP control (Table S1) |
| `07_structure_lambda.py` | collinearity of temperature with structure; lambda_GC; structure-conditioned RDA |
| `08_contrast_under_structure.py` | SNP versus TE contrast under matched conditioning (Table S2) |
| `09_chr9_characterisation.py` | internal linkage and leading axis of the chr9 candidate region |
| `10_chr9_annotation.py` | gene density of the chr9 region (RefGen_v4) and overlap with known inversions |
| `11_teosinte_te_introgression.py` | TE frequencies versus the maize-ancestry gradient |
| `12_lfmm2_teosinte.R`, `run_lfmm2.sh` | LFMM2 (LEA) of SNPs and TEs on temperature, K = 1 to 5 (Table S3) |
| `13_maize_gea.py` | maize replication: contrast, per-locus FDR, full versus fully observed, leave-one-out, lambda_GC (Table S6) |
| `14_maize_noisematch.py` | maize noise-matched SNP control (Table S4) |
| `15_maize_snp_peaks.py` | maize significant SNPs by chromosome and inside inv4m |
| `16_lfmm2_maize.R` | maize LFMM2 (Table S5) |
| `17_gene_distance.py` | TE frequency by distance to the nearest gene (RefGen_v5), Figure 4 |
| `18_figures_teosinte.py`, `19_manhattan_teosinte.py` | Figures 1 to 3, S1 and S2 |
| `20_figures_maize.py` | Figures S3 to S6 |
| `21_table1_teosinte.py`, `22_table2_maize.py` | Tables 1 and 2 |

The gene annotations and the inversion table are downloaded into `data/annotation/` on first use.
The LFMM2 scripts require R with the LEA package and the SNP matrices.

### Data (data/)

- `te_hierarchy.txt`: the TE hierarchy used by PoPoolationTE2 (1,477 library sequences, 23
  families, 8 orders; family = order_superfamily).
- `te_teosinte/hilo.teinsertions`, `te_maize/maize.teinsertions`: PoPoolationTE2 insertion calls
  (sample index, chromosome, position, family, order, support, strand, frequency);
  `pop_order.txt` in each directory maps the sample index to the population pool.
- `te_teosinte/te_freq_matrix.tsv`, `te_maize/te_freq_matrix.tsv`: polymorphic TE insertions
  (2,871 teosinte, 3,352 maize) with per-population frequencies, number of populations in which
  the insertion was detected (`nDetected`) and among-population variance; `te_freq_matrix_all.tsv`
  holds all clustered loci.
- `te_teosinte/gea_results.tsv`: per-locus statistics of the 2,871 teosinte TE insertions from
  `scripts/teosinte/11_gea.R` (RDA1 loading; LFMM2 p-value on climate PC1 with K = 2; candidate flags).
- `metadata/`: `climate_table_hilo.tsv` (WorldClim bio1 to bio19), `elevation_srtm30m.tsv`,
  `maize_ancestry_perpop.tsv`, `mexicana_ancestry_maize_perpop.tsv`,
  `sample_sheet_hilo_mexicana.tsv`, `sample_sheet_hilo_maize.tsv`.
- `tables/`: `Table1_population_metadata.tsv`, `Table2_maize_population_metadata.tsv` and
  `TableS1.tsv` to `TableS6.tsv` (captions in `tables/supplementary_tables.md`).

## Software

conda environment (`scripts/teosinte/environment.yml`): SRA Toolkit 3.4.1, fastp 1.3.5,
FastQC 0.12.1, MultiQC 1.35, BWA 0.7.19, SAMtools 1.23.1, BEDTools 2.31.1, PoPoolationTE2
1.10.04, EDTA 1.8.5 (TE library), R 4.3.3 (LEA 3.14, terra, geodata, data.table), Python 3
(NumPy, SciPy, Matplotlib).

## Citation

If you use this code or the TE-insertion resources, please cite the associated manuscript and the
underlying data (Calfee et al. 2021, PLOS Genetics, doi:10.1371/journal.pgen.1009810).

## License

The code and the processed resources in this repository are released under the MIT License (see `LICENSE`).
