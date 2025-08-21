# Selection and adaptive introgression scans for the PIBv1 manuscript

## Scripts and analyses by Daniela Tejada Martinez, PhD (@dtejadam)

## Background

We performed selection scans on 17 target populations from Near Oceania using two different selection statistics:

1. The normalized version (PBSn1) of the Population Branch Statistic (PBS, [Yi et al. 2010](https://doi.org/10.1126/science.1190371))
introduced in [Malaspinas et al. 2016](https://doi.org/10.1038/nature18299)
2. Cross-Population Extended Haplotype Homozygosity (XP-EHH, [Sabeti et al. 2007](https://doi.org/10.1038/nature06250))

PBSn1 was evaluated in 20 SNV sliding windows with a 5 SNV step using the [Hudson et al. 1992](https://doi.org/10.1093/genetics/132.2.583),
while XP-EHH was evaluated on a per-SNV basis. We then normalized results by converting statistic values to the corresponding quantile rank
in each empirical distribution, and then combined across selection statistics by calculating a [Fisher's combined score](https://en.wikipedia.org/wiki/Fisher%27s_method) (FCS)
from one minus the quantile rank of each statistic. (N.B. There are some important technical caveats for interpreting the resulting FCS.)
PBS and XP-EHH were calculated using the Python package [scikit-allel](https://pypi.org/project/scikit-allel/), and FCS was evaluated with
custom R code.

We then intersected genomic windows that were outliers in the FCS distribution with high frequency archaic core haplotypes (see the supplement for more
details) in order to establish a set of genomic regions with reliable evidence for adaptive introgression.

These genomic regions were then annotated with the closest protein-coding gene in the GENCODE build 38 annotation to identify likely functional
effects.

## Pipeline

### Dependencies

This pipeline consists of a series of bash scripts for SLURM submissions, as well as some Python and R scripts run by the bash scripts.
Be sure to have the following dependencies installed before running the pipeline:

- [BEDtools](https://github.com/arq5x/bedtools2)
- Python 3
  - NumPy
  - SciPy
  - Pandas
  - scikit-allel
- R
  - tidyverse
  - data.table
  - ggrepel
  - qqman

These dependencies can generally be installed with conda using the following commands:

```bash
conda create -n selectionscans -c conda-forge -c bioconda "python==3" \
   numpy scipy pandas "scikit-allel" \
   bedtools r \
   "r-tidyverse" "r-data.table" "r-ggrepel" "r-qqman"
```

Be sure that the conda `selectionscans` environment is activated in each script for the above approach to work.

Alternatively, you can do the same thing with your base conda environment, or temporarily adjust your `.bash_profile`
to activate the `selectionscans` environment, or don't bother with conda and just make sure you have Python 3, R, and
BEDtools with the appropriate packages installed in your base environment.

### Prerequisite files and directory structure

Aside from these dependencies, certain input files and a specific directory structure are necessary for the pipeline to work properly.
This Github repository is intentionally structured in such a way, so all you need to do is add the appropriate files in the corresponding
directories.  These key files include:

- `key_files/gencode.v38lift37.annotation.codingGENES_sorted.bed`: A BED of protein-coding genes generated from the GENCODE build 38 annotation
lifted to b37 in GFF3 format using [gencodeGFFtoBED.awk](https://github.com/YourePrettyGood/HumanPopGenScripts/RefPrep/gencodeGFFtoBED.awk)
- `key_files/gmap_files/chr*.b37.gmap`: Genetic map files in b37 reference coordinates with three tab-separated columns and a header 
(Physical position `pos` in bp, chromosome `chr`, and genetic position `cM` in centiMorgans)
- `key_files/header_file.txt`: 
- `key_files/hs37d5_wChr.genome`: A file derived from the first two columns of the .fai index of the hs37d5 reference assembly, but with "chr" prepended to
each chromosome (thus consisting of two columns: chromosome and chromosome length, but no header)
- `key_files/PIBv1_Sprime_putative_adaptive_introgressed_corehaps_OCNnoVanuatu_wGeneList_lumpedpopsrenamed.bed`: A BED-like file of high-frequency archaic
introgressed core haplotypes that gets intersected with the selection scan outlier genomic regions
- `vcf_files/*.vcf.gz`: One bgzipped VCF per autosome with accompanying tabix indices
- `run_gaia_selection/Oceanic_pops_FCS/chr_info.tsv`: A file derived from the .fai index of the hs37d5 reference assembly, but truncated to only the autosomes
- `run_gaia_selection/Oceanic_pops_PBS/pops_oceania/*.txt`: Files listing sample IDs corresponding to each population (one ID per line), where `pop1_`
indicates a target population, `pop2_` indicates the reference population, and `pop3_` indicates the outgroup population(s)
- `run_gaia_selection/Oceanic_pops_XPEHH/pops_oceania/*.txt`: Files listing sample IDs corresponding to each population (one ID per line), where `pop1_`
indicates a target population and `pop2_` indicates the reference population (`pop3_` is ignored for XP-EHH)

It is worth noting that the scripts are currently written assuming the VCFs are based on the hs37d5/b37 human reference genome build.
This affects the content of `key_files/hs37d5_wChr.genome`, `key_files/gmap_files/*`, and `run_gaia_selection/Oceanic_pops_FCS/chr_info.tsv`.
Furthermore, this affects the BED file of protein-coding genes `key_files/gencode.v38lift37.annotation.codingGENES_sorted.bed`, which is
needed for the annotation steps of the pipeline.

### High-level workflow

The selection scan pipeline consists of three main stages:

1. Running PBS on each population
2. Running XP-EHH on each population (can be run in parallel with stage 1)
3. Combining results of stages 1 and 2 into windowed FCS

At the end of each stage, genomic regions are filtered at different thresholds for outliers and these outlier regions are annotated
with the closest protein-coding gene. Furthermore, for each annotated overlapping or nearby protein-coding gene, the genomic window
with the strongest selection signal is annotated (see the "topgene" column).

The scripts for these three stages can be found in the three subdirectories of `run_gaia_selection`:

1. 
  a. `Oceanic_pops_PBS`
  b. `Oceanic_pops_XPEHH`
2. `Oceanic_pops_FCS`

### Quick start

Once you've installed the dependencies and prepared the key files (including VCFs), you should adjust the SLURM submission scripts
in each subdirectory of `run_gaia_selection` (i.e. the `*.sh` scripts) appropriately so that the SLURM array indices match the size
of your data. Some scripts don't use SLURM arrays, so leave those be. Others have a comment indicating whether the SLURM array should
be per-population (e.g. for the 17 populations in the paper, `#SBATCH -a 0-16`) or per-population AND per-autosome (e.g. for the 17
populations in the paper and 22 autosomes, `#SBATCH 0-373`, since 17*22=374 and we're using 0-based task indices so that Python is happy).
The scripts are currently set the way they were for the analyses in the paper.

Within each stage, scripts should be run in the following order.

> [!TIP] Stages 1a and 1b can be run simultaneously, as they are independent of each other. Stage 2 must be run after
> both Stage 1a and 1b are complete.

Stage 1a (PBS):
1. `Oceanic_pops_PBS/run_pbs.sh`
2. `Oceanic_pops_PBS/run_cat_populations.sh` (bash script, not a SLURM submission script)
3. `Oceanic_pops_PBS/run_bedtools_closest_sort_2024.sh`
4. `Oceanic_pops_PBS/run_topGENES_quantiles_final_2024.sh`
5. `Oceanic_pops_PBS/run_quantile_files.sh`
6. `Oceanic_pops_PBS/run_bedtools_intersect_introgressed_tracks.sh` (optional)

Stage 1b (XP-EHH):
1. `Oceanic_pops_XPEHH/run_xpehh.sh`
2. `Oceanic_pops_XPEHH/run_cat_populations.sh` (bash script, not a SLURM submission script)
3. `Oceanic_pops_XPEHH/run_bedtools_closest_sort_2024.sh`
4. `Oceanic_pops_XPEHH/run_topGENES_quantiles_final_2024.sh`
5. `Oceanic_pops_XPEHH/run_quantile_files.sh`
6. `Oceanic_pops_XPEHH/run_bedtools_intersect_introgressed_tracks.sh` (optional)

Stage 2 (FCS):
1. `Oceanic_pops_FCS/run_final-stats_HMP_pbs-xpehh_2024.sh`
2. `Oceanic_pops_FCS/run_FCS.sh`
3. `Oceanic_pops_FCS/run_plot_manhattan_Fcs_loop.sh`

> [!WARNING] Making the final Manhattan plots may take a long time, since each plot is generated sequentially, and there
> are a lot of points to render genome-wide.

### Manual pipeline sanity checks


