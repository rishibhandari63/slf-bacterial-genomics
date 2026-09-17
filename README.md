# Comparative genomics of culturable bacteria from the spotted lanternfly

Analysis code for six bacterial genomes isolated from *Lycorma delicatula*
(spotted lanternfly), covering taxonomic placement, functional annotation, and
identification of genus-specific gene content in two candidate novel species.

Written for the USDA-ARS SCINet **Ceres** cluster (SLURM). Paths and module
names are Ceres-specific; see [Adapting this to another system](#adapting-this-to-another-system).

---

## The six isolates

| ID | Assignment | Size (Mb) | G+C (%) | dDDH to closest type strain |
|---|---|---|---|---|
| Sample1 | *Acinetobacter schindleri* | 3.28 | 42.5 | 73.0% |
| Sample2 | *Chryseobacterium takakiae* | 4.83 | 37.1 | 90.9% |
| Sample3 | *Sphingomonas yabuuchiae* | 4.33 | 66.0 | 73.6% |
| Sample4 | *Massilia* sp. nov. | 5.71 | 65.9 | 51.3% |
| Sample5 | *Glutamicibacter mishrai* | 3.54 | 59.5 | 88.4% |
| Sample8 | *Agromyces* sp. nov. | 4.28 | 71.3 | 33.2% |

There is no Sample6 or Sample7: eight colonies were picked and six carried
through.

---

## Pipeline

Each stage was chosen in response to the result of the one before it, so the
numbering is the run order.

### `00_setup/` — environments and databases

| Script | Purpose |
|---|---|
| `01_download_gtdbtk_db.slurm` | GTDB R220 reference data |
| `04_install_fastani_env.sh` | fastANI conda env |
| `08_install_pyani_env.sh` | pyani conda env |
| `11_install_functional_envs.sh` | eggNOG-mapper, dbCAN, antiSMASH, MGE/defense envs |
| `12_download_functional_dbs.slurm` | all seven functional databases |

### `01_taxonomy/` — what are these organisms?

| Script | Purpose |
|---|---|
| `02_checkm2.slurm` | completeness and contamination |
| `03_gtdbtk_classify.slurm` | GTDB-Tk `classify_wf` against R220 |
| `05_download_refs_and_fastani.slurm` | fetch comparators, fastANI screen |
| `06_gtdbtk_denovo.slurm` | bac120 de novo tree, 42 genomes |
| `09_run_pyani_anim.slurm` | alignment-based ANI (run with `-m ANIb`) |
| `10_patch_pyani_zerodiv.py` | patches a divide-by-zero in pyani |
| `07_plot_ani_heatmap.py` | ANI heatmap |
| `16_TYGS_dDDH_instructions.md` | dDDH via the TYGS web server (not scriptable) |

### `02_annotation/` — what do they encode?

| Script | Purpose |
|---|---|
| `13a_prokka_samples.slurm` | Prokka over the six isolates |
| `13_funct_annotation.slurm` | eggNOG-mapper + run_dbcan |
| `17_summarise_functional.py` | COG, KEGG and CAZyme summary tables |

### `03_specialised/` — biosynthesis, mobility, defense

| Script | Purpose |
|---|---|
| `14_antismash.slurm` | biosynthetic gene clusters |
| `18_summarise_antismash.py` | BGC tables from the region GenBanks |
| `15_mge_defense.slurm` | geNomad, MOB-suite, AMRFinderPlus, DefenseFinder |
| `19_summarise_mge_defense.py` | MGE / AMR / defense summary tables |

### `04_genus_specific/` — what makes the two novel species novel?

| Script | Purpose |
|---|---|
| `roary.sh` | genus pan-genome and core alignment (`-i 80 -cd 95 -g 600000`) |
| `20_roary_unique_genes.py` | joins Roary output to annotation |
| `21_pangenome_diagnostics.py` | detects clustering failure — **run this before trusting Roary** |
| `22a_preflight.sh` | checks paths and self-inclusion before the real run |
| `22_unique_genes_diamond.slurm` | the analysis that worked: direct homology search |
| `23_genus_placement.slurm` | POCP / AAI against type species of candidate genera |
| `24_pocp_aai.py` | POCP and AAI calculator |

### `05_figures/`

matplotlib scripts. `build_fig1.py` composites four panels into a single
7.5 × 8 in figure, embedding externally-produced PDFs (iTOL tree, pyani
heatmap) as vector while redrawing all labelling at figure-appropriate size.

---

## Things that went wrong, and why they matter

These cost real time. They are documented in-script, and worth reading before
adapting any of this.

**Roary cannot cluster across a genus at default settings.** At `-i 95` it
recovered a core genome of 26 genes from 153 *Agromyces* genomes and called 51%
of one proteome "unique" — an over-splitting artefact, not a result.
`21_pangenome_diagnostics.py` detects this: if most clusters sit in ≤2 genomes
and no single genome is an outlier, the clustering has failed rather than the
input. For gene-content questions we replaced it with direct DIAMOND search
(`22_unique_genes_diamond.slurm`); for the core alignment we lowered `-i` to 80
and raised `-g` to 600000.

**Prokka reassigns locus tags on every run.** Two Prokka runs of the same genome
give identical gene calls under different identifiers, so joining DIAMOND output
to an eggNOG table from a *different* Prokka run silently returns 0% matches.
Script 22 checks for this and says so.

**An isolate inside its own reference set produces a silent zero.** Every protein
gets a perfect self-hit and nothing looks unique. Script 22 carries four
independent guards plus a post-hoc detector that aborts if ≥80% of proteins have
a full-length 100% identity hit.

**`set -o pipefail` plus an early-exiting pipeline consumer gives SIGPIPE.**
`zcat file | awk '{print; exit}'` kills a job in three seconds with exit 141 and
no error message. Use `awk` without the early exit, and `find -print -quit`
instead of `find | head -1`.

**`[ x -ne y ] && echo ...` as a script's last statement inverts its exit
status.** When the counts match — the good case — the test is false, the script
exits 1, and SLURM reports a successful job as FAILED. Use a full `if` block.

**COG category V is not an anti-phage proxy.** It lumps restriction-modification
components with ABC-type efflux transporters. The isolate with the highest
category V representation had zero defense systems by DefenseFinder.

---

## Adapting this to another system

Every script sets `PROJECT=/project/kanglab/rishi.bhandari` near the top; change
that first. Beyond paths you will need to replace the `module load` lines (Ceres
provides GTDB-Tk, CheckM2 and Roary as modules; everything else is installed
into conda envs by `00_setup/`), and the SLURM account, partition names and
resource requests.

Scripts are deliberately not parameterised beyond `PROJECT` and, in script 22,
`SAMPLE` — which is set at submit time so the file never needs editing:

```bash
sbatch --export=ALL,SAMPLE=Sample4 22_unique_genes_diamond.slurm
```

---

## Software

Versions used are recorded in each script header. Principal tools: CheckM2,
GTDB-Tk 2.6.1 (GTDB R220), fastANI, skani, pyani, FastTree, TYGS/GGDC 4.0,
Prokka, Roary 3.13.0, eggNOG-mapper 2.1.15 (eggNOG 5.0), DIAMOND, run_dbcan 5,
antiSMASH 8.0.4, geNomad, MOB-suite 3.1.9, AMRFinderPlus 3.12.8,
DefenseFinder 3.0.0, RAxML. Assemblies were produced by Plasmidsaurus using
Filtlong, Autocycler, Flye, hifiasm, Plassembler, dnaapler and Medaka.

## Data availability

Genome assemblies: *[add BioProject / BioSample / GenBank accessions]*.
Reference genome accessions used in each comparison are listed in
`reference_genomes_used.txt`, written by script 22.

## Citation

*[Add once the manuscript is submitted.]*

## License

MIT — see `LICENSE`.
