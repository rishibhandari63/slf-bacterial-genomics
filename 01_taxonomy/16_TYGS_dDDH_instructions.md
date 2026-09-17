# dDDH via TYGS (web) — step-by-step

TYGS (Type Strain Genome Server, https://tygs.dsmz.de/) computes digital
DNA-DNA hybridization (dDDH, via GGDC) between your genomes and the type
strains of related species, plus a genome-based tree — exactly the evidence
reviewers expect next to ANI for a novel-species claim. It's a free web
service run by DSMZ; there is nothing to install.

## Submitting

1. Get your 6 genome FASTAs onto your local machine (from Ceres:
   `scp rishi.bhandari@ceres.scinet.usda.gov:/project/<project>/phylo/ANI/Sample*.fasta .`).
2. Go to https://tygs.dsmz.de/ → "Start your TYGS analysis".
3. Upload all 6 FASTA files in ONE submission (they'll be analyzed
   together and against type strains automatically — you do NOT need to
   upload the reference genomes yourself; TYGS picks type strains from its
   own database).
4. Enter your email address; submit. Runs typically take a few hours to a
   day depending on queue. Results arrive as an email link.

## Reading the results

- **dDDH (d4 / formula 2) is the value to report** — it's the
  recommended formula for draft genomes (it's independent of genome
  length and robust to incomplete assemblies).
- Species boundary: **dDDH < 70%** against every described type strain =
  supports a novel species. Also note the accompanying 95% confidence
  interval TYGS reports.
- Subspecies boundary: 79% is the conventional cutoff, occasionally
  worth mentioning for borderline cases.
- TYGS also returns a GBDP tree with branch support — you can cite it as
  corroborating the GTDB-Tk/core-genome trees, but keep your RAxML
  core-genome tree as the primary phylogeny figure (it's built from far
  more signal).
- Check the "Pairwise comparisons" table for each sample: if TYGS's
  closest type strain differs from GTDB-Tk's closest reference, say so
  and explain (usually a genome available in one database but not the
  other) — reviewers notice mismatches between tables.

## Reporting in the Methods

Suggested sentence (fill versions/dates from the TYGS result page footer,
which lists the exact GGDC/TYGS versions used for your run):

> Digital DNA–DNA hybridization (dDDH) values between the isolates and
> the type strains of related species were calculated with the Type
> (Strain) Genome Server (TYGS; Meier-Kolthoff & Göker, 2019), using the
> Genome BLAST Distance Phylogeny method and recommended formula d4;
> species and subspecies boundaries were assessed using the standard 70%
> and 79% dDDH thresholds, respectively.

References to add:

- Meier-Kolthoff, J.P., & Göker, M. (2019). TYGS is an automated
  high-throughput platform for state-of-the-art genome-based taxonomy.
  Nature Communications, 10, 2182.
- Meier-Kolthoff, J.P., Auch, A.F., Klenk, H.-P., & Göker, M. (2013).
  Genome sequence-based species delimitation with confidence intervals
  and improved distance functions. BMC Bioinformatics, 14, 60.
  (the GGDC/dDDH method itself)

Verify both citations against the publisher record before submission.
