# GATK-SV-inspired Nextflow pipeline

This directory holds a from-scratch Nextflow reimplementation of the parts of
[GATK-SV](../README.md) we actually need, rather than a line-by-line port of
the WDL. This document is the design plan; implementation is tracked
separately.

## Why not port the WDL directly

We looked at porting `wdl/` mechanically to Nextflow and concluded it's a
large, multi-month undertaking dominated by plumbing that doesn't apply to
us:

- 116 WDL files (~36K lines), 310 tasks, 109 scatter blocks, 4-5 levels of
  nested sub-workflow calls, hundreds of pass-through runtime-override
  params.
- A genuine non-parallel join (`CombineBatches`) that requires every batch in
  a cohort to be gathered before downstream steps can run.
- Cromwell/GCS-specific behaviour (`localization_optional` streaming reads,
  `glob()`-based dynamic shard counts) with no clean Nextflow equivalent.
- The pipeline is architected for a fixed, known-upfront cohort (biobank
  joint-calling) or a single static reference panel (single-sample mode).
  Neither matches our situation: heterogeneous projects, data arriving
  incrementally, sometimes as singletons and sometimes as trios.

Instead, we're building a smaller pipeline inspired by GATK-SV's core idea —
run an ensemble of SV/CNV callers, then harmonise their output into one
genotyped, annotated, seqR-importable VCF — using:

- SV/CNV callers from nf-core where available (otherwise thin wrappers
  around the same tools GATK-SV uses: Manta, Wham, Scramble; depth callers
  cn.MOPS / GATK-gCNV).
- GATK's `SVCluster` as the harmonisation/clustering engine (the same tool
  GATK-SV itself uses at every scale, from within-batch clustering to
  cross-batch cohort merging).
- GATK's `SVAnnotate` for functional/AF annotation.

## How GATK-SV's harmonisation actually works

Two distinct mechanisms do the "combine multiple tools into one VCF" job:

1. **Clustering** (`SVCluster`): a deterministic geometric merge of records
   across input VCFs (from different callers and/or different samples) that
   overlap in position/size/type/strand within tunable thresholds. Not
   statistical, nothing to train — the same operation whether merging 2
   samples or 2,000. This is what makes clustering incremental-friendly: a
   new sample's raw calls can be clustered directly against an
   already-clustered site list without redoing the whole thing.

2. **Statistical genotyping/filtering**: needs cohort context to calibrate
   (variant metric cutoffs, GQ recalibration model). In GATK-SV these are
   *trained artifacts* (`rf_cutoffs`, `gq_recalibrator_model_file`) that are
   produced once and then reused as plain file inputs — not retrained on
   every run. GATK-SV's own single-sample mode relies on exactly this: a
   frozen reference panel bundle (evidence matrices, gCNV model, clustered
   site VCF, trained cutoffs) that a new case sample is merged against.

We're adopting the same split.

## Design: panel-based, two-tier pipeline

Split into two Nextflow pipelines sharing modules:

### 1. `build-panel`

Run once per project initially, and re-run manually later when we judge
there are enough new samples to be worth it (see "Panel refresh policy"
below).

- Input: all currently-available samples **for one project** (BAMs/CRAMs).
  Panels are per-project, not global — different projects use different
  sequencing protocols and are subject to different ethics approvals, so
  their samples must not be pooled.
- Steps:
  - Run SV/CNV callers per sample (Manta/Wham/Scramble-equivalent; gCNV
    training in cohort mode; cn.MOPS).
  - Collect PE/SR/SD/coverage evidence per sample, merge into panel-wide
    matrices.
  - `SVCluster` across the panel to produce a merged site VCF.
  - Train/derive genotyping cutoffs (see "Where to simplify" below — this
    step may end up much lighter than GATK-SV's random-forest + ML GQ
    recalibration).
- Output: a versioned **panel bundle** (see "Panel bundle contents"), stored
  per project, e.g. `panels/<project>/vN/`.

### 2. `genotype-new-sample`

Run per arrival — one sample, one trio, or a batch of 10-30 singletons — the
fast incremental path.

- Input: new sample(s) for a project + that project's pinned panel bundle
  version.
- Steps:
  - Run the same per-sample caller/evidence steps as above.
  - Run gCNV in case mode against the panel's pre-trained model tar (no
    retraining).
  - Concatenate the new sample's evidence onto the panel's frozen matrices
    (file/array concatenation, not recomputation).
  - `SVCluster` the new sample's raw calls against the panel's merged site
    VCF, producing a VCF with both old and new genotypes anchored to the
    same site backbone.
  - Genotype using the panel's frozen cutoffs.
  - `SVAnnotate` + allele-frequency annotation.
  - Track case-only novel calls that don't match any panel site (mirrors
    GATK-SV single-sample's `non_genotyped_unique_depth_calls` output) —
    input to the panel refresh decision.
- Output: a VCF scoped to this arrival's samples, genotyped against the
  shared site backbone, ready for seqR import. Because every arrival
  clusters against the same frozen site list, site IDs/coordinates stay
  consistent across arrivals from the same panel version.

### Panel bundle contents

Per project, per version, produced by `build-panel` and consumed by
`genotype-new-sample`:

| Artifact | Replaces (vs. recomputing from scratch) |
|---|---|
| Merged PE/SR/SD/coverage matrices | Batch-wide evidence recomputation |
| gCNV model tar (+ ploidy model) | Retraining the depth model |
| Panel CNV beds (cn.MOPS) | Panel's own depth calls, merged with new sample's |
| Merged/clustered site VCF | Re-clustering the whole panel from scratch |
| Panel per-caller raw VCF tars | Inputs `SVCluster` needs alongside the merged site VCF |
| Trained genotyping cutoffs | Re-deriving metric cutoffs per run |
| Sample list, PED | Cohort context (background rates, sex-aware logic) |

### Trios

No attempt to replicate GATK-SV's simultaneous cross-sample evidence
pooling. Each trio member is clustered/genotyped independently against the
project's panel (same as a singleton), then inheritance/de novo analysis is
done as a downstream step over the resulting joint VCF — consistent
genotypes across family members in one VCF is what seqR's trio-based
filtering needs; it doesn't require GATK-SV-style joint calling to get
there.

### Panel refresh policy

Simple, manual, for now:

- Initial `build-panel` run per project uses **all currently available
  samples** for that project.
- Re-run `build-panel` **manually**, per project, when we judge there are
  enough new samples accumulated since the last panel version to be worth
  the cost of a rebuild. No automatic trigger (e.g. count/time threshold) is
  implemented yet — the case-only novel-call tracking from
  `genotype-new-sample` is the intended future input to that decision, but
  isn't wired up to anything automatic.
- Old panel versions are kept; a project's in-flight `genotype-new-sample`
  runs pin a specific panel version so a refresh doesn't retroactively
  change coordinates/site IDs for already-processed samples.

## Inputs

### Alignment format

All GATK-SV callers and evidence-collection tools consume aligned, indexed
BAM or CRAM (hg38) directly — nothing takes FASTQ as a pipeline input. (The
one place `fastq` appears in the upstream WDL is internal to the Scramble
task, which uses `samtools fastq` to extract soft-clipped reads from the BAM
in-memory for a BWA realignment step on DRAGEN-aligned inputs — not an
external input.) Some tools are two hops from the BAM: gCNV/cn.MOPS consume
read-depth counts (`CollectReadCounts` output), and Scramble consumes both
the BAM and the coverage counts plus Manta's VCF (for realignment around
candidate MEI sites), so it runs downstream of coverage collection and
Manta rather than independently from raw reads. None of this affects the
sample sheet, which only needs to supply the BAM/CRAM per sample — these
intermediate hops are internal to the pipeline.

### Sample sheet

Tab-separated (TSV), one row per sample:

```text
family_id	sample_id	bam	bai
FAM001	SAMPLE_A	/data/project1/SAMPLE_A.bam	/data/project1/SAMPLE_A.bam.bai
FAM001	SAMPLE_B	/data/project1/SAMPLE_B.bam	/data/project1/SAMPLE_B.bam.bai
FAM001	SAMPLE_C	/data/project1/SAMPLE_C.cram	/data/project1/SAMPLE_C.cram.crai
SINGLETON01	SAMPLE_D	/data/project1/SAMPLE_D.bam	/data/project1/SAMPLE_D.bam.bai
```

- `bam`: path to the alignment file, `.bam` or `.cram` — format inferred
  from the extension where a tool needs to branch on it (e.g. Wham has
  separate BAM/CRAM code paths upstream). One column regardless of format,
  matching how GATK-SV itself treats `bam_or_cram_file` as a single
  parameter.
- `bai`: path to the matching index (`.bai` or `.crai`), required
  explicitly — not auto-derived from the alignment path, so indexes staged
  somewhere non-standard don't fail silently.
- `family_id`: redundant with the pedigree file but kept here too, so
  channel construction can group samples by family without a join, and so
  schema validation can catch a sample-sheet/pedigree mismatch early.

### Pedigree file

A separate, standard 6-column PED file per project, adjacent to the sample
sheet:

```text
#family_id  individual_id  paternal_id  maternal_id  sex  phenotype
FAM001      SAMPLE_A       0            0            1    1
FAM001      SAMPLE_B       0            0            2    1
FAM001      SAMPLE_C       SAMPLE_A     SAMPLE_B     1    2
SINGLETON01 SAMPLE_D       0            0            2    1
```

- `sex`: 0=unknown, 1=male, 2=female. Sex chromosome aneuploidies should be
  entered as 0 (GATK-SV convention — downstream ploidy-aware genotyping
  reads this).
- `phenotype`: 0=unknown, 1=unaffected, 2=affected — identifies the proband
  for trio de novo analysis.
- Singletons get a family of one, with `paternal_id`/`maternal_id` = 0.
- One pedigree file per project, matching the panel being per-project.

### Sample and family ID constraints

Carried over from GATK-SV, which imposes these to avoid parsing errors in
shell-based tooling downstream (e.g. `grep`). Validate at sample-sheet/PED
ingestion time rather than letting these surface as opaque failures deep in
the pipeline:

- Alphanumeric and underscores only — no dashes, whitespace, or other
  special characters.
- Unique within the project; not a substring of another ID in the same
  project.
- Should not be purely numeric, and should not contain `chr`, `name`,
  `DEL`, `DUP`, `CPX`, or `CHROM` as a substring.
- Applies to both `sample_id` and `family_id`, in both the sample sheet and
  the pedigree file.

## Where to simplify relative to upstream GATK-SV

Candidates to prototype both ways before committing:

- **Genotyping cutoffs**: GATK-SV's `FilterBatch` random forest + separate
  ML-based GQ recalibration (`TrainGqRecalibrator`) exist to squeeze out
  false positives at biobank scale across heterogeneous batches. For
  smaller, more homogeneous per-project panels, `SVCluster`'s
  support-count/concordance fields plus simple hard filters may be
  sufficient.
- **No first-class "batch" concept.** GATK-SV's WDL threads batch-array
  inputs (order-matched arrays of per-batch files) through most of the
  pipeline to support N-way batch fan-in. We only need panel-version +
  incoming-sample-set; that plumbing goes away entirely.

## Composability

Individual pipeline pieces should be runnable standalone as well as as part
of `build-panel` / `genotype-new-sample`, starting with the SV/CNV callers.
Each caller (and each evidence-collection step) is its own subworkflow with
its own entry point, taking a sample sheet + reference/resource params and
emitting per-sample VCFs or evidence files. `build-panel` and
`genotype-new-sample` are themselves just compositions of these
subworkflows plus the harmonisation steps (`SVCluster`, genotyping,
`SVAnnotate`) — nothing in a caller subworkflow should assume it's being
run as part of either.

This also gives us a natural place to validate each caller against GATK-SV's
own output on the same samples before trusting it in the full pipeline.

### Directory layout

Following [nf-core pipeline convention](https://nf-co.re/docs/contributing/modules)
(as seen in e.g. [sarek](https://github.com/nf-core/sarek)): thin per-tool
wrappers live in `modules/`, pipeline-specific composition lives in
`subworkflows/local/`, and each standalone-runnable entry point lives in
`workflows/`, invoked from a top-level `main.nf` via `-entry`.

```text
nextflow/
├── main.nf                        # entry-point dispatcher (-entry <name>)
├── nextflow.config
├── nextflow_schema.json           # top-level schema; per-workflow schemas alongside each workflow
├── assets/
│   └── schema_samplesheet.json    # nf-schema validation for the TSV sample sheet
├── conf/
│   ├── modules.config             # per-module resource/publishDir/args config
│   └── panels.config              # per-project panel bundle paths/versions
├── modules/
│   ├── nf-core/                   # unmodified nf-core modules (installed via nf-core CLI)
│   │   ├── manta/germline/
│   │   ├── whamg/
│   │   ├── scramble/clusteridentifier/
│   │   ├── scramble/clusteranalysis/
│   │   ├── gatk4/collectreadcounts/
│   │   ├── gatk4/collectsvevidence/
│   │   ├── gatk4/preprocessintervals/
│   │   ├── gatk4/filterintervals/
│   │   ├── gatk4/annotateintervals/
│   │   ├── gatk4/determinegermlinecontigploidy/
│   │   ├── gatk4/germlinecnvcaller/
│   │   ├── gatk4/postprocessgermlinecnvcalls/
│   │   ├── gatk4/svcluster/
│   │   └── gatk4/svannotate/
│   └── local/                     # no nf-core module exists
│       └── cnmops/
├── subworkflows/
│   └── local/
│       ├── bam_call_manta/            # per-caller: BAM/CRAM -> raw caller VCF
│       ├── bam_call_wham/
│       ├── bam_call_scramble/         # needs coverage counts + manta VCF as well as BAM
│       ├── bam_collect_evidence/      # PE/SR/SD (CollectSVEvidence)
│       ├── bam_collect_counts/        # read-depth counts (CollectReadCounts)
│       ├── counts_call_gcnv_cohort/   # gCNV training, build-panel only
│       ├── counts_call_gcnv_case/     # gCNV scoring against a panel model, genotype-new-sample only
│       ├── counts_call_cnmops/
│       ├── vcfs_cluster_svcluster/    # SVCluster, used both within-panel and panel+case
│       ├── vcf_genotype/              # apply cutoffs -> genotyped VCF
│       ├── vcf_annotate_svannotate/
│       └── panel_bundle_io/           # read/write the versioned panel bundle (see Panel bundle contents)
└── workflows/
    ├── call_manta.nf              # standalone entry: sample sheet -> Manta VCFs
    ├── call_wham.nf               # standalone entry: sample sheet -> Wham VCFs
    ├── call_scramble.nf           # standalone entry: sample sheet -> Scramble VCFs
    ├── call_cnmops.nf             # standalone entry: sample sheet -> cn.MOPS calls
    ├── call_gcnv.nf               # standalone entry: sample sheet -> gCNV calls (cohort mode)
    ├── build_panel.nf             # full build-panel pipeline
    └── genotype_new_sample.nf     # full genotype-new-sample pipeline
```

Each `workflows/call_*.nf` is a thin wrapper: read the sample sheet, call the
matching `subworkflows/local/bam_call_*` subworkflow, publish the VCFs. This
is the first slice to implement — it exercises the sample sheet schema, the
module wrappers, and per-caller resource config, without needing the panel
bundle, harmonisation, or genotyping logic yet. `build_panel.nf` and
`genotype_new_sample.nf` come later and compose the same subworkflows.

## nf-core modules

Checked against what GATK-SV's WDL actually calls
([SV/CNV callers](https://broadinstitute.github.io/gatk-sv/docs/gs/sv_callers)):
Manta, Wham, and Scramble are the raw PE/SR-based callers (MELT is legacy/
licensed, superseded by Scramble); GATK-gCNV and cn.MOPS are the depth-based
CNV callers. Canvas and Delly are **not** part of GATK-SV's toolchain — they
were an earlier draft based on general nf-core SV-caller availability rather
than what GATK-SV specifically uses, and have been dropped below in favour
of Scramble and the gCNV chain.

### Callers

- [`manta/germline`](https://nf-co.re/modules/manta_germline/)
- [`whamg`](https://nf-co.re/modules/whamg/)
- [`scramble/clusteridentifier`](https://nf-co.re/modules/scramble_clusteridentifier/) +
  [`scramble/clusteranalysis`](https://nf-co.re/modules/scramble_clusteranalysis/) —
  maps directly onto GATK-SV's own two-step Scramble task
  (`ScramblePart1`/cluster identification → `ScramblePart2`/cluster
  analysis). Need to confirm these wrap the same Scramble version/behaviour
  GATK-SV pins, in particular the DRAGEN-alignment soft-clip realignment
  path (see [Alignment format](#alignment-format)).
- cn.MOPS — **no nf-core module exists** (it's a Bioconductor R package;
  nf-core doesn't wrap it). Needs a thin custom process, same container
  strategy as GATK-SV's own `dockerfiles/cnmops`.

### GATK (evidence collection, depth-based CNV, harmonisation)

nf-core's `gatk4/` namespace already covers nearly the full GATK-side
toolchain GATK-SV relies on, as separate modules:

- [`gatk4/collectreadcounts`](https://nf-co.re/modules/gatk4_collectreadcounts/),
  [`gatk4/collectsvevidence`](https://nf-co.re/modules/gatk4_collectsvevidence/) —
  PE/SR/SD evidence + read-depth counts (GATK-SV's `CollectSVEvidence` /
  `CollectCoverage`).
- [`gatk4/preprocessintervals`](https://nf-co.re/modules/gatk4_preprocessintervals/),
  [`gatk4/filterintervals`](https://nf-co.re/modules/gatk4_filterintervals/),
  [`gatk4/annotateintervals`](https://nf-co.re/modules/gatk4_annotateintervals/) —
  interval prep for gCNV.
- [`gatk4/determinegermlinecontigploidy`](https://nf-co.re/modules/gatk4_determinegermlinecontigploidy/),
  [`gatk4/germlinecnvcaller`](https://nf-co.re/modules/gatk4_germlinecnvcaller/),
  [`gatk4/postprocessgermlinecnvcalls`](https://nf-co.re/modules/gatk4_postprocessgermlinecnvcalls/) —
  the gCNV chain, in both cohort mode (`build-panel`) and case mode
  (`genotype-new-sample`), matching GATK-SV's own `TrainGCNV` /
  `GermlineCNVCohort` vs. `GermlineCNVCase` split.
- [`gatk4/svcluster`](https://nf-co.re/modules/gatk4_svcluster/) — the
  harmonisation/clustering engine (see
  [How GATK-SV's harmonisation actually works](#how-gatk-svs-harmonisation-actually-works)).
- [`gatk4/svannotate`](https://nf-co.re/modules/gatk4_svannotate/) —
  functional/AF annotation.

## Local development / validation

No containers or real BAM/CRAM data are needed to validate pipeline wiring
(channel shapes, module inputs/outputs, sample sheet parsing). Nextflow's
`-stub-run` executes each process's `stub:` block instead of its real
`script:`, which nf-core modules ship by default.

One gotcha: `-stub-run` still executes any `eval(...)`-typed process output
for real (this is how Nextflow captures tool version strings), even though
the rest of the process uses its stub. `manta/germline` declares
`eval("configManta.py --version")` for its `versions_manta` output, so a
stub run fails with `configManta.py: command not found` unless something
by that name is on `PATH`. Options, in order of preference:

- Don't consume the `versions_*` output from local subworkflows until
  running for real (what `bam_call_manta` currently does).
- Put a throwaway shim on `PATH` for local stub runs only, e.g.:

  ```sh
  mkdir -p /tmp/fake-bin
  printf '#!/usr/bin/env bash\necho "fake version"\n' > /tmp/fake-bin/configManta.py
  chmod +x /tmp/fake-bin/configManta.py
  PATH="/tmp/fake-bin:$PATH" nextflow run workflows/call_manta.nf -stub-run ...
  ```

Expect the same issue with other nf-core modules that use `eval()` for
version capture.

Example, using the fixtures in `tests/data/`:

```sh
nextflow run workflows/call_manta.nf \
  --input tests/data/samplesheet.tsv \
  --fasta tests/data/ref.fa \
  --fasta_fai tests/data/ref.fa.fai \
  --outdir /tmp/results \
  -stub-run
```

Real runs (with actual containers and data) happen on HPC via Slurm +
Apptainer, not in this environment — use the `apptainer` profile in
`nextflow.config`.

### A wiring pitfall worth remembering

`MANTA_GERMLINE` (like most nf-core modules) declares its reference genome
inputs (`fasta`, `fasta_fai`) as ordinary channel inputs, not `value`
channels. A per-sample channel combined positionally with a single-emission
reference channel exhausts that reference after the first sample, and
Nextflow silently drops every sample after it — no error, just missing
output. Fix: `.collect()` the reference channels before passing them to the
process, turning them into a broadcastable value reused for every item in
the per-sample channel (`bam_call_manta/main.nf` does this for `fasta`/
`fasta_fai`). Apply the same pattern to every other caller subworkflow.

## HPC (Slurm) config {#hpc-slurm-config}

Real runs happen on HPC via Slurm + Apptainer, not in this environment.
Cluster-specific values (partition/queue names, account string, Apptainer
bind paths) don't belong in the committed `nextflow.config` — that file
stays portable across every environment (local stub runs, this cluster,
anyone else's cluster). They go in a separate, gitignored config layered
on top at runtime with `-c`:

1. Copy the template: `cp conf/hpc.config.example conf/hpc.config`
   (`conf/hpc.config` is gitignored — see the repo's `.gitignore` — so your
   real queue/account values never get committed).
2. Fill in the `TODO`s: partition name(s), Slurm account (`clusterOptions`),
   any Apptainer bind paths your job scripts need, and the reference genome
   paths (`params.fasta`/`params.fasta_fai`).
3. Run with both the portable engine profile and the site-specific config
   layered together:

   ```sh
   nextflow run workflows/call_manta.nf -profile apptainer \
     -c conf/hpc.config --input samplesheet.tsv --outdir /path/to/results
   ```

`-profile apptainer` (in `nextflow.config`) enables the Apptainer engine
itself — that's portable and stays committed. `conf/hpc.config` adds only
what's specific to this cluster on top of it: `process.executor = 'slurm'`,
`process.queue`, `process.clusterOptions`, an `executor {}` block tuning
Slurm job submission/polling rate, and reference genome paths under
`params {}`. `--fasta`/`--fasta_fai` are ordinary pipeline params (declared
in `nextflow.config`, consumed by every caller subworkflow, not specific to
Manta) — since the reference genome is fixed per cluster, setting it once
in `conf/hpc.config` is preferable to retyping the path on every
invocation; a `--fasta` flag on the command line would still override it
for a one-off run against a different reference. See the comments in
`conf/hpc.config.example` for the full shape, modelled on
[nf-core/configs](https://github.com/nf-core/configs)' institutional
profile convention.

If different processes need different partitions (e.g. a high-memory queue
for gCNV once that's implemented), add `withLabel:` overrides inside the
`process {}` block in `conf/hpc.config` rather than hardcoding a queue name
into any module or subworkflow — resource/placement config belongs in
config files, not in pipeline code.

## Status

`bam_call_manta` (wrapping nf-core's `manta/germline`) and its standalone
entry point `workflows/call_manta.nf` are implemented and validated via
`-stub-run` against a mixed BAM/CRAM sample sheet (see
`tests/data/samplesheet.tsv`). This is the reference pattern for the
remaining callers/evidence-collection subworkflows in the
[directory layout](#directory-layout): `MANTA_GERMLINE`'s wrapper is the
template for `BAM_CALL_WHAM`, `BAM_CALL_SCRAMBLE`, etc. — vendor the
nf-core module with `nf-core modules install`, wrap it in a
`subworkflows/local/` composition matching the channel/`.collect()`
patterns above, add a thin `workflows/call_*.nf` entry point, validate with
`-stub-run`.

Not yet started: cn.MOPS (needs a local module, no nf-core equivalent), the
gCNV chain, `SVCluster`/`SVAnnotate` harmonisation, genotyping, panel
bundle I/O, and the full `build_panel.nf` / `genotype_new_sample.nf`
compositions.
