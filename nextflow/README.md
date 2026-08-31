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

Tab-separated (TSV), one row per sample. Column order doesn't matter —
`nf-schema` matches by header name, not position — but we write it as
`family_id, ped_id, sample_id, bam, bai`:

```text
family_id	ped_id	sample_id	bam	bai
FAM001	PED_A	23PS00385-RDN0206	/data/project1/23PS00385-RDN0206.bam	/data/project1/23PS00385-RDN0206.bam.bai
FAM001		SAMPLE_B	/data/project1/SAMPLE_B.bam	/data/project1/SAMPLE_B.bam.bai
SINGLETON01		SAMPLE_C	/data/project1/SAMPLE_C.cram	/data/project1/SAMPLE_C.cram.crai
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
- `ped_id` (optional): the pedigree file's `individual_id` for this sample,
  if it differs from `sample_id`. Defaults to `sample_id` when blank/absent
  (`SAMPLE_B`, `SAMPLE_C` above). GATK-SV itself assumes `sample_id` and
  the PED file's `individual_id` are the same string — our real data
  doesn't always agree: `sample_id` is often a lab/external reference
  (`23PS00385-RDN0206`), while the pedigree file may use a different
  internal family/relationship-position convention (`PED_A`, e.g.
  `RDN0206-00`). Every step that keys off the pedigree file (ploidy
  tables; the VCF sample-column rewrites in
  `modules/local/{manta,wham}/fix_sample_id`, which write `ped_id`, not
  `sample_id`, into the VCF's internal sample name) uses `ped_id`.
  `sample_id` remains the identifier used for output filenames/display.
  This exact mismatch caused a real `KeyError` in GATK-SV's
  `format_svtk_vcf_for_gatk.py` on a real run before `ped_id` existed —
  the ploidy table was keyed by the PED file's `individual_id`, the VCF's
  sample column had been rewritten to `sample_id`, and the two didn't
  match.

### Pedigree file

A separate, standard 6-column PED file per project, adjacent to the sample
sheet:

```text
#family_id  individual_id  paternal_id  maternal_id  sex  phenotype
FAM001      PED_A          0            0            1    1
FAM001      SAMPLE_B       0            0            2    1
SINGLETON01 SAMPLE_C       0            0            2    1
```

- `individual_id`: matches the sample sheet's `ped_id` column when present
  (`PED_A` above, for the same sample whose `sample_id` is
  `23PS00385-RDN0206` in the sample sheet example), otherwise `sample_id`
  itself (`SAMPLE_B`, `SAMPLE_C`).
- `sex`: 0=unknown, 1=male, 2=female. Sex chromosome aneuploidies should be
  entered as 0 (GATK-SV convention — downstream ploidy-aware genotyping
  reads this).
- `phenotype`: 0=unknown, 1=unaffected, 2=affected — identifies the proband
  for trio de novo analysis.
- Singletons get a family of one, with `paternal_id`/`maternal_id` = 0.
- One pedigree file per project, matching the panel being per-project.

### Sample and family ID constraints

Carried over from GATK-SV, which imposes these to avoid parsing errors in
shell-based tooling downstream — some of it does unquoted `grep`/`awk`
matching directly against `sample_id` (e.g.
`wdl/GatherBatchEvidence.wdl:543`, `wdl/SingleSampleFiltering.wdl`).
Validate at sample-sheet/PED ingestion time rather than letting these
surface as opaque failures deep in the pipeline:

- `family_id`: alphanumeric and underscores only — no dashes, whitespace,
  or other special characters. Matches GATK-SV's rule as-is.
- `sample_id`: alphanumeric, underscores, and **hyphens** — our real sample
  IDs contain hyphens (e.g. `21PS00884-RDN0043`), so
  `assets/schema_samplesheet.json` relaxes GATK-SV's rule for this field
  specifically. This reintroduces the class of bug the stricter rule exists
  to prevent, in any place we end up calling GATK-SV-derived shell logic
  verbatim rather than GATK's compiled tools (`SVCluster`/`SVAnnotate`
  take IDs as proper arguments and aren't at risk). Before vendoring any
  GATK-SV shell script or WDL task logic, check it for unquoted sample-ID
  interpolation into `grep`/`awk`/test expressions first.
- Both fields: unique within the project; not a substring of another ID in
  the same project. Should not be purely numeric, and should not contain
  `chr`, `name`, `DEL`, `DUP`, `CPX`, or `CHROM` as a substring.
- Applies in both the sample sheet and the pedigree file.

### Static reference resources

Several `params.*` are static, sample-independent hg38 resource files (see
[Panel bundle contents](#panel-bundle-contents) for the broader
static-vs-per-sample split). They're GATK-SV's own resources — same files
their WDL pipeline uses — published on public GCS buckets, so a plain
`curl`/`wget` over HTTPS works with no GCP credentials needed:

| Param | Source (GATK-SV's `resources_hg38.json` key) | URL |
| --- | --- | --- |
| `primary_contigs_list` | `primary_contigs_list` | `https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/sv-resources/resources/v1/primary_contigs.list` |
| `primary_contigs_fai` | `primary_contigs_fai` | `https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/sv-resources/resources/v1/contig.fai` |
| `pesr_exclude_intervals` (+ `_tbi`) | `pesr_exclude_list` (+ `_index`) | `https://storage.googleapis.com/gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/PESR.encode.peri_all.repeats.delly.hg38.blacklist.sorted.bed.gz` (+ `.tbi`) |
| `manta_region_bed` (+ `_tbi`) | `manta_region_bed` (+ `_index`) | `https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/sv-resources/resources/v1/primary_contigs_plus_mito.bed.gz` (+ `.tbi`) |
| `reference_dict` | `reference_dict` | `https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict` |
| `clustering_config_part1` | `clustering_config_part1` | `https://storage.googleapis.com/gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/clustering_config.part_one.tsv` |
| `clustering_config_part2` | `clustering_config_part2` | `https://storage.googleapis.com/gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/clustering_config.part_two.tsv` |
| `stratification_config_part1` | `stratification_config_part1` | `https://storage.googleapis.com/gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/stratify_config.part_one.tsv` |
| `stratification_config_part2` | `stratification_config_part2` | `https://storage.googleapis.com/gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/stratify_config.part_two.tsv` |
| `clustering_track_sr` | `clustering_tracks[0]` | `https://storage.googleapis.com/gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/hg38.SimpRep.sorted.pad_100.merged.bed` |
| `clustering_track_sd` | `clustering_tracks[1]` | `https://storage.googleapis.com/gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/hg38.SegDup.sorted.merged.bed` |
| `clustering_track_rm` | `clustering_tracks[2]` | `https://storage.googleapis.com/gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/hg38.RM.sorted.merged.bed` |

`primary_contigs_list` (plain list) and `primary_contigs_fai` (`.fai` —
contig *and length*) are both needed, for different tools: Wham's `-c`
flag and `vcfs_cluster_svcluster`'s ploidy-table/format-conversion scripts
tolerate either format, but `svtk standardize --contigs` (used by
`SVTK_STANDARDIZE`, see [Status](#status)) needs lengths to build its
output VCF header, so it specifically needs the `.fai` form.
`min_svsize` (default `50`, GATK-SV's own default) has no file to
fetch — it's a plain integer param, also consumed by `svtk standardize`.
The `clustering_config_part*`/`stratification_config_part*`/
`clustering_track_*` params are stage-2-only (`GroupedSVCluster`'s
context-aware re-clustering — see [Harmonisation stage
2](#harmonisation-stage-2-cross-caller-merge)); everything else is used
starting from stage 1.

#### `scripts/fetch_static_resources.sh`

Fetches every resource in the table above into a target directory and
prints a ready-to-use Nextflow config snippet (a `params { ... }` block)
on stdout, pointing at the downloaded paths:

```sh
scripts/fetch_static_resources.sh /path/to/static-resources > conf/hpc.config
```

Idempotent — re-running skips any file already present in the target
directory by name, so it's safe to re-run after an interrupted fetch, or
to pick up newly-added resources (e.g. after this pipeline grows a new
stage that needs one) without re-downloading everything. Doesn't fetch
the reference genome itself (`fasta`/`fasta_fai`) — those aren't
GATK-SV-specific and you likely already have a cluster-wide copy; set
those yourself in the printed config (or merge the script's output into
an existing `conf/hpc.config` rather than overwriting one that already
has them — the script only emits the `params {}` block, it doesn't read
or preserve anything already in the target file).

Use GATK-SV's own files rather than regenerating equivalents, to match
their validated behavior exactly (e.g. `manta_region_bed` is specifically
primary contigs *plus mito*, not just `primary_contigs_list` reformatted
as a BED).

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
│       ├── utils_input_channels/      # shared: sample sheet -> [meta,bam,bai] + reference channels
│       ├── bam_call_manta/            # per-caller: BAM/CRAM -> raw caller VCF
│       ├── bam_call_wham/
│       ├── bam_call_scramble/         # needs coverage counts + manta VCF as well as BAM
│       ├── bam_collect_evidence/      # PE/SR/SD (CollectSVEvidence)
│       ├── bam_collect_counts/        # read-depth counts (CollectReadCounts)
│       ├── counts_call_gcnv_cohort/   # gCNV training, build-panel only
│       ├── counts_call_gcnv_case/     # gCNV scoring against a panel model, genotype-new-sample only
│       ├── counts_call_cnmops/
│       ├── vcfs_cluster_svcluster/    # SVCluster, used both within-panel and panel+case (harmonisation stage 1)
│       ├── vcfs_combine_batches/      # cross-caller merge (harmonisation stage 2)
│       ├── vcf_genotype/              # apply cutoffs -> genotyped VCF
│       ├── vcf_annotate_svannotate/
│       └── panel_bundle_io/           # read/write the versioned panel bundle (see Panel bundle contents)
└── workflows/
    ├── call_manta.nf              # standalone entry: sample sheet -> Manta VCFs
    ├── call_wham.nf               # standalone entry: sample sheet -> Wham VCFs
    ├── call_scramble.nf           # standalone entry: sample sheet -> Scramble VCFs
    ├── call_cnmops.nf             # standalone entry: sample sheet -> cn.MOPS calls
    ├── call_gcnv.nf               # standalone entry: sample sheet -> gCNV calls (cohort mode)
    ├── call_all.nf                # standalone entry: sample sheet -> VCFs from every implemented caller
    ├── cluster_manta.nf           # standalone entry: sample sheet -> clustered Manta site VCF (harmonisation stage 1)
    ├── cluster_wham.nf            # standalone entry: sample sheet -> clustered Wham site VCF (harmonisation stage 1)
    ├── combine_batches.nf         # standalone entry: sample sheet -> cross-caller merged site VCF (harmonisation stage 2)
    ├── build_panel.nf             # full build-panel pipeline
    └── genotype_new_sample.nf     # full genotype-new-sample pipeline
```

Each `workflows/call_*.nf` is a thin wrapper: read the sample sheet, call the
matching `subworkflows/local/bam_call_*` subworkflow, publish the VCFs. This
is the first slice to implement — it exercises the sample sheet schema, the
module wrappers, and per-caller resource config, without needing the panel
bundle, harmonisation, or genotyping logic yet.

`call_all.nf` sits between the single-caller entries and `build_panel.nf` /
`genotype_new_sample.nf`: it parses the sample sheet once (via
`utils_input_channels`) and runs every implemented caller subworkflow
against the same samples, so you get output from all callers without
invoking each `call_*.nf` separately — but it still has no panel bundle,
harmonisation, or genotyping. Each caller subworkflow it calls is still
independently runnable through its own `call_*.nf`; nothing in `call_all.nf`
is caller-specific logic, it's just composition. Add a new `include` +
call block there as each caller subworkflow lands. `build_panel.nf` and
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

**`-stub-run` validates wiring, not everything.** Three real bugs shipped
past `-stub-run` and were only caught by an HPC run (or a real local
`-profile docker` run reproducing it):

- A vendored script called by bare name that was never actually on `PATH`
  (`script:` only, never exercised by any `stub:` block).
- A `docker://`-prefixed container reference that broke `-profile docker`,
  which `-stub-run` never noticed because it never invokes the container
  engine at all regardless of profile.
- A genuine data-correctness bug that needed *real* data flowing through
  *real* tool logic to surface at all: Manta writes the BAM's `SM` tag as
  the VCF's sample column, not `sample_id`, and nothing failed until a
  downstream step (`SVCluster`'s ploidy-table lookup) actually needed the
  two to match. No amount of stub testing — or even a correct, running
  container — would have caught this; it only showed up once a real
  sample's VCF met a real ploidy table with a real, non-matching sample
  name. See [Status](#status) for the fix.

See the [wiring pitfalls](#a-wiring-pitfall-worth-remembering) below for
the first two. For any process that shells out to a vendored script or
declares its own `container` directive, do at least one real (non-stub)
run — a local `-profile docker` run is usually enough, doesn't need HPC
access — before considering it validated. But even that isn't sufficient
proof of *correctness*: when wrapping a bare tool GATK-SV's own WDL task
does more with (renames, reheaders, converts), check the WDL task's full
command block, not just whether the tool runs — a process that runs
successfully and produces a plausible-looking VCF can still be silently
wrong in a way that only breaks two or three steps later.

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
`.map { meta, f -> f }.collect()` (strip meta, then collect) if the target
process wants a bare path rather than a `[meta, path]` tuple — collecting
first wraps the whole tuple in a one-element list instead and breaks the
downstream `.map` (`bam_call_wham/main.nf` needed this variant, since
`WHAMG` takes bare `path(fasta)`/`path(fasta_fai)`, unlike Manta).

Second: `file()` behaves differently inside a `conf/modules.config`
`ext.args` closure than in ordinary workflow script code. There, `file`
resolves to `nextflow.script.ScriptBinding.file()` (a zero-arg accessor)
instead of the top-level `file(path)` factory function, so a bare
`file(params.some_path)` call inside `ext.args = { ... }` fails with
`No signature of method: ...file() is applicable for argument types: ()
values: []` — a genuine namespace collision specific to config-closure
scope. Fix: call the fully-qualified `nextflow.Nextflow.file(...)` instead
(see `conf/modules.config`'s `WHAMG` block, which reads
`params.primary_contigs_list`'s contents this way to build Wham's `-c`
argument). A `lib/*.groovy` helper class was tried first, on the
expectation that config closures could call into it, but classes there
weren't visible from `conf/`-included config closures in this Nextflow
version (`No such variable`) — inlining in the config file directly was
what actually worked.

Third: a `container` directive value is engine-specific syntax, not one
universal reference. `docker://<image>` is Apptainer/Singularity syntax
for "pull from Docker Hub"; plain `docker` does not understand that
scheme and fails with `docker: invalid reference format`. Our own
vendored-script modules (`ploidy_table_from_ped`,
`format_svtk_vcf_for_gatk`, `format_gatk_vcf_for_svtk`) originally
hardcoded `container 'docker://drtomc/gatk-sv-nf-sv-scripts:0.1.0'`,
which worked under `-profile apptainer` but broke `-profile docker`
outright — undetected until a real (non-stub) local run, since
`-stub-run` never invokes the container engine at all. Fixed by using
the bare image reference (`drtomc/gatk-sv-nf-sv-scripts:0.1.0`, no
scheme), which both engines accept. For a container that genuinely needs
different references per engine (e.g. a Docker Hub image plus a
separately-hosted Singularity build, as most nf-core modules do), use
the `workflow.containerEngine in ['singularity', 'apptainer'] ? ... : ...`
conditional instead — see any nf-core-installed module, or
`modules/local/genome_file/main.nf`, for the pattern.

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

Implemented and validated via `-stub-run` against a mixed BAM/CRAM sample
sheet (see `tests/data/samplesheet.tsv`):

- `bam_call_manta` (wraps nf-core's `manta/germline`) and its standalone
  entry point `workflows/call_manta.nf`. Confirmed working against real
  BAM/CRAM data on HPC (Slurm + Apptainer). GATK-SV restricts Manta's own
  calling to primary contigs + mito via `--callRegions`
  (`wdl/Manta.wdl:137`, the `manta_region_bed` resource, wired into
  `MANTA_GERMLINE`'s `target_bed`). Raw output is then standardized with
  `SVTK_STANDARDIZE` — see below.
  - **Not implemented**: GATK-SV also runs Manta's raw output through
    `convertInversion.py` (converts Manta's inversion-signaling `BND` pairs
    into proper `INV` records) before standardization. Skipped for now;
    revisit once inversions specifically need validating.
- `bam_call_wham` (wraps nf-core's `whamg`) and its standalone entry point
  `workflows/call_wham.nf`. Contig restriction via `-c primary_contigs_list`
  is wired through `ext.args` in `conf/modules.config` (see the pitfalls
  above for why this needed `nextflow.Nextflow.file()` rather than a bare
  `file()` call). Raw output is then standardized with `SVTK_STANDARDIZE`,
  same as Manta.
  - **Deliberately not implemented**: GATK-SV also scatters Wham calls
    over an `include_bed_file` region whitelist and concatenates, to bound
    runtime and avoid regions Wham struggles with. We currently run
    `whamg` genome-wide in one shot instead. Revisit if runtime on real
    data warrants it.
- `SVTK_STANDARDIZE` (`modules/local/svtk_standardize`), used by both
  `bam_call_manta` and `bam_call_wham` — wraps GATK-SV's own `svtk
  standardize` (`src/svtk/`, their `StandardizeVCFs` task,
  `wdl/PESRPreprocessing.wdl`). In one step: rewrites the VCF's sample
  column to `ped_id`, restricts to primary contigs, applies minimum-size
  filtering, and sets the INFO fields GATK's `SVCluster` requires
  (`SVTYPE`/`CHR2`/`END`/`STRANDS`/`SVLEN`/`ALGORITHMS`).
  - **Supersedes three earlier point-fixes**
    (`modules/local/{manta,wham}/fix_sample_id`'s `bcftools reheader`
    rename, and `modules/local/vcf_primary_contigs_only`'s post-hoc
    `bcftools` filter, all removed) that were each independently
    reconstructing part of what `svtk standardize` already does correctly.
    Discovered incomplete via three separate real failures on real data:
    a `KeyError` on sample name (the rename target was wrong —
    `sample_id` instead of `ped_id`), a `KeyError` on an ALT contig
    (`chr1_KI270706v1_random` — Manta's own `--callRegions` restriction
    didn't fully prevent it), and finally a GATK `IllegalArgumentException:
    Expected ALGORITHMS field` — none of the three point-fixes set
    `ALGORITHMS` at all, since that wasn't the bug any of them were
    written to fix.
  - **Not vendored like `bin/*.py`**: `svtk` is a real Python package
    (Cython extension, `pybedtools`/`bedtools` dependency) that GATK-SV
    itself builds from a multi-stage Dockerfile
    (`dockerfiles/sv-pipeline/Dockerfile`) on top of their own base
    images — not something to replicate from scratch. Runs in GATK-SV's
    own published `sv_pipeline_docker` image
    (`us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline`) instead, which
    is publicly pullable (verified — no GCP credentials needed).
  - **Wham needs its raw, un-renamed output**, specifically: `svtk`'s own
    Wham standardizer (`src/svtk/svtk/standardize/std_wham.py`) determines
    genotypes by checking whether each *raw* sample name (from the input
    VCF's own sample column) appears in the `TAGS` INFO field, then maps
    to `--sample-names` (`ped_id`) as it writes output — it does this
    matching itself. Feeding it already-renamed input (what the old
    `WHAM_FIX_SAMPLE_ID` did) would only work by coincidence; verified
    directly that feeding it Wham's true raw output produces the correct
    genotype and sample rename.
  - Needs a second static resource beyond `primary_contigs_list`:
    `primary_contigs_fai` (`.fai` format — contig *and length*, not just
    names — GATK-SV's `primary_contigs_fai` resource), since `svtk
    standardize --contigs` needs lengths to build the standardized VCF
    header. `primary_contigs_list` (plain list) is still needed separately
    for Wham's `-c` flag and `vcfs_cluster_svcluster`'s ploidy-table/
    format-conversion scripts, which happen to tolerate either format.
  - Verified directly (not just wired/stub-tested): ran the actual `svtk
    standardize` command against hand-built VCFs reproducing both real
    failures (an ALT-contig Manta-shaped record, a `TAGS`-based
    Wham-shaped record) in the real `sv_pipeline_docker` container before
    committing this change.
- `workflows/call_all.nf`, running every implemented caller (Manta, Wham)
  against the same sample sheet in one invocation — see
  [directory layout](#directory-layout) for what this is and isn't.
  Sample-sheet parsing and reference-channel construction, previously
  duplicated in `call_manta.nf` and `call_wham.nf`, are now factored into
  `subworkflows/local/utils_input_channels`, shared by all three entry
  points.

This is the reference pattern for the remaining callers/evidence-collection
subworkflows in the [directory layout](#directory-layout): vendor the
nf-core module with `nf-core modules install`, wrap it in a
`subworkflows/local/` composition (adding local pre/post-processing steps
where GATK-SV's own task does more than the bare tool invocation — check
the corresponding `wdl/*.wdl` file, don't assume the nf-core module is a
drop-in match), add a thin `workflows/call_*.nf` entry point, validate with
`-stub-run`.

### Harmonisation stage 1: per-caller, cross-sample clustering

`subworkflows/local/vcfs_cluster_svcluster` and its two standalone entry
points, `workflows/cluster_manta.nf` and `workflows/cluster_wham.nf`,
implement stage 1 of harmonisation (see
[How GATK-SV's harmonisation actually works](#how-gatk-svs-harmonisation-actually-works)):
clustering one caller's SV calls across all samples in a run into
representative sites, mirroring GATK-SV's `ClusterPESR` workflow
(`wdl/PESRClustering.wdl`). This is deliberately *not* stage 2
(cross-caller merging, `CombineBatches`/`MergeBatchSites`) — that's not
implemented yet. `vcfs_cluster_svcluster` is caller-agnostic (takes
`caller` as a plain string param, e.g. `'manta'`/`'wham'`); adding
`cluster_wham.nf` alongside `cluster_manta.nf` needed no changes to the
subworkflow itself, only a near-identical entry point swapping
`BAM_CALL_MANTA` for `BAM_CALL_WHAM` (and dropping `manta_region_bed`,
Manta-only) — confirms the subworkflow's genericity holds for a second
caller, not just designed-for-one-and-hoped.

Both entry points confirmed working end-to-end against real BAM/CRAM data
on HPC (Slurm + Apptainer), not just `-stub-run`.

Steps, matching `ClusterPESR`:

1. Ploidy table from the pedigree file (`PLOIDY_TABLE_FROM_PED`, wraps
   GATK-SV's `ploidy_table_from_ped.py`).
2. Per-sample format conversion (`FORMAT_SVTK_VCF_FOR_GATK`, wraps
   GATK-SV's `format_svtk_vcf_for_gatk.py`) + interval-exclusion filtering
   (`VCF_ENDS_BED` → `bedtools/intersect` → `EXCLUDED_VARIANT_IDS` →
   `EXCLUDE_VARIANTS_BY_ID`, the same `bcftools query | awk | sort |
   bedtools intersect | cut | sort | uniq | bcftools view` shape GATK-SV
   uses inline in both `PreparePESRVcfs` and `ExcludeIntervalsByEndpoints`
   — one local module chain, reused for both).
3. `SVCluster` (`gatk4/svcluster`) across all samples' prepared VCFs in
   one call.
4. Post-clustering interval exclusion (same module chain as step 2,
   applied a second time — mirrors `ExcludeIntervalsByEndpoints`).
5. Format back to svtk (`FORMAT_GATK_VCF_FOR_SVTK`, wraps GATK-SV's
   `format_gatk_vcf_for_svtk.py`).

Deliberately not implemented: GATK-SV scatters step 3 by contig then
concatenates, for parallelism at cohort scale. We run one `SVCluster` call
across the whole genome instead — revisit if runtime on real data warrants
it (same simplification already made for Wham's region scatter). GATK-SV's
`min_size` filtering (drop SVs below a minimum length) is threaded through
as a param but not yet applied — only interval exclusion is, so far.

#### Vendored GATK-SV scripts

Three of GATK-SV's own Python scripts (`ploidy_table_from_ped.py`,
`format_svtk_vcf_for_gatk.py`, `format_gatk_vcf_for_svtk.py`) are copied
verbatim into `bin/`. See `bin/README.md` for exactly which upstream
commit they were vendored from, and for two known gaps:

**These scripts are passed to their processes as explicit `path` inputs**
(e.g. `PLOIDY_TABLE_FROM_PED(ped, contigs, file(ploidy_script))`, invoked
as `python3 ${script}`), not called by bare name relying on Nextflow's
`bin/`-auto-PATH mechanism. That mechanism ties `PATH` to
`scriptFile.main.parent.resolve('bin')` — the directory of the *script
passed to `nextflow run`* — unconditionally, not the launch/current
directory. Since every entry point lives in `workflows/`, not the
repository root where `bin/` actually is, `bin/` auto-PATH silently never
applies to any `workflows/*.nf` script. This produced a real failure on an
HPC run (`exit 127: ploidy_table_from_ped.py: command not found`) that
`-stub-run` never caught, because a process's `stub:` block never calls
the vendored script by name — only `script:` does, and stub validation
never executes `script:`. The fix threads each script's path down from
the entry point (where `pipelineRoot` is already computed correctly for
`assets/schema_samplesheet.json`, for the same reason) through the
subworkflow to the process. **Lesson for future modules that shell out to
a vendored script**: don't rely on bin/-auto-PATH from a `workflows/*.nf`
entry point; pass the script as an explicit input instead, and validate
with a real (non-stub, e.g. `-profile docker`) run at least once, since
stub runs won't exercise the actual invocation.

- **pysam version**: the two format-conversion scripts import `pysam`.
  GATK-SV's own container pins `pysam==0.15.4` specifically (built from
  source) because newer pysam/htslib reject VCF records with `END < POS`,
  which occurs for `BND`/`CTX` (breakend/translocation) records. Our
  container uses current pysam instead, for build simplicity.
  `format_svtk_vcf_for_gatk.py`'s `_parse_bnd_ends()` already works around
  part of this by manually text-parsing `BND`/`CTX` records rather than
  trusting pysam's own parsed `END` value — but whether current pysam's
  VCF *iteration itself* fails outright on such records (independent of
  what the script does with them afterward) is untested. **Test this
  empirically against real Manta/Wham output** (which does produce `BND`
  records) before trusting this path's `BND` handling.
- These scripts, plus `tabix`, run in our own image,
  `drtomc/gatk-sv-nf-sv-scripts:0.1.0`, built from
  `dockerfiles/sv-scripts/Dockerfile` and pushed to Docker Hub — see that
  Dockerfile for build/push instructions if it needs updating.

#### A local module, not a repurposed nf-core one

`EXCLUDE_VARIANTS_BY_ID` (`bcftools view -i 'ID!=@file' | tabix`) is a
small local module rather than nf-core's `bcftools/view`. The `ID!=@file`
bcftools expression needs the exclusion file's staged path known at
script-render time; nf-core's `bcftools/view` only exposes
`--regions-file`/`--targets-file`/`--samples-file`, none of which map onto
an arbitrary `-i`/`-e` expression referencing an external file. Generic
nf-core modules cover the *tool*, not every *invocation shape* GATK-SV
uses that tool with — check the actual flags a module's `script:` block
constructs before assuming it fits, not just its name.

#### A patched nf-core module

`modules/nf-core/gatk4/svcluster/main.nf`'s `stub:` block is patched
locally (diverges from `github.com/nf-core/modules`): upstream's stub only
creates `*.vcf.gz`, but the module's own `output:` block also declares
`clustered_vcf_index` (`*.vcf.gz.tbi`), so `-stub-run` failed with
"Missing output file(s)" until fixed. Only the `stub:` block is
touched — `script:` (the real invocation) is unmodified upstream code.

#### Writing files from a subworkflow — do it in a process, not inline Groovy

An earlier version of `vcfs_cluster_svcluster` built the bedtools genome
file (`cut -f1,2 fasta.fai`) with inline Groovy
(`fai.resolveSibling("genome.file"); out.text = ...`) directly in the
subworkflow body. This silently wrote a `genome.file` next to whatever
`fasta_fai` pointed at — for the test fixtures, `tests/data/`, polluting
the source tree; on a real run, potentially a read-only or shared
reference directory. `-stub-run` never caught it because the closure ran
regardless of stub/real mode — inline Groovy in a subworkflow body isn't
process-scoped, so it runs at graph-construction time either way, unlike a
process's `script:`/`stub:` blocks. Fixed by moving it into a proper
`GENOME_FILE` process (`modules/local/genome_file`), whose output lands in
the task's own work directory like everything else. Prefer a process over
inline Groovy file-writing whenever output is more than an in-memory value
passed to the next channel step.

#### `whamg` needs much more memory than its nf-core label gives it

A real HPC run of `cluster_wham.nf` (`~50x`-coverage WGS BAM) had `whamg`
OOM-killed by Slurm (`.command.log`: `slurmstepd: error: Detected 1
oom_kill event`) partway through — its own log showed it had reached
"processed 1800Mb of the genome" out of ~3.1Gb, after first "Loading
discordant reads into forest" (an in-memory index built across the whole
genome, not something proportional to a chunk of input the way most tasks
are). The `whamg` nf-core module uses `label 'process_medium'`, which
`nextflow.config` sets to 12GB — shared with several other, much lighter
processes, and not enough here. GATK-SV's own WDL sizing formula for this
task (`wdl/Whamg.wdl`, `mem_per_bam_size * bam_size + mem_bam_offset`)
defaults to a flat ~5GB (`mem_bam_offset=4.9`) when its Picard-metrics
inputs (`pf_reads_improper_pairs`/`pct_exc_total`) aren't supplied — which
we don't currently wire up — so it's both lower than what already failed
and not a formula worth chasing for this. Fixed with a flat, generous
per-process override (`memory = 32.GB`, plus explicit `cpus`/`time`) in
`conf/modules.config`'s `withName: 'WHAMG'` block, rather than
reverse-engineering GATK-SV's metrics-dependent formula. Revisit with real
per-sample memory usage data once more samples have run, if 32GB turns out
to be too high (waste) or still too low (rare, larger-input samples).

#### The OOM above should have failed the task outright, and didn't — a `pipefail` gap

Nextflow reported the OOM-killed `whamg` run above as a **successful**
task. The actual failure only surfaced several steps later, in an
unrelated downstream process, as a confusing `pysam` error ("invalid
file ... is it VCF/BCF format?") on what turned out to be a truncated
28-byte gzip stream. Root cause: `whamg`'s own `script:` block
(`modules/nf-core/whamg/main.nf`) pipes its output through `bgzip`:

```sh
whamg ... | bgzip --threads ... --stdout > out.vcf.gz
```

Nextflow's default task shell is `['/bin/bash', '-ue']` — no `pipefail`.
Bash's default pipeline exit status is the *last* command's only, so when
Slurm SIGKILLed `whamg` mid-stream, `bgzip` just saw its stdin close early
and exited 0 on the truncated data — bash then reported the whole pipeline
(and therefore the whole task) as exit 0, hiding the OOM completely from
Nextflow's success/failure tracking. **This is not a Slurm/HPC-profile
weakness** — the same silent failure would happen with any executor,
local runs included; it's purely a shell-strictness gap in how the task
script is invoked, independent of *why* `whamg` died (OOM here, but a
segfault or any other signal would be swallowed identically). Any other
module using a `cmd1 | cmd2 > out` pattern has the same latent risk.

Fixed globally, not per-module: `nextflow.config`'s top-level `process {}`
block now sets `shell = ['/bin/bash', '-ueo', 'pipefail']`, so a killed
first stage in *any* pipe now fails the task immediately with a non-zero
exit code, the way Nextflow's own error handling expects. Chosen over
patching individual modules' `script:` blocks (which would mean
diverging every affected nf-core module from upstream, the exact
one-off-patch cost the [patched `gatk4/svcluster`
stub](#a-patched-nf-core-module) already documents) — a shell default
is the natural place to fix a shell-semantics gap. Re-validated all three
`-stub-run` entry points (`call_manta`, `cluster_wham`,
`combine_batches`) afterward to confirm `-ueo pipefail` doesn't break any
existing `script:`/`stub:` block (none rely on tolerating an unset
variable or a failing pipeline stage).

### Harmonisation stage 2: cross-caller merge

`subworkflows/local/vcfs_combine_batches` and its standalone entry point
`workflows/combine_batches.nf` implement stage 2 of harmonisation:
merging multiple callers' already-clustered site VCFs (Manta's and
Wham's, from stage 1) into one cross-caller site list. Mirrors GATK-SV's
`CombineBatches` workflow (`wdl/CombineBatches.wdl`). `combine_batches.nf`
composes the full pipeline in one entry point — sample sheet in, one
merged VCF out — running `BAM_CALL_MANTA`/`BAM_CALL_WHAM` and both
callers' `VCFS_CLUSTER_SVCLUSTER` calls itself, same pattern as
`call_all.nf`.

Steps, matching `CombineBatches`:

1. Naive join across callers (`GATK4_SVCLUSTER_JOIN`, near-zero overlap
   thresholds — only merges exact-position duplicates).
2. Real clustering pass (`GATK4_SVCLUSTER_SITES`, realistic overlap
   thresholds, plus a `--variant-prefix` keyed by `cohort_name`).
3. Two rounds of context-aware re-clustering
   (`GATK4_GROUPEDSVCLUSTER_PART1`/`PART2`, wraps GATK's `GroupedSVCluster`
   walker — see below), stratified by genomic context tracks (simple
   repeats, segmental duplications, RepeatMasker).
4. Format back to svtk (`FORMAT_GATK_VCF_FOR_SVTK`, reused from stage 1).

This subworkflow builds its *own* ploidy table
(`PLOIDY_TABLE_FROM_PED_COMBINE_BATCHES`, aliased — see below), with
`retain_female_chr_y=true`, separate from stage 1's own ploidy table
(`retain_female_chr_y` left at its default `false`). That flag isn't
something the vendored `ploidy_table_from_ped.py` script itself
supports — it's WDL-task-level post-processing
(`wdl/TasksClusterBatch.wdl`'s `CreatePloidyTableFromPed`:
`sed -e 's/\t0/\t1/g'`, rewriting every `0` ploidy value to `1`, which for
females only affects chrY) that our `PLOIDY_TABLE_FROM_PED` module now
supports via `task.ext.retain_female_chr_y`.

Deliberately not implemented: GATK-SV's cross-*batch* SR evidence flag
reconciliation (`ExtractSRVariantLists`/`CombineSRBothsidePass`/
`SetSRVariantFlags`) is skipped — it exists to combine
`BOTHSIDES_SUPPORT`/`HIGH_SR_BACKGROUND` flags across multiple *batches*
of the same caller category (GATK-SV's own cohort-mode "batch" concept),
which has no analog here: we have one VCF per *caller*, not per *batch*,
and this pipeline doesn't have a first-class batch concept at all (see
[Where to simplify](#where-to-simplify-relative-to-upstream-gatk-sv)).
Also not scattered by contig — same simplification as stage 1.

#### A new GATK subcommand, no nf-core module

`modules/local/gatk4/grouped_sv_cluster` wraps GATK's `GroupedSVCluster`
walker (marked **BETA — WORK IN PROGRESS** by GATK itself) — no nf-core
module exists for it. Runs from the same container as `gatk4/svcluster`
(`community.wave.seqera.io/library/gatk4-main_gcnvkernel`) — GATK-SV's own
WDL runs both `SVCluster` and `GroupedSVCluster` from the same
`gatk_docker` image with no separate image param, confirming they're
different walkers in the same GATK jar, not separate tools. Verified
directly: pulled the image and ran `gatk GroupedSVCluster --help` to
confirm the tool and every flag name (`--track-intervals`, `--track-name`,
`--stratify-overlap-fraction`, etc.) match what GATK-SV's WDL task uses,
before writing the module.

New static resources, same public-bucket pattern as everything else:
`clustering_config_part1`/`_part2`, `stratification_config_part1`/`_part2`
(TSV clustering/stratification configs), and three context-track BEDs
(`clustering_track_sr`/`_sd`/`_rm` — simple repeats, segmental
duplications, RepeatMasker; exposed as three separate params rather than
one ordered list, so nothing can silently get the track/name pairing out
of order).

#### A wiring bug this surfaced: process-name-derived output filenames can collide across pipeline stages

`GATK4_SVCLUSTER`'s (and `GATK4_GROUPEDSVCLUSTER`'s) output filename is
derived from `meta.id` (`prefix = task.ext.prefix ?: "${meta.id}"`). An
earlier version of this subworkflow reused the same `meta` (built once
from `cohort_name`) across all five stages — so stage *N*'s output and
stage *N+1*'s staged input both wanted the filename `cohort.vcf.gz`. Under
`-stub-run` this produced a confusing "Missing output file(s) `*.vcf.gz`"
error (the stub's `echo ... > cohort.vcf.gz` was overwriting its own
staged input symlink in place, not producing a distinct new file the glob
matcher could find) — and it would have broken real execution the same
way, not just the stub. Fixed by giving each stage's `meta.id` a distinct,
stage-specific suffix (matching GATK-SV's own `output_prefix` convention,
e.g. `"~{cohort_name}.combine_batches.~{contig}.join_vcfs"`), renaming
back to a plain `cohort_name` only for the final published output.
**Lesson for any future multi-stage chain reusing the same process
multiple times**: check whether the process derives its output filename
from `meta.id`, and if so, make sure `meta.id` actually changes between
producer and consumer stages — reusing one static `meta` across an entire
chain is a trap specifically because it looks correct (every individual
process call is well-formed) right up until two stages collide on a
filename.

### Not yet started

Scramble (needs coverage counts + Manta's VCF as additional inputs, not
just the BAM — more involved than Wham), cn.MOPS (needs a local module, no
nf-core equivalent), the gCNV chain, `MergeBatchSites` (GATK-SV's
own site-list-only shortcut for cohorts that skip GATK-gCNV — not
implemented; `CombineBatches`, above, is), `SVAnnotate`, genotyping, panel
bundle I/O, and the full `build_panel.nf` / `genotype_new_sample.nf`
compositions.
