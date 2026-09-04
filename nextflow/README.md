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
| `sd_locs_vcf` (+ `_idx`) | `sd_locs_vcf` | `https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf` (+ `.idx`) |
| `preprocessed_intervals` | `preprocessed_intervals` | `https://storage.googleapis.com/gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/preprocessed_intervals.interval_list` |

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
2](#harmonisation-stage-2-cross-caller-merge)); `sd_locs_vcf`/
`preprocessed_intervals` are used by `bam_collect_evidence` (see
[Status](#status)); everything else is used starting from stage 1.
`sd_locs_vcf` is the one exception to "everything here is bgzip/`.tbi`" —
it's a plain, uncompressed VCF with a GATK-style `.idx` sibling index, not
`.vcf.gz`/`.tbi` like every other VCF this pipeline handles.

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

#### A real bug: the ploidy table's own output glob matched its intermediate file too

A real HPC run of `combine_batches.nf` failed with `SVCluster`'s `A USER
ERROR has occurred: Illegal argument value: Positional arguments were
provided ',tmp.ploidy.tsv}' but no positional argument is defined for this
tool` — the actual `--ploidy-table` invocation showed two filenames
(`ped.ploidy.tsv tmp.ploidy.tsv`) where `SVCluster` expects exactly one.
Root cause: `PLOIDY_TABLE_FROM_PED`'s `output:` block declares `path("*.
ploidy.tsv")`, and — only when `retain_female_chr_y=true`, i.e. only at
this subworkflow's `PLOIDY_TABLE_FROM_PED_COMBINE_BATCHES` call site — the
`script:` block's own intermediate file was itself named `tmp.ploidy.tsv`,
which matches that same glob alongside the real, sed-processed
`${prefix}.ploidy.tsv`. So the `ploidy_table` output channel silently
emitted *two* files instead of one; `GATK4_SVCLUSTER`'s `path
ploidy_table` (a single-value input) staged both, and Groovy's
space-joined interpolation of the two into `--ploidy-table x y` is what
`SVCluster` then choked on. Invisible under `-stub-run`, whose `stub:`
block only ever `touch`es one file — the glob never had two real
candidates to collide on. Fixed by renaming the intermediate file to
`raw.ploidy_table.tsv` (doesn't match `*.ploidy.tsv`), then verified
directly with a small standalone real (non-stub, `-profile docker`)
workflow isolating just this process with `retain_female_chr_y=true`
against a female sample in `tests/data/cohort.ped` — confirmed exactly one
file is now emitted, and that the female chrY ploidy rewrite (`0` → `1`)
still applies correctly. **Lesson for any process whose `output:` block
uses a glob**: check every intermediate filename the `script:` block
creates, not just the declared final output name, for accidental matches —
especially under a config-driven branch (`task.ext.*`) that only some call
sites exercise, since a glob collision that's specific to one branch is
easy to miss when validating with the default/most-common configuration.

#### A real bug: `ECN` missing because our stage 1 round-trips through svtk format, GATK-SV's doesn't

A real HPC run of `combine_batches.nf` got past the ploidy-table fix above
and failed one step later, in `GATK4_SVCLUSTER_SITES`, with `SVCluster`'s
`java.lang.IllegalArgumentException: Genotype missing required field ECN`.
`ECN` (Expected Copy Number) is a per-genotype `FORMAT` field GATK's
clustering walkers require to compute sample overlap — added by
`FORMAT_SVTK_VCF_FOR_GATK`/`format_svtk_vcf_for_gatk.py` from the ploidy
table, and explicitly *stripped* by `FORMAT_GATK_VCF_FOR_SVTK`/
`format_gatk_vcf_for_svtk.py` (`remove_formats.add('ECN')`) when
converting back to svtk format for publishing.

Root cause: `vcfs_cluster_svcluster` (stage 1) converts to GATK format,
clusters, then converts *back* to svtk format at its own published
output — that's `caller_vcfs`, this subworkflow's input — so `ECN` is
already gone by the time `combine_batches.nf` hands it to stage 2.
Checking GATK-SV's own WDL for how they avoid this: `CombineBatches.wdl`
takes `pesr_vcfs`/`depth_vcfs` as plain inputs with no visible
GATK-formatting call of its own (its `import "FormatVcfForGatk.wdl"` is
unused dead code in the current file) — because in their real pipeline
those inputs are `GenotypeBatch.genotyped_pesr_vcf`, and `GenotypeBatch`'s
own `GenotypeSVs` GATK task re-adds `ECN` internally as part of
genotyping. In other words: **GATK-SV's real ClusterBatch → CombineBatches
path always has a genotyping stage in between that happens to restore
`ECN`**; it's not that their stage 2 tolerates svtk-format input. We don't
implement genotyping yet, so `caller_vcfs` here really is stage 1's final,
`ECN`-stripped output, and stage 2 has to redo the svtk→GATK conversion
itself before any of `CombineBatches`' own clustering steps can run.

Fixed by adding a step 0 to `vcfs_combine_batches`: re-run
`FORMAT_SVTK_VCF_FOR_GATK` (aliased `FORMAT_SVTK_VCF_FOR_GATK_COMBINE_BATCHES`,
reusing the same vendored script and stage-2 ploidy table) on each
caller's VCF before `GATK4_SVCLUSTER_JOIN`. Confirmed safe to reuse
as-is on a multi-sample, cross-sample-clustered VCF (not just the
per-sample VCFs it's used on at stage 1): the script's `convert()`
function iterates `record.samples.items()` generically, with no
single-sample assumption anywhere in it. `FORMAT_SVTK_VCF_FOR_GATK` itself
emits no `.tbi` (stage 1 gets one for free from a downstream
interval-exclusion step's own `tabix` call, which stage 2 doesn't have),
so a small new local module, `modules/local/tabix` (`TABIX`), was added
to index explicitly — preferred over nf-core's `tabix/tabix` (checked
first; it's deprecated and now hard-asserts false, pointing at
`htslib/bgziptabix` instead) or `htslib/bgziptabix` itself (a general
compress-and/or-index module with file-type-sniffing branches that don't
apply when the input is always already bgzipped, as it is here) — a
plain `tabix <file>` local module matches the existing small-local-module
style already used for `EXCLUDE_VARIANTS_BY_ID`/`GENOME_FILE`. Invisible
under `-stub-run`: stub VCFs carry no real `FORMAT` fields to begin with,
so there's nothing for GATK to find missing. **Lesson for any subworkflow
that consumes another subworkflow's *published* (svtk-format) output as
input to more GATK-walker processing**: check whether the producer's
publish step stripped anything the consumer's GATK tools need back — the
producer and consumer being individually correct doesn't mean the round
trip between them is lossless.

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

### Per-sample evidence collection (`bam_collect_evidence`)

`subworkflows/local/bam_collect_evidence` and its standalone entry point
`workflows/collect_evidence.nf` collect raw per-sample PE (discordant
pairs)/SR (split reads)/SD (site depth at known SNP sites)/RD (binned
coverage) evidence from a BAM/CRAM — the sibling of `bam_call_manta`/
`bam_call_wham` at the same tier (one BAM in, per-sample output out), not
downstream of clustering/harmonisation. Mirrors the evidence-collection
half of GATK-SV's `GatherSampleEvidence` workflow
(`wdl/GatherSampleEvidence.wdl`), narrowed to just PE/SR/SD/RD collection
(not the SV calling it also bundles — already covered here by
`bam_call_manta`/`bam_call_wham`).

This is deliberately **not** wired into anything downstream yet — no
merge step, no genotyping. It exists because a future genotyping stage
(GATK's `TrainSVGenotyping`/`GenotypeSVs` walkers) needs this evidence,
merged panel-wide, as an input; collecting it per-sample now is the first
of several unbuilt pieces on that path (see [Not yet
started](#not-yet-started)). See the subworkflow's own top-of-file note
for the full reasoning chain (traced through GATK-SV's WDL: their real
`ClusterBatch` → `GenotypeBatch` → `CombineBatches` ordering, and why
`ECN` — see the fix above — is a symptom of us not yet having a
genotyping stage in between).

Two GATK walkers:

1. `CollectSVEvidence` (`GATK_COLLECT_SV_EVIDENCE`) → PE/SR/SD in one BAM
   pass. Reproduces `wdl/CollectSVEvidence.wdl`'s `RunCollectSVEvidence`
   task.
2. `CollectReadCounts` (`GATK_COLLECT_READ_COUNTS`) → RD (binned
   coverage), over `preprocessed_intervals`. Reproduces
   `wdl/CollectCoverage.wdl`'s `CollectCounts` task.

#### Neither walker's output can be compressed/indexed in the same container

Both walkers live in the same `gatk4-main_gcnvkernel` container already
used by `gatk4/svcluster`/`grouped_sv_cluster` — verified directly
(`gatk CollectSVEvidence --help` / `gatk CollectReadCounts --help`) rather
than assumed, avoiding a pull of GATK-SV's own separate `gatk_docker`
image just for these two walkers. But that container has neither `bgzip`
nor `tabix` (verified directly: `which bgzip tabix` found nothing, only
`sed`) — unlike GATK-SV's own WDL tasks, which call `tabix`/`bgzip`
inline in the same command block as the `gatk` invocation, assuming a
richer image. A Nextflow process has exactly one container, so this had
to become three small modules, not one:

- `CollectSVEvidence`'s PE/SR/SD output needs only *indexing* — GATK's
  walker writes it already bgzipped via htsjdk's own feature-file writer
  (the same mechanism GATK's VCF writers use for `.vcf.gz`), confirmed by
  `wdl/CollectSVEvidence.wdl` itself calling `tabix` directly on the
  walker's raw output with no separate `bgzip` step first. New module
  `modules/local/tabix_sv_evidence` (`TABIX_SV_EVIDENCE`) — a `tabix -f
  -0 -s1 -b2 -e2` call, using the same `bcftools_htslib` container as
  `modules/local/tabix`. Not folded into that existing module: PE/SR/SD
  files aren't VCFs, so they need those explicit format flags — plain
  `tabix <file>` (auto-detect, `modules/local/tabix`'s one real use case)
  doesn't apply.
- `CollectReadCounts`'s RD output (`--format TSV`) is genuinely plain
  text, needing an actual `bgzip` call GATK-SV's own WDL makes explicitly
  after a `sed` rewrite (of the `@RG` header's `SM:` tag, to `ped_id` —
  see [Sample and family ID
  constraints](#sample-and-family-id-constraints) for why `ped_id`, not
  `sample_id`, is used for every PED-keyed value in this pipeline). New
  module `modules/local/bgzip` (`BGZIP`) for the compression step; `sed`
  itself is present in the GATK container so that part stays inline.

`modules/local/gatk4/collect_sv_evidence`'s own top-of-file note has the
same reasoning, closer to the code.

#### `sd_locs_vcf` is the one static resource that isn't bgzip/`.tbi`

GATK-SV's `sd_locs_vcf` resource
(`Homo_sapiens_assembly38.dbsnp138.vcf`) is a plain, uncompressed VCF
with a GATK-style `.idx` sibling index — every other VCF this pipeline
handles is bgzipped with a `.tbi`. `CollectSVEvidence`'s
`-F`/`--site-depth-locs-vcf` flag takes only the VCF path (confirmed via
`--help`; no separate index flag), so `sd_locs_vcf_idx` is threaded
through as a plain sibling `path` input for GATK's own index
auto-discovery to find, not passed on the command line itself.

#### `CollectReadCounts` OOM'd — our 3GB default was a mis-transcription of GATK-SV's own task

A real HPC run of `collect_evidence.nf` (~50x WGS BAM) had
`GATK_COLLECT_READ_COUNTS` complete a full 91-minute whole-genome
traversal (729M reads processed) and then hit `java.lang.OutOfMemoryError:
Java heap space` in `CollectReadCounts.onTraversalSuccess` — building the
final binned-count list in memory before writing output, a real cost
proportional to genome-wide interval count, not per-chunk. The module was
falling back to its own hardcoded 3GB default (no `task.memory` set for
`process_single`, same fallback-logging pattern as `GATK4_SVCLUSTER`/
`GroupedSVCluster`).

That 3GB came from misreading `wdl/CollectCoverage.wdl`'s `CollectCounts`
task: its `RuntimeAttr default_attr` sets `mem_gb: 3.0`, which looks like
the task's real default — but that field is never actually consulted for
the JVM heap size. The task's `command <<<...>>>` block computes
`machine_mem_gb` from a *separate* top-level `mem_gb` input with its own
independent fallback: `Float machine_mem_gb = select_first([mem_gb,
12.0])`. GATK-SV's real effective default is **12GB**, 4x what we'd
copied. Fixed in `conf/modules.config` with a `withName:
'GATK_COLLECT_READ_COUNTS'` override, `memory = 12.GB` — matching
GATK-SV's own validated value directly, unlike `WHAMG`'s fixed 32GB
override above (picked from one observed failure point, no equivalent
upstream number to match). `CollectSVEvidence`'s own `RuntimeAttr` doesn't
have this trap — its `mem_gb: 3.75` *is* what its memory calculation
actually resolves to, no hidden second input shadowing it — so it was
left untouched; only `CollectReadCounts` had this bug. **Lesson for
reading any GATK-SV WDL task's `RuntimeAttr default_attr` block**: check
that `mem_gb` (or whichever field) is actually threaded into the memory
calculation used by `runtime { memory: ... }`, not shadowed by a
same-named-but-different top-level input with its own default — a
`RuntimeAttr` object existing in the task doesn't guarantee it's live.

### Panel-wide evidence merging (`vcfs_merge_evidence`, `vcfs_merge_read_counts`)

Two subworkflows, composed by the standalone entry point
`workflows/merge_evidence.nf`, merge every sample's per-sample evidence
(from `bam_collect_evidence`) into one panel-wide matrix each — the last
unbuilt prerequisite for genotyping (`TrainSVGenotyping`/`GenotypeSVs`,
still not implemented — see [Not yet started](#not-yet-started)), which
needs `--rd-file`/`--discordant-pairs-file`/`--split-reads-file` inputs at
exactly this panel-wide shape.

Split into two subworkflows, not one, because they're genuinely different
in kind:

- `vcfs_merge_evidence` merges PE, SR, and BAF, mirroring GATK-SV's
  `BatchEvidenceMerging` workflow (`wdl/BatchEvidenceMerging.wdl`). Two
  GATK walkers: `PrintSVEvidence` (`GATK_PRINT_SV_EVIDENCE`, aliased
  twice — once for PE, once for SR) merges same-type files across every
  sample into one file; `SiteDepthtoBAF` (`GATK_SITE_DEPTH_TO_BAF`)
  converts+merges every sample's SD file into one BAF file in the same
  step, using `sd_locs_vcf`. GATK-SV's `BatchEvidenceMerging` accepts
  *either* pre-existing BAF files or SD files (via `SDtoBAF`) — only the
  SD path is implemented, since `bam_collect_evidence` only ever produces
  SD, never raw BAF.
- `vcfs_merge_read_counts` merges RD (binned coverage), mirroring
  GATK-SV's `MakeBincovMatrix` workflow (`wdl/MakeBincovMatrix.wdl`).
  Plain shell (`awk`/`paste`/`bgzip`), not a GATK walker — three new local
  modules (`BINCOV_SET_BINS`, `BINCOV_MAKE_COLUMNS`, `BINCOV_ZPASTE`),
  same "small chained shell modules, `bcftools_htslib` container" pattern
  as `TABIX`/`BGZIP`/`EXCLUDE_VARIANTS_BY_ID` elsewhere in this pipeline.
  Genuinely different shape from every other subworkflow so far: one
  sample's counts are picked as the reference bin grid
  (`BINCOV_SET_BINS`, via `.first()` on the collected channel), *then*
  every sample (including the reference one) is reshaped against that
  grid and pasted together column-wise — a "pick one, then fan out"
  structure none of this pipeline's other merge steps have needed.

Deliberately simpler than the WDL in one respect, for both subworkflows:
GATK-SV's `MergeEvidence`/`SDtoBAF` tasks have an optional
`rename_samples` step (rewriting each file's sample column to the batch's
canonical sample list) because their per-sample evidence files are keyed
by whatever `sample_id` was passed to `GatherSampleEvidence`, which may
not match the batch's own naming. This pipeline doesn't have that
mismatch to correct: `bam_collect_evidence`'s `GATK_COLLECT_SV_EVIDENCE`
already writes every PE/SR/SD file with `--sample-name` set to
`meta.ped_id` at collection time — the same PED-keyed identity every
other lookup in this pipeline already uses — so there's nothing left to
rename by the time evidence reaches the merge step. A real, if quiet,
payoff from the `ped_id` design decision made much earlier in this
project (see [Sample and family ID
constraints](#sample-and-family-id-constraints)).

#### `env` process outputs need a quoted string, not a bare identifier

`BINCOV_SET_BINS` needed to emit a shell variable (`BINSIZE`, computed at
runtime from the data) as a process output, for
`BINCOV_MAKE_COLUMNS`/`BINCOV_ZPASTE` downstream to consume as a `val`.
The first attempt, `env BINSIZE, emit: binsize` (bare identifier, by
analogy with `path`/`val`'s own syntax), failed at *script compilation*
with `` `BINSIZE` is not defined `` — not a runtime error, so
`-stub-run` caught it immediately, before any real execution. The correct
syntax is `env 'BINSIZE', emit: binsize` (a quoted string) — confirmed
empirically with an isolated test process, since this distinction isn't
obvious from output-block syntax that otherwise looks uniform
(`path("...")`, `val(...)`, `env(...)` would all parse, but only the
quoted-string form actually works). Also needs a real `export` in the
`script:`/`stub:` body (`export BINSIZE=...`, not a bare `BINSIZE=...`
shell assignment) for Nextflow's output capture to see it at all.

#### `.combine()` flattens list-valued channel items instead of nesting them

`vcfs_merge_read_counts` needed to pair one collected list (every
sample's bincov column, `List<Path>`) with one single value (the shared
bin-locations file) into a single process call's input tuple. The
straightforward `columns_channel.combine(bin_locs_channel)` produced a
flat list — `[colA, colB, locsfile]` — not the nested `[[colA, colB],
locsfile]` a `.map { cols, locs -> ... }` destructuring assumed, causing
`Invalid method invocation call with arguments: [...] (java.util.LinkedList)
on _closure_ type`. Confirmed via an isolated test
(`channel.of([1,2,3]).combine(channel.of('x'))` emits `[1,2,3,'x']`, not
`[[1,2,3],'x']`) — `.combine()` flattens *any* list-valued item against
whatever it's paired with, rather than treating a list as one atomic
tuple element. Fixed by wrapping the list in an extra list first
(`.map { cols -> [ cols ] }`) before combining, which keeps it
distinguishable from the scalar side afterward. **Lesson for pairing a
collected list with a single broadcast value via `.combine()`**: verify
the actual emitted shape (e.g. with a throwaway `.view()`) rather than
assuming standard list/tuple semantics — this behavior isn't obvious from
the operator's name and wasn't caught by `-stub-run`'s wiring check until
the actual channel shapes met at runtime (a compile-time-clean, then
runtime `Invalid method invocation` error, one step later than the `env`
syntax issue above).

#### `zcat | head -c 1` under global `pipefail` SIGPIPEs `zcat` — a real HPC failure

A real HPC run of `merge_evidence.nf` failed with `BINCOV_SET_BINS`
exiting 141 (SIGPIPE, 128+13) on `zcat -f ... | sed ... | ... > tmp_locs`
— the whole pipeline, not any one command's own logic, and with empty
`Command output`. Root cause: an *earlier* line in the same script,
`firstchar=$(zcat -f ${count_file} | head -c 1)`, used `head -c 1` to
read just the first byte — a completely standard "peek at the start of a
stream" idiom. `head -c 1` exits as soon as it has its one byte, closing
its end of the pipe while `zcat` may still be writing; `zcat` then gets
`SIGPIPE`. Under a normal (non-strict) shell this is silent and harmless
— it's what `head` is *for*. But this pipeline's global `pipefail`
(`nextflow.config`, deliberately set — see "The OOM above should have
failed the task outright, and didn't" earlier, in the harmonisation
stage 1 section) means a command substitution or pipeline containing that SIGPIPE
reports exit 141 as *its own* result — and because this was inside a
`firstchar=$(...)` assignment under `set -e` (part of the same global
`pipefail` line), that 141 killed the whole script, not just the
substitution. The exact same pattern (`sed -n 'Np'`, `head`, anything
that deliberately stops reading early) anywhere downstream of a live pipe
has the identical risk — confirmed a second, not-yet-exercised instance
of the same bug in `BINCOV_MAKE_COLUMNS` (same `zcat | head -c 1` line)
while fixing the first.

Fixed in both modules by decompressing to a plain file first (`zcat -f
${count_file} > tmp_raw`), then reading from that file with `head`/`sed`
afterward — a `head`/`sed -n` reading from a *file*, not a pipe, has no
upstream process left to `SIGPIPE` when it stops early. Verified for real
(not just `-stub-run`, which never exercises `script:` and so never had
a chance to catch this): a standalone test workflow running
`BINCOV_SET_BINS` and `BINCOV_MAKE_COLUMNS` against a synthetic
5000-line, multi-KB gzipped counts file, under the same `-ueo pipefail`
shell setting as production, both completed successfully with correct
output. **Lesson for writing shell inside a `script:` block on this
pipeline specifically**: any `head -n`/`head -c`/`sed -n 'Np'`/similar
early-terminating command must read from a file, not sit downstream of a
still-live pipe — global `pipefail` makes an otherwise-idiomatic and
harmless pattern into a real, silent-until-you-hit-it failure mode. When
porting shell from a GATK-SV WDL task (most of which run under plain
`set -euo pipefail` too, so this risk already existed there in theory,
just apparently never triggered), check every pipe for an early-exiting
stage, not just the ones that look unusual.

### Path to genotyping: the real prerequisite chain

`merge_evidence.nf` validated cleanly on a real HPC run (both `-profile
docker` locally and a real-data HPC run), closing out panel-wide evidence
merging. The obvious next step is genotyping
(`TrainSVGenotyping`/`GenotypeSVs`, GATK walkers in `GenotypeBatch.wdl`),
which consumes the merged PE/SR/BAF/RD evidence directly. But tracing
`GenotypeBatch.wdl`'s own inputs back through GATK-SV's WDL shows a longer
dependency chain than that framing suggests — `rf_cutoffs` (the random
forest genotyping cutoffs `TrainSVGenotyping` needs) is not a static
resource, it's computed per-cohort by `FilterBatchSites`/`AdjudicateSV`
(`svtk adjudicate`), which itself consumes `evidence_metrics` — a whole
separate metrics-aggregation subworkflow, `GenerateBatchMetrics`, not yet
built here. Decided to port this faithfully rather than stub or simplify
`rf_cutoffs` early, so the real build order is three pieces in sequence,
not one:

1. **`GenerateBatchMetrics`** (done — see task breakdown below):
   per-variant evidence metrics from the sites VCF + merged PE/SR/BAF/RD,
   producing the flat `evidence_metrics.tsv` that `AdjudicateSV` needs.
2. **`FilterBatchSites`/`AdjudicateSV`** (done — see its own section
   below): runs `svtk adjudicate` against those metrics to emit
   `rf_cutoffs` (and `scores`). Originally scoped to just `AdjudicateSV`
   itself (`FilterAnnotateVcf`/`PlotSVCountsPerSample` seemed unnecessary
   for genotyping) — reversed once `GenotypeBatch.wdl`'s real `vcf` input
   was traced through `wdl/GATKSVPipelineBatch.wdl`'s own wiring and
   turned out to be `FilterAnnotateVcf`'s own sites-filtered output
   (merged across callers), not `GenerateBatchMetrics`'/`CombineBatches`'
   raw sites VCF — so `FilterAnnotateVcf` ended up built after all (see
   its own section below). `PlotSVCountsPerSample` (outlier-sample QC
   plotting) remains out of scope — genuinely not needed by anything
   downstream.
3. **`GenotypeBatch`** (`TrainSVGenotyping`/`GenotypeSVs`): the genotyping
   step itself, now unblocked. **Done** — see its own section below.

#### `GenerateBatchMetrics` task breakdown

Traced from `wdl/GenerateBatchMetrics.wdl`. Chain: concat all-caller VCFs
→ scatter into shards → per shard (`FormatVcf` svtk→GATK conversion →
`SVRegionOverlap` segdup/repeat-mask annotation → `AggregateSVEvidence`
PE/SR/BAF-based metrics → `AggregateDepthEvidence` RD-based metrics) →
concat shards back → `AggregateTests` (`aggregate.py`) produces the final
`evidence_metrics.tsv`.

| Step | Status |
| --- | --- |
| `CreatePloidyTableFromPed` | existing (`ploidy_table_from_ped`) — reuse |
| `ConcatVcfs` (input + output) | **done** — `CONCAT_VCFS` module (`bcftools concat`), aliased as `CONCAT_INPUT_VCFS`/`CONCAT_OUTPUT_VCFS` with per-alias `ext.args` (`--allow-overlaps`/`--naive`, `conf/modules.config`) — this did NOT already exist anywhere in this pipeline despite an earlier (wrong) note here saying "reuse"; genuinely new. Validated for real: a 2-shard round-trip (scatter then concat) reproduced the original 2-record VCF exactly, byte-for-byte record content |
| `FormatVcf` (svtk→GATK) | existing (`FORMAT_SVTK_VCF_FOR_GATK`) — reuse per-shard |
| `median_coverage` | **done** — `MEDIAN_COVERAGE` module, ports `MedianCov.wdl`'s `WGD/bin/medianCoverage.R` (vendored verbatim as `bin/medianCoverage.R`); needed here (by `AggregateSVEvidence`/`AggregateDepthEvidence`), one layer earlier than genotyping itself needs it. Validated for real (not just `-stub-run`): a synthetic bincov matrix through a standalone test workflow, output checked by hand against the expected transpose |
| `ScatterVcf` | **done** — `SCATTER_VCF` module, wraps `bcftools +scatter`. Validated for real: a 2-record VCF split at `-n 1` (2 shards, records not duplicated/dropped, checked by hand) and consolidated at `-n 10` (1 shard); also validated the WDL's own empty-VCF placeholder-shard fallback against a genuinely empty (0-record) VCF, confirming `bcftools +scatter` itself produces zero output files in that case (checked empirically, not assumed from the WDL's comment) |
| `SVRegionOverlap` | **done** — `GATK_SVREGIONOVERLAP` module; static resources `segdups`/`rmsk` already in `resources_hg38.json`, wired as two fixed named tracks (not a general list — this pipeline only ever has these two). Validated for real against a synthetic VCF + synthetic segdup/rmsk BED tracks; output overlap fractions/endpoint counts checked by hand |
| `AggregateSVEvidence` | **done** — `GATK_AGGREGATESVEVIDENCE` module; consumes merged PE/SR/BAF + `median_coverage`. Validated for real end-to-end, including genuinely-collected (not hand-faked) PE/SR/BAF evidence from a real synthetic BAM run through actual `CollectSVEvidence`/`SiteDepthtoBAF` |
| `AggregateDepthEvidence` | **done** — `GATK_AGGREGATEDEPTHEVIDENCE` module; consumes merged RD + `median_coverage`. Validated for real |
| `aggregate.py` → `AggregateTests` | **done** — `AGGREGATE_TESTS` module, vendors `src/sv-pipeline/02_evidence_assessment/02e_metric_aggregation/scripts/aggregate.py` into `bin/aggregate.py` (adapted, not verbatim — see `bin/README.md`'s own note: inlines the one `svtk.utils` function it needed, `get_called_samples`, to avoid a `pybedtools`/`bedtools` container dependency for code the script never otherwise calls). Validated for real against a synthetic annotated VCF, output metrics checked by hand (vf, rmsk, RDQ, poor_region_cov all correct) |
| `vcfs_generate_batch_metrics` subworkflow + `workflows/generate_batch_metrics.nf` | **done** — composes all of the above (plus `bam_call_manta`/`bam_call_wham`/`vcfs_cluster_svcluster` for stage-1 per-caller clustering and `bam_collect_evidence`/`vcfs_merge_evidence`/`vcfs_merge_read_counts`/`median_coverage` for panel-wide evidence, same "one entry point runs the full chain" pattern as `combine_batches.nf`/`merge_evidence.nf`). **Deliberately does not run `vcfs_combine_batches`**: confirmed against `wdl/GATKSVPipelinePhase1.wdl`'s own wiring that `GenerateBatchMetrics` takes each caller's stage-1 clustered VCF directly (`ClusterBatch.clustered_*_vcf`), not `CombineBatches`' later cross-caller merged output — an easy wrong assumption from the workflow's name alone, checked against the WDL's actual call site rather than guessed. Validated for real end-to-end (real per-shard evidence annotation, not just wiring) and via full-sample-sheet `-stub-run` (78 processes, `BAM_CALL_MANTA` through `AGGREGATE_TESTS`) |

Built `median_coverage` and vendored `aggregate.py` first since they were
the two independent, no-interdependency pieces — no wiring to any other
new module required to validate either in isolation. Both containers are
new (`dockerfiles/median-coverage`: `r-base` + CRAN `optparse`;
`dockerfiles/aggregate-tests`: `python:3.12-slim` + `pysam`/`numpy`/
`pandas`), built and validated locally, and now pushed to Docker Hub
(`drtomc/gatk-sv-nf-median-coverage:0.1.0`,
`drtomc/gatk-sv-nf-aggregate-tests:0.1.0`) — available on HPC/any other
machine, not just this one's local Docker daemon.

#### A significant finding: five GATK-SV walkers only exist in GATK-SV's own GATK fork, not bioconda

Before writing `SVRegionOverlap`/`AggregateSVEvidence`/`AggregateDepthEvidence`,
tried the obvious thing: reuse `gatk4-main_gcnvkernel`, the same container
already used by `GATK4_SVCLUSTER`/`GATK_PRINT_SV_EVIDENCE`/every other
`gatk4/*` module here. `gatk --list` on that image has no
`SVRegionOverlap`/`AggregateSVEvidence`/`AggregateDepthEvidence` at all —
confirmed directly, not assumed. Pulling GATK-SV's own `gatk_docker`
(`inputs/values/dockers.json`: `us.gcr.io/broad-dsde-methods/gatk-sv/gatk:mw-gatk-sv-53d5c2d`)
and running `gatk --list` there shows all three, plus (relevant later)
`TrainSVGenotyping`/`GenotypeSVs`. The two images are genuinely different
GATK builds — `gatk4-main_gcnvkernel` is bioconda's `gatk4-main=4.7.0.0`;
GATK-SV's own image is a patched fork, `4.6.2.0-92-g6379d28-SNAPSHOT` —
not just a version difference but a different codebase, since these five
walkers are exclusive to GATK-SV's fork and aren't in stock/bioconda GATK
at all.

Decided to point these three new modules' `container` directly at
GATK-SV's own GCR image, same as upstream WDL's `gatk_docker`, rather
than trying to get a smaller substitute built. This breaks from every
other `gatk4/*` module's small/Wave-hosted-image pattern (a genuine
inconsistency, not an oversight) — accepted because the alternative
(waiting on/building a from-source container for GATK-SV's own fork) is
much more work for uncertain payoff. **This same constraint will apply to
`TrainSVGenotyping`/`GenotypeSVs`** (`GenotypeBatch`, still queued after
this) — plan for the same container choice there, don't re-discover this.

The `container` line itself is a plain image reference (not `docker://`,
not engine-conditional) — same reasoning as
`modules/local/ploidy_table_from_ped`'s own note: Apptainer/Singularity
accept a bare `registry/path:tag` reference (defaulting to
Docker-Hub-style pulling for any unscoped reference, not just Docker Hub
itself), and plain `docker` requires the bare form anyway, so one
reference works for both engines.

No `environment.yml` added for these three modules, unlike every other
`gatk4/*` module here — those exist for nf-core-template-convention
reasons but are dead weight in practice (no `conda` profile is configured
anywhere in `nextflow.config`/`conf/`), and adding one here specifically
would be actively misleading: there is no bioconda package that provides
these walkers to pin a version against.

All three modules validated for real, not just wiring: a synthetic BAM
(hand-built via `samtools`/a literal SAM file, since `tests/data/`'s own
`SAMPLE_A.bam` is a 0-byte `-stub-run`-only placeholder, not a real BAM —
confirmed by trying to actually run `merge_evidence.nf` non-stub against
it, which failed for exactly that reason) run through the real
`CollectSVEvidence`/`SiteDepthtoBAF` walkers to get genuinely-shaped
PE/SR/SD/BAF evidence, rather than hand-faking the tabix-indexed evidence
file schemas and risking a false-confidence test. One invocation mistake
caught and fixed this way: `CollectSVEvidence`'s `-F`/site-depth-locs flag
was initially guessed from `--help`'s wording
(`--site-depth-locs-vcf`, which doesn't exist) instead of checked against
`wdl/CollectSVEvidence.wdl`'s own actual invocation (`-F`) — corrected
before it went in any module, but a reminder to verify a real invocation
against the WDL's own command line, not `--help` text, when in doubt.

Next: the `vcfs_generate_batch_metrics` subworkflow and
`workflows/generate_batch_metrics.nf` entry point, wiring everything
above together -- the last piece of `GenerateBatchMetrics`.

#### Wiring `vcfs_generate_batch_metrics`: three real bugs, all caught by real (non-stub) validation

`GenerateBatchMetrics` is now fully built: `vcfs_generate_batch_metrics`
composes every module above, `workflows/generate_batch_metrics.nf` composes
that with `bam_call_manta`/`bam_call_wham`/`vcfs_cluster_svcluster` (stage
1) and `bam_collect_evidence`/`vcfs_merge_evidence`/
`vcfs_merge_read_counts`/`median_coverage` (panel-wide evidence) -- same
"one entry point runs the full chain" pattern as `combine_batches.nf`/
`merge_evidence.nf`. Validated for real end-to-end (a genuine 2-record
VCF through the whole chain to `evidence_metrics.tsv`, output checked by
hand) and via a full-sample-sheet `-stub-run` (78 processes). None of the
three bugs below were visible from writing the subworkflow code alone --
each one only appeared once real data was pushed all the way through:

- **`GATK_AGGREGATESVEVIDENCE`'s `-V`/`-O` resolved to the same file.**
  Every per-shard stage
  (`GATK_SVREGIONOVERLAP`/`GATK_AGGREGATESVEVIDENCE`/
  `GATK_AGGREGATEDEPTHEVIDENCE`) derives its own output filename from
  `meta.id`, but the shard's `meta.id` was carried unchanged through the
  whole chain -- so `AggregateSVEvidence` read and overwrote its own input
  file in place, silently succeeding with **0 variants processed** rather
  than erroring (empty output looks like valid output at a glance). Same
  root cause `vcfs_combine_batches` already documents for its own
  multi-stage `GATK4_SVCLUSTER`/`GATK4_GROUPEDSVCLUSTER` chain -- reading
  that note first didn't prevent writing the same bug again here, since
  the actual renaming step was missed, not misunderstood. Fixed by giving
  each stage its own `meta.id` suffix (`.svregionoverlap`,
  `.aggregatesvevidence`).
- **`CONCAT_OUTPUT_VCFS` failed on real shard output: "Unsorted positions
  ... 400 followed by 100."** WDL/Cromwell's `scatter` always preserves
  index order in its output array; Nextflow's `.collect()` on a scattered
  channel gathers emissions in whatever order the parallel shard
  processes actually finish, not shard order. Fixed by sorting
  shards by filename before the final `--naive` concat (safe here since
  `SCATTER_VCF`'s own shard names are zero-padded and lexically sortable)
  -- the WDL's own `ConcatVcfs` task offers exactly this as its
  `sort_vcf_list` option, unused at this particular WDL call site since
  Cromwell doesn't need it there.
- **The sort fix's first attempt silently did nothing.** Sorted on
  `it[0].toString()` (the shard's full absolute path, e.g.
  `.../work/c0/c41846.../shard_000000.vcf.gz`) instead of `it[0].name`
  (just the filename) -- this sorts by the **work-directory hash
  prefix**, which is random per run, not by the shard number the sort was
  actually meant to fix. Confirmed by testing the sort in isolation before
  concluding it was correct, then re-testing against the real subworkflow
  and finding the output order genuinely unchanged -- the isolated test
  used files with meaningful literal paths (`/tmp/a...`, `/tmp/z...`), so
  it passed while hiding the exact real-path failure mode.

#### `MEDIAN_COVERAGE` OOM'd for real on a real HPC run — no memory default existed at all

A real HPC run of `genotype_batch.nf` (a few tens of samples) had
`MEDIAN_COVERAGE` killed outright (exit 137, SIGKILL) partway through
`Rscript medianCoverage.R`, right after it finished loading the
`optparse` package — easy to misread as a missing-R-package problem from
the log alone, but exit 137 with no R error message is the standard
signature of an out-of-memory kill, not a package/environment issue.
Root cause: unlike most other memory-hungry modules in this pipeline
(`WHAMG`, `GATK_COLLECT_READ_COUNTS`, both already have their own
documented fixes above), `MEDIAN_COVERAGE` had **no memory default at
all** — `label 'process_single'` alone, and `process_single` itself has
no `memory` set in `nextflow.config` (only `process_medium` does) — so it
ran with whatever minimal implicit default Nextflow/the executor falls
back to. `medianCoverage.R`'s own top-of-file comment already documents
the risk ("Note: loads entire coverage matrix into memory... may pose a
problem for large matrices"), so this was a predictable gap in hindsight,
just not one caught by this module's own earlier validation (a tiny
synthetic bincov matrix, far too small to hit any real memory ceiling).
Fixed with an explicit `memory = 24.GB` in `conf/modules.config` —
deliberately *not* GATK-SV's own `wdl/MedianCov.wdl` figure (a flat 80GB,
"~0.5Gb per sample"), which is calibrated for their biobank-scale cohorts
(hundreds-thousands of samples); this pipeline's own panel scale (a few
tens of samples) doesn't need that much, so 24GB was chosen as generous
headroom over the raw per-sample matrix size without provisioning for a
cohort two orders of magnitude larger than this pipeline ever runs.
**Revisit with real observed memory usage if this OOMs again** — 24GB is
a reasoned estimate, not a value derived from an actual successful run's
peak usage the way `GATK_COLLECT_READ_COUNTS`'s 12GB fix was.

### `FilterBatchSites`/`AdjudicateSV`: random-forest cutoff derivation

`ADJUDICATE_SV` module, wrapping GATK-SV's `svtk adjudicate` (vendored,
adapted — see `bin/README.md`'s own note: calls the vendored
`svtk.adjudicate` logic directly rather than installing `svtk` itself, to
avoid a `pybedtools`/`bedtools`/Cython container dependency for logic that
doesn't need it). Takes `GenerateBatchMetrics`'s `evidence_metrics.tsv`,
emits `scores`/`cutoffs`/`RF_intermediate_files.tar.gz` — `cutoffs` is
`rf_cutoffs`, the direct input `TrainSVGenotyping` needs.
`PlotSVCountsPerSample` (outlier-sample QC plotting) remains out of scope
— see [Path to genotyping](#path-to-genotyping-the-real-prerequisite-chain)
above. `FilterAnnotateVcf`, this WDL's third task, is built separately
below (its own section) once `GenotypeBatch`'s real `vcf` input turned
out to need it after all.

Validated for real against a synthetic `evidence_metrics.tsv` and via
`-stub-run`. Output cutoffs file format double-checked column-by-column
against `wdl/GenotypeBatch.wdl`'s own consumption
(`awk -F '\t' '{if ($5=="PEQ") print $2 }'` — column 5 is this file's own
`metric` column, column 2 is `cutoff`; confirmed by grepping the real
output, not assumed from column names alone).

#### Getting realistic synthetic training data took several iterations — a real lesson about this specific random forest, not a code bug

Building a `evidence_metrics.tsv` fixture large enough to exercise the
random forest hit a genuine dead end worth recording: a naive synthetic
dataset with two cleanly-separated "strong signal" / "weak signal"
clusters (no overlap between them) reliably made `adjudicate_PESR` (the
7th and last of `adjudicate_SV`'s internal passes) fail with "No Fail
variants included in training set" — even though the earlier six passes
(`BAF1`, `SR1`, `RD`, `PE`) all succeeded on the exact same two-cluster
shape. Traced (not guessed) to `random_forest.py`'s `learn_cutoffs()`:
for `name in ("PE_prob", ...)` it takes `min(learn_cutoff_dist(...),
learn_cutoff_fdr(...))`, and `learn_cutoff_fdr`'s 5%-false-positive-rate
ROC threshold degenerates to a value extremely close to zero when the
training data has literally zero distributional overlap between Pass and
Fail — meaning even a "weak" synthetic row nearly always exceeds the
learned cutoff and gets marked as passing, regardless of how small its
raw feature value is. Confirmed directly: calling `RandomForest.run()` by
hand on the same two-cluster data gave the textbook-correct answer (weak
rows correctly failed); only `adjudicate_PE`'s own `name == "PE_prob"`
branch (invoking the FDR path) produced the wrong-looking cutoff. This
is not a bug in the vendored code — real GATK-SV metrics have continuous,
naturally-overlapping distributions this FDR logic is built for; a
synthetic fixture needs the same property (each row's metrics as noisy
functions of one continuous "true positive-ness" signal, not two discrete
clusters) to exercise this classifier meaningfully. **Lesson for any
future synthetic test data for this module**: don't hand-pick two
non-overlapping value ranges for "should pass" vs "should fail" — sample
every metric as continuous noise around a single underlying quality
signal instead, or the FDR-based cutoff logic specifically won't behave
like it would on real data.

### `FilterAnnotateVcf`: RF-score sites filtering

`FILTER_ANNOTATE_VCF` module + `filter_batch_sites` subworkflow. Built
after `AdjudicateSV`, once tracing `GenotypeBatch.wdl`'s real `vcf` input
through `wdl/GATKSVPipelineBatch.wdl`'s own wiring showed it's
`FilterAnnotateVcf`'s own sites-filtered output (merged across callers via
`MergePesrDepthVcfs`, `allow_overlaps=true`), not `GenerateBatchMetrics`'
raw sites VCF and not `CombineBatches`' output either — an easy wrong
assumption to make from the workflow names alone, the same kind of
WDL-call-site-not-workflow-name check that mattered for
`vcfs_generate_batch_metrics` earlier.

Reproduces `FilterBatchSites.wdl`'s `FilterAnnotateVcf` task exactly:
score-based filtering (`score >= 0.5` from `AdjudicateSV`'s `scores`),
INV/BND/INS-scoring-passers reclassified as `BND`, an inline Python
`END2`/`CHR2` backfill, then two vendored scripts
(`rewrite_SR_coords.py`/`annotate_RF_evidence.py`) for SR-based breakpoint
coordinate correction and per-record evidence-type annotation. `svtype`
per record. `filter_batch_sites` runs this once per caller (currently
Manta + Wham, both PESR-type — this pipeline has no depth caller yet, see
[Not yet started](#not-yet-started)), then merges across callers.

New container (`dockerfiles/filter-annotate-vcf`): needs both
`bcftools`/`bgzip` (the task's own inline shell pipeline) *and*
`pysam`/`numpy`/`pandas` (the two vendored scripts) in the same
container, since a single Nextflow process has exactly one container —
neither this pipeline's existing `bcftools_htslib` image nor its
`aggregate-tests` image (has the Python packages, no `bcftools`) covers
both, so this combines both rather than splitting one WDL task's single
inline pipeline across two Nextflow processes. Built, validated for real,
and pushed (`drtomc/gatk-sv-nf-filter-annotate-vcf:0.1.0`).

Validated for real end-to-end against a synthetic sites VCF + real
`AdjudicateSV` output (`scores`/`cutoffs` from the earlier synthetic
`evidence_metrics.tsv` run) — confirmed by hand that a `BND`-rescored
record was correctly reclassified, had `END2`/`CHR2` backfilled, had its
coordinates rewritten using its own `SR1POS`/`SR2POS` metrics (matching
the exact 0-based conversion the script performs), and was annotated with
the right `EVIDENCE` classes matching its own per-evidence-type scores.
One test-data-only gap found and fixed: the synthetic VCF's header
initially didn't declare `END2`/`CHR2` INFO fields, which the real
pipeline's earlier stages already carry — not a module bug.

### `GenotypeBatch`: genotyping itself

`genotype_batch` subworkflow + `workflows/genotype_batch.nf` entry point
— the actual goal of the whole prerequisite chain. Reproduces
`wdl/GenotypeBatch.wdl` in full: `TrainSVGenotyping` (once, whole sites
VCF, using `rf_cutoffs`' `PEQ`/`SRQ` values) → per contig (`PrintSVEvidence`
×3, subsetting panel-wide RD/PE/SR to that contig, then tabix-indexed →
`GenotypeSVs`, using the trained RD/PE/SR cutoff tables) → concat contig
shards back → `SeparateDepthPesr` (split by `INFO/ALGORITHMS`) →
`GenerateRegenoCoverageMedians` (`RD_MCR` extraction, vendored
`extract_format_table.py` — for `RegenotypeCNVs`, itself not built here,
kept only for output parity with upstream).

Five new modules: `GATK_TRAINSVGENOTYPING`, `GATK_GENOTYPESVS`, a new
`GATK_PRINT_SV_EVIDENCE_CONTIG` (distinct from the existing
`GATK_PRINT_SV_EVIDENCE`: that one merges *multiple per-sample* files
with `-F <list-file>`; this one subsets *one already-merged* file to a
single contig with a bare `--evidence-file`, matching
`GenotypeSVs`' own per-contig `PrintSVEvidence` calls exactly),
`SEPARATE_DEPTH_PESR`, `GENERATE_REGENO_COVERAGE_MEDIANS`. Same
GATK-SV-fork-only container as the three `GenerateBatchMetrics` walker
modules for `TrainSVGenotyping`/`GenotypeSVs` (confirmed directly, same as
before: neither exists in `gatk4-main_gcnvkernel`).

#### Two real Nextflow wiring bugs, caught by real (non-stub) validation

- **`.combine()` flattens a wrapped single-tuple value against a
  multi-element scattered tuple, not just a plain list.** The
  wrap-collect-unwrap pattern (`.map { x -> [[x]] }.collect().map {
  it[0] }`) that correctly broadcasts a *pair* (`[file, index]`) across a
  scatter elsewhere in this pipeline still gets re-flattened by
  `.combine()` when the scattered side itself already has multiple
  positional elements — found for real via an actual Groovy closure arity
  mismatch at runtime (`Invalid method invocation`), not assumed: a
  2-argument closure (`contig, pair ->`) needed to become a 3-argument one
  (`contig, f, tbi ->`) once the broadcast value's own `[f, tbi]` pair got
  flattened into the combined tuple's own positional slots. Same root
  lesson `vcfs_merge_read_counts` already documents for `.combine()`, but
  the exact flattening depth depends on both sides' shapes, not just the
  broadcast side's — verify the actual emitted shape with `.view()` each
  time, don't assume the earlier fix generalizes unchanged.
- **`GATK_PRINT_SV_EVIDENCE_CONTIG`'s output codec detection needs a
  filename prefix, and the caller must give each per-contig,
  per-evidence-type call a distinct one.** Same "process-name-derived
  output filename collides with input" pitfall documented earlier for
  `vcfs_combine_batches`/`vcfs_generate_batch_metrics`, plus a second,
  narrower gotcha specific to this module: GATK's own evidence-file codec
  detection requires the output filename to both have a real prefix
  segment (a bare `rd.txt.gz` fails "no suitable codecs found") and end in
  `.rd.txt.gz`/`.pe.txt.gz`/`.sr.txt.gz` specifically (case-insensitive).
  `conf/modules.config` sets `ext.prefix = { "${meta.id}.rd" }` (etc.) for
  each of the three `GATK_PRINT_SV_EVIDENCE_CONTIG` aliases accordingly.

#### An honest limitation: `TrainSVGenotyping`'s PE/SR evidence pass was never exercised to full completion with real data

`TrainSVGenotyping`'s RD/depth training pass succeeded for real
("Training on 1 CNV sites", "Training completed") against a genuine
synthetic VCF + median-coverage + RD evidence file, confirming the
module's own invocation (every flag, the `PEQ`/`SRQ` extraction from
`rf_cutoffs`, file staging) is correct. Its PE evidence pass, however,
consistently failed with `IllegalStateException: No discordant pair
counts after first pass` (`DiscordantPairEvidenceGenotyper
.finalizeFirstPass`) against every synthetic PE evidence file tried —
including a two-sample, four-discordant-pair file, and regardless of
whether the sites VCF had any PESR-typed (non-depth) variant at all. Not
traced to a specific root cause (unlike the `AdjudicateSV` FDR-cutoff
case, which was tracked down and understood) — likely some internal
minimum this walker wants that a hand-built, few-record synthetic PE file
doesn't reach, but this was **not confirmed**, and is left as a real,
open gap rather than papered over. Consequence: `GATK_GENOTYPESVS`'s own
invocation was validated with real data only up to the point it consumes
RD/PE/SR cutoff tables (a clean, precisely-located `rd_table.tsv`-format
error when fed an empty placeholder table) — never with a real,
successfully-trained table, since `TrainSVGenotyping` itself never fully
succeeded end-to-end with synthetic data. The full chain (all 116
processes, `TrainSVGenotyping` through `GenerateRegenoCoverageMedians`)
**is** confirmed correct via `-stub-run`, and every module's own
real-data invocation up to this specific evidence-volume wall is
independently confirmed correct — but a real, successful, non-stub
`TrainSVGenotyping`+`GenotypeSVs` pass together, with genuine PE/SR
evidence, remains unverified. **Test this for real against actual HPC
data** before trusting this path blindly, same caution this README
already gives every other real-vs-stub gap.

### Not yet started

Scramble (needs coverage counts + Manta's VCF as additional inputs, not
just the BAM — more involved than Wham), cn.MOPS (needs a local module, no
nf-core equivalent), the gCNV chain, `MergeBatchSites` (GATK-SV's own
site-list-only shortcut for cohorts that skip GATK-gCNV — not
implemented; `CombineBatches`, above, is), `SVAnnotate`, panel bundle I/O,
`RegenotypeCNVs` (depth-call re-genotyping QC; `GenotypeBatch`'s own
`GenerateRegenoCoverageMedians` output exists only for this, unconsumed
for now), and the full `build_panel.nf` / `genotype_new_sample.nf`
compositions. The whole genotyping prerequisite chain
(`GenerateBatchMetrics` → `FilterBatchSites`/`AdjudicateSV` →
`GenotypeBatch`) is now built — see [Path to
genotyping](#path-to-genotyping-the-real-prerequisite-chain) above for
the one open real-data validation gap (`TrainSVGenotyping`'s PE evidence
pass, see `GenotypeBatch`'s own section). `FilterBatchSites.wdl`'s
`PlotSVCountsPerSample` (outlier-sample QC plotting) remains deliberately
unbuilt — genuinely not needed by anything downstream.
