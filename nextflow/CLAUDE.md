# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

An independent, in-progress Nextflow reimplementation of the parent repo's GATK-SV pipeline (see `../CLAUDE.md` for the original WDL pipeline this is inspired by). Not a port — a simpler, panel-based reinterpretation for a specific situation: heterogeneous per-project cohorts (singleton batches of 10-30, or trios arriving one at a time), where the goal is a harmonized multi-caller SV VCF importable into seqR, not biobank-scale joint calling.

**`README.md` in this directory is the primary design/status doc** — read it before making non-trivial changes. It covers the two-tier `build-panel`/`genotype-new-sample` design, sample sheet and pedigree formats, static reference resources, deliberate simplifications relative to upstream GATK-SV, and — in its `Status` section — a detailed, still-growing log of every real bug found on actual HPC runs, with root cause and fix. That bug log is not just history: several of the pitfalls it documents (see below) are easy to reintroduce in new code written the same way.

## Commands

### Validate wiring (no containers, no real data)
```shell
nextflow run workflows/<name>.nf -profile stub -stub-run \
  --input tests/data/samplesheet.tsv --fasta tests/data/ref.fa --fasta_fai tests/data/ref.fa.fai \
  [... other required params, see the workflow's own header comment for its full example] \
  --outdir /tmp/results
```
Every `workflows/*.nf` entry point has a full example invocation in its own top-of-file comment — check that first rather than guessing params. Test fixtures live in `tests/data/` (tiny synthetic BAMs/CRAMs and placeholder static-resource files).

**`-stub-run` validates wiring only, not correctness.** It executes each process's `stub:` block instead of `script:`, so it never invokes a container, never runs real tool logic, and stub VCFs carry no real FORMAT/INFO fields — several real bugs (a vendored script never on `PATH`, a `docker://`-prefixed image breaking `-profile docker`, a missing VCF field only a real tool run would populate) shipped past `-stub-run` and were only caught on a real HPC run. Treat a clean `-stub-run` as necessary, not sufficient.

One recurring gotcha: `-stub-run` still executes any `eval(...)`-typed process output for real (e.g. Manta's module declares `eval("configManta.py --version")`), so a stub run can fail on `command not found` for a tool that's otherwise never invoked. Workaround when validating a pipeline that includes Manta:
```shell
mkdir -p /tmp/fake-bin
printf '#!/usr/bin/env bash\necho "fake version"\n' > /tmp/fake-bin/configManta.py
chmod +x /tmp/fake-bin/configManta.py
PATH="/tmp/fake-bin:$PATH" nextflow run workflows/call_manta.nf -stub-run ...
```

### Real local validation (no HPC needed)
For any process that shells out to a vendored script or declares its own `container`, do at least one real (non-stub) run before considering it validated — a local `-profile docker` run is usually enough:
```shell
nextflow run workflows/<name>.nf -profile docker --input ... --outdir /tmp/results
```

### Real HPC runs
```shell
nextflow run workflows/<name>.nf -profile apptainer -c conf/hpc.config --input ... --outdir results
```
`conf/hpc.config` is gitignored (site-specific Slurm queue/account/bind-paths + static resource paths); copy `conf/hpc.config.example` to create it. `scripts/fetch_static_resources.sh <target-dir>` downloads every GATK-SV static hg38 resource this pipeline needs and prints a ready-to-merge `params {}` config block — idempotent, safe to re-run.

### Adding an nf-core module
```shell
nf-core modules install <tool>/<subtool>
```
Then patch its `container` line if needed for plain `docker` (see pitfall below), and add any per-process resource/args overrides in `conf/modules.config` (never edit `nextflow.config` for per-module tuning).

## Architecture

- **Layout**: `modules/nf-core/` (unmodified vendored nf-core modules), `modules/local/` (hand-written — either no nf-core module exists, or a real GATK walker with no nf-core wrapper), `subworkflows/local/` (pipeline-specific composition), `workflows/` (standalone, independently-runnable entry points — no `main.nf` dispatcher; run each `workflows/*.nf` directly via `nextflow run`). Every subworkflow is composable: a later entry point (e.g. `combine_batches.nf`) runs earlier subworkflows (`bam_call_manta`, `vcfs_cluster_svcluster`) itself rather than assuming their output already exists on disk.
- **Two-stage harmonisation**: stage 1 (`vcfs_cluster_svcluster`, per caller) clusters one caller's calls across all samples into representative sites (mirrors GATK-SV's `ClusterPESR`). Stage 2 (`vcfs_combine_batches`) merges multiple callers' stage-1 output into one cross-caller site list (mirrors `CombineBatches`). Both use GATK's `SVCluster`/`GroupedSVCluster` walkers.
- **`ped_id` vs `sample_id`**: the sample sheet has both. `sample_id` (→ `meta.id`) is the lab/external reference name, used for filenames/display. `ped_id` is the pedigree file's `individual_id`, used for *every* PED-keyed value — ploidy table lookups, VCF sample-column rewrites, `--sample-name` on evidence-collection walkers. The two commonly differ in real data; using the wrong one produces real `KeyError`s downstream that `-stub-run` cannot catch (stub data has no real sample names to mismatch). When adding a process that consumes or writes a sample identifier, check whether it's PED-keyed and use `meta.ped_id`, not `meta.id`, if so.
- **svtk-format vs GATK-format VCFs**: callers (Manta, Wham) produce svtk-format VCFs; GATK's clustering/genotyping walkers need a different format (notably the `ECN` FORMAT field, added by `FORMAT_SVTK_VCF_FOR_GATK` and *stripped* by `FORMAT_GATK_VCF_FOR_SVTK`). A subworkflow that feeds one stage's *published* (svtk-format) output into another GATK-walker-based stage must re-run the svtk→GATK conversion first — the producer and consumer both being individually correct doesn't mean the round trip between them is lossless. This exact gap caused a real `Genotype missing required field ECN` failure; see the README's Status section for the full trace through GATK-SV's own WDL.
- **Vendored Python scripts** (`bin/*.py`, copied from GATK-SV's `src/sv-pipeline/scripts/`, see `bin/README.md` for provenance): always pass these as explicit `path` process inputs invoked via `python3 ${script}`, never rely on Nextflow's `bin/`-auto-PATH mechanism — that mechanism ties PATH to the directory of the script passed to `nextflow run`, which is `workflows/`, not `nextflow/` where `bin/` actually lives, so it silently never fires for any entry point here.
- **Container references must be conditional on engine**: `container "${workflow.containerEngine in ['singularity','apptainer'] && !task.ext.singularity_pull_docker_container ? 'https://.../data' : 'community.wave.seqera.io/...'}"` — the nf-core convention. A bare `docker://<image>` reference (Apptainer/Singularity-only syntax) breaks plain `docker` outright; `-stub-run` never catches this since it never invokes any container engine regardless of profile.
- **GATK walker containers are minimal**: the shared `gatk4-main_gcnvkernel` image (used by `GATK4_SVCLUSTER`, `GroupedSVCluster`, and every `modules/local/gatk4/*` module) has `gatk` and coreutils only — no `bgzip`, no `tabix`. Any new GATK-walker module whose output needs compressing or indexing needs a *separate* small chained module in a container that has those tools (`modules/local/tabix`, `modules/local/tabix_sv_evidence`, `modules/local/bgzip` all exist for this reason) — verify with `docker run --rm <image> which <tool>` before assuming a tool is present, don't assume from a WDL task's inline `tabix`/`bgzip` call that the equivalent tool is available in whatever container you're using here.
- **Nextflow gotchas hit repeatedly enough to check for explicitly**:
  - A process `output:` block using a glob (e.g. `path("*.foo.tsv")`) will silently match an *intermediate* file the script itself created, not just the intended final output, if the intermediate's name happens to match too — check every filename the `script:` block creates, not just the declared final name.
  - A process whose output filename is derived from `meta.id` (via `task.ext.prefix ?: "${meta.id}"`) will collide with itself across multiple stages of one subworkflow if `meta.id` isn't given a distinct value per stage — give each stage in a multi-step chain its own `meta.id` suffix.
  - `env` process outputs need a quoted string (`env 'NAME', emit: x`), not a bare identifier (`env NAME` fails at script compile time) — and the script needs a real `export NAME=...`, not a bare shell assignment.
  - `.combine()` flattens list-valued channel items against whatever they're paired with, rather than nesting them as one tuple element — `channel.of([1,2]).combine(channel.of('x'))` emits `[1,2,'x']`, not `[[1,2],'x']`. Wrap a list in an extra list first (`.map { l -> [l] }`) if you need it to stay distinguishable after a `combine`.
  - `nextflow.config`'s top-level `process {}` block sets `shell = ['/bin/bash', '-ueo', 'pipefail']` deliberately (not the Nextflow default `-ue`, no pipefail) — a killed first stage of any `cmd1 | cmd2 > out` pipeline would otherwise report the *second* command's exit code as the task's result, hiding e.g. an OOM kill as a silent success. Don't override this per-process without a specific reason.
