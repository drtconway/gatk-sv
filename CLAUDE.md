# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

GATK-SV is a cloud-native structural variant (SV) discovery pipeline for Illumina short-read WGS data, written in WDL and orchestrated via Cromwell (run through Cromshell or Terra). It joint-genotypes SVs across multiple callers, harmonizing their output into one VCF.

This repository also contains an independent, in-progress **Nextflow reimplementation** under `nextflow/` — a simpler, panel-based reinterpretation of this pipeline for a specific use case (heterogeneous per-project cohorts, incremental sample arrival). See `nextflow/CLAUDE.md` for that subtree; it has its own build/run/test conventions, unrelated to the ones below.

## Repository structure

- `wdl/` — WDLs implementing the pipeline. One module per WDL file (e.g. `ClusterBatch.wdl`, `GatherBatchEvidence.wdl`, `FilterBatch.wdl`), composed into top-level entry points: `GATKSVPipelineBatch.wdl` (cohort/batch mode), `GATKSVPipelinePhase1.wdl`, `GATKSVPipelineSingleSample.wdl`.
- `src/` — pipeline scripts and packages invoked from WDL tasks, run inside Docker containers:
  - `sv-pipeline/` — the bulk of the Python/R scripts used throughout the pipeline
  - `svtk/` — Python toolkit for SV VCF standardization/parsing
  - `sv_utils/` — Python package with its own `setup.py` and `tests/` (pytest)
  - `svqc/`, `svtest/` — metric QC and summary-metric generation
  - `RdTest/`, `WGD/`, `denovo/`, `str/`, `stripy/` — task-specific script groups
- `inputs/` — `templates/` (input JSON templates) + `values/` (values used to populate them, e.g. `resources_hg38.json` for static reference resources, `dockers.json` for pinned image tags).
- `dockerfiles/` — Dockerfile sources for pipeline images; see `scripts/docker/README.md` for build instructions.
- `scripts/` — `inputs/` (generate input JSONs from templates), `test/` (WDL validation), `docker/`, `notebooks/`.
- `website/` — pipeline documentation site (published at https://broadinstitute.github.io/gatk-sv/); that site is the canonical source for how to actually *run* the pipeline.

## Commands

### Python lint
```shell
tox -e lint       # flake8, per tox.ini's ignore/exclude list
```

### WDL validation
CI (`testwdls.yaml`) runs on any push/PR touching `wdl/**`, `inputs/**`, or the scripts below:
```shell
scripts/inputs/build_default_inputs.sh -d .                       # generate input JSONs from templates
python scripts/test/miniwdl_validation.py --imports-dir wdl wdl/*.wdl   # syntax check via miniwdl
scripts/test/validate.sh -t -d . -j womtool.jar                    # validate JSON inputs against WDL, via womtool
```
`validate.sh -t` also validates Terra cohort-mode input JSONs, not just test inputs. `womtool.jar` is downloaded from the [Cromwell releases page](https://github.com/broadinstitute/cromwell/releases).

### Python package tests
`src/sv_utils` is a proper Python package (`setup.py`) with its own `tests/` directory — run with pytest from within `src/sv_utils`.

## Architecture notes

- **WDL module composition**: the top-level pipelines (`GATKSVPipelineBatch.wdl` etc.) import and chain the per-module WDLs in `wdl/`. To understand what a stage actually does, read the module WDL directly (e.g. `ClusterBatch.wdl` for cross-sample clustering, `GatherBatchEvidence.wdl` for PE/SR/BAF/RD evidence collection, `FilterBatch.wdl` for genotyping-cutoff derivation) rather than inferring from the top-level pipeline alone — the top-level files are mostly wiring.
- **Docker-per-task**: nearly every WDL task specifies its own `docker` image (pinned in `inputs/values/dockers.json`), and most tasks shell out to a script in `src/` rather than embedding logic inline in WDL. When changing task behavior, the actual logic is almost always in `src/sv-pipeline/scripts/` or one of the other `src/` packages, not in the `.wdl` file's `command <<<...>>>` block itself (which is often just argument marshalling).
- **Resource pinning via `inputs/values/`**: static reference resources (hg38 contigs, exclude lists, clustering configs) and Docker image tags are centralized in `inputs/values/resources_hg38.json` and `dockers.json`, not hardcoded per-WDL. `scripts/inputs/build_default_inputs.sh` renders `inputs/templates/` against these values to produce runnable input JSONs.
- **`RuntimeAttr` sizing pattern**: WDL tasks generally define a `RuntimeAttr default_attr` object plus an optional `RuntimeAttr? runtime_attr_override` input, then `select_first([runtime_attr_override, default_attr])`. When reading a task's actual memory/CPU/disk sizing, check whether the field you care about (e.g. `mem_gb`) is genuinely threaded into the `runtime { ... }` block's calculation — some tasks compute their real resource values from a *separate* top-level input with its own independent default, leaving the `default_attr` field looking authoritative but actually unused (seen firsthand in `CollectCoverage.wdl`'s `CollectCounts` task, where `default_attr.mem_gb: 3.0` is dead and the real default is `12.0` from a different input).
