# Vendored GATK-SV scripts

Python scripts copied verbatim from GATK-SV's own `src/sv-pipeline/scripts/`,
vendored (not called from GATK-SV's `sv_pipeline_docker` image) so we
control the container these run in rather than depending on GATK-SV's
specific image build. See `nextflow/README.md` for the harmonisation design
these support.

| Script | Vendored from | Upstream commit |
|---|---|---|
| `bin/ploidy_table_from_ped.py` | `src/sv-pipeline/scripts/ploidy_table_from_ped.py` | `483973d6db3a8a9ba8fa72c5756a8ab53e7877f4` (2025-10-02) |
| `bin/format_svtk_vcf_for_gatk.py` | `src/sv-pipeline/scripts/format_svtk_vcf_for_gatk.py` | `483973d6db3a8a9ba8fa72c5756a8ab53e7877f4` (2025-10-02) |
| `bin/format_gatk_vcf_for_svtk.py` | `src/sv-pipeline/scripts/format_gatk_vcf_for_svtk.py` | `483973d6db3a8a9ba8fa72c5756a8ab53e7877f4` (2025-10-02) |
| `bin/medianCoverage.R` | `src/WGD/bin/medianCoverage.R` | `f091af0be6836446a9118f5db7cc3de9bde44fa6` (2025-09-02) |
| `bin/aggregate.py` | `src/sv-pipeline/02_evidence_assessment/02e_metric_aggregation/scripts/aggregate.py` | `2d5d8a44802548a184d7e57bcc967d2e611e8349` (2025-12-18) — **adapted, not verbatim**, see below |
| `bin/adjudicate_sv.py` | `src/svtk/svtk/adjudicate/adjudicate_sv.py` | `2d5d8a44802548a184d7e57bcc967d2e611e8349` (2025-12-18) — **adapted, not verbatim**, see below |
| `bin/random_forest.py` | `src/svtk/svtk/adjudicate/random_forest.py` | `9c6fbf562f750002253a37873d132db030657f7c` (2022-06-10) |
| `bin/labelers.py` | `src/svtk/svtk/adjudicate/labelers.py` | `9c6fbf562f750002253a37873d132db030657f7c` (2022-06-10) |
| `bin/adjudicate.py` | `src/svtk/svtk/cli/adjudicate.py` | `9f5620c1e3a6c633bfe3eae3d24f973cea8e8d1e` (2021-11-04) — **adapted, not verbatim**, see below |
| `bin/extract_format_table.py` | `src/sv-pipeline/scripts/extract_format_table.py` | `797b760468d9dba80610c51aeaa9b71506363ab6` (2026-02-06) |
| `bin/rewrite_SR_coords.py` | `src/sv-pipeline/03_variant_filtering/scripts/rewrite_SR_coords.py` | `2d5d8a44802548a184d7e57bcc967d2e611e8349` (2025-12-18) |
| `bin/annotate_RF_evidence.py` | `src/sv-pipeline/03_variant_filtering/scripts/annotate_RF_evidence.py` | `9f5620c1e3a6c633bfe3eae3d24f973cea8e8d1e` (2021-11-04) |

Each commit above is the file's own last-touching commit (`git log -1 --
-- <path>`), not necessarily a shared commit across every row — re-sync
per file, not by pinning one commit for the whole table.

Since these are copies, not symlinks or a submodule, they will silently
drift from upstream if GATK-SV fixes a bug in them. Re-sync by diffing
against the commit above and re-copying if GATK-SV's `main` moves past it
in a way that matters (check `git log <upstream commit>..main --
src/sv-pipeline/scripts/<script>` from the repo root).

## `medianCoverage.R`: first vendored R script

Only third-party dependency is the CRAN `optparse` package (base R
otherwise) — no bioconductor or compiled genomics library needed, unlike
GATK-SV's own `sv_pipeline_qc_docker` (the large, general-purpose image
this task runs in upstream). Given its own container
(`dockerfiles/median-coverage`, `r-base` + `optparse`), rather than
pulling that image just for this one script — same reasoning as the
Python scripts above having their own purpose-built container instead of
GATK-SV's `sv-pipeline-virtual-env`.

## `aggregate.py`: adapted, not verbatim

Upstream's `aggregate.py` imports `svtk.utils` for one function,
`get_called_samples` (~15 lines, pysam-only logic — walks a variant
record's genotypes and, for `CNV` records, non-diploid copy states).
Importing `svtk.utils` at all pulls in `pybedtools` as a module-level
import side effect, which in turn needs the `bedtools` binary — a real
container dependency for code this script itself never calls.
`get_called_samples`'s body (plus its one constant, `NULL_GT`) is inlined
directly into the vendored copy here instead, and the `import svtk.utils
as svu` line (and `svu.` prefix at its one call site) removed — everything
else is unchanged from upstream. This keeps the container to
`pysam`+`numpy`+`pandas`, consistent with this project's vendoring
rationale generally (control the container, don't replicate upstream's
image), rather than pulling in `svtk` as a whole dependency (with its own
`setup.py`, a Cython extension, and `pybedtools`/`bedtools`) for one small
function. **If re-syncing from upstream**, diff against
`src/sv-pipeline/02_evidence_assessment/02e_metric_aggregation/scripts/aggregate.py`
directly — the `svu.get_called_samples(variant)` call site will show as a
diff against this file's inlined `get_called_samples(variant)`; that's
expected, not drift to fix.

## `adjudicate.py`/`adjudicate_sv.py`: adapted, not verbatim (svtk adjudicate)

Wraps GATK-SV's `svtk adjudicate` (`wdl/FilterBatchSites.wdl`'s
`AdjudicateSV` task — the random-forest cutoff derivation genotyping needs,
see `nextflow/README.md`'s own "Path to genotyping" section). The actual
logic (`svtk/adjudicate/{adjudicate_sv,random_forest,labelers}.py`) has no
`pysam`/`pybedtools` dependency at all — only `pandas`/`numpy`/
`scikit-learn` — but invoking it via the real `svtk adjudicate` CLI would
still require installing the whole `svtk` package (its own `setup.py`,
Cython extension, `pybedtools`/`bedtools`), because `scripts/svtk`'s own
dispatcher (`import svtk.cli as cli`) eagerly imports every subcommand
module at startup, several of which do need those. Vendored the three
logic files verbatim (`adjudicate_sv.py`, `random_forest.py`,
`labelers.py`) plus a thin CLI wrapper (`adjudicate.py`, adapted from
`svtk/cli/adjudicate.py`'s own ~15 lines) instead — same rationale as
`aggregate.py`'s own adaptation, avoiding an unrelated dependency for code
that doesn't need it.

`adjudicate_sv.py`'s only adaptation: `from svtk.adjudicate import
rf_classify, labelers` became `from random_forest import rf_classify` +
`import labelers` (flat imports, since these are now sibling files in the
same directory rather than a real Python package) — everything else
unchanged. `random_forest.py`/`labelers.py` needed no changes at all (no
cross-module imports beyond stdlib/numpy/pandas/sklearn).
`adjudicate.py`'s own CLI wrapper mirrors `svtk/cli/adjudicate.py`
directly, with the same import adaptation.

**If re-syncing `adjudicate_sv.py` from upstream**, diff against
`src/svtk/svtk/adjudicate/adjudicate_sv.py` directly — the `from
svtk.adjudicate import rf_classify, labelers` line will show as a diff
against this file's flat imports; that's expected, not drift to fix.

## Known gap: pysam version

`format_svtk_vcf_for_gatk.py` and `format_gatk_vcf_for_svtk.py` both import
`pysam`. GATK-SV's own container pins `pysam==0.15.4` specifically (see
`dockerfiles/sv-pipeline-virtual-env/Dockerfile`), built from source,
because newer pysam/htslib reject VCF records with `END < POS` — which
occurs for `BND`/`CTX` (breakend/translocation) records. Our container
(`environment.yml` in this directory) uses a current pysam instead, for
build simplicity, rather than replicating that old from-source build.

`format_svtk_vcf_for_gatk.py`'s `_parse_bnd_ends()` already works around
part of this by manually text-parsing `BND`/`CTX` records' `END` field
directly from the VCF file rather than trusting pysam's own parsed value —
so the *value* used for those records should be correct regardless of
pysam version. What's untested is whether current pysam's VCF *iteration*
itself raises before that manual parsing code path is ever reached (i.e.
whether opening/iterating a VCF containing an "invalid" END fails outright
on a modern pysam/htslib, independent of what the script does with that
record afterward).

**Test this empirically** against real Manta/Wham output (which does
produce `BND` records) before trusting this path's `BND` handling. If
current pysam does fail on `BND`/`CTX` records, options are: pin an older
pysam here too, or patch the vendored copy to open the VCF in a way that
tolerates it (e.g. reading start-to-finish as text where necessary, as
`_parse_bnd_ends` already partly does).
