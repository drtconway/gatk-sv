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

Since these are copies, not symlinks or a submodule, they will silently
drift from upstream if GATK-SV fixes a bug in them. Re-sync by diffing
against the commit above and re-copying if GATK-SV's `main` moves past it
in a way that matters (check `git log <upstream commit>..main --
src/sv-pipeline/scripts/<script>` from the repo root).

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
