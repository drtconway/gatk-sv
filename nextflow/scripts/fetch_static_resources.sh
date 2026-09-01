#!/usr/bin/env bash
# Fetch GATK-SV's static, sample-independent hg38 resource files (see
# README.md#static-reference-resources) into a target directory, and print
# a Nextflow config snippet on stdout that sets the corresponding params to
# the downloaded paths.
#
# Usage:
#   scripts/fetch_static_resources.sh <target-dir> > conf/hpc.config
#
# Idempotent: skips any file that already exists in <target-dir> (by name),
# so re-running after adding this pipeline's reference genome/fasta.fai to
# the same conf/hpc.config, or after a partial/interrupted fetch, only
# downloads what's missing. Doesn't fetch the reference genome itself
# (fasta/fasta_fai) -- those aren't GATK-SV-specific and you likely already
# have a copy on your cluster; set them yourself in the printed config, or
# merge this script's output into an existing conf/hpc.config.
#
# All URLs are public GCS buckets (no GCP credentials needed) -- same files
# GATK-SV's own WDL pipeline uses, from inputs/values/resources_hg38.json.

set -ueo pipefail

if [ $# -ne 1 ]; then
    echo "Usage: $0 <target-dir>" >&2
    exit 1
fi

target_dir="$1"
mkdir -p "$target_dir"
target_dir="$(cd "$target_dir" && pwd)"

# name -> URL. One resource per line: <param-name-fragment> <url>
# (param-name-fragment is used to derive both the local filename and the
# params.* key -- see the emit loop below.)
resources="
primary_contigs_list https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/sv-resources/resources/v1/primary_contigs.list
primary_contigs_fai https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/sv-resources/resources/v1/contig.fai
reference_dict https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict
pesr_exclude_intervals https://storage.googleapis.com/gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/PESR.encode.peri_all.repeats.delly.hg38.blacklist.sorted.bed.gz
pesr_exclude_intervals_tbi https://storage.googleapis.com/gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/PESR.encode.peri_all.repeats.delly.hg38.blacklist.sorted.bed.gz.tbi
manta_region_bed https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/sv-resources/resources/v1/primary_contigs_plus_mito.bed.gz
manta_region_bed_tbi https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/sv-resources/resources/v1/primary_contigs_plus_mito.bed.gz.tbi
clustering_config_part1 https://storage.googleapis.com/gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/clustering_config.part_one.tsv
clustering_config_part2 https://storage.googleapis.com/gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/clustering_config.part_two.tsv
stratification_config_part1 https://storage.googleapis.com/gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/stratify_config.part_one.tsv
stratification_config_part2 https://storage.googleapis.com/gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/stratify_config.part_two.tsv
clustering_track_sr https://storage.googleapis.com/gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/hg38.SimpRep.sorted.pad_100.merged.bed
clustering_track_sd https://storage.googleapis.com/gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/hg38.SegDup.sorted.merged.bed
clustering_track_rm https://storage.googleapis.com/gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/hg38.RM.sorted.merged.bed
sd_locs_vcf https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
sd_locs_vcf_idx https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx
preprocessed_intervals https://storage.googleapis.com/gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/preprocessed_intervals.interval_list
"

declare -a param_lines=()

while read -r name url; do
    [ -z "$name" ] && continue
    filename="$(basename "$url")"
    dest="$target_dir/$filename"

    if [ -f "$dest" ]; then
        echo "==> $name: already present ($dest), skipping" >&2
    else
        echo "==> $name: fetching $url" >&2
        curl -fsSL "$url" -o "$dest.tmp"
        mv "$dest.tmp" "$dest"
    fi

    param_lines+=("    $name = '$dest'")
done <<< "$resources"

echo "==> Done. hg38 static resources are in $target_dir" >&2
echo "==> min_svsize (plain integer, no file) already defaults to 50 in nextflow.config -- no param needed here unless overriding." >&2

cat <<EOF
// Static, sample-independent hg38 resource files fetched by
// scripts/fetch_static_resources.sh -- see
// nextflow/README.md#static-reference-resources. Regenerate by re-running
// that script against the same target directory (idempotent: only
// downloads files that aren't already there).
//
// Not included here: fasta/fasta_fai (your own reference genome copy, not
// GATK-SV-specific) -- set those yourself alongside this block.
params {
$(IFS=$'\n'; echo "${param_lines[*]}")
}
EOF
