#!/usr/bin/env bash
#
# First-pass QC of a genotyped SV VCF (e.g. genotype_batch's
# <cohort>.genotyped.pesr.vcf.gz): site counts, sizes and frequencies by
# SVTYPE; per-sample counts with upstream's IQR outlier rule; and, for
# every complete trio in the PED, Mendelian consistency, de novo rate and
# transmission rate on autosomes, per trio and then pooled by SVTYPE, size,
# AF, caller and RF evidence, plus a GQ-filter sweep (trio_strata.tsv).
# Usage:
#   scripts/qc_genotyped_vcf.sh <vcf> <ped> <outdir> [n_iqr]
#
# Needs bcftools and awk (mawk is fine) on PATH, e.g. in the pipeline's own
# bcftools image:
#   apptainer exec -B "$PWD" \
#     docker://community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a \
#     bash scripts/qc_genotyped_vcf.sh outputs/genotype_batch/cohort.genotyped.pesr.vcf.gz \
#       pedigree.ped qc/pesr
#
# Writes TSVs to <outdir> and prints summary.txt. Outliers use upstream's
# rule (src/sv-pipeline/03_variant_filtering/scripts/
# get_outliers_from_svcounts.helper.R): a sample fails if its count for any
# SVTYPE is outside [Q1 - n*IQR, Q3 + n*IQR]; upstream's default n is 6.
#
set -euo pipefail

if [[ $# -lt 3 ]]; then
    sed -n '9,11p' "$0" >&2
    exit 1
fi

vcf=$1
ped=$2
out=$3
n_iqr=${4:-6}
mkdir -p "$out"

# ---- Sites ----------------------------------------------------------------

bcftools query -f '%CHROM\t%POS\t%ID\t%INFO/SVTYPE\t%INFO/SVLEN\t%FILTER\t%INFO/AC\t%INFO/AN\n' "$vcf" \
    > "$out/sites.tsv"

# Per SVTYPE: sites, sites with any alt GT, singletons, median |SVLEN|.
sort -t $'\t' -k4,4 "$out/sites.tsv" | awk -F '\t' -v OFS='\t' '
    function flush() { if (t != "") print t, n, poly, single }
    BEGIN { print "svtype", "sites", "polymorphic", "singletons" }
    $4 != t { flush(); t = $4; n = poly = single = 0 }
    { n++; if ($7 != "." && $7 > 0) poly++; if ($7 == 1) single++ }
    END { flush() }
' > "$out/sites_by_svtype.tsv.tmp"

# Median |SVLEN| per SVTYPE, via a numeric sort within each type.
awk -F '\t' -v OFS='\t' '$5 != "." { l = $5 < 0 ? -$5 : $5; print $4, l }' "$out/sites.tsv" \
    | sort -t $'\t' -k1,1 -k2,2n \
    | awk -F '\t' -v OFS='\t' '
        function flush() { if (t != "") print t, (n % 2 ? v[(n + 1) / 2] : (v[n / 2] + v[n / 2 + 1]) / 2) }
        $1 != t { flush(); t = $1; n = 0; delete v }
        { v[++n] = $2 }
        END { flush() }
    ' > "$out/median_svlen.tmp"

awk -F '\t' -v OFS='\t' '
    FILENAME == ARGV[1] { med[$1] = $2; next }
    FNR == 1 { print $0, "median_svlen"; next }
    { print $0, ($1 in med ? med[$1] : ".") }
' "$out/median_svlen.tmp" "$out/sites_by_svtype.tsv.tmp" > "$out/sites_by_svtype.tsv"
rm -f "$out/sites_by_svtype.tsv.tmp" "$out/median_svlen.tmp"

# Size bins per SVTYPE (|SVLEN|; "." for types without one, e.g. BND).
awk -F '\t' -v OFS='\t' '
    BEGIN {
        nb = split("<100 100-500 500-1k 1k-5k 5k-50k 50k-1M >=1M", bins, " ")
        split("100 500 1000 5000 50000 1000000", edge, " ")
    }
    {
        t = $4; types[t] = 1
        if ($5 == ".") { c[t, "."]++; dot[t] = 1; next }
        l = $5 < 0 ? -$5 : $5
        for (i = 1; i < nb && l >= edge[i]; i++) ;
        c[t, bins[i]]++
    }
    END {
        printf "svtype"; for (i = 1; i <= nb; i++) printf "\t%s", bins[i]; print "\tno_svlen"
        for (t in types) {
            printf "%s", t
            for (i = 1; i <= nb; i++) printf "\t%d", c[t, bins[i]]
            printf "\t%d\n", c[t, "."]
        }
    }
' "$out/sites.tsv" | { read -r h; echo "$h"; sort; } > "$out/size_bins.tsv"

# Allele-frequency bins per SVTYPE (from AC/AN).
awk -F '\t' -v OFS='\t' '
    BEGIN { nb = split("AC=0 AC=1 AC=2-5 AF<1% AF=1-5% AF=5-50% AF>=50%", bins, " ") }
    {
        t = $4; types[t] = 1
        ac = $7; an = $8
        if (ac == "." || an == "." || an == 0) { c[t, "noAN"]++; next }
        af = ac / an
        b = ac == 0 ? 1 : ac == 1 ? 2 : ac <= 5 ? 3 : af < 0.01 ? 4 : af < 0.05 ? 5 : af < 0.5 ? 6 : 7
        c[t, bins[b]]++
    }
    END {
        printf "svtype"; for (i = 1; i <= nb; i++) printf "\t%s", bins[i]; print "\tno_AN"
        for (t in types) {
            printf "%s", t
            for (i = 1; i <= nb; i++) printf "\t%d", c[t, bins[i]]
            printf "\t%d\n", c[t, "noAN"]
        }
    }
' "$out/sites.tsv" | { read -r h; echo "$h"; sort; } > "$out/af_bins.tsv"

awk -F '\t' '{ c[$6]++ } END { print "filter\tsites"; for (f in c) print f "\t" c[f] }' "$out/sites.tsv" \
    | { read -r h; echo "$h"; sort -t $'\t' -k2,2nr; } > "$out/filters.tsv"

# ---- Per sample -----------------------------------------------------------

bcftools query -l "$vcf" > "$out/samples.txt"

# Non-ref GTs per sample per SVTYPE (long format, as upstream's svcounts),
# plus no-call GTs per sample. Samples with zero of a type get an explicit 0.
bcftools query -f '%INFO/SVTYPE[\t%SAMPLE=%GT]\n' "$vcf" | awk -F '\t' -v OFS='\t' -v samples="$out/samples.txt" '
    BEGIN { while ((getline s < samples) > 0) order[++ns] = s }
    {
        t = $1; types[t] = 1
        for (i = 2; i <= NF; i++) {
            eq = index($i, "="); s = substr($i, 1, eq - 1); gt = substr($i, eq + 1)
            if (gt ~ /\./) { miss[s]++; continue }
            if (gt ~ /[1-9]/) alt[s, t]++
        }
    }
    END {
        for (i = 1; i <= ns; i++) {
            s = order[i]
            for (t in types) print s, t, alt[s, t] + 0
            print s, "NOCALL", miss[s] + 0
        }
    }
' | { echo -e "sample\tsvtype\tcount"; sort -t $'\t' -k2,2 -k1,1; } > "$out/per_sample_counts.tsv"

# Upstream IQR rule, per SVTYPE (NOCALL excluded), with R's default (type 7)
# quantiles.
tail -n +2 "$out/per_sample_counts.tsv" | awk -F '\t' '$2 != "NOCALL"' \
    | sort -t $'\t' -k2,2 -k3,3n \
    | awk -F '\t' -v OFS='\t' -v k="$n_iqr" '
        function q(p,   h, lo) { h = (n - 1) * p + 1; lo = int(h); return lo >= n ? v[n] : v[lo] + (h - lo) * (v[lo + 1] - v[lo]) }
        function flush(   i, q1, q3, lo, hi) {
            if (t == "") return
            q1 = q(0.25); q3 = q(0.75); lo = q1 - k * (q3 - q1); hi = q3 + k * (q3 - q1)
            print t, n, q(0.5), q1, q3, lo, hi > cut
            # A class most samples lack (q3 = 0, e.g. CTX) has a zero-width
            # cutoff, so any carrier would be flagged: skip it.
            if (q3 == 0) return
            for (i = 1; i <= n; i++) if (v[i] < lo || v[i] > hi) print s[i], t, v[i], lo, hi
        }
        BEGIN { cut = "/dev/stderr"; print "sample", "svtype", "count", "low_cutoff", "high_cutoff" }
        $2 != t { flush(); t = $2; n = 0; delete v; delete s }
        { n++; s[n] = $1; v[n] = $3 }
        END { flush() }
    ' 2> "$out/per_sample_cutoffs.tsv.tmp" > "$out/outliers.tsv"
{ echo -e "svtype\tsamples\tmedian\tq1\tq3\tlow_cutoff\thigh_cutoff"; cat "$out/per_sample_cutoffs.tsv.tmp"; } \
    > "$out/per_sample_cutoffs.tsv"
rm -f "$out/per_sample_cutoffs.tsv.tmp"

# ---- Trios ----------------------------------------------------------------

# Complete trios: PED rows whose child, father and mother are all in the VCF.
awk -F '[\t ]+' -v samples="$out/samples.txt" '
    BEGIN { while ((getline s < samples) > 0) in_vcf[s] = 1 }
    /^#/ || NF < 4 { next }
    ($2 in in_vcf) && ($3 in in_vcf) && ($4 in in_vcf) { print $2 "\t" $3 "\t" $4 }
' "$ped" > "$out/trios.txt"

# Shared by both trio passes. Autosomes only: chrX/chrY need sex-aware
# rules. GTs are reduced to alt-allele dosage (0/1/2); anything else
# (no-call, haploid, multiallelic CN) counts as missing.
trio_awk_lib='
    function dose(g) {
        if (g == "0/0" || g == "0|0") return 0
        if (g == "0/1" || g == "1/0" || g == "0|1" || g == "1|0") return 1
        if (g == "1/1" || g == "1|1") return 2
        return -1
    }
    # Can child dosage c come from one allele of f and one of m?
    function ok(c, f, m,   a, b) {
        for (a = (f == 2); a <= (f > 0); a++)
            for (b = (m == 2); b <= (m > 0); b++)
                if (a + b == c) return 1
        return 0
    }
    function read_trios(file,   line, p) {
        while ((getline line < file) > 0) { split(line, p, "\t"); nt++; C[nt] = p[1]; F[nt] = p[2]; M[nt] = p[3] }
    }
    function pct(a, b) { return b ? sprintf("%.1f", 100 * a / b) : "." }
'

if [[ -s "$out/trios.txt" ]]; then
    members=$(tr '\t' '\n' < "$out/trios.txt" | sort -u | paste -sd, -)
    bcftools query -s "$members" -f '%CHROM\t%INFO/SVTYPE[\t%SAMPLE=%GT]\n' "$vcf" \
        | awk -F '\t' -v OFS='\t' -v trios="$out/trios.txt" "$trio_awk_lib"'
            function add(k, t, field) { cnt[k, t, field]++; cnt[k, "ALL", field]++; types[t] = 1 }
            BEGIN { read_trios(trios) }
            $1 ~ /^(chr)?[XYM]/ { next }
            {
                t = $2
                delete g
                for (i = 3; i <= NF; i++) { eq = index($i, "="); g[substr($i, 1, eq - 1)] = substr($i, eq + 1) }
                for (k = 1; k <= nt; k++) {
                    c = dose(g[C[k]]); f = dose(g[F[k]]); m = dose(g[M[k]])
                    gc = g[C[k]]; gf = g[F[k]]; gm = g[M[k]]
                    anyalt = (gc ~ /[1-9]/) || (gf ~ /[1-9]/) || (gm ~ /[1-9]/)
                    if (!anyalt) continue
                    if (c < 0 || f < 0 || m < 0) { add(k, t, "missing"); continue }
                    add(k, t, "informative")
                    if (!ok(c, f, m)) add(k, t, "inconsistent")
                    if (c > 0 && f == 0 && m == 0) add(k, t, "de_novo")
                    # One het parent, other hom-ref: expect ~50% transmission.
                    if ((f == 1 && m == 0) || (f == 0 && m == 1)) {
                        add(k, t, "het_parent")
                        if (c > 0) add(k, t, "transmitted")
                    }
                }
            }
            END {
                print "child", "father", "mother", "svtype", "informative", "inconsistent", "pct_inconsistent",
                      "de_novo", "pct_de_novo", "missing", "het_parent", "transmitted", "pct_transmitted"
                types["ALL"] = 1
                for (k = 1; k <= nt; k++)
                    for (t in types) {
                        inf = cnt[k, t, "informative"] + 0; hp = cnt[k, t, "het_parent"] + 0
                        printf "%s\t%s\t%s\t%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%d\t%s\n", C[k], F[k], M[k], t,
                            inf, cnt[k, t, "inconsistent"], pct(cnt[k, t, "inconsistent"], inf),
                            cnt[k, t, "de_novo"], pct(cnt[k, t, "de_novo"], inf),
                            cnt[k, t, "missing"], hp, cnt[k, t, "transmitted"], pct(cnt[k, t, "transmitted"], hp)
                    }
            }
        ' | { read -r h; echo "$h"; sort -t $'\t' -k1,1 -k4,4; } > "$out/trios.tsv"
else
    echo -e "child\tfather\tmother\tsvtype" > "$out/trios.tsv"
fi

# ---- Trio error strata ----------------------------------------------------

# The same counts pooled over all trios, split by site properties (SVTYPE,
# size, AF, caller, RF evidence), then a GQ sweep: at each gq_min, any trio
# member's GT with GQ < gq_min is treated as a no-call (as a per-genotype
# GQ filter would do). child_alt is the number of child non-ref calls left,
# a proxy for what each filter costs in sensitivity.
hdr=$(bcftools view -h "$vcf")
info_or_dot() { if grep -q "^##INFO=<ID=$1," <<< "$hdr"; then echo "%INFO/$1"; else echo "."; fi }
has_gq=0; grep -q '^##FORMAT=<ID=GQ,' <<< "$hdr" && has_gq=1
gq_thresholds="1 10 20 30 50 100 200 300 500 700 900 999"

if [[ -s "$out/trios.txt" ]]; then
    fmt="%CHROM\t%INFO/SVTYPE\t%INFO/SVLEN\t$(info_or_dot AC)\t$(info_or_dot AN)\t$(info_or_dot ALGORITHMS)\t$(info_or_dot EVIDENCE)"
    fmt+="[\t%SAMPLE=%GT:$([[ $has_gq == 1 ]] && echo '%GQ' || echo '.')]\n"
    bcftools query -s "$members" -f "$fmt" "$vcf" \
        | awk -F '\t' -v OFS='\t' -v trios="$out/trios.txt" -v has_gq="$has_gq" -v thresholds="$gq_thresholds" "$trio_awk_lib"'
            function sizebin(l) {
                if (l == ".") return "no_svlen"
                if (l < 0) l = -l
                return l < 100 ? "<100" : l < 500 ? "100-500" : l < 1000 ? "500-1k" : l < 5000 ? "1k-5k" : \
                       l < 50000 ? "5k-50k" : l < 1000000 ? "50k-1M" : ">=1M"
            }
            function afbin(ac, an,   af) {
                if (ac == "." || an == "." || an == 0) return "no_AN"
                af = ac / an
                return ac <= 1 ? "AC<=1" : ac <= 5 ? "AC=2-5" : af < 0.01 ? "AF<1%" : af < 0.05 ? "AF=1-5%" : \
                       af < 0.5 ? "AF=5-50%" : "AF>=50%"
            }
            function tally(key, c, f, m) {
                if (!(key in seen)) { seen[key] = 1; order[++nk] = key }
                if (c < 0 || f < 0 || m < 0) { n[key, "missing"]++; return }
                n[key, "informative"]++
                if (c > 0) n[key, "child_alt"]++
                if (!ok(c, f, m)) n[key, "inconsistent"]++
                if (c > 0 && f == 0 && m == 0) n[key, "de_novo"]++
                if ((f == 1 && m == 0) || (f == 0 && m == 1)) {
                    n[key, "het_parent"]++
                    if (c > 0) n[key, "transmitted"]++
                }
            }
            function gtgq(s, which,   v, colon) {
                v = g[s]; colon = index(v, ":")
                return which == "gt" ? substr(v, 1, colon - 1) : substr(v, colon + 1)
            }
            BEGIN { read_trios(trios); nth = split(thresholds, th, " ") }
            $1 ~ /^(chr)?[XYM]/ { next }
            {
                t = $2; sz = sizebin($3); af = afbin($4, $5); alg = $6; ev = $7
                delete g
                for (i = 8; i <= NF; i++) { eq = index($i, "="); g[substr($i, 1, eq - 1)] = substr($i, eq + 1) }
                for (k = 1; k <= nt; k++) {
                    gc = gtgq(C[k], "gt"); gf = gtgq(F[k], "gt"); gm = gtgq(M[k], "gt")
                    if (!((gc ~ /[1-9]/) || (gf ~ /[1-9]/) || (gm ~ /[1-9]/))) continue
                    c = dose(gc); f = dose(gf); m = dose(gm)
                    tally("svtype\t" t "\t0", c, f, m)
                    tally("svtype\tALL\t0", c, f, m)
                    tally("size\t" t " " sz "\t0", c, f, m)
                    tally("af\t" af "\t0", c, f, m)
                    tally("algorithms\t" alg "\t0", c, f, m)
                    tally("svtype_algorithms\t" t " " alg "\t0", c, f, m)
                    tally("evidence\t" ev "\t0", c, f, m)
                    if (!has_gq) continue
                    # GQ "." fails every threshold.
                    qc = gtgq(C[k], "gq"); qf = gtgq(F[k], "gq"); qm = gtgq(M[k], "gq")
                    for (j = 1; j <= nth; j++) {
                        T = th[j]
                        cc = (qc == "." || qc + 0 < T) ? -1 : c
                        ff = (qf == "." || qf + 0 < T) ? -1 : f
                        mm = (qm == "." || qm + 0 < T) ? -1 : m
                        tally("gq_min\t" t "\t" T, cc, ff, mm)
                        tally("gq_min\tALL\t" T, cc, ff, mm)
                    }
                }
            }
            END {
                print "dimension", "stratum", "gq_min", "informative", "child_alt", "inconsistent", "pct_inconsistent",
                      "de_novo", "pct_de_novo", "missing", "het_parent", "transmitted", "pct_transmitted"
                for (i = 1; i <= nk; i++) {
                    key = order[i]; inf = n[key, "informative"] + 0; hp = n[key, "het_parent"] + 0
                    print key, inf, n[key, "child_alt"] + 0, n[key, "inconsistent"] + 0, pct(n[key, "inconsistent"], inf),
                          n[key, "de_novo"] + 0, pct(n[key, "de_novo"], inf), n[key, "missing"] + 0,
                          hp, n[key, "transmitted"] + 0, pct(n[key, "transmitted"], hp)
                }
            }
        ' | { read -r h; echo "$h"; sort -t $'\t' -k1,1 -k2,2 -k3,3n; } > "$out/trio_strata.tsv"
else
    echo -e "dimension\tstratum\tgq_min" > "$out/trio_strata.tsv"
fi

# ---- Summary --------------------------------------------------------------

{
    echo "QC of $vcf"
    echo "$(wc -l < "$out/samples.txt" | tr -d ' ') samples, $(wc -l < "$out/sites.tsv" | tr -d ' ') sites, $(wc -l < "$out/trios.txt" | tr -d ' ') complete trios"
    echo
    echo "== Sites by SVTYPE";            cat "$out/sites_by_svtype.tsv"; echo
    echo "== Sizes (|SVLEN|)";            cat "$out/size_bins.tsv"; echo
    echo "== Allele frequency";           cat "$out/af_bins.tsv"; echo
    echo "== FILTER";                     cat "$out/filters.tsv"; echo
    echo "== Per-sample non-ref counts (cutoffs at ${n_iqr}x IQR)"; cat "$out/per_sample_cutoffs.tsv"; echo
    echo "== Outlier samples ($(($(wc -l < "$out/outliers.tsv") - 1)) sample/SVTYPE pairs)"
    head -n 41 "$out/outliers.tsv"
    [[ $(wc -l < "$out/outliers.tsv") -gt 41 ]] && echo "... (see outliers.tsv)"
    echo
    echo "== Trios, autosomes, all SVTYPEs (per-SVTYPE rows in trios.tsv)"
    awk -F '\t' 'NR == 1 || $4 == "ALL"' "$out/trios.tsv"
    echo
    echo "== Trio errors pooled over all trios, by stratum (trio_strata.tsv)"
    for dim in svtype size af algorithms svtype_algorithms evidence; do
        echo "-- $dim"
        awk -F '\t' -v d="$dim" 'NR == 1 || ($1 == d && $3 == 0)' "$out/trio_strata.tsv" | cut -f 2,4-
    done
    echo "-- GQ sweep (GTs with GQ < gq_min treated as no-calls)"
    [[ $has_gq == 1 ]] || echo "(no FORMAT/GQ in this VCF)"
    awk -F '\t' 'NR == 1 || $1 == "gq_min" || ($1 == "svtype" && $3 == 0)' "$out/trio_strata.tsv" \
        | cut -f 2- | { read -r h; echo "$h"; sort -t $'\t' -k1,1 -k2,2n; }
} | tee "$out/summary.txt"
