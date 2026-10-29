//
// Reduce a PED file to what sex-aware steps need: individual ID and sex,
// with family, parents and phenotype zeroed and rows sorted by ID.
//
// The table is built as a value at launch (pedSexes() below) and this
// task only writes it, so the task hash depends on the ID+sex content
// alone, not on the PED file's path/size/mtime. Edits that change no one's
// sex (parentage, phenotype, family IDs, row order) then leave this
// task, its output path and everything downstream cached (README: "PED
// edits no longer rerun the cohort"). Any step that reads the PED for sex
// only should take this output, not the PED.
//
// Usage: PED_SEXES(ped.map { meta, f -> [ meta, pedSexes(f) ] })
//

// Same rows the earlier awk version wrote (`!/^#/ && NF >= 5`, sorted by
// ID in C order), byte for byte, so its consumers' caches carry over.
def pedSexes(ped) {
    ped.readLines()
        .findAll { line -> !line.startsWith('#') && line.split('\t', -1).size() >= 5 }
        .collect { line -> def t = line.split('\t', -1); [ t[1], t[4] ] }
        .sort { a, b -> a[0] <=> b[0] }
        .collect { id, sex -> "0\t${id}\t0\t0\t${sex}\t0\n" }
        .join('')
}

process PED_SEXES {
    tag "${meta.id}"

    input:
    tuple val(meta), val(sexes)

    output:
    tuple val(meta), path("${task.ext.prefix ?: meta.id}.sexes.ped"), emit: ped

    when:
    task.ext.when == null || task.ext.when

    exec:
    def prefix = task.ext.prefix ?: "${meta.id}"
    task.workDir.resolve("${prefix}.sexes.ped").text = sexes
}
