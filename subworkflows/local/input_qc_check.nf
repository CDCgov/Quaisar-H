//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_QC_CHECK } from '../../modules/local/samplesheet_qc_check'

workflow INPUT_QC_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_QC_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_qc_channels(it) }
        .set {ch_samples }

    // Create channels to match upstream processes
    // FASTP_TRIMD.out.reads -> tuple val(meta), path('*.trim.fastq.gz'), optional:true, emit: reads
    ch_samples.map{meta, reads, fastp_pass_json, fastp_failed_json, spades, mlst, quast, amrfinderplus -> 
        [ [id:meta.id, single_end:true], reads]
    }.set { ch_reads }

    // GATHERING_READ_QC_STATS: tuple val(meta), path(fastp_trimd_json), path(fastp_singles_json)
    ch_samples.map{meta, reads, fastp_pass_json, fastp_failed_json, spades, mlst, quast, amrfinderplus -> 
        [ [id:meta.id, single_end:true], fastp_pass_json, fastp_failed_json]
    }.set { ch_fastp_json }

    // SPADES_WF.out.spades_ch -> SPADES.out.scaffolds.map{meta, scaffolds -> [ [id:meta.id, single_end:true], scaffolds]}
    ch_samples.map{meta, reads, fastp_pass_json, fastp_failed_json, spades, mlst, quast, amrfinderplus -> 
        [ [id:meta.id, single_end:true], spades]
    }.set { ch_spades }

    // MLST.out.tsv -> tuple val(meta), path("*.tsv"), emit: tsv
    ch_samples.map{meta, reads, fastp_pass_json, fastp_failed_json, spades, mlst, quast, amrfinderplus -> 
        [ [id:meta.id, single_end:true], mlst]
    }.set { ch_mlst }

    // QUAST.out.report_tsv -> tuple val(meta), path("*.tsv"), emit: tsv
    ch_samples.map{meta, reads, fastp_pass_json, fastp_failed_json, spades, mlst, quast, amrfinderplus -> 
        [ [id:meta.id, single_end:true], quast]
    }.set { ch_quast }
    
    // AMRFINDERPLUS_RUN.out.report -> tuple val(meta), path("${meta.id}_all_genes.tsv"), emit: report
    ch_samples.map{meta, reads, fastp_pass_json, fastp_failed_json, spades, mlst, quast, amrfinderplus -> 
        [ [id:meta.id, single_end:true], amrfinderplus]
    }.set { ch_amrfinderplus }

    emit:
    reads             = ch_reads
    fastp_json        = ch_fastp_json
    spades            = ch_spades
    mlst              = ch_mlst
    quast             = ch_quast
    amrfinderplus     = ch_amrfinderplus
    valid_samplesheet = SAMPLESHEET_QC_CHECK.out.csv
    versions          = SAMPLESHEET_QC_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_qc_channels(LinkedHashMap row) {
    def meta = [:]
    meta.id           = row.sample
    meta.single_end   = row.single_end.toBoolean()
    missing_input = false

    def array = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }

    if (!meta.single_end) {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
    }

    // Check remaining files
    if (!file(row.fastp_pass_json).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Fastp passed reads JSON file does not exist!\n${row.fastp_pass_json}"
    }

    if (!file(row.fastp_failed_json).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Fastp failed reads JSON file does not exist!\n${row.fastp_failed_json}"
    }

    if (!file(row.spades).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> SPAdes assembly file does not exist!\n${row.spades}"
    }

    if (!file(row.mlst).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> MLST TSV report file does not exist!\n${row.mlst}"
    }

    if (!file(row.quast).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> QUAST TSV report file does not exist!\n${row.quast}"
    }

    if (!file(row.amrfinderplus).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> AMRFinder+ report file does not exist!\n${row.amrfinderplus}"
    }

    if (meta.single_end) {
        array = [ meta, [ file(row.fastq_1) ], file(row.fastp_pass_json), file(row.fastp_failed_json), file(row.spades), file(row.mlst), file(row.quast), file(row.amrfinderplus) ]
    } else {
        array = [ meta, [ file(row.fastq_1), file(row.fastq_2) ], file(row.fastp_pass_json), file(row.fastp_failed_json), file(row.spades), file(row.mlst), file(row.quast), file(row.amrfinderplus) ]
    }

    return array
}
