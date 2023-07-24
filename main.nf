#!/usr/bin/env nextflow
/*
========================================================================================
    CDCgov/phoenix
========================================================================================
    Github : https://github.com/CDCgov/phoenix
    Slack  : https://staph-b-dev.slack.com/channels/phoenix-dev
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

WorkflowMain.initialise(workflow, params, log)

//Check coverage is above its threshold
if (params.coverage < 30) { exit 1, 'The minimum coverage allowed for QA/QC purposes is 30 and is the default. Please choose a value >=30.' }
//Check path of kraken2db
if (params.kraken2db == null) { exit 1, 'Input path to kraken2db not specified!' }

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { PHOENIX_EXTERNAL       } from './workflows/phoenix'
include { PHOENIX_QC_EXTERNAL    } from './workflows/phoenix_qc'
include { PHOENIX_EXQC           } from './workflows/cdc_phoenix'
include { SCAFFOLDS_EXTERNAL     } from './workflows/scaffolds'
include { SCAFFOLDS_EXQC         } from './workflows/cdc_scaffolds'
include { SRA_PREP               } from './workflows/sra_prep'

//
// WORKFLOW: Run main cdcgov/phoenix analysis pipeline
//
workflow PHOENIX {
    // Validate input parameters
    // Check input path parameters to see if they exist
    def checkPathParamList = [ params.input, params.multiqc_config, params.kraken2db] //removed , params.fasta to stop issue w/connecting to aws and igenomes not used
    for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

    // Check mandatory parameters

    //input on command line
    if (params.input) { ch_input = file(params.input) } else { exit 1, 'For -entry PHOENIX: Input samplesheet not specified!' }
    ch_versions = Channel.empty() // Used to collect the software versions

    main:
        PHOENIX_EXTERNAL ( ch_input, ch_versions )
    emit:
        scaffolds        = PHOENIX_EXTERNAL.out.scaffolds
        trimmed_reads    = PHOENIX_EXTERNAL.out.trimmed_reads
        mlst             = PHOENIX_EXTERNAL.out.mlst
        amrfinder_report = PHOENIX_EXTERNAL.out.amrfinder_report
        gamma_ar         = PHOENIX_EXTERNAL.out.gamma_ar
        summary_report   = PHOENIX_EXTERNAL.out.summary_report
}

//
// WORKFLOW: Run QC steps from main cdcgov/phoenix analysis pipeline
//
workflow PHOENIX_QC {
    main:
        PHOENIX_QC_EXTERNAL ()
}

//
// WORKFLOW: Run internal version of cdcgov/phoenix analysis pipeline that includes BUSCO, SRST2 and KRAKEN_ASMBLED
//
workflow CDC_PHOENIX {
    // Validate input parameters
    // Check input path parameters to see if they exist
    def checkPathParamList = [ params.input, params.multiqc_config, params.kraken2db]
    for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

    // Check mandatory parameters

    //input on command line
    if (params.input) { ch_input = file(params.input) } else { exit 1, 'For -entry CDC_PHOENIX: Input samplesheet not specified!' }
    ch_versions = Channel.empty() // Used to collect the software versions

    main:
        PHOENIX_EXQC ( ch_input, ch_versions )

    emit:
        scaffolds        = PHOENIX_EXQC.out.scaffolds
        trimmed_reads    = PHOENIX_EXQC.out.trimmed_reads
        mlst             = PHOENIX_EXQC.out.mlst
        amrfinder_report = PHOENIX_EXQC.out.amrfinder_report
        gamma_ar         = PHOENIX_EXQC.out.gamma_ar
        summary_report   = PHOENIX_EXQC.out.summary_report
}

/*
========================================================================================
    RUN SRA WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Run internal version of phx based on sample SRA names, these will be pulled from NCBI for you. 
//
workflow SRA {
    // Validate input parameters
    // Check input path parameters to see if they exist
    def checkPathParamList = [ params.input_sra, params.multiqc_config, params.kraken2db ]
    for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

    // Check mandatory parameters

    //input on command line
    if (params.input_sra) { ch_input = file(params.input_sra) } else { exit 1, 'For -entry SRA: Input samplesheet not specified! Make sure to use --input_sra NOT --input' }

    main:
        // pull data and create samplesheet for it.
        SRA_PREP ( ch_input )
        // pass samplesheet to PHOENIX
        PHOENIX_EXTERNAL ( SRA_PREP.out.samplesheet, SRA_PREP.out.versions )

    emit:
        scaffolds        = PHOENIX_EXTERNAL.out.scaffolds
        trimmed_reads    = PHOENIX_EXTERNAL.out.trimmed_reads
        mlst             = PHOENIX_EXTERNAL.out.mlst
        amrfinder_report = PHOENIX_EXTERNAL.out.amrfinder_report
        gamma_ar         = PHOENIX_EXTERNAL.out.gamma_ar
        summary_report   = PHOENIX_EXTERNAL.out.summary_report
}

//
// WORKFLOW: Run cdc version of phx based on sample SRA names, the fastq will be pulled from NCBI for you. 
//

workflow CDC_SRA {
    // Validate input parameters
    // Check input path parameters to see if they exist
    def checkPathParamList = [ params.input_sra, params.multiqc_config, params.kraken2db]
    for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

    // Check mandatory parameters

    //input on command line
    if (params.input_sra) { ch_input = file(params.input_sra) } else { exit 1, 'For -entry CDC_SRA: Input samplesheet not specified! Make sure to use --input_sra NOT --input' }

    main:
        // pull data and create samplesheet for it.
        SRA_PREP ( ch_input )
        // pass samplesheet to PHOENIX
        PHOENIX_EXQC ( SRA_PREP.out.samplesheet, SRA_PREP.out.versions )

    emit:
        scaffolds        = PHOENIX_EXQC.out.scaffolds
        trimmed_reads    = PHOENIX_EXQC.out.trimmed_reads
        mlst             = PHOENIX_EXQC.out.mlst
        amrfinder_report = PHOENIX_EXQC.out.amrfinder_report
        gamma_ar         = PHOENIX_EXQC.out.gamma_ar
        summary_report   = PHOENIX_EXQC.out.summary_report
}

/*
========================================================================================
    RUN SCAFFOLD WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Entry point to analyze scaffold file(s) and run everything after Spades
//
workflow SCAFFOLDS {
    // Validate input parameters
    // Check input path parameters to see if they exist
    if (params.input != null ) {  // if a samplesheet is passed
        if (params.indir != null ) { //if samplesheet is passed and an input directory exit
            exit 1, 'For -entry SCAFFOLDS: You need EITHER an input samplesheet or a directory! Just pick one.' 
        } else { // if only samplesheet is passed check to make sure input is an actual file
            def checkPathParamList = [ params.input, params.multiqc_config, params.kraken2db ]
            for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }
            ch_input_indir = null //keep input directory null if not passed
            // get full path for input and make channel
            if (params.input) { ch_input = file(params.input) }
        }
    } else {
        if (params.indir != null ) { // if no samplesheet is passed, but an input directory is given
            ch_input = null //keep samplesheet input null if not passed
            def checkPathParamList = [ params.indir, params.multiqc_config, params.kraken2db ]
            for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }
            ch_input_indir = Channel.fromPath(params.indir, relative: true)
        } else { // if no samplesheet is passed and no input directory is given
            exit 1, 'For -entry SCAFFOLDS: You need EITHER an input samplesheet or a directory!' 
        }
    }

    main:
        SCAFFOLDS_EXTERNAL ( ch_input, ch_input_indir )

    emit:
        scaffolds        = SCAFFOLDS_EXTERNAL.out.scaffolds
        mlst             = SCAFFOLDS_EXTERNAL.out.mlst
        amrfinder_report = SCAFFOLDS_EXTERNAL.out.amrfinder_report
        gamma_ar         = SCAFFOLDS_EXTERNAL.out.gamma_ar
        summary_report   = SCAFFOLDS_EXTERNAL.out.summary_report
}

//
// WORKFLOW: Entry point to analyze scaffold file(s) and run everything after Spades
//
workflow CDC_SCAFFOLDS {
    // Validate input parameters
    // Check input path parameters to see if they exist
    if (params.input != null ) {  // if a samplesheet is passed
        // allow input to be relative
        //input_samplesheet_path = Channel.fromPath(params.input, relative: true)
        if (params.indir != null ) { //if samplesheet is passed and an input directory exit
            exit 1, 'For -entry CDC_SCAFFOLDS: You need EITHER an input samplesheet or a directory! Just pick one.' 
        } else { // if only samplesheet is passed check to make sure input is an actual file
            def checkPathParamList = [ params.input, params.multiqc_config, params.kraken2db ]
            for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }
            ch_input_indir = null //keep input directory null if not passed
            // get full path for input and make channel
            if (params.input) { ch_input = file(params.input) }
        }
    } else {
        if (params.indir != null ) { // if no samplesheet is passed, but an input directory is given
            ch_input = null //keep samplesheet input null if not passed
            def checkPathParamList = [ params.indir, params.multiqc_config, params.kraken2db ]
            for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }
            ch_input_indir = Channel.fromPath(params.indir, relative: true)
        } else { // if no samplesheet is passed and no input directory is given
            exit 1, 'For -entry CDC_SCAFFOLDS: You need EITHER an input samplesheet or a directory!' 
        }
    }

    main:
        SCAFFOLDS_EXQC ( ch_input, ch_input_indir )

    emit:
        scaffolds        = SCAFFOLDS_EXQC.out.scaffolds
        mlst             = SCAFFOLDS_EXQC.out.mlst
        amrfinder_report = SCAFFOLDS_EXQC.out.amrfinder_report
        gamma_ar         = SCAFFOLDS_EXQC.out.gamma_ar
        summary_report   = SCAFFOLDS_EXQC.out.summary_report
}

/*
========================================================================================
    THE END
========================================================================================
*/
