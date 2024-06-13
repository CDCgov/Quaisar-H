process CLIA_GRIPHIN {
    label 'process_low'
    container 'quay.io/jvhagey/phoenix:base_v2.1.0'

    input:
    path(pipeline_stats_files)
    path(fairy_summaries)
    path(original_samplesheet)
    path(db)
    path(outdir) // output directory used as prefix for the summary file
    val(coverage)
    path(spades_outcome_files)

    output:
    path("*_GRiPHin_Summary.xlsx"),    emit: griphin_report
    path("*_GRiPHin_Summary.tsv"),     emit: griphin_tsv_report
    path("Phoenix_Summary.tsv"),       emit: phoenix_tsv_report
    path("Directory_samplesheet.csv"), emit: converted_samplesheet
    path("versions.yml"),              emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    if (params.ica==false) { ica = "" } 
    else if (params.ica==true) { ica = "python ${params.bin_dir}" }
    else { error "Please set params.ica to either \"true\" if running on ICA or \"false\" for all other methods." }
    // define variables
    def container = task.container.toString() - "quay.io/jvhagey/phoenix:"
    """
    full_path=\$(readlink -f ${outdir})

    ${ica}CLIA_GRiPHin.py -d \$full_path -a $db --output ${outdir} --coverage ${coverage}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
       python: \$(python --version | sed 's/Python //g')
       griphin.py: \$(${ica}GRiPHin.py --version)
       phoenix_base_container: ${container}
    END_VERSIONS
    """
}