process CHECK_MLST_WITH_SRST2 {
    tag "$meta.id"
    label 'process_single'
    container "quay.io/jvhagey/phoenix:base_v2.1.0"

    input:
    tuple val(meta), path(mlst_file), path(srst2_file), path(taxonomy_file), val(status), path(local_dbases)

    output:
    tuple val(meta), path("*_combined.tsv"), emit: checked_MLSTs
    tuple val(meta), path("*_status.txt"),   emit: status
    path("versions.yml"),                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    if (params.ica==false) {
        ica = ""
    } else if (params.ica==true) {
        ica = "python ${workflow.launchDir}/bin/"
    } else {
        error "Please set params.ica to either \"true\" if running on ICA or \"false\" for all other methods."
    }
    def container = task.container.toString() - "quay.io/jvhagey/phoenix:"
    """
    if [[ "${status[0]}" == "True" ]]; then
        ${ica}fix_MLST2.py --input $mlst_file --srst2 $srst2_file --taxonomy $taxonomy_file --mlst_database $local_dbases
    elif [[ "${status[0]}" == "False" ]]; then
        ${ica}fix_MLST2.py --input $mlst_file --taxonomy $taxonomy_file --mlst_database $local_dbases
    else 
        echo "Something went very wrong, please open an issue on Github for the PHoeNIx developers to address."
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}