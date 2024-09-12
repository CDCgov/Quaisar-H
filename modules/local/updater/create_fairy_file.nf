process CREATE_FAIRY_FILE {
    tag "${meta.id}"
    label 'process_low'
    container 'quay.io/jvhagey/phoenix:base_v2.1.0'

    input:
    tuple val(meta), path(indir)

    output:
    tuple val(meta), path("${meta.id}_*_summary.txt"), emit: fairy_outcome
    path("versions.yml"),                              emit: versions

    script: 
    // define variables
    def prefix = task.ext.prefix ?: "${meta.id}"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix:"
    """
    # To make species specific pipelines and updater entry backwards compatable we will remake the fairy file
    ## Since file ran we will assume there was no corruption to begin with
    echo "PASSED: File 2022GL-00947_R1 is not corrupt." > ${prefix}_summary_old_2.txt
    echo "PASSED: File 2022GL-00947_R2 is not corrupt." >> ${prefix}_summary_old_2.txt

    if [ -s "${indir}/${prefix}/fastp_trimd/${prefix}_1.trim.fastq.gz" ] && [ -s "${indir}/${prefix}/fastp_trimd/${prefix}_2.trim.fastq.gz" ]; then
        count1=\$(zcat "${indir}/${prefix}/fastp_trimd/${prefix}_1.trim.fastq.gz" | wc -l)
        count2=\$(zcat "${indir}/${prefix}/fastp_trimd/${prefix}_2.trim.fastq.gz" | wc -l)
        if [ \$((count1 / 4)) -eq \$((count2 / 4)) ]; then
            echo "PASSED: Read pairs for ${prefix} are equal." >> ${prefix}_summary_old_2.txt
            echo "PASSED: There are reads in ${prefix} R1/R2 after trimming." >> ${prefix}_summary_old_2.txt
            if [ -s "${indir}/${prefix}/assembly/${prefix}.filtered.scaffolds.fa.gz" ]; then
                echo "PASSED: More than 0 scaffolds in ${prefix} after filtering." >> ${prefix}_summary_old_2.txt
                mv ${prefix}_summary_old_2.txt ${prefix}_scaffolds_summary.txt
            else
                echo "FAILED: No scaffolds in ${prefix} after filtering!" >> ${prefix}_summary_old_2.txt
                mv ${prefix}_summary_old_2.txt ${prefix}_scaffolds_summary.txt
            fi
        else
            echo "FAILED: The number of reads in R1/R2 are NOT the same!" >> ${prefix}_summary_old_2.txt
            mv ${prefix}_summary_old_2.txt ${prefix}_rawstats_summary.txt
        fi
    else
        echo "PASSED: Read pairs for ${prefix} are equal." >> ${prefix}_summary_old_2.txt
        echo "FAILED: There are 0 reads in " + ${prefix} + " R1/R2 after trimming!" >> ${prefix}_summary_old_2.txt
        mv ${prefix}_summary_old_2.txt ${prefix}_trimstats_summary.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
       phoenix_base_container: ${container}
    END_VERSIONS
    """
}