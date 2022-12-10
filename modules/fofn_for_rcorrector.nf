process FOFN_FOR_RCORRECTOR {
    tag "Combining TrimGalore output for Rcorrector step"
    label "local_process"

    input:
    val filenames

    output:
    path "fq_all.txt"

    script:
    """
    echo ${filenames} | tr "," "\n" | tr -d " ][" | grep '.gz' | sort > fq_all.txt
    """
}
