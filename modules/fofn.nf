process FOFN {
    tag "Make a FOFN from a Nextflow tuple with sample names and file names"
    label "local_process"

    input:
    val filenames
    val pattern

    output:
    path "fofn.txt"

    script:
    """
    echo ${filenames} | tr "," "\n" | tr -d " ][" | grep "${pattern}" | sort > fofn.txt.txt
    """
}
