process MV_IN_1DIR {
    tag "Move (e.g.) FASTQ files into a single dir"
    label "local_process"

    input:
    path files

    output:
    path "files_all"

    script:
    """
    mkdir -p files_all
    mv -v ${files} files_all
    """
}
