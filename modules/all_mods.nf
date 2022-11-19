// Processes
process SUBSET_FASTQ {
    tag "Subsample FASTQ files for $sample_id"
    publishDir "${params.outdir}/subset_fastq", mode: 'copy'

    input:
    tuple val(sample_id), path(fq_pair)
    val(n_reads)

    output:
    tuple val(sample_id), path("out/*fastq.gz"), emit: fq 
    path "out/logs/slurm*log", emit: log
    path "out/logs/version.txt", emit: version

    script:
    """
    fqsub.sh --R1_in ${fq_pair[0]} --outdir out --n_reads ${n_reads}
    
    cp .command.log out/logs/slurm-${sample_id}.log
    """
}

process FASTQC {
    tag "FastQC to QC FASTQ files for $sample_id"
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(fq_pair)

    output:
    path "*zip", emit: zip
    path "*html", emit: html
    path "logs/slurm*log", emit: log

    script:
    """
    fastqc $fq_pair
    
    mkdir -p logs
    cp .command.log logs/slurm-${sample_id}.log
    """
}

process MULTIQC_RAW {
    tag "MultiQC for FastQC output for raw FASTQ files"
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path(all_logs)

    output:
    path "*html", emit: report
    path "logs/slurm*log", emit: log

    script:
    """
    multiqc --filename multiqc_raw.html --interactive $all_logs
    
    mkdir -p logs
    cp .command.log logs/slurm-multiqc_raw.log
    """
}

process MULTIQC_TRIM {
    tag "MultiQC for trimmed output"
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path(all_logs)

    output:
    path "*html", emit: report
    path "logs/slurm*log", emit: log

    script:
    """
    multiqc --filename multiqc_trimmed.html --interactive $all_logs
    
    mkdir -p logs
    cp .command.log logs/slurm-multiqc_trimmed.log
    """
}

process MULTIQC_PRE {
    tag "MultiQC for other pre-assembly steps"
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path(all_logs)

    output:
    path "*html", emit: report
    path "logs/slurm*log", emit: log

    script:
    """
    multiqc --filename multiqc_preprocess.html --interactive $all_logs

    mkdir -p logs
    cp .command.log logs/slurm-multiqc_preprocess.log
    """
}

process TRIMGALORE {
    tag "TrimGalore to remove adapters and quality-trim for $sample_id"
    publishDir "${params.outdir}/trimgalore", mode: 'copy'

    input:
    tuple val(sample_id), path(fq_pair)

    output:
    path "trimmed/*fastq.gz", emit: fq_trimmed
    path "logs/*_trimming_report.txt", emit: trim_report
    path "logs/slurm*log", emit: log
    path "logs/version.txt", emit: version
    path "fastqc/*zip", emit: fqc_zip
    path "fastqc/*html", emit: fqc_html

    script:
    """
    if [[ ${params.trim_nextseq} = true ]]; then
        nextseq_arg="--nextseq"
    else
        nextseq_arg=""
    fi

    trimgalore.sh --R1 ${fq_pair[0]} --outdir . --quality 5 --length 36 \${nextseq_arg}
    cp .command.log logs/slurm-${sample_id}.log
    """
}

process CAT_FILENAMES {
    tag "Combining TrimGalore output for Rcorrector step"
    label "local_process"

    input:
    val filename

    output:
    path "fq_all.txt"

    script:
    """
    echo $filename | tr "," "\n" | tr -d " ][" | grep '.gz' | sort > fq_all.txt
    """
}

process RCORRECTOR {
    tag "Rcorrector read correction"
    publishDir "${params.outdir}/rcorrector", mode: 'copy'

    input:
    val fq_list

    output:
    path "*fq.gz", emit: fq
    path "logs/slurm*log", emit: log
    path "fofn*txt", emit: fofn

    script:
    """
    echo "$fq_list" | grep . > fofn.txt
    subset_id=\$(head -n 1 fofn.txt | xargs -I{} basename {} .fastq.gz)

    rcorrector.sh -I fofn.txt -o .
    
    mv fofn.txt fofn_\$subset_id.txt
    cp .command.log logs/slurm-rcorrector.log
    """
}

process RCORRFILTER {
    tag "Rcorrfilter removal of unfixable reads for $sample_id"
    publishDir "${params.outdir}/rcorrfilter", mode: 'copy'
    
    input:
    tuple val(sample_id), path(fq_pair)

    output:
    tuple val(sample_id), path("*fastq.gz"), emit: fq 
    path "logs/slurm*log", emit: log

    script:
    """
    rcorrfilter.sh -i ${fq_pair[0]} -o .
    
    cp .command.log logs/slurm-${sample_id}.log
    """
}

process SORTMERNA {
    tag "SortMeRNA rRNA removal for $sample_id"
    publishDir "${params.outdir}/sortmerna", mode: 'copy'

    input:
    tuple val(sample_id), path(fq_pair)

    output:
    tuple val(sample_id), path("unmapped/*fastq.gz"), emit: fq_unmapped
    path "mapped/*fastq.gz", emit: fq_mapped 
    path "logs/*log", emit: log
    path "logs/version.txt", emit: version
    
    script:
    """
    sortmerna.sh --R1 ${fq_pair[0]} --outdir .
    
    cp .command.log logs/slurm-${sample_id}.log
    """
}

process KRAKEN {
    tag "Kraken contamination detection for $sample_id"
    publishDir "${params.outdir}/kraken", mode: 'copy'

    input:
    tuple val(sample_id), path(fq_pair)

    output:
    tuple val(sample_id), path("unclassified/*fastq.gz"), emit: fq
    path "unclassified/*fastq.gz", emit: fq_list
    path "classified/*fastq.gz", emit: fq_classified
    tuple val(sample_id), path("*main.txt"), emit: main
    path "*report.txt", emit: report
    path "logs/slurm*log", emit: log
    path "logs/version.txt", emit: version
    
    script:
    """
    kraken.sh \
        --infile ${fq_pair[0]} \
        --outdir . \
        --db-dir $params.krakendb \
        --classified-out \
        --unclassified-out \
        --confidence 0.5
    
    cp .command.log logs/slurm-${sample_id}.log
    """
}

process KRONA {
    tag "Krona plots of Kraken results for $sample_id"
    publishDir "${params.outdir}/krona", mode: 'copy'

    input:
    tuple val(sample_id), path(kraken_res)

    output:
    path "*html"
    path "logs/slurm*log"
    
    script:
    """
    krona.sh --infile ${kraken_res} --outfile ./krona_${sample_id}.html
    
    cp .command.log logs/slurm-${sample_id}.log
    """
}

process ORNA {
    tag "ORNA read normalization for $sample_id"
    publishDir "${params.outdir}/orna", mode: 'copy'

    input:
    tuple val(sample_id), path(fq_pair)

    output:
    path "orna_out/*fastq.gz", emit: fq
    path "orna_out/logs/slurm*log", emit: log
    path "orna_out/logs/version.txt", emit: version
    
    script:
    """
    orna.sh --R1_in ${fq_pair[0]} --outdir orna_out
    
    cp .command.log orna_out/logs/slurm-${sample_id}.log
    """
}

process JOIN_KRAKEN {
    tag "Combining Kraken output for assembly processes"
    label "local_process"

    input:
    path kraken_fastqs

    output:
    path "kraken_fastqs_all"

    script:
    """
    mkdir kraken_fastqs_all
    mv -v $kraken_fastqs kraken_fastqs_all
    """
}

process JOIN_ORNA {
    tag "Combining ORNA output for assembly processes"
    label "local_process"

    input:
    path orna_fastqs

    output:
    path "orna_fastqs_all"

    script:
    """
    mkdir orna_fastqs_all
    mv -v $orna_fastqs orna_fastqs_all
    """
}

process INDEX_GENOME {
    tag "Index a reference genome"
    publishDir "${params.outdir}/index_genome", mode: 'copy'

    input:
    path ref_fasta

    output:
    path "index_dir", emit: index
    path "index_dir/logs/version.txt", emit: version
    
    script:
    """
    star_index.sh -i ${ref_fasta} -o index_dir
    
    cp .command.log index_dir/logs/slurm.log
    
    """
}

process MAP2GENOME {
    tag "Map reads to a reference genome"
    publishDir "${params.outdir}/map2genome", mode: 'copy'

    input:
    path ref_fasta
    tuple val(sample_id), path(fq_pair)

    output:
    path "mapped/*bam", emit: bam
    path "logs/slurm*log", emit: log
    path "logs/version.txt", emit: version
    
    script:
    """
    star_align.sh -i ${fq_pair[0]} -r ${ref_fasta} -o .
    
    cp .command.log logs/slurm.log
    """
}

process TRINITY {
    tag "Transcriptome assembly with Trinity"
    publishDir "${params.outdir}/trinity", mode: 'copy'

    input:
    path dir_with_all_fqs
    each norm

    output:
    path "trinity_out/Trinity.fasta", emit: assembly
    path "trinity_out/logs/slurm*log", emit: log
    path "trinity_out/logs/version.txt", emit: version
    
    script:
    """
    trinity.sh \
        -i ${dir_with_all_fqs} \
        -o trinity_out \
        --min_contig_length ${params.min_contig_length} \
        --normalize ${norm} \
        --strandedness ${params.strandedness}
    
    mv -v Trinity.fasta trinity_norm${norm}.fasta

    cp .command.log trinity_out/logs/slurm.log
    """
}

process TRINITY_GUIDED {
    tag "Transcriptome assembly with Trinity (genome-guided)"
    publishDir "${params.outdir}/trinity_guided", mode: 'copy'

    input:
    path bam

    output:
    path "trinity_out/trinity_gg.fasta", emit: assembly
    path "trinity_out/logs/slurm*log", emit: log
    path "trinity_out/logs/version.txt", emit: version
    
    script:
    """
    trinity.sh \
        -i ${bam} \
        -o trinity_out \
        --genome_guided \
        --genome_guided_max_intron 250000 \
        --strandedness ${params.strandedness}
    
    mv -v Trinity.fasta trinity_gg.fasta

    cp .command.log trinity_out/logs/slurm-${sample_id}.log
    """
}

process TRANSABYSS {
    tag "Transcriptome assembly with Trans-ABySS"
    publishDir "${params.outdir}/transabyss", mode: 'copy'

    input:
    path dir_with_all_fqs
    each kmer_size
    val norm

    output:
    path "*final.fa", emit: assembly
    path "logs/slurm*log", emit: log
    path "logs/version.txt", emit: version
    
    script:
    """
    transabyss.sh \
        -i ${dir_with_all_fqs} \
        -o . \
        --id transabyss_norm${norm}_k${kmer_size} \
        --min_contig_length ${params.min_contig_length} \
        --kmer_size ${kmer_size} \
        --strandedness ${params.strandedness}
    
    cp .command.log logs/slurm.log
    """
}

process SPADES {
    tag "Transcriptome assembly with SPADES"
    publishDir "${params.outdir}/spades", mode: 'copy'

    input:
    path dir_with_all_fqs
    each kmer_size
    val norm

    output:
    path "spades_*.fasta", emit: assembly
    path "logs/slurm*log", emit: log
    path "logs/version.txt", emit: version
    
    script:
    """
    spades.sh \
        --indir ${dir_with_all_fqs} \
        --outdir . \
        --mode rna \
        --kmer_size ${kmer_size} \
        --strandedness ${params.strandedness}
    
    mv transcripts.fasta spades_norm${norm}_k${kmer_size}.fasta

    cp .command.log logs/slurm.log
    """
}

process CONCAT_ASSEMBLIES {
    tag "Concatenate all assemblies into a single FASTA file"
    publishDir "${params.outdir}/concat_assemblies", mode: 'copy'

    input:
    x

    output:
    path "concat_assembly.fasta", emit: assembly
    path "logs/slurm*log", emit: log
    
    script:
    """
    x
    """
}

process EVIGENE {
    tag "Merge all assemblies into a single FASTA file"
    publishDir "${params.outdir}/merged_assembly", mode: 'copy'

    input:
    path concat_assembly

    output:
    path "final/evigene_all.fasta", emit: assembly_all
    path "final/evigene_primarytrans", emit: assembly_primarytrans
    path "final/okayset", emit: okayset
    path "logs/slurm*log", emit: log
    
    script:
    """
    evigene.sh -i ${concat_assembly} -o .

    cp .command.log logs/slurm.log
    """
}

process BUSCO {
    tag "Evaluate an assembly with Busco"
    publishDir "${params.outdir}/busco", mode: 'copy'

    input:
    path assembly

    output:
    path "XX", emit: out
    path "logs/slurm*log", emit: log
    
    //TODO Need to save with assembly ID - runs for each assembly

    script:
    """
    busco.sh -i ${assembly} -d ${params.busco_db} -o .

    cp .command.log logs/slurm.log
    """
}

process RNAQUAST {
    tag "Evaluate an assembly with RNAQuast"
    publishDir "${params.outdir}/rnaquast", mode: 'copy'

    input:
    path assembly

    output:
    path "XX", emit: out
    path "logs/slurm*log", emit: log
    
    //TODO Run for all assemblies at once?

    script:
    """
    rnaquast.sh -i ${assembly} -o .

    cp .command.log logs/slurm.log
    """
}

process DETONATE {
    tag "Evaluate an assembly with Detonate"
    publishDir "${params.outdir}/detonate", mode: 'copy'

    input:
    path assembly
    path dir_with_all_fqs
    //TODO Make sure to get FASTQs that were used for the assembly...

    output:
    path "XX", emit: out
    path "logs/slurm*log", emit: log
    
    script:
    """
    detonate.sh -i ${assembly} -I ${dir_with_all_fqs} -o .

    cp .command.log logs/slurm.log
    """
}

process ENTAP_CONFIG {
    tag "Configure EnTap"
    publishDir "${params.outdir}/entap", mode: 'copy'

    input:
    path refseq_db
    path uniprot_db
    path entap_config

    output:
    path "XX", emit: prepped_dbs
    path "logs/slurm*log", emit: log
    
    script:
    """
    dbs="$refseq_db $uniprot_db"

    entap_config.sh -d \${dbs} -c ${config} -o .

    cp .command.log logs/slurm.log
    """
}

process ENTAP {
    tag "Annotate an assembly with EnTap"
    publishDir "${params.outdir}/entap", mode: 'copy'

    input:
    path assembly
    path refseq_db
    path uniprot_db
    path entap_config

    output:
    path "XX", emit: out
    path "logs/slurm*log", emit: log
    
    script:
    """
    dbs="$refseq_db $uniprot_db"

    entap.sh -i ${assembly} -d \${dbs} -c ${config} -o .

    cp .command.log logs/slurm.log
    """
}

process KALLISTO_INDEX {
    tag "Index the final transcriptome assembly with Kallisto"
    publishDir "${params.outdir}/kallisto_index", mode: 'copy'

    input:
    path merged_assembly

    output:
    path "assembly.idx", emit: index
    path "logs/slurm*log", emit: log
    
    script:
    """
    kallisto_index.sh -i ${merged_assembly} -o assembly.idx

    cp .command.log logs/slurm.log
    """
}

process KALLISTO {
    tag "Index the final transcriptome assembly with Kallisto"
    publishDir "${params.outdir}/kallisto_index", mode: 'copy'

    input:
    tuple val(sample_id), path(fq_pair)
    path assembly_index

    output:
    path "assembly.idx", emit: index
    path "logs/slurm*log", emit: log
    
    script:
    """
    kallisto_quant.sh -i ${fq_pair[0]} -r ${assembly_index} -o .

    cp .command.log logs/slurm.log
    """
}
