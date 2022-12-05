// Processes
process SUBSET_FASTQ {
    tag "Subsample FASTQ files for $sample_id"
    publishDir "${params.outdir}/subset_fastq", mode: 'copy'

    input:
    tuple val(sample_id), path(fq_pair)
    val n_reads

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
    path all_logs

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
    path all_logs

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
    path all_logs

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
    
    mv fofn.txt fofn_\${subset_id}.txt
    cp .command.log logs/slurm-rcorrector_\${subset_id}.log
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
    tag "Map reads to a reference genome for $sample_id"
    publishDir "${params.outdir}/map2genome", mode: 'copy'

    input:
    tuple val(sample_id), path(fq_pair)
    path ref_index

    output:
    path "bam/*bam", emit: bam
    path "logs/slurm*log", emit: log
    path "star_logs/*", emit: star_log
    path "logs/version.txt", emit: version
    
    script:
    """
    star_align.sh -i ${fq_pair[0]} -r ${ref_index} -o .
    
    cp .command.log logs/slurm-${sample_id}.log
    """
}

process MERGE_BAM {
    tag "Merge BAM files for genome-guided assembly"
    publishDir "${params.outdir}/map2genome", mode: 'copy'

    input:
    path bamfiles

    output:
    path "merged.bam", emit: bam
    path "logs/slurm*log", emit: log
    path "logs/*_version.txt", emit: version
    
    script:
    """
    merge_bam.sh -o merged.bam ${bamfiles}
    
    cp .command.log logs/slurm-merge.log
    """
}

process TRINITY {
    tag "Assembly with Trinity - normalization $norm"
    publishDir "${params.outdir}/trinity", mode: 'copy'

    input:
    path dir_with_all_fqs
    each norm

    output:
    path "trinity_out/trinity_norm*fasta", emit: assembly
    path "trinity_out/trinity_norm*gene_trans_map", emit: gene2trans
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
    
    mv -v trinity_out.Trinity.fasta trinity_out/trinity_norm${norm}.fasta
    mv -v trinity_out.Trinity.fasta.gene_trans_map trinity_out/trinity_norm${norm}.gene_trans_map

    cp .command.log trinity_out/logs/slurm.log
    """
}

process TRINITY_GUIDED {
    tag "Assembly with Trinity (genome-guided)"
    publishDir "${params.outdir}/trinity_guided", mode: 'copy'

    input:
    path bam

    output:
    path "trinity_out/trinity_gg.fasta", emit: assembly
    path "trinity_out/trinity_gg.gene_trans_map", emit: gene2trans
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
    
    mv -v trinity_out/Trinity-GG.fasta trinity_out/trinity_gg.fasta
    mv -v trinity_out/Trinity-GG.fasta.gene_trans_map trinity_out/trinity_gg.gene_trans_map

    cp .command.log trinity_out/logs/slurm.log
    """
}

process TRANSABYSS {
    tag "Assembly with Trans-ABySS - kmer $kmer_size - normalization $norm"
    publishDir "${params.outdir}/transabyss", mode: 'copy'

    input:
    path dir_with_all_fqs
    each kmer_size
    val norm

    output:
    path "transabyss_norm*.fasta", emit: assembly
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
    
    mv -v transabyss_norm${norm}_k${kmer_size}-final.fa transabyss_norm${norm}_k${kmer_size}.fasta

    cp .command.log logs/slurm.log
    """
}

process SPADES {
    tag "Assembly with SPADES - kmer $kmer_size - normalization $norm"
    publishDir "${params.outdir}/spades", mode: 'copy'

    input:
    path dir_with_all_fqs
    each kmer_size
    val norm

    output:
    path "spades_norm*.fasta", emit: assembly
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
    publishDir "${params.outdir}/concat_assembly", mode: 'copy'

    input:
    path assemblies

    output:
    path "concat_assembly.fasta", emit: assembly
    path "logs/slurm*log", emit: log
    
    script:
    """
    concat_assemblies.sh -o concat_assembly.fasta ${assemblies}

    cp .command.log logs/slurm.log
    """
}

process EVIGENE {
    tag "Merge all assemblies into a single FASTA file"
    publishDir "${params.outdir}/evigene", mode: 'copy'

    input:
    path concat_assembly

    output:
    path "out/final/evigene_all.fasta", emit: all
    path "out/final/evigene_primarytrans.fasta", emit: primarytrans
    path "out/okayset", emit: okayset
    path "out/logs/slurm*log", emit: log
    
    script:
    """
    evigene.sh -i ${concat_assembly} -o out

    cp .command.log out/logs/slurm.log
    """
}

process BUSCO {
    tag "Evaluate an assembly with Busco"
    publishDir "${params.outdir}/busco", mode: 'copy'

    input:
    path assembly
    val busco_db

    output:
    path "*", emit: out
    path "logs/slurm*log", emit: log
    
    script:
    """
    busco.sh -i ${assembly} -d ${busco_db} -o .

    cp .command.log logs/slurm.log
    """
}

process RNAQUAST {
    tag "Evaluate an assembly with RNAQuast"
    publishDir "${params.outdir}/rnaquast", mode: 'copy'

    input:
    path assembly

    output:
    path "*report*", emit: out
    path "logs/slurm*log", emit: log
    
    script:
    """
    rnaquast.sh \
        --assemblies ${assembly} \
        --outdir . \
        --strandedness ${params.strandedness}

    cp .command.log logs/slurm.log
    """
}

process DETONATE {
    tag "Evaluate an assembly with Detonate"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path assembly
    path dir_with_all_fqs

    output:
    path "detonate", emit: out
    path "detonate/logs/slurm.log", emit: log
    path "detonate/logs/version.txt", emit: version
    
    script:
    """
    detonate.sh \
        --assembly ${assembly} \
        --fq-dir ${dir_with_all_fqs} \
        --outdir detonate

    cp .command.log detonate/logs/slurm.log
    """
}

process MAP2TRANSCRIPTOME {
    tag "Map the reads back to the transcriptome"

    input:
    path assembly
    path dir_with_all_fqs

    output:
    path "map2trans.bam", emit: bam
    path "logs/slurm.log", emit: log
    path "logs/version.txt", emit: version
    
    script:
    """
    bowtie2.sh \
        --assembly ${assembly} \
        --fq-dir ${dir_with_all_fqs} \
        --bam map2trans.bam

    cp .command.log logs/slurm.log
    """
}

process DOWNLOAD_REFSEQ {
    tag "Download an NCBI RefSeq database"
    publishDir "${params.outdir}/dbs/refseq", mode: 'copy'

    input:
    val refseq_type

    output:
    path "refseq_*.fasta"
    
    script:
    """
    wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/${refseq_type}/*protein*faa.gz
    cat complete* | gunzip -c > refseq_${refseq_type}.fasta
    """
}

process DOWNLOAD_NR {
    tag "Download the NCBI NR database"
    publishDir "${params.outdir}/dbs/nr", mode: 'copy'

    output:
    path "nr_database.fasta"

    script:
    """
    wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nr.*.tar.gz
    
    for archive in nr.*tar.gz; do
        tar -xvzf \$archive
    done

    cat nr.* > nr_database.fasta
    """
}

process ENTAP_CONFIG {
    tag "Configure EnTap"
    publishDir "${params.outdir}/entap", mode: 'copy'

    input:
    path protein_dbs

    output:
    path "db_dir", emit: db_dir
    path "entap_config.ini", emit: config
    
    script:
    """
    entap_config.sh \
        --config_out entap_config.ini \
        --db_dir db_dir \
        --taxon ${params.entap_taxon} \
        --contam ${params.entap_contam} \
        ${protein_dbs}

    cp .command.log db_dir/logs/slurm.log
    """
}

process ENTAP {
    tag "Annotate an assembly with EnTap"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path assembly
    path db_dir
    path config
    path bam

    output:
    path "entap", emit: out
    path "logs/slurm*log", emit: log
    
    script:
    """
    entap.sh \
        --assembly ${assembly} \
        --db-dir ${db_dir} \
        --config ${config} \
        --bam ${bam} \
        --outdir entap

    cp .command.log entap/logs/slurm.log
    """
}

process ENTAP_PROCESS {
    tag "Post-processing of EnTap output"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path entap_out
    path assembly_1trans
    path assembly_alltrans

    output:
    path "entap_ed", emit: out
    path "entap_ed/assembly.fasta", emit: assembly
    path "logs/slurm*log", emit: log
    
    script:
    """
    entap_process.sh \
        --entap-dir ${entap_out} \
        --in-1trans ${assembly_1trans} \
        --in-alltrans ${assembly_alltrans} \
        --outdir entap_ed

    cp .command.log entap_ed/logs/slurm.log
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
    tag "Quantify read counts with Kallisto for $sample_id"
    publishDir "${params.outdir}/kallisto", mode: 'copy'

    input:
    tuple val(sample_id), path(fq_pair)
    path assembly_index

    output:
    path "${sample_id}/abundance*h5", emit: counts_h5
    path "${sample_id}/abundance*tsv", emit: counts_tsv
    path "${sample_id}/run_info.json", emit: info
    path "${sample_id}/logs/slurm*log", emit: log
    
    script:
    """
    kallisto_quant.sh -i ${fq_pair[0]} -r ${assembly_index} -o ${sample_id}

    cp .command.log ${sample_id}/logs/slurm.log
    """
}
