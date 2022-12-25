#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=nf_tram
#SBATCH --output=slurm-nf_tram-%j.out

# ==============================================================================
#                                FUNCTIONS
# ==============================================================================
# Help function
Print_help() {
    echo
    echo "================================================================================"
    echo "                        $0"
    echo "             Run the TRanscriptome AsseMbly (tram) Nextflow workflow"
    echo "================================================================================"
    echo "REQUIRED OPTIONS:"
    echo "  -i/--indir      <dir>   Input dir with FASTQ files"
    echo "  -o/--outdir     <dir>   Output directory for workflow results"
    echo "  --busco_db      <str>   BUSCO database name (see https://busco.ezlab.org/list_of_lineages.html)"
    echo "  --taxon         <str>   Taxon name for EnTAP, using format 'homo_sapiens' (underscores, no spaces)"
    echo "  --contam        <str>   Comma-separated list of contaminant taxa for EnTAP (e.g., 'viruses,bacteria')"
    echo
    echo "OTHER INPUT DATA OPTIONS:"
    echo "  -q/--fq_pattern <str>   FASTQ file pattern (in single quotes)                   [default: '*_R{1,2}*.fastq.gz']"
    echo "  --ref_fasta     <file>  Genome FASTA file for Trinity ref-guided assembly       [default: none]"
    echo "  --subset_fq     <int>   Subset (subsample) FASTQ files to <int> reads           [default: use all reads]"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --trim_nextseq          Use NextSeq/NovaSeq polyG-trimming TrimGalore option    [default: don't use]"
    echo "  --kraken_db_url <URL>   URL to a URL to Kraken database/index from https://benlangmead.github.io/aws-indexes/k2 [default: https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20220926.tar.gz]"
    echo "  --more_args     <str>   Quoted string with additional arguments to pass to 'nextflow run'"
    echo
    echo "NEXTFLOW-RELATED OPTIONS:"
    echo "  -n/--nf_file    <file>  Nextflow workflow definition file                       [default: 'workflows/nf-transcriptome-assembly/main.nf']"
    echo "  -no-resume              Don't attempt to resume workflow run, but start over    [default: resume]"
    echo "  -profile        <str>   Profile from any of the config files to use             [default: 'conda,normal']"
    echo "  --container_dir <dir>   Singularity container dir                               [default: '/fs/project/PAS0471/containers']"
    echo "                            - This is where any containers used in the workflow will be downloaded to"
    echo "  -work-dir       <dir>   Scratch (work) dir for the workflow                     [default: Nextflow default = 'work']"
    echo "                            - This is where the workflow results will be stored before final results are copied to the specified output dir"
    echo "                            - This should preferable be a dir in OSC's scratch dir rather than in the main project dir"
    echo "  -config         <file>  Additional config file                                  [default: none]"
    echo "                            - Settings in this file will override default settings"
    echo "                            - Note that the mcic-scripts OSC config file will always be included, too"
    echo "                              (https://github.com/mcic-osu/mcic-scripts/blob/main/nextflow/osc.config)"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --dryrun                Dry run: don't execute commands, only parse arguments and report"
    echo "  --debug                 Run the script in debug mode (print all code)"
    echo "  -h                      Print this help message and exit"
    echo "  --help                  Print the Nextflow workflow's help and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i data/fastq -o results/tram"
}

# Load the software
Load_software() {
    # Load OSC's Conda module
    module load miniconda3/4.12.0-py39

    # Activate the Nextflow Conda environment
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate; done
    source activate /fs/project/PAS0471/jelmer/conda/nextflow

    ## Singularity container dir - any downloaded containers will be stored here
    export NXF_SINGULARITY_CACHEDIR="$container_dir"
    # Limit memory for Nextflow main process - see https://www.nextflow.io/blog/2021/5_tips_for_hpc_users.html
    export NXF_OPTS='-Xms1g -Xmx4g'
}

# Print help for the focal program
Print_help_workflow() {
    Load_software
    nextflow run "$nf_file" --help
}

# Print SLURM job resource usage info
Resource_usage() {
    echo
    ${e}sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%60,Elapsed,CPUTime,MaxVMSize | \
        grep -Ev "ba|ex"
    echo
}

# Resource usage information
Time() {
    /usr/bin/time -f \
        '\n# Ran the command:\n%C \n\n# Run stats by /usr/bin/time:\nTime: %E   CPU: %P    Max mem: %M K    Exit status: %x \n' \
        "$@"
}   

# Exit upon error with a message
Die() {
    error_message=${1}
    error_args=${2-none}
    
    echo
    echo "====================================================================="
    printf "$0: ERROR: %s\n" "$error_message" >&2
    echo -e "\nFor help, run this script with the '-h' / '--help' option"
    if [[ "$error_args" != "none" ]]; then
        echo -e "\nArguments passed to the script:"
        echo "$error_args"
    fi
    echo -e "\nEXITING..." >&2
    echo "====================================================================="
    echo
    exit 1
}


# ==============================================================================
#                     CONSTANTS AND DEFAULTS
# ==============================================================================
# URL to OSC Nextflow config file
OSC_CONFIG_URL=https://raw.githubusercontent.com/mcic-osu/mcic-scripts/main/nextflow/osc.config
osc_config=mcic-scripts/nextflow/osc.config  # Will be downloaded if not present here

# Option defaults
fq_pattern='*_R{1,2}*.fastq.gz'
container_dir=/fs/project/PAS0471/containers
nf_file="workflows/nf-transcriptome-assembly/main.nf"
profile="conda,normal"
resume=true && resume_arg="-resume"
trim_nextseq=false && nextseq_arg=""
kraken_db_url="https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20220926.tar.gz"

debug=false
dryrun=false && e=""
slurm=true


# ==============================================================================
#                     PARSE COMMAND-LINE OPTIONS
# ==============================================================================
# Placeholder defaults
indir=""
outdir=""
busco_db=""
taxon=""
contam=""
config_file="" && config_arg=""
more_args=""
work_dir="" && work_dir_arg=""
subset_fq="" && subset_arg=""
ref_fasta="" && ref_arg=""

# Parse command-line options
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --indir )          shift && indir=$1 ;;
        -q | --fq_pattern )     shift && fq_pattern=$1 ;;
        -o | --outdir )         shift && outdir=$1 ;;
        -n | --nf_file )        shift && nf_file=$1 ;;
        --busco_db )            shift && busco_db=$1 ;;
        --taxon )               shift && taxon=$1 ;;
        --contam )              shift && contam=$1 ;;
        --ref_fasta )           shift && ref_fasta=$1 ;;
        --kraken_db_url )       shift && kraken_db_url=$1 ;;
        --trim_nextseq )        trim_nextseq=true ;;
        --subset_fq )           shift && subset_fq=$1 ;;
        --container_dir )       shift && container_dir=$1 ;;
        --more_args )           shift && more_args=$1 ;;
        -config )               shift && config_file=$1 ;;
        -profile )              shift && profile=$1 ;;
        -work-dir )             shift && work_dir=$1 ;;
        -no-resume )            resume=false ;;
        --debug )               debug=true ;;
        --dryrun )              dryrun=true && e="echo ";;
        -h )                    Print_help; exit ;;
        --help )                Print_help_workflow; exit ;;
        * )                     Die "Invalid option $1" "$all_args" ;;
    esac
    shift
done


# ==============================================================================
#                          OTHER SETUP
# ==============================================================================
[[ "$debug" = true ]] && set -o xtrace

# Check if this is a SLURM job
[[ -z "$SLURM_JOB_ID" ]] && slurm=false

# Load Conda environment
[[ "$dryrun" = false ]] && Load_software

# Bash strict settings
set -ueo pipefail

# Check input
[[ "$indir" = "" ]] && Die "Please specify an input dir with -i" "$all_args"
[[ "$outdir" = "" ]] && Die "Please specify an input dir with -o" "$all_args"
[[ "$busco_db" = "" ]] && Die "Please specify an BUSCO db name with --busco_db" "$all_args"
[[ "$taxon" = "" ]] && Die "Please specify a taxon name with --taxon" "$all_args"
[[ "$contam" = "" ]] && Die "Please specify a list of contaminant taxa --contam" "$all_args"
[[ ! -d "$indir" ]] && Die "Input dir $indir does not exist"

# Get the OSC config file
if [[ ! -f "$osc_config" ]]; then
    osc_config="$outdir"/$(basename "$OSC_CONFIG_URL")
    if [[ ! -f $(basename "$OSC_CONFIG_URL") ]]; then
        wget -q -O "$osc_config" "$OSC_CONFIG_URL"
    fi
fi

# Define trace output dir
trace_dir="$outdir"/pipeline_info

# Build the config argument
config_arg="-c $osc_config"
[[ "$config_file" != "" ]] && config_arg="$config_arg -c ${config_file/,/ -c }"

# Build other Nextflow arguments
[[ "$resume" = false ]] && resume_arg=""
[[ "$subset_fq" != "" ]] && subset_arg="--subset_fq $subset_fq"
[[ "$trim_nextseq" = true ]] && nextseq_arg="--trim_nextseq"
[[ "$ref_fasta" != "" ]] && ref_arg="--ref_fasta $ref_fasta"
[[ "$work_dir" != "" ]] && work_dir_arg="-work-dir $work_dir"

# Report
echo
echo "=========================================================================="
echo "                         STARTING SCRIPT NF_TRAM.SH"
date
echo "=========================================================================="
echo "INPUT/OUTPUT DATA OPTIONS:"
echo "  Input dir:                       $indir"
echo "  FASTQ pattern:                   $fq_pattern"
echo "  Output dir:                      $outdir"
[[ "$ref_fasta" != "" ]] && echo "  Reference genome FASTA:          $ref_fasta"
[[ "$subset_fq" != "" ]] && echo "  Nr. of reads to subset FASTQ to: $subset_fq"
echo
echo "RUN SETTINGS:"
echo "  BUSCO database:                  $busco_db"
echo "  NextSeq/NovaSeq polyG-trimming:  $trim_nextseq"
echo "  Kraken database URL:             $kraken_db_url"
echo "  Taxon name:                      $taxon"
echo "  Contaminant taxa:                $contam"
[[ "$more_args" != "" ]] && echo "  Additional arguments:            $more_args"
echo
echo "NEXTFLOW-RELATED OPTIONS:"
echo "  Resume previous run:             $resume"
echo "  Nextflow workflow file:          $nf_file"
echo "  Work (scratch) dir:              $work_dir"
echo "  Container dir:                   $container_dir"
echo "  Config 'profile':                $profile"
[[ "$config_file" != "" ]] && echo "  Additional config file:          $config_file"
[[ $dryrun = true ]] && echo -e "\nTHIS IS A DRY-RUN"
echo "=========================================================================="
echo

# ==============================================================================
#                               RUN
# ==============================================================================
# Make necessary dirs
${e}mkdir -pv "$work_dir" "$container_dir" "$outdir"/logs "$trace_dir"

# Remove old trace files
[[ -f "$trace_dir"/report.html ]] && ${e}rm "$trace_dir"/report.html
[[ -f "$trace_dir"/trace.txt ]] && ${e}rm "$trace_dir"/trace.txt
[[ -f "$trace_dir"/timeline.html ]] && ${e}rm "$trace_dir"/timeline.html
[[ -f "$trace_dir"/dag.png ]] && ${e}rm "$trace_dir"/dag.png

# Run the workflow
echo -e "# Starting the workflow...\n"

[[ "$dryrun" = false ]] && set -o xtrace

${e}Time nextflow run \
        "$nf_file" \
        --reads "$indir/$fq_pattern" \
        --outdir "$outdir" \
        --busco_db "$busco_db" \
        --kraken_db_url "$kraken_db_url" \
        --entap_taxon "$taxon" \
        --entap_contam "$contam" \
        -ansi-log false \
        -with-report "$trace_dir"/report.html \
        -with-trace "$trace_dir"/trace.txt \
        -with-timeline "$trace_dir"/timeline.html \
        -with-dag "$trace_dir"/dag.png \
        -profile "$profile" \
        $nextseq_arg \
        $work_dir_arg \
        $ref_arg \
        $subset_arg \
        $config_arg \
        $resume_arg \
        $more_args

[[ "$debug" = false ]] && set +o xtrace


# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo
echo "========================================================================="
if [[ "$dryrun" = false ]]; then
    [[ "$slurm" = true ]] && Resource_usage
fi
echo "# Done with script"
date
