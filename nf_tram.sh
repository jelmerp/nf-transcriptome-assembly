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
## Help function
Print_help() {
    echo
    echo "================================================================================"
    echo "                        $0"
    echo "             Run the Nextflow TRanscriptome AsseMbly (tram) pipeline"
    echo "================================================================================"
    echo "REQUIRED OPTIONS:"
    echo "  -i/--indir      <dir>   Input dir with FASTQ files"
    echo "  -o/--outdir     <dir>   Output directory for workflow results"
    echo "  --busco-db      <str>   BUSCO database name (see https://busco.ezlab.org/list_of_lineages.html)"
    echo "  --taxon         <str>   #TODO"
    echo "  --contam        <str>   #TODO"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  -q/--fq-pattern <str>   FASTQ file pattern (in single quotes)                   [default: '*_R{1,2}*.fastq.gz']"
    echo "  --ref-fasta     <file>  Genome FASTA file for Trinity ref-guided assembly       [default: none]"
    echo "  --trim-nextseq          Use NextSeq/NovaSeq polyG-trimming TrimGalore option    [default: don't use]"
    echo "  --subset-fastq  <int>   Subset (subsample) FASTQ files to <int> reads           [default: use all reads]"
    echo "  --more-args     <str>   Quoted string with additional arguments to pass to 'nextflow run'"
    echo
    echo "NEXTFLOW-RELATED OPTIONS:"
    echo "  -n/--nf-file    <file>  Nextflow workflow definition file                       [default: 'workflows/nf-transcriptome-assembly/main.nf']"
    echo "  -no-resume              Don't attempt to resume workflow run, but start over    [default: resume]"
    echo "  -profile        <str>   Profile from any of the config files to use             [default: 'conda,normal']"
    echo "  --container-dir <dir>   Singularity container dir                               [default: '/fs/project/PAS0471/containers']"
    echo "                            - This is where any containers used in the workflow will be downloaded to"
    echo "  -work-dir       <dir>   Scratch (work) dir for the workflow                     [default: '/fs/scratch/PAS0471/$USER/nf_tram']"
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
    echo "  -h/--help               Print this help message and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i data/fastq -o results/tram"
}

## Load the software
Load_software() {
    ## Load OSC's Conda module
    module load miniconda3/4.12.0-py39

    ## Activate the Nextflow Conda environment
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate; done
    source activate /fs/project/PAS0471/jelmer/conda/nextflow

    ### Singularity container dir - any downloaded containers will be stored here
    export NXF_SINGULARITY_CACHEDIR="$container_dir"
    ## Limit memory for Nextflow main process - see https://www.nextflow.io/blog/2021/5_tips_for_hpc_users.html
    export NXF_OPTS='-Xms1g -Xmx4g'
}

## Print SLURM job resource usage info
Resource_usage() {
    echo
    ${e}sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%60,Elapsed,CPUTime,MaxVMSize | \
        grep -Ev "ba|ex"
    echo
}

## Resource usage information
Time() {
    /usr/bin/time -f \
        '\n# Ran the command:\n%C \n\n# Run stats by /usr/bin/time:\nTime: %E   CPU: %P    Max mem: %M K    Exit status: %x \n' \
        "$@"
}   

## Exit upon error with a message
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
## URL to OSC Nextflow config file
OSC_CONFIG_URL=https://raw.githubusercontent.com/mcic-osu/mcic-scripts/main/nextflow/osc.config
osc_config=mcic-scripts/nextflow/osc.config  # Will be downloaded if not present here

## Option defaults
fq_pattern='*_R{1,2}*.fastq.gz'
container_dir=/fs/project/PAS0471/containers
work_dir=/fs/scratch/PAS0471/$USER/nf_tram
nf_file="workflows/nf-transcriptome-assembly/main.nf"
profile="conda,normal"
resume=true && resume_arg="-resume"
trim_nextseq=false && nextseq_arg=""

debug=false
dryrun=false && e=""
slurm=true


# ==============================================================================
#                     PARSE COMMAND-LINE OPTIONS
# ==============================================================================
## Placeholder defaults
indir=""
outdir=""
busco_db=""
taxon=""
contam=""
config_file="" && config_arg=""
more_args=""
work_dir_arg=""
subset_fastq="" && subset_arg=""
ref_fasta="" && ref_arg=""

## Parse command-line options
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --indir )          shift && indir=$1 ;;
        -q | --fq-pattern )     shift && fq_pattern=$1 ;;
        -o | --outdir )         shift && outdir=$1 ;;
        -n | --nf_file )        shift && nf_file=$1 ;;
        --busco-db )            shift && busco_db=$1 ;;
        --taxon )               shift && taxon=$1 ;;
        --contam )              shift && contam=$1 ;;
        --ref-fasta )           shift && ref_fasta=$1 ;;
        --trim-nextseq )        trim_nextseq=true ;;
        --subset-fastq )        shift && subset_fastq=$1 ;;
        --container-dir )       shift && container_dir=$1 ;;
        --more-args )           shift && more_args=$1 ;;
        -config )               shift && config_file=$1 ;;
        -profile )              shift && profile=$1 ;;
        -work-dir )             shift && work_dir_arg=$1 ;;
        -no-resume )            resume=false ;;
        --debug )               debug=true ;;
        --dryrun )              dryrun=true && e="echo ";;
        -h | --help )           Print_help; exit ;;
        * )                     Die "Invalid option $1" "$all_args" ;;
    esac
    shift
done


# ==============================================================================
#                          OTHER SETUP
# ==============================================================================
[[ "$debug" = true ]] && set -o xtrace

## Check if this is a SLURM job
[[ -z "$SLURM_JOB_ID" ]] && slurm=false

## Load Conda environment
[[ "$dryrun" = false ]] && Load_software

## Bash strict settings
set -ueo pipefail

## Check input
[[ "$indir" = "" ]] && Die "Please specify an input dir with -i" "$all_args"
[[ "$outdir" = "" ]] && Die "Please specify an input dir with -o" "$all_args"
[[ "$busco_db" = "" ]] && Die "Please specify an BUSCO db name with --busco_db" "$all_args"
[[ ! -d "$indir" ]] && Die "Input dir $indir does not exist"

## Get the OSC config file
if [[ ! -f "$osc_config" ]]; then
    [[ ! -f $(basename "$OSC_CONFIG_URL") ]] && wget -q "$OSC_CONFIG_URL"
    osc_config=$(basename "$OSC_CONFIG_URL")
fi

## Define trace output dir
trace_dir="$outdir"/pipeline_info

## Build the config argument
config_arg="-c $osc_config"
if [[ "$config_file" != "" ]]; then
    config_arg="$config_arg -c ${config_file/,/ -c }"
fi

## Build other Nextflow arguments
[[ "$resume" = false ]] && resume_arg=""
[[ "$subset_fastq" != "" ]] && subset_arg="--subset_fastq $subset_fastq"
[[ "$trim_nextseq" = true ]] && nextseq_arg="--trim_nextseq"
[[ "$ref_fasta" != "" ]] && ref_arg="--ref_fasta $ref_fasta"

## Work dir
if [[ "$work_dir_arg" = "" ]]; then
    ## If using the default, add a run ID using the outdir
    work_dir=$work_dir/$(basename "$outdir")
else
    ## If a work_dir was provided as an arg, use that one as-is
    work_dir="$work_dir_arg"
fi

## Report
echo
echo "=========================================================================="
echo "                         STARTING SCRIPT NF_TRAM.SH"
date
echo "=========================================================================="
echo "RUN OPTIONS:"
echo "Input dir:                       $indir"
echo "FASTQ pattern:                   $fq_pattern"
echo "Output dir:                      $outdir"
echo "BUSCO database:                  $busco_db"
echo "NextSeq/NovaSeq polyG-trimming:  $trim_nextseq"
echo "Taxon name:                      $taxon"
echo "Contaminant taxa:                $contam"
[[ "$ref_fasta" != "" ]] && echo "Reference genome FASTA:          $ref_fasta"
[[ "$subset_fastq" != "" ]] && echo "Nr. of reads to subset FASTQ to: $subset_fastq"
[[ "$more_args" != "" ]] && echo "Additional arguments:            $more_args"
echo
echo "NEXTFLOW-RELATED OPTIONS:"
echo "Resume previous run:             $resume"
echo "Nextflow workflow file:          $nf_file"
echo "Work (scratch) dir:              $work_dir"
echo "Container dir:                   $container_dir"
echo "Config 'profile':                $profile"
[[ "$config_file" != "" ]] && echo "Additional config file:          $config_file"
[[ $dryrun = true ]] && echo -e "\nTHIS IS A DRY-RUN"
echo "=========================================================================="
echo

# ==============================================================================
#                               RUN
# ==============================================================================
## Make necessary dirs
${e}mkdir -p "$work_dir" "$container_dir" "$outdir"/logs "$trace_dir"

## Remove old trace files
[[ -f "$trace_dir"/report.html ]] && ${e}rm "$trace_dir"/report.html
[[ -f "$trace_dir"/trace.txt ]] && ${e}rm "$trace_dir"/trace.txt
[[ -f "$trace_dir"/timeline.html ]] && ${e}rm "$trace_dir"/timeline.html
[[ -f "$trace_dir"/dag.png ]] && ${e}rm "$trace_dir"/dag.png

## Run the workflow
echo -e "# Starting the workflow...\n"

[[ "$dryrun" = false ]] && set -o xtrace

${e}Time nextflow run \
        "$nf_file" \
        --reads "$indir/$fq_pattern" \
        --outdir "$outdir" \
        --busco_db "$busco_db" \
        --taxon "$taxon" \
        --taxon "$contam" \
        -work-dir "$work_dir" \
        -ansi-log false \
        -with-report "$trace_dir"/report.html \
        -with-trace "$trace_dir"/trace.txt \
        -with-timeline "$trace_dir"/timeline.html \
        -with-dag "$trace_dir"/dag.png \
        -profile "$profile" \
        $nextseq_arg \
        $ref_arg \
        $subset_arg \
        $config_arg \
        $resume_arg \
        $more_args

[[ "$debug" = false ]] && set +o xtrace

#TODO - Don't need report files, with config?


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
