#!/usr/bin/env bash

set -eu -o pipefail

# ACTIVATE 'TSO500ichorCNA' CONDA ENVIRONMENT!

# EXECUTE FROM sy096!

#: Title       : From sortedBAM to PDFs (requires sorted BAMs with .bai files)
#: Date        : 07-07-2025
#: Author      : Benjamin Spiegl based on work of Raúl Mejía and Isaac Lazzeri; maintained also by Matvii Mykhailichenko
#: Version     : 1.1
#: Description :
#: Arguments   : $1 path to reference genome (It has to be already indexed. In PATH/TO/REFGENOME/ there should also be the index file fai or bai extensions) (PATH/TO/REFGENOME/ucsc.hg19.fasta)
#:             : $2 path to HMM copy (irrelevant - might be deleted)
#:             : $3 path to ichorCNA installation (parent folder)
#: Notes       : install the software environment using the conda/TSO500ichorCNA.yml file!

DEBUGGING=false

N_CORES_MAX=24  # also scales the requested amount of RAM -> this number * 2 = requested RAM in GB
MB_PER_SCAFFOLD=2000
MB_FOR_INSERT_SIZES=15000
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
if [ "$DEBUGGING" = "true" ]; then
  echo "this is the value of the variable script_dir in 'TSO500_ichorCNA_from_offANDontarget_reads_PARALLEL.sh': ${script_dir}"  # for debugging
fi
SINGLE_SAMPLE_TSO500_OFFTARGET_PIPELINE_SCRIPT="${script_dir}/TSO500_ichorCNA_singleSample_calledTarget.sh"
MARK_DUPLICATES_ABSOLUTE_PATH="${script_dir}/mark_duplicates_and_insert_sizes_for_TSO500_ichorCNA.py"  # MUST be inside of same parent folder as this script!! (written by Benjamin Spiegl -> use up to 24 cores! (N_CORES_MAX variable))
if [ ! -f "${MARK_DUPLICATES_ABSOLUTE_PATH}" ]
then
  echo "ERROR - parallelized mark duplicates Python3 script '${MARK_DUPLICATES_ABSOLUTE_PATH}' not found in current working directory '$(pwd)'! Exiting.."
  exit 2
fi
CUSTOM_TSO500_PON_PATH="${script_dir}/accessory_files/PoN_of_10controls_fromTSO_1000kb.txt_median.rds"
if [ ! -f "${CUSTOM_TSO500_PON_PATH}" ]
then
  echo "ERROR - TSO500 liquid panel of normals (PoN) file '${CUSTOM_TSO500_PON_PATH}' not found! Exiting.."
  exit 2
fi

#--------------------------------------------------------------------------------------------------
#*************************************FUNCTIONS****************************************************
#--------------------------------------------------------------------------------------------------

function apply_to_each_bam(){
  # use sbatch calls instead of serially processing samples!
  #       sbatch --nodelist=sy096 --ntasks=1 --cpus-per-task=24 --mem=240000 [MB]
  #arg1 = single-sample pipeline bash script path
  #arg2 = input directory path (where BAM files are)
  #arg3 = /PATH/TO/ichorCNAdir
  #arg4 = number of cores to use
  #echo "---------------------------------"
  #echo "parameters received:"
  #echo $1
  #echo $2
  #echo $3
  #echo $4
  #echo "---------------------------------"
  single_sample_ichor_script="${1}"
  if [ ! -f "${single_sample_ichor_script}" ]
  then
    echo "ERROR - single sample ichorCNA computation script '${single_sample_ichor_script}' either is not accessible or does not exist! Exiting.."
    exit 2
  fi
  input_directory_path="${2}"
  if [ ! -d "${input_directory_path}" ]
  then
    echo "ERROR - path to alleged input directory '${input_directory_path}' does not exist! Exiting.."
    exit 2
  fi
  ichor_installation_path="${3}"
  if [ ! -d "${ichor_installation_path}" ]
  then
    echo "ERROR - provided path to ichorCNA directory directory '${ichor_installation_path}' either is inaccessible or does not exist! Exiting.."
    exit 2
  fi
  use_scaff_mem_mb=$(expr $4 \* $MB_PER_SCAFFOLD)
  use_total_mem_mb=$(expr $use_scaff_mem_mb \+ $MB_FOR_INSERT_SIZES)

  echo "will use ${use_total_mem_mb} MB in total per sample!"
  for sample_name in "${NAMES[@]}"; do  # WARNING - NAMES is a global array!
    #echo 'This is $i: ' "$i"
    #echo "would call now:
    echo "INFO - putting ichor job for sample '${sample_name}' into the slurm queue.."
    #arg1 = base_unique_sample_name
    #arg2 = input directory path
    #arg3 = /PATH/TO/ichorCNAdir
    #arg4 = number of processes to use
    #arg5 = path to parallelized mark duplicates python script
    timestamp=$(date +"%Y-%m-%d_%H-%M-%S")
    slurm_log_file="${input_directory_path}/${sample_name}_${timestamp}_SLURM.log"
    slurm_error_log="${input_directory_path}/${sample_name}_${timestamp}_SLURM.err"
    sbatch  --partition=cpu --reservation=HumanGenetik --ntasks=1 --cpus-per-task=$4 --mem=$use_total_mem_mb --output="${slurm_log_file}" --error="${slurm_error_log}" "${single_sample_ichor_script}" "${sample_name}" "${input_directory_path}" "${ichor_installation_path}" $4 "${MARK_DUPLICATES_ABSOLUTE_PATH}" "${CUSTOM_TSO500_PON_PATH}"
    # bash "${single_sample_ichor_script}" "${sample_name}" "${input_directory_path}" "${ichor_installation_path}" $4 "${MARK_DUPLICATES_ABSOLUTE_PATH}" "${CUSTOM_TSO500_PON_PATH}"  # bash call for non-slurm-managed environments
  done
}

#--------------------------------------------------------------------------------------------------
#***************************************MAIN*******************************************************
#--------------------------------------------------------------------------------------------------
POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -idir|--input_dir)
    INPUT_DIR="$2"
    shift # past argument
    shift # past value
    ;;
    -phmm|--path_to_hmm)
    PATH_TO_HMM="$2"
    shift # past argument
    shift # past value
    ;;
    -pichor|--path_to_ichorCNA)
    PATH_TO_ICHOR="$2"
    shift # past argument
    shift # past value
    ;;
    -h|--help)
    echo ""
    echo "------------"
    echo "Description:"
    echo "------------"
    echo "	"
    echo "	"
    echo "	"
    echo ""
    echo "------"
    echo "Usage:"
    echo "------"
    echo "	bash TSO500_ichorCNA_from_offANDontarget_reads_PARALLEL.sh -idir path_to_input_dir -phmm path_to_hmm"
    echo "	-pichor path_to_ichorcna"
    echo ""
    echo "----------"
    echo "Arguments:"
    echo "----------"
    echo "	-idir or --input_dir: path to the directory containing paired ended fastq files to be analyzed"
    echo "	-phmm or --path_to_hmm: path to the hmm"
    echo "	-gref or --reference_genome: path to the genome. Need to be indexed"
    echo "	-pichor or --path_to_ichorCNA: path to the ichorCNA directory cloned by the ichorCNA repo"
    exit 1
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done

#adding HMM and anaconda to the path (if anaconda is not needed just uncomment it)
export PATH=$PATH:$PATH_TO_HMM  #/usr/local/bin/hmmcopy_utils/bin  #"/PATH/TO/HMMCOPY"

echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "+                You are running parallelized version of the TSO500 off-target read-based ichorCNA pipeline            +"
echo "+  (adapted from https://github.com/broadinstitute/ichorCNA by Benjamin Spiegl, Raul Mejia Pedroza and Isaac Lazzeri)  +"
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo ""
echo " --> Specified input directory: $INPUT_DIR"
echo " --> Specified ichorCNA path was: $PATH_TO_ICHOR"
echo ""

cd "${INPUT_DIR}" || (echo "ERROR - failed to cd into '${INPUT_DIR}'. Exiting.."; exit 1)

# extract unique names of the samples + make sure the index files exist! (skip directories named after bam files)
NAMES=()  # empty array to be filled! used by function 'apply_to_each_bam'
shopt -s extglob  # extglob is not enabled by default in Bash scripts! used by pattern below:
for bam_file in !(*.markdup).bam
do
  if [ -f "$bam_file" ]
  then
    expected_index_file="$bam_file".bai
    if [ ! -f "${expected_index_file}" ]
    then
      echo "WARNING - index file for BAM file '$bam_file' not found in input directory! Will create the index file now.."
      samtools index --bai --output "${expected_index_file}" "$bam_file"
    fi
    NAMES+=("${bam_file%.bam}")  # used by function 'apply_to_each_bam'
  else
    echo "Skipping directory '$bam_file'.."
  fi
done

echo "INFO - BAM files of these samples will be processed now in parallel by the TSO500 off-target reads ichorCNA pipeline:"
printf "  + %s\n" "${NAMES[@]}"


#execute the pipline for each sample
apply_to_each_bam "${SINGLE_SAMPLE_TSO500_OFFTARGET_PIPELINE_SCRIPT}" "${INPUT_DIR}" "${PATH_TO_ICHOR}" ${N_CORES_MAX}
#arg1 = single sample ichorCNA TSO500 offtarget reads script path
#arg2 = input directory path
#arg3 = /PATH/TO/ichorCNAdir
#arg4 = number of processes to use
timestamp_now=$(date +"%d.%m.%Y, at %H:%M:%S")
echo "Done enqueuing sbatch jobs for TSO500 liquid samples for ichorCNA analysis at: ${timestamp_now}"
