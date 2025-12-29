#!/usr/bin/env bash

set -e -o pipefail

hostname

# ACTIVATE CONDA ENVIRONMENT
eval "$(conda shell.bash hook)"
conda activate TSO500ichorCNA

# CHECK COMANDLINE ARGUMENTS:
#arg1 = base_unique_sample_name
#arg2 = input directory path
#arg3 = /PATH/TO/ichorCNAdir
#arg4 = number of processes to use
#arg5 = path to parallelized mark duplicates python script
#arg6 = path to TSO500 liquid-specific PoN file
# Check that exactly 5 arguments are provided
if [[ $# -ne 6 ]]; then
    echo "Usage: $0 <base_unique_sample_name> <input_dir> <path_to_ichorCNA_installation> <threads> <mark_duplicates_python_script_path> <PoN_file_path>"
    exit 1
fi

DEBUGGING=false

# check passed arguments:
CUSTOM_TSO500_PON_PATH="${6}"
if [ ! -f "${CUSTOM_TSO500_PON_PATH}" ]
then
  echo "ERROR - custom PoN file not found/not accessible at '${CUSTOM_TSO500_PON_PATH}'. Exiting .."
  exit 1
fi
REF_BUILD_CENTROMERE_GAP_TABLE="${3}"/inst/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt
if [ ! -f "${REF_BUILD_CENTROMERE_GAP_TABLE}" ]
then
  echo "ERROR - reference build centromere gap table not found/not accessible at '${REF_BUILD_CENTROMERE_GAP_TABLE}'. Exiting .."
  exit 1
fi
REF_GC_WIG_FILE="${3}/inst/extdata/gc_hg19_1000kb.wig"
if [ ! -f "${REF_GC_WIG_FILE}" ]
then
  echo "ERROR - reference build GC content WIG file not found/not accessible at '${REF_GC_WIG_FILE}'. Exiting .."
  exit 1
fi
REF_MAPPING_WIG_FILE="${3}/inst/extdata/map_hg19_1000kb.wig"
if [ ! -f "${REF_MAPPING_WIG_FILE}" ]
then
  echo "ERROR - reference build mapping WIG file not found/not accessible at '${REF_MAPPING_WIG_FILE}'. Exiting .."
  exit 1
fi
ICHOR_CNA_SCRIPT_PATH="${3}/scripts/runIchorCNA.R"
if [ ! -f "${ICHOR_CNA_SCRIPT_PATH}" ]
then
  echo "ERROR - ichorCNA R script not found/not accessible at '${ICHOR_CNA_SCRIPT_PATH}'. Exiting .."
  exit 1
fi

#--------------------------------------------------------------------------------------------------
#*************************************FUNCTIONS****************************************************
#--------------------------------------------------------------------------------------------------

function pipeline(){
  #arg1 = base_unique_sample_name
  sample_name="${1}"
  #arg2 = input directory path
  input_directory_path="${2}"
  #arg3 = /PATH/TO/ichorCNAdir
  ichor_installation_path="${3}"
  #arg4 = number of processes to use
  proc_number=${4}
  #arg5 = path to parallelized mark duplicates python script
  markdup_script_path="${5}"

  # markduplicates script uses now MarkDuplicatesWithMateCigar with patterned flowcell-specific parameters
  # for parallelized duplicates marking!
  # (python3.13 on sy096 in TSO500ichorCNA conda environment)

  cd "${input_directory_path}" || (echo "ERROR - failed to cd into '${input_directory_path}'. Exiting.."; exit 1)
  user_temp_dir="${TEMP}"
  if [ ! -d "${user_temp_dir}" ]
  then
    # try to use home directory
    user_temp_dir=~
    if [ ! -d "${user_temp_dir}" ]
    then
      # use alternative on isilon (works only on MedBioNode cluster)
      user_temp_dir='/home/isilon/HumGenAnalysis'
      if [ -d "${user_temp_dir}" ]
      then
        user_temp_dir="${user_temp_dir}/Oncoservice/ichorCNA_analyses/tmp"
        mkdir -p "${user_temp_dir}"
      fi
    else
      user_temp_dir="${user_temp_dir}/TSO500ichor_tmp"
      mkdir -p "${user_temp_dir}"
    fi
  fi
  # final check if variable is properly set
  if [ ! -d "${user_temp_dir}" ]; then
    echo "UNEXPECTED ERROR - user_temp_dir could not be set. Exiting..."
    exit 1
  fi
  echo "INFO - will use temporary directory '${user_temp_dir}'"
  markdup_temp_dir="${user_temp_dir}/markDuplicates-${sample_name}"  # THIS MUST BE SAMPLE-SPECIFIC!
  mkdir -p "${markdup_temp_dir}"
  trap "echo 'cleaning up temporary directory for marking duplicates..' && rm -rdf ${markdup_temp_dir}; conda deactivate" SIGINT SIGTERM
    # Skip markduplicates if already exists
  if [ -f "${input_directory_path}/${sample_name}.markdup.bam" ]; then
    echo "STEP1: MarkDuplicates already completed for ${sample_name}, skipping..."
  else
    echo "***************************************************************"
    echo "  STEP1: Running MarkDuplicates for sample: ${sample_name}     "
    echo "***************************************************************"
    python3 "${markdup_script_path}" --input-bam "${input_directory_path}/${sample_name}.bam" --output-dir "${input_directory_path}" --processes ${proc_number} --temp-dir "${markdup_temp_dir}"
  fi

  # delete temporary directory for sample if it still exists after marking duplicates with picard
  if [ -d "${markdup_temp_dir}" ]
  then
    rm -rf "${markdup_temp_dir}"
  fi
  ichor "${sample_name}" "${input_directory_path}" "${ichor_installation_path}"
}

function ichor(){
  #arg1 = base_unique_sample_nam
  #arg2 = parent path to input BAM file
  #arg3 = /PATH/TO/ichorCNAdir
  if [ -z "$1" ]; then
    echo "Error: Missing sample ID (first argument)."
    exit 1
  fi
  # Step 2: run readCounter
  if [ -f "${2}/${1}.wig" ]; then
    echo "STEP2: Already done"
  else
    echo "***************************************************************"
    echo "  STEP2: Running readCounter now for sample: ${1}              "
    echo "***************************************************************"
    readCounter --window 1000000 --quality 20 \
    --chromosome "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY" \
    "${2}/${1}.markdup.bam" > "${2}/${1}.wig"
  fi
  # Step 3: actually run ichorCNA now!
  if [ -d "${2}/Results_ichorCNA_${1}" ]; then
    echo "STEP3: readCounter was already run."
  else
    echo "***************************************************************"
    echo "  STEP3: Running ichorCNA now for sample: ${1}                 "
    echo "***************************************************************"
    mkdir -p "${2}/Results_ichorCNA_${1}"  # output directory (absolute path)
    Rscript "${ICHOR_CNA_SCRIPT_PATH}" \
    --WIG "${2}/${1}.wig" \
    --gcWig "${REF_GC_WIG_FILE}" \
    --mapWig "${REF_MAPPING_WIG_FILE}" \
    --ploidy "c(2,3,4)" \
    --normal "c(0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95)" \
    --maxCN 7 \
    --id "${1}" \
    --estimateNormal TRUE \
    --estimatePloidy TRUE \
    --includeHOMD FALSE \
    --chrs "c(1:22,\"X\")" \
    --chrTrain "c(1:18)" \
    --centromere "${REF_BUILD_CENTROMERE_GAP_TABLE}" \
    --normalPanel "${CUSTOM_TSO500_PON_PATH}" \
    --chrTrain "c(1:18)" \
    --txnE 0.999 \
    --txnStrength 100000 \
    --scStates "c(1,3)" \
    --estimateScPrevalence TRUE \
    --maxFracGenomeSubclone 0.5 \
    --maxFracCNASubclone 0.7 \
    --minSegmentBins 50 \
    --altFracThreshold 0.7 \
    --lambdaScaleHyperParam 3 \
    --libdir "${3}" \
    --outDir "${2}/Results_ichorCNA_${1}"
  fi
}

#--------------------------------------------------------------------------------------------------
#***************************************MAIN*******************************************************
#--------------------------------------------------------------------------------------------------

sample_name="${1}"
# log node we are on for debugging purposes
echo "INFO - running analysis for sample '${sample_name}' on node '$(hostname)'"

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
number_threads=${4}
markduplicates_script="${5}"
if [ ! -f "${markduplicates_script}" ]
then
  echo "ERROR - parallelized mark duplicates Python3 script '${markduplicates_script}' not found in current working directory '$(pwd)'! Exiting.."
  exit 2
else
  if [ "$DEBUGGING" = "true" ]; then
    echo "this is the value of the variable markduplicates_script in 'TSO500_ichorCNA_singleSample_calledTarget.sh': ${markduplicates_script}"  # for debugging
  fi
fi

pipeline "${sample_name}" "${input_directory_path}" "${ichor_installation_path}" $number_threads "${markduplicates_script}"
#arg1 = base_unique_sample_name
#arg2 = input directory path
#arg3 = /PATH/TO/ichorCNAdir
#arg4 = number of processes to use
#arg5 = path to parallelized mark duplicates python script
timestamp_now=$(date +"%d.%m.%Y, at %H:%M:%S")
echo "Done processing TSO500 liquid sample '${sample_name}' for ichorCNA analysis at: ${timestamp_now}"
conda deactivate
#--------------------------------------------------------------------------------------------------