#!/usr/bin/bash

set -e -o pipefail  # exit as soon as a command fails

# info - conda env is activated in the called target script 'drv_TSO500_offtarget_ichorCNA_parallel.sh'
#        the sbatch inside the parallelization runs the script on any CPU partition node

while [[ $# -gt 0 ]]; do
    case "$1" in
        -d|--indir)
            if [[ -z "$2" || "$2" == -* ]]; then
                echo "Error: --indir requires a non-empty argument."
                exit 1
            fi
            input_path="$2"
            shift 2
            ;;
        -p|--ichorpath)
            if [[ -z "$2" || "$2" == -* ]]; then
                echo "Error: --ichorpath requires a non-empty argument."
                exit 1
            fi
            path_to_ichor="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 --indir input --ichorpath path_to_ichor_folder"
            exit 0
            ;;
        *)
            # default commands
            echo "Unknown parameter passed: $1"
            exit 1
            ;;
    esac
    # shift  # DEACTIVATED because:
    #        For flags that take an argument, you need to shift twice — once for the flag itself,
    #        once for its argument. Your script shifts only once inside the case and once outside,
    #        which is correct if and only if the argument exists and is in $2. But if $2 is missing
    #        (i.e., the flag is the last argument or missing its value), your script will silently
    #        assign an empty value or even crash. It’s safer to check if the argument exists and shift
    #        correctly.
done

# check if parameters are available
if [[ -z "$input_path" || -z "$path_to_ichor" ]]; then
    echo "❌ Error: --indir and --ichorpath are required."  # e.g., /home/gpfs/o_spiegl/installations/ichorCNA
    echo "Usage: $0 -d/--indir <file> -p/--ichorpath <file>"
    exit 1
fi

echo "Input directory supposedly containing TSO500 ctDNA sample BAM files inside: $input_path"
echo "Directory contianing the ichorCNA installation: $path_to_ichor"

# check parameters
if [ ! -d "${input_path}" ]
then
  echo "ERROR - input directory '${input_path}' does not exist"
  exit 1
fi
if [ ! -d "${path_to_ichor}" ]
then
  echo "ERROR - expected ichorCNA installation does not exist in directory '${path_to_ichor}'"
  exit 1
fi

script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
path_hmm="/usr/bin"
# Raul's reference genome was not available

# RUN THE THANG!
# (analysis of everything with BAM in its name from input directory; sample-by-sample)
bash "${script_dir}"/TSO500_ichorCNA_from_offANDontarget_reads_PARALLEL.sh -idir "${input_path}" -phmm "${path_hmm}" -pichor "${path_to_ichor}" 2>&1 | tee "${input_path}"/TSO500_liquid_offtarget_ichorCNA.log
