#!/bin/bash

#########################################################################################################
###      This script will be called via a cron task to execute control_launcher.py (with QOCT-RA)     ###
#########################################################################################################

####################################
#       Script configuration       #
####################################

# Load your Python distribution (including YAML and Jinja2)

source "/CECI/home/ulb/cqp/niacobel/CHAINS/load_modules.sh"

####################################
#         Preparation step         #
####################################

# Command line arguments

cluster_name=$1
out_dir=$2

# Define the path towards CHAINS' directory (parent directory from the directory of this script)

script_dir=$(dirname "$(readlink -f "$0")")
chains_path=$(dirname "$(readlink -f "${script_dir}")")

# Define the path towards CHAINS' configuration file

chains_config="${chains_path}/configs/chains_config.yml"

# Create the function allowing us to read YAML files (this is done through Python, see https://stackoverflow.com/a/47791935/14608112)

yaml() {
    python3 -c "import yaml;print(yaml.safe_load(open('$1'))$2)"
}

# Define the directory containing the QCHEM output files

WATCH_DIR=$(yaml "${chains_config}" "['output_qchem']")

# Define the extension of the QCHEM output files

OUT_FILEPATH="${WATCH_DIR}/*.out"

# Define the path towards the CONTROL LAUNCHER directory

control_dir="${chains_path}/control_launcher"

# Define the directory where the log files will be stored

control_logs="${WATCH_DIR}/control_logs"

####################################
#         Start of execution       #
####################################

# Exit immediately if there's no file to process

if [ $(ls ${OUT_FILEPATH} 2>/dev/null | wc -l) -eq 0 ]; then
  exit

# Otherwise execute control_launcher.py for each file present in the WATCH_DIR directory

else

  file_list=$(ls ${OUT_FILEPATH} 2>/dev/null)
  mkdir -p "${control_logs}"

  for filepath in ${file_list}
  do
    filename="$(basename -- "${filepath}")"
    MOL_NAME=${filename%.*}
    mkdir -p "${out_dir}"
    python "${control_dir}/control_launcher.py" -p qoctra -s "${filepath}" -cf "${chains_path}/configs/qoctra/opc.yml" -o "${out_dir}" -cl "${cluster_name}" -ow -as > "${control_logs}/$(date +"%Y%m%d_%H%M%S")_${MOL_NAME}_OPC.log"
  done

  echo -e "$(date +"%Y-%m-%d %T")\tINFO - Successfully processed:\n${file_list}"

fi
