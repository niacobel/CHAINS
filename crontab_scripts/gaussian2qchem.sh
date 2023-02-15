#!/bin/bash

#########################################################################################################
###              This script will be called via a cron task to execute abin_launcher.py               ###
###                  (with the qchem profile, the qchem and qchem_TZVP config files)                  ###
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

# Define the directory containing the GAUSSIAN optimized geometry files

WATCH_DIR=$(yaml "${chains_config}" "['output_gaussian']")

# Define the extension of the GAUSSIAN optimized geometry files

XYZ_FILEPATH="${WATCH_DIR}/*.xyz"

# Define the path towards the ABIN LAUNCHER directory

abin_dir="${chains_path}/abin_launcher"

# Define the directory where the log files will be stored

abin_logs="${WATCH_DIR}/abin_logs"

####################################
#         Start of execution       #
####################################

# Exit immediately if there's no file to process

if [ $(ls ${XYZ_FILEPATH} 2>/dev/null | wc -l) -eq 0 ]; then
  exit

# Otherwise execute abin_launcher.py for each file present in the WATCH_DIR directory

else

  file_list=$(ls ${XYZ_FILEPATH} 2>/dev/null)
  mkdir -p "${abin_logs}"

  for filepath in ${file_list}
  do
    filename="$(basename -- "${filepath}")"
    MOL_NAME=${filename%.*}
    mkdir -p "${out_dir}"
    python "${abin_dir}/abin_launcher.py" -p qchem -m "${filepath}" -cf "${chains_path}/configs/qchem/SVP.yml" -o "${out_dir}" -cl "${cluster_name}" -ow -kc > "${abin_logs}/$(date +"%Y%m%d_%H%M%S")_${MOL_NAME}.log"

  done

  echo -e "$(date +"%Y-%m-%d %T")\tINFO - Successfully processed:\n${file_list}"

fi
