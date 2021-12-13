#!/bin/bash

#########################################################################################################
###    This script will be called via a cron task to execute the various results treatment scripts    ###
#########################################################################################################

####################################
#       Script configuration       #
####################################

# Load your Python distribution (including YAML and Jinja2)

source "/CECI/home/ulb/cqp/niacobel/CHAINS/load_modules.sh"

####################################
#         Preparation step         #
####################################

# Define the path towards CHAINS' directory (parent directory from the directory of this script)

script_dir=$(dirname "$(readlink -f "$0")")
chains_path=$(dirname "$(readlink -f "${script_dir}")")

# Define the path towards CHAINS' configuration file

chains_config="${chains_path}/configs/chains_config.yml"

# Create the function allowing us to read YAML files (this is done through Python, see https://stackoverflow.com/a/47791935/14608112)

yaml() {
    python3 -c "import yaml;print(yaml.safe_load(open('$1'))$2)"
}

# Define the directory containing the finish files

WATCH_DIR=$(yaml "${chains_config}" "['output_qoctra']")

# Define the extension of the finish files

END_FILEPATH="${WATCH_DIR}/*.end"

# Define the directory where the results are stored

results_dir=$(yaml "${chains_config}" "['results_dir']")

# Define the path towards the RESULTS TREATMENT directory

treatment_dir="${chains_path}/results_treatment"

# Define the directory where the log files will be stored

treatment_logs="${WATCH_DIR}/treatment_logs"

# Define the directory where the treated results will be stored

treated_dir=$(yaml "${chains_config}" "['treated_dir']")

####################################
#         Start of execution       #
####################################

# Exit immediately if there's no file to process

if [ $(ls ${END_FILEPATH} 2>/dev/null | wc -l) -eq 0 ]; then
  exit

# Otherwise execute the various results treatment scripts for each file present in the WATCH_DIR directory

else

  file_list=$(ls ${END_FILEPATH} 2>/dev/null)
  mkdir -p "${treatment_logs}"
  mkdir -p "${WATCH_DIR}/treated"
  
  for filepath in ${file_list}
  do
    filename="$(basename -- "${filepath}")"
    MOL_NAME=${filename%_opc*.end}
    logfile="${treatment_logs}/$(date +"%Y%m%d_%H%M%S")_${filename%.end}.log"
    mkdir -p "${treated_dir}/pulses_plots"
    python "${treatment_dir}/plot_pulses.py" -s "${results_dir}/${MOL_NAME}" -o "${treated_dir}/pulses_plots" -qt 0.1 > ${logfile}
    python "${treatment_dir}/results_treatment.py" -s "${results_dir}/${MOL_NAME}" -o "${treated_dir}/comp_results.yml" >> ${logfile}
    mv "${filepath}" "${WATCH_DIR}/treated"
  done

  #python "${treatment_dir}/plot_results.py" -i "${treated_dir}/comp_results.yml" -o "${treated_dir}" > "${treatment_logs}/$(date +"%Y%m%d_%H%M%S")_plot_results.log"

  echo -e "$(date +"%Y-%m-%d %T")\tINFO - Successfully processed:\n${file_list}"

fi
