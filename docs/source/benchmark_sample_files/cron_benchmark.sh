#!/bin/bash

#########################################################################################################
###                This script will be called via a cron task to execute benchmark.py                 ###
#########################################################################################################

####################################
#       Script configuration       #
####################################

# Define the path to the directory where the benchmark CSV files will be created

benchmark_path="/home/users/n/i/niacobel/abin_docs_sample/benchmark"

# Load your Python distribution

module --force purge
module load releases/2018b
module load Python/3.6.6-foss-2018b

####################################
#         Preparation step         #
####################################

# Command line arguments

prefix=$1

# Define the timestamp

timestamp=$(date +"%Y%m%d_%H%M%S")

# Pretty print for log messages

log_msg () {
  echo -e "$(date +"%Y-%m-%d %T")\t$1"
}

# Define the tmp file we want to scan

WATCH_FILE="${benchmark_path}/${prefix}_tmp.csv"

# Define the folder where the tmp file will be archived

ARCHIVE="${benchmark_path}/archive"

# Define the folder where the log files will be stored

BENCH_LOGS="${benchmark_path}/bench_logs"

# Define the path towards the benchmark.py script (same directory as this script)

path_script=$(dirname "$(readlink -f "$0")")

####################################
#         Start of execution       #
####################################

# Exit immediately if there's no file to process

if [ ! -f "${WATCH_FILE}" ]; then
  exit

# Otherwise execute benchmark.py

else

  # Archive the original tmp file

  filename="$(basename -- ${WATCH_FILE})"
  mkdir -p ${ARCHIVE}
  mv ${WATCH_FILE} ${ARCHIVE}/${filename%.*}_${timestamp}.csv

  # Execute benchmark.py

  mkdir -p ${BENCH_LOGS}
  python ${path_script}/benchmark.py --tmp ${ARCHIVE}/${filename%.*}_${timestamp}.csv --final ${benchmark_path}/${prefix}_final.csv > ${BENCH_LOGS}/${prefix}_${timestamp}.log

  log_msg "INFO - Processed new lines in ${WATCH_FILE}"

fi