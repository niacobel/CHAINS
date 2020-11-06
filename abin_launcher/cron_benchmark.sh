#!/bin/bash

#########################################################################################################
###                This script will be called via a cron task to execute benchmark.py                 ###
#########################################################################################################

####################################
#       Script configuration       #
####################################

# Define the name of the directory where the benchmark CSV files will be created

benchmark_dir="/CECI/home/ulb/cqp/niacobel/BENCHMARK"

# Load your Python distribution

source /CECI/home/ulb/cqp/niacobel/CHAINS/load_modules.sh

####################################
#         Preparation step         #
####################################

# Command line arguments

PROGRAM=$1
CLUSTER_NAME=$2

# Define the timestamp

timestamp=$(date +"%Y%m%d_%H%M%S")

# Pretty print for log messages

log_msg () {
  echo -e "$(date +"%Y-%m-%d %T")\t$1"
}

# Define the tmp file we want to scan

WATCH_FILE="${benchmark_dir}/${PROGRAM}_${CLUSTER_NAME}_tmp.csv"

# Define the folder where the tmp file will be archived

ARCHIVE="${benchmark_dir}/archive"

# Define the folder where the log files will be stored

BENCH_LOGS="${benchmark_dir}/bench_logs"

# Define the path towards the benchmark.py script

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
  python ${path_script}/benchmark.py --tmp ${ARCHIVE}/${filename%.*}_${timestamp}.csv --final ${benchmark_dir}/${PROGRAM}_${CLUSTER_NAME}.csv > ${BENCH_LOGS}/${PROGRAM}_${CLUSTER_NAME}_${timestamp}.log

  log_msg "INFO - Processed new lines in ${WATCH_FILE}"

fi