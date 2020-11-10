#!/bin/bash

#########################################################################################################
###                This script will be called via a cron task to execute benchmark.py                 ###
#########################################################################################################

####################################
#       Script configuration       #
####################################

# Load your Python distribution

source /CECI/home/ulb/cqp/niacobel/CHAINS/load_modules.sh

####################################
#         Preparation step         #
####################################

# Command line arguments

prefix=$1
benchmark_path=$2

# Define the timestamp

timestamp=$(date +"%Y%m%d_%H%M%S")

# Define the temporary CSV file we want to scan

WATCH_FILE="${benchmark_path}/${prefix}_tmp.csv"

# Define the directory where the temporary file will be archived

archive="${benchmark_path}/archive"

# Define the directory where the log files will be stored

bench_logs="${benchmark_path}/bench_logs"

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

  # Archive the temporary file

  filename="$(basename -- ${WATCH_FILE})"
  mkdir -p ${archive}
  mv ${WATCH_FILE} ${archive}/${filename%.*}_${timestamp}.csv

  # Execute benchmark.py

  mkdir -p ${bench_logs}
  python ${path_script}/benchmark.py --tmp ${archive}/${filename%.*}_${timestamp}.csv --final ${benchmark_path}/${prefix}_final.csv --prob ${benchmark_path}/${prefix}_prob.csv > ${bench_logs}/${prefix}_${timestamp}.log

  echo -e "$(date +"%Y-%m-%d %T")\tINFO - Processed new lines in ${WATCH_FILE}"

fi