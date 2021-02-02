################################################################################################################################################
##                                                   Gaussian Output Quality Control Script                                                   ##
##                                                                                                                                            ##
##                            This script scans a given gaussian log file looking for possible errors and warnings                            ##
##                                           that might interfere with the rest of the calculations                                           ##
################################################################################################################################################

import sys

with open(sys.argv[1], 'r') as output:
  output_content = output.read().splitlines()

print("\nChecking GAUSSIAN log file...")

# Check if the second to last line corresponds to what we should have
if "Normal termination of Gaussian" not in output_content[-1]:
  print("\nERROR: GAUSSIAN's log file does not mention having terminated normally. Check the log file content before proceeding further.")
  print("Aborting ...")
  exit(1)

print("No errors detected in the log file.")