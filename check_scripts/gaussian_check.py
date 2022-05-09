#!/usr/bin/env python3
 
################################################################################################################################################
##                                                   Gaussian Output Quality Control Script                                                   ##
##                                                                                                                                            ##
##                            This script scans a given gaussian log file looking for possible errors and warnings                            ##
##                                           that might interfere with the rest of the calculations                                           ##
##                                                                                                                                            ##
##      Possible exit status: [0] No errors detected - [1] Errors detected, termination required - [2] Errors handled, relaunch required      ##
################################################################################################################################################

import os
import sys
import time

# =================================================================== #
# =================================================================== #
#                        FUNCTIONS DEFINITIONS                        #
# =================================================================== #
# =================================================================== #

def update_file(file_path:str,file_content:str,archive=False):

  """Updates a file (given by path) with the new content (given by content). If archive is set to True, the old file will be renamed by adding the current time as a suffix to its basename, before the extension."""

  if not os.path.isfile(file_path):
    print("ERROR: The update_file function has been called with a path argument that does not lead to a file.")
    print("Aborting ...")
    sys.exit(1)

  print("Updating the input file ... ", end="")
  if archive:
    file_basename = os.path.splitext(os.path.basename(file_path))[0]
    file_ext = os.path.splitext(os.path.basename(file_path))[1]
    file_dirname = os.path.dirname(file_path)
    timestr = time.strftime("%Y_%m_%d-%H_%M_%S")
    os.rename(file_path,os.path.join(file_dirname,file_basename + "_" + timestr + file_ext))
  else:
    os.remove(file_path)

  with open(file_path, "w", encoding='utf-8') as new_file:
    new_file.write(file_content)
  print("%12s" % "[DONE]")

# =================================================================== #
# =================================================================== #
#                            MAIN FUNCTION                            #
# =================================================================== #
# =================================================================== #

def main():

  # Preparation step
  # ================

  # Important files

  output_file = sys.argv[1]
  input_file = os.path.splitext(output_file)[0] + ".com"

  with open(output_file, 'r') as output:
    output_content = output.read().splitlines()

  with open(input_file, 'r') as inp:
    input_content = inp.read().splitlines()

  # Copy the input content to allow altering it later during iteration

  input_copy = input_content.copy()

  # Start checking the log file
  # ===========================

  print("\nChecking GAUSSIAN log file...")

  # Normal termination of Gaussian
  if "Normal termination of Gaussian" in output_content[-1]:
    print("Normal termination of Gaussian.")
    sys.exit(0)

  # Inaccurate quadrature in CalDSu (resolution provided by https://wongzit.github.io/gaussian-common-errors-and-solutions/#l502-l1002)
  elif any("Inaccurate quadrature in CalDSu".lower() in line.lower() for line in output_content[-100:]):
    print("\nDetected the following error: Inaccurate quadrature in CalDSu.")

    # Handle the error
    for line in input_copy:
      if line.startswith("#"):
        if "Int=ultrafine" in line and "SCF=novaracc" in line:
          print("\nERROR: The keywords 'Int=ultrafine' and 'SCF=novaracc' have already been added. Check the log file content before proceeding further.")
          print("Aborting ...")
          sys.exit(1)
        else:
          idx = input_copy.index(line)
          input_content[idx] += " Int=ultrafine SCF=novaracc"

    # Update the input file
    update_file(input_file,"\n".join(input_content),archive=True)

    # Send the exit code 2 (the error has been handled and the calculation needs to be relaunched)
    sys.exit(2)
    
  # Unknown error
  else:
    print("\nERROR: GAUSSIAN's log file does not mention having terminated normally. Check the log file content before proceeding further.")
    print("Aborting ...")
    sys.exit(1)

# =================================================================== #
# =================================================================== #
#                          CALL MAIN FUNCTION                         #
# =================================================================== #
# =================================================================== #

if __name__ == "__main__":
  main()
