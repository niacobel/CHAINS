#!/usr/bin/env python3

########################################################################################################################################################
##                                                       Constraints variation results treatment                                                      ##
##                                                                                                                                                    ##
##                     This script compiles the results of the constraints variation, then identifies the most interesting combos                     ##
##                          and copies their associated pulse and PCP files from the GLOBALSCRATCH to the target directory.                           ##
##                                                                                                                                                    ##
##                                      Extended documentation is available at https://chains-ulb.readthedocs.io/                                     ##
########################################################################################################################################################

import argparse
import csv
import os
import re
import shutil
from statistics import mean

# =================================================================== #
# =================================================================== #
#                        EXCEPTIONS DEFINITIONS                       #
# =================================================================== #
# =================================================================== #

class Error(Exception):
    """Base class for exceptions in this script."""
    pass

class CV_Error(Error):
    """Exception raised for errors specific to certain instructions in this script.

    Attributes
    ----------
    message : str
        Proper error message explaining the error.
    """

    def __init__(self, message):
        self.message = message

# =================================================================== #
# =================================================================== #
#                        FUNCTIONS DEFINITIONS                        #
# =================================================================== #
# =================================================================== #

def check_abspath(path:str,context:str,type="either"):
    """Checks if a path towards a file or directory exists and is of the correct type. If it is, the function returns its absolute version.

    Parameters
    ----------
    path : str
        The path towards the file or directory you want to test.
    context : str
        Message to show on screen to give more information in case of an exception (e.g. the role of the directory or file that was checked, where the checked path was given, etc.).
    type : str, optional
        The type of element for which you would like to test the path ('file', 'directory' or 'either').
        By default, checks if the path leads to either a file or a directory (type = 'either').
    
    Returns
    -------
    abspath : str
        Normalized absolute version of the path.

    Raises
    ------
    ValueError
        If the specified type when calling the function is not 'file', 'directory' or 'either'.
    AD_Error
        If the type does not match what is given in the path, or if the path does not exist.
    """

    # Check "type" argument

    if type not in ["file","directory","either"]:
      raise ValueError ("The specified type for which the check_abspath function has been called is not one of 'file', 'directory' or 'either'")

    # Prepare to print a helpful error message in case of problem with the given path

    msg = "\nSomething went wrong when checking the path " + path + "\nContext: " + context + "\n"

    # Check path
    
    if not os.path.exists(path):
      raise CV_Error (msg + "ERROR: %s does not seem to exist." % path)
    elif type == "file":
      if not os.path.isfile(path):
        raise CV_Error (msg + "ERROR: %s is not a file" % path)
    elif type == "directory":
      if not os.path.isdir(path):
        raise CV_Error (msg + "ERROR: %s is not a directory" % path)
    elif type == "either":
      if not os.path.isdir(path) and not os.path.isfile(path):
        raise CV_Error (msg + "ERROR: %s is neither a file nor a directory" % path)

    # If everything went well, get the normalized absolute version of the path
    
    abspath = os.path.abspath(path)

    return abspath

# =================================================================== #
# =================================================================== #
#                       COMMAND LINE ARGUMENTS                        #
# =================================================================== #
# =================================================================== #

# Define the arguments needed for the script (here they are defined as named arguments rather than positional arguments, check https://stackoverflow.com/questions/24180527/argparse-required-arguments-listed-under-optional-arguments for more info).

parser = argparse.ArgumentParser(add_help=False, description="This script compiles the results of the constraints variation, then identifies the most interesting combos and copies their associated pulse and PCP files from the GLOBALSCRATCH to the target directory.")

required = parser.add_argument_group('Required arguments')
required.add_argument("-r","--results_dir", type=str, help="Path to the directory containing the results of the constraints variation.", required=True)
required.add_argument('-o', '--out_dir', type=str, help="Path to the directory where the compiled results and the most interesting combos files will be stored.", required=True)

optional = parser.add_argument_group('Optional arguments')
optional.add_argument('-h','--help',action='help',default=argparse.SUPPRESS,help='Show this help message and exit')

# =================================================================== #
# =================================================================== #
#                            MAIN FUNCTION                            #
# =================================================================== #
# =================================================================== #

def main(): 

  # =================================================================== #
  # =================================================================== #
  #                           PREPARATION STEP                          #
  # =================================================================== #
  # =================================================================== #
  
  try:
    
    # Get the size of the terminal in order to have a prettier output, if you need something more robust, go check http://granitosaurus.rocks/getting-terminal-size.html

    columns, rows = shutil.get_terminal_size()

    # Output Header

    print("".center(columns,"*"))
    print("")
    print("ALPHA/DURATION PARAMETERS SEARCH RESULTS TREATMENT BEGINS NOW".center(columns))
    print("")
    print("".center(columns,"*"))

    # ========================================================= #
    # Read command line arguments                               #
    # ========================================================= #

    args = parser.parse_args()

    # Required arguments

    results_dir = args.results_dir           # Path to the directory containing the results of the alpha-duration parameters search
    out_dir = args.out_dir                   # Path to the directory where the compiled results and the best parameters combo files will be stored

    # ========================================================= #
    # Check arguments                                           #
    # ========================================================= #

    # Check the existence of the directories passed as arguments

    results_dir = check_abspath(results_dir,"Command line argument -r / --results_dir","directory")
    print ("{:<40} {:<100}".format('\nResults directory:',results_dir))

    out_dir = check_abspath(out_dir,"Command line argument -o / --out_dir","directory")
    print ("{:<40} {:<100}".format('\nOutput directory:',out_dir))

  # ========================================================= #
  # Exception handling for the preparation step               #
  # ========================================================= #

  except CV_Error as error:
    print("")
    print(error)
    exit(-1)

  # =================================================================== #
  # =================================================================== #
  #                         RESULTS COMPILATION                         #
  # =================================================================== #
  # =================================================================== #

  try:

    print ("{:<40}".format('\nCompiling the results ...'), end="")

    # Initialize the compiled results list of dictionaries

    comp_results = []

    # Look for directories in the results directory (see https://stackoverflow.com/questions/800197/how-to-get-all-of-the-immediate-subdirectories-in-python for reference).

    dir_list_all = [dir.name for dir in os.scandir(results_dir) if dir.is_dir()]

    # Only keep the 'xxx_alphaXXX_durXXX'-type directories

    dir_list = []

    for dirname in dir_list_all:
 
      # Define the 'xxx_fluX_windX' regex and apply it to the dirname

      pattern = re.compile(r"^\d+_(?P<suffix>flu\d+_wind\d+)$")
      matching_dir = pattern.match(dirname)

      # If it matches the regex, collect the data

      if matching_dir is not None:

        dir_list.append(dirname)
        suffix = matching_dir.group("suffix")

        param_file_path = os.path.join(results_dir, dirname, "param_%s.nml" % suffix)

        with open(param_file_path, 'r') as param_file:
          param_content = param_file.read().splitlines()

        fluence_pattern = re.compile(r'^\s+fluence\s+=\s+(?P<fluence>\d+\.\d+d[+-]\d+)')
        window_pattern = re.compile(r'^\s+spectral_filter_fwhm\s+=\s+(?P<window>\d+\.\d+d[+-]\d+)')

        for line in param_content:

          if fluence_pattern.match(line): 
            fluence = float((fluence_pattern.match(line).group('fluence')).replace('d','e'))
         
          elif window_pattern.match(line): 
            window = float((window_pattern.match(line).group('window')).replace('d','e'))

        # Prepare to store information specific to this calculation

        data = {
            "Directory" : dirname,
            "Fluence" : fluence,
            "Window (cm-1)" : window
        }

        # Iterate over all the PCP directories

        pcp_dirs = [dir.name for dir in os.scandir(os.path.join(results_dir, dirname)) if dir.is_dir() and (dir.name).startswith("PCP_")]
        orientations = [] # List of all possible orientations

        for pcp_dir in pcp_dirs:

          pcp_file_path = os.path.join(results_dir, dirname, pcp_dir, "obj.res")
          orientation = pcp_dir.partition("PCP_")[2]

          with open(pcp_file_path, 'r') as pcp_f:
            file_content = pcp_f.read()

          # Define the expression patterns for the lines of the iterations file
          # For example "      0     1  1sec |Proba_moy  0.693654D-04 |Fidelity(U)  0.912611D-01 |Chp  0.531396D-04 -0.531399D-04 |Aire -0.202724D-03 |Fluence  0.119552D-03 |Recou(i)  0.693654D-04 |Tr_dist(i) -0.384547D-15 |Tr(rho)(i)  0.100000D+01 |Tr(rho^2)(i)  0.983481D+00 |Projector  0.100000D+01"
          rx_pcp_line = re.compile(r"^\s+\d+\s+\d+\s+\d+(?:sec|min)\s\|Proba_moy\s+\d\.\d+D[+-]\d+\s\|Fidelity\(U\)\s+(?P<fidelity>\d\.\d+D[+-]\d+)\s\|Chp\s+\d\.\d+D[+-]\d+\s+-?\d\.\d+D[+-]\d+\s\|Aire\s+-?\d\.\d+D[+-]\d+\s\|Fluence\s+\d\.\d+D[+-]\d+\s\|Recou\(i\)\s+(?P<overlap>\d\.\d+D[+-]\d+)\s\|Tr_dist\(i\)\s+-?\d\.\d+D[+-]\d+\s\|Tr\(rho\)\(i\)\s+\d\.\d+D[+-]\d+\s\|Tr\(rho\^2\)\(i\)\s+\d\.\d+D[+-]\d+\s\|Projector\s+(?P<projector>\d\.\d+D[+-]\d+)")

          # Get the values

          line_data = rx_pcp_line.match(file_content)
          if line_data is not None:
            for key in ['projector','overlap','fidelity']:

              raw_value = line_data.group(key)
              value = float(re.compile(r'(\d*\.\d*)[dD]([-+]?\d+)').sub(r'\1E\2', raw_value)) # Replace the possible d/D from Fortran double precision float format with an "E", understandable by Python)
              value = min (value, 1)
              
              data_key = orientation + "_" + key.capitalize()
              data.update({
                data_key : value
              })
            
            orientations.append(orientation) # Will be reused later

          else:
            raise CV_Error ("ERROR: Unable to get the values from the file %s" % pcp_file_path) 

        # Store information specific to this calculation

        comp_results.append(data)

    if dir_list == []:
      raise CV_Error ("ERROR: Can't find any 'xxx_alphaXXX_durXXX' directory in %s" % results_dir)

    print("[ DONE ]")

  # ========================================================= #
  # Exception handling for the results compilation            #
  # ========================================================= #

  except CV_Error as error:
    print("")
    print(error)
    exit(-1)

  # =================================================================== #
  # =================================================================== #
  #                    MOST INTERESTING COMBINATION                     #
  # =================================================================== #
  # =================================================================== #

  print ("{:<40}".format('\nIdentifying most interesting combos ...'), end="")

  # Identify the best orientation

  highest_mean = -float('inf')
  for orientation in orientations:
    proj_mean = mean([data[orientation + "_Projector"] for data in comp_results])
    if proj_mean > highest_mean:
      highest_mean = proj_mean
      best_ori = orientation

  #! Most interesting combinations: the ones just higher than 90%, 75%, 50% and 25% ?

  print("[ DONE ]")

  # =================================================================== #
  # =================================================================== #
  #                  FILES MANIPULATION AND GENERATION                  #
  # =================================================================== #
  # =================================================================== #

  # Write the compiled results into a CSV file

  print ("{:<40}".format('\nCreating CSV file ...'), end="")

  csv_header = [key for key in comp_results[0].keys() if key != 'Directory']

  with open(os.path.join(out_dir,'convar_comp_results.csv'), 'w', newline='', encoding='utf-8') as csvfile:

    csv_writer = csv.DictWriter(csvfile, fieldnames=csv_header, delimiter=';', quoting=csv.QUOTE_MINIMAL)
    csv_writer.writeheader()

    for data in comp_results:
      del data['Directory']
      csv_writer.writerow(data)  

  print("[ DONE ]")

  print("")
  print("".center(columns,"*"))
  print("")
  print("END OF EXECUTION".center(columns))
  print("")
  print("".center(columns,"*"))

# =================================================================== #
# =================================================================== #
#                          CALL MAIN FUNCTION                         #
# =================================================================== #
# =================================================================== #

if __name__ == "__main__":
  main()   