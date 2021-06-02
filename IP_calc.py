#!/usr/bin/env python3

################################################################################################################################################
##                                               Ionization potentials calculator for molecules                                               ##
##                                                                                                                                            ##
##                                 This script computes the ionization potentials of molecules by extracting                                  ##
##                         information from a given source file then writes the resulting values in a given CSV file.                         ##
################################################################################################################################################

import argparse
import csv
import os
import re
import shutil
from tempfile import NamedTemporaryFile

# =================================================================== #
# =================================================================== #
#                        EXCEPTIONS DEFINITIONS                       #
# =================================================================== #
# =================================================================== #

class Error(Exception):
    """Base class for exceptions in this script."""
    pass

class IP_Error(Error):
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
    IP_Error
        If the type does not match what is given in the path, or if the path does not exist.
    """

    # Check "type" argument

    if type not in ["file","directory","either"]:
      raise ValueError ("The specified type for which the check_abspath function has been called is not one of 'file', 'directory' or 'either'")

    # Prepare to print a helpful error message in case of problem with the given path

    msg = "\nSomething went wrong when checking the path " + path + "\nContext: " + context + "\n"

    # Check path
    
    if not os.path.exists(path):
      raise IP_Error (msg + "ERROR: %s does not seem to exist." % path)
    elif type == "file":
      if not os.path.isfile(path):
        raise IP_Error (msg + "ERROR: %s is not a file" % path)
    elif type == "directory":
      if not os.path.isdir(path):
        raise IP_Error (msg + "ERROR: %s is not a directory" % path)
    elif type == "either":
      if not os.path.isdir(path) and not os.path.isfile(path):
        raise IP_Error (msg + "ERROR: %s is neither a file nor a directory" % path)

    # If everything went well, get the normalized absolute version of the path
    
    abspath = os.path.abspath(path)

    return abspath

def gaussian_dft(source_content:list):
    """Parses the content of a Gaussian DFT calculation log file, looking to extract the HOMO energy and the SCF energy of the molecule. If a cation single point calculation was performed, the function also extracts the SCF energy of the cation. If a cation geometry optimization was performed, the function extracts both the initial and final SCF energies of the cation. The calculation of the cation must directly follow the calculation of the neutral molecule in the log file. The charge of the cation must be equal to the neutral charge plus one.

    Parameters
    ----------
    source_content : list
        Content of the Gaussian log file. Each element of the list is a line of the file.
    
    Returns
    -------
    energies : dict
        The extracted information of the source file. It contains two mandatory keys and two optional keys with their associated values: ``homo``, ``neutral``, ``cation`` and ``cation_opt`` where
        
        - ``homo`` is the energy of the HOMO (in Hartree)
        - ``neutral`` is the energy of the neutral molecule (in Hartree)
        - ``cation`` is the energy of the cation in the neutral geometry (in Hartree) - only if a cation calculation was performed (single point or geometry optimization)
        - ``cation_opt``is the energy of the cation in its optimized geometry (in Hartree) - only if a cation geometry optimization was performed

    Raises
    ------
    IP_Error
        If some of the needed values are missing or unknown.

    """

    # ========================================================= #
    #                     Neutral molecule                      #
    # ========================================================= #

    homo = False
    neutral = False

    # Define the expression patterns for the lines containing information about the neutral molecule.

    neutral_rx = {

      # Look for the "Alpha  occ. eigenvalues --   -0.35890  -0.35890  -0.35890  -0.28056  -0.28056" type of line (and store the last value of the line)
      'homo_value': re.compile(r'^\s*Alpha  occ\. eigenvalues -- \s*(?:-?\d+\.\d+\s*)*\s*(?P<value>-?\d+\.\d+)$'),

      # Look for the "SCF Done:  E(RB3LYP) =  -1454.81079575     A.U. after   11 cycles" type of line (and store the energy value)
      'scf_value': re.compile(r'^\s*SCF Done:\s*E\(\w*\)\s*=\s*(?P<value>-?\d+\.\d+)\s*A\.U\. after\s*\d+\s*cycles$'),

      # Look for the "Charge =  0 Multiplicity = 1" type of line (and store the charge value)
      'charge': re.compile(r'^\s*Charge\s*=\s*(?P<charge>\d+)\s*Multiplicity\s*=\s*\d+\s*$'),

      # No need to go further than the "Normal termination of Gaussian" line (which signals the end of the neutral calculation)
      'end': re.compile(r'^\s*Normal termination of Gaussian.*$')

    }

    # Parse the source file to get the information

    for line in source_content:

      # Define when the section ends

      if neutral_rx['end'].match(line):
        break  

      # Store the charge of the neutral molecule to compare it with the cation charge later

      elif neutral_rx['charge'].match(line):
        neutral_charge = int(neutral_rx['charge'].match(line).group('charge'))

      # If the line matches our orbitals energy pattern, store the last energy value of the line (this value will be overwritten everytime a new line of this type is encountered)

      elif neutral_rx['homo_value'].match(line):
        homo = float(neutral_rx['homo_value'].match(line).group('value'))

      # If the line matches our SCF energy pattern, store the energy value of the line (this value will be overwritten everytime a new line of this type is encountered)

      elif neutral_rx['scf_value'].match(line):
        neutral = float(neutral_rx['scf_value'].match(line).group('value'))

    # Raise an exception if the HOMO energy has not been found

    if not homo:
      raise IP_Error ("ERROR: Unable to find the HOMO energy in the source file")

    # Raise an exception if the SCF energy of the neutral molecule has not been found

    if not neutral:
      raise IP_Error ("ERROR: Unable to find the SCF energy of the neutral molecule in the source file")

    # ========================================================= #
    #               SCF energy(ies) of the cation               #
    # ========================================================= #

    cation = False
    cation_opt = False
    section_found = False
    first = True

    # Define the expression patterns for the lines containing information about the SCF energy.

    cation_rx = {

      # Look for the "Normal termination of Gaussian" line (which signals the end of the neutral calculation but also the end of the cation calculation)
      'limit': re.compile(r'^\s*Normal termination of Gaussian'),

      # Look for the "Charge =  1 Multiplicity = 2" type of line (and store the charge value)
      'charge': re.compile(r'^\s*Charge\s*=\s*(?P<charge>\d+)\s*Multiplicity\s*=\s*\d+\s*$'),

      # Look for the "SCF Done:  E(RB3LYP) =  -1454.81079575     A.U. after   11 cycles" type of line (and store the energy value)
      'value': re.compile(r'^\s*SCF Done:\s*E\(\w*\)\s*=\s*(?P<value>-?\d+\.\d+)\s*A\.U\. after\s*\d+\s*cycles$')

    }

    # Parse the source file to get the information

    for line in source_content:

      # Define when the section begins and ends

      if not section_found:
        if cation_rx['limit'].match(line): # The line matches the limit pattern
          section_found = True

      elif section_found and cation_rx['limit'].match(line):
        # If the section was already found and a second occurrence of the limit pattern is found, it is time to leave
        break

      # Check the charge of the cation

      elif cation_rx['charge'].match(line):

        cation_charge = int(cation_rx['charge'].match(line).group('charge'))
        if cation_charge and (cation_charge - neutral_charge != 1):
          raise IP_Error("ERROR: The cation charge (%s) does not correspond to the charge of the neutral molecule (%s) plus one. Ionization potentials cannot be defined." % (cation_charge,neutral_charge))

      # If the line matches our SCF energy pattern, store the energy value of the line 

      elif cation_rx['value'].match(line):

        if first:
          cation = float(cation_rx['value'].match(line).group('value'))     # Only store the first occurrence of SCF energy in the cation variable
          first = False
        else:        
          cation_opt = float(cation_rx['value'].match(line).group('value')) # This value will be overwritten everytime a new line of this type is encountered (except the first value)

    # ========================================================= #
    #                     Return the values                     #
    # ========================================================= #   

    energies = {
      "homo": homo,
      "neutral": neutral
    }

    if cation:
      energies.update({
        "cation": cation
      })

    if cation_opt:
      energies.update({
        "cation_opt": cation_opt
      })

    print("")
    print(''.center(50, '-'))
    print("{:>30} {:<20}".format("HOMO energy: ", "{:.5f}".format(homo)))
    print("{:>30} {:<20}".format("Neutral molecule energy: ", "{:.5f}".format(neutral)))
    print("{:>30} {:<20}".format("Cation energy: ", "{:.5f}".format(cation) if cation else "N/A"))
    print("{:>30} {:<20}".format("Cation energy (optimized): ", "{:.5f}".format(cation_opt) if cation_opt else "N/A"))
    print(''.center(50, '-'))

    return energies

# =================================================================== #
# =================================================================== #
#                       COMMAND LINE ARGUMENTS                        #
# =================================================================== #
# =================================================================== #

# Define the arguments needed for the script (here they are defined as named arguments rather than positional arguments, check https://stackoverflow.com/questions/24180527/argparse-required-arguments-listed-under-optional-arguments for more info).

parser = argparse.ArgumentParser(add_help=False, description="This script computes the ionization potentials of molecules by extracting information from a given source file then writes the resulting values in a given CSV file.")

required = parser.add_argument_group('Required arguments')
required.add_argument("-s","--source", type=str, help="Path to the source file containing all the necessary information that needs to be processed.", required=True)
required.add_argument("-o", "--output", type=str, help="Path towards the CSV file that will contain the resulting values, extension must be .csv.", required=True)
#required.add_argument("-p", "--parsing", type=str, help="Name of the parsing function that will extract the necessary information from the source file.", required=True)

optional = parser.add_argument_group('Optional arguments')
optional.add_argument('-h','--help',action='help',default=argparse.SUPPRESS,help='Show this help message and exit')

# =================================================================== #
# =================================================================== #
#                            MAIN FUNCTION                            #
# =================================================================== #
# =================================================================== #

# We encapsulate all the instructions in a main function, this is done so that the documentation (or other scripts) can import this file without immediately executing it (see https://realpython.com/python-main-function/ for details)
def main(): 

  # =================================================================== #
  # =================================================================== #
  #                           PREPARATION STEP                          #
  # =================================================================== #
  # =================================================================== #

  # For more information on try/except structures, see https://www.tutorialsteacher.com/python/exception-handling-in-python
  try:

    # Get the size of the terminal in order to have a prettier output, if you need something more robust, go check http://granitosaurus.rocks/getting-terminal-size.html

    columns, rows = shutil.get_terminal_size()

    # Output Header

    print("".center(columns,"*"))
    print("")
    print("EXECUTION OF THE IONIZATION POTENTIALS CALCULATOR BEGINS NOW".center(columns))
    print("")
    print("".center(columns,"*"))

    # ========================================================= #
    # Read command line arguments                               #
    # ========================================================= #

    args = parser.parse_args()

    # Required arguments

    source = args.source                     # Source file containing all the necessary information
    output = args.output                     # CSV file that will contain the resulting values
    parsing_fct = "gaussian_dft"             # Parsing function that will extract the necessary information from the source file (if you decide to add it as a command line argument, replace "gaussian_dft" by args.parsing)
  
    # ========================================================= #
    # Check arguments                                           #
    # ========================================================= #

    # Check the existence of the source file, then get its name and the name of the directory where it is located

    source = check_abspath(source,"Command line argument -s / --source","file")
    print ("{:<40} {:<100}".format('\nSource file:',source))

    source_filename = os.path.basename(source)
    source_name = str(source_filename.split('.')[0]) # Getting rid of the format extension to get the name of the source

    # If the output file already exists, make sure it's a file with the .csv extension

    if os.path.exists(output):

      if not os.path.isfile(output):
        raise IP_Error ("ERROR: The command line argument -o / --output (%s) points towards a nonfile object." % output)

      if os.path.isfile(output) and os.path.splitext(output)[-1].lower() != (".csv"):
        raise IP_Error ("ERROR: %s is probably not a CSV file (the extension must be '.csv')." % output)

    print ("{:<40} {:<100}".format('\nOutput CSV file:',output))

    #!TODO Check the existence of the parsing function that will fetch the relevant information from the source file

#    if not callable(parsing_fct):
#      raise IP_Error ("ERROR: There is no parsing function named %s defined in this script." % parsing_fct)

    print ("{:<40} {:<100}".format('\nParsing function:',parsing_fct))

  # ========================================================= #
  # Exception handling for the preparation step               #
  # ========================================================= #

  except IP_Error as error:
    print("")
    print(error)
    exit(-1)

  # =================================================================== #
  # =================================================================== #
  #                        SOURCE FILE TREATMENT                        #
  # =================================================================== #
  # =================================================================== #

  # For more information on try/except structures, see https://www.tutorialsteacher.com/python/exception-handling-in-python
  try:

    # ========================================================= #
    # Load the source file                                      #
    # ========================================================= #

    print ("{:<40}".format('\nLoading the source file ...'), end="")
    with open(source, 'r') as source_file:
      source_content = source_file.read().splitlines()
    print("[ DONE ]")

    # Cleaning up the source file from surrounding spaces and blank lines

    source_content = list(map(str.strip, source_content))   # Remove leading & trailing blank/spaces
    source_content = list(filter(None, source_content))     # Remove blank lines/no char

    # ========================================================= #
    # Call the parsing function                                 #
    # ========================================================= #

    print ("\nParsing the source file ...")

    # Call the parsing function (fixed for now, but might change in the future)

    energies = eval(parsing_fct)(source_content)

    # Check the energies dictionary

    required_keys = frozenset({"homo", "neutral"})

    if not isinstance(energies, dict):
      raise IP_Error ('ERROR: The "energies" variable returned by the %s parsing function is not a dictionary.' % parsing_fct) 

    for key in required_keys:
      if key not in energies:
        raise IP_Error ('ERROR: There is no defined "%s" key in the energies dictionary returned by the %s parsing function.' % (key, parsing_fct))  

    print("\nThe source file has been succesfully parsed.")

  # ========================================================= #
  # Exception handling for the source file treatment          #
  # ========================================================= #

  except IP_Error as error:
    print("")
    print(error)
    exit(-1)

  # =================================================================== #
  # =================================================================== #
  #                  IONIZATION POTENTIALS CALCULATION                  #
  # =================================================================== #
  # =================================================================== #

  # ========================================================= #
  # Compute the ionization potentials (IP)                    #
  # ========================================================= #

  print ("{:<40}".format('\nComputing ionization potentials ...'), end="")

  # Compute the IP value according to Koopmans' theorem (opposite of the HOMO energy)

  ip_koopmans = - float(energies['homo'])

  # If cation was available, compute the vertical IP value

  if energies.get('cation'):
    ip_vertical = float(energies['cation']) - float(energies['neutral'])
  else:
    ip_vertical = "N/A"

  # If cation_opt was available, compute the adiabatic IP value

  if energies.get('cation_opt'):
    ip_adiabatic = float(energies['cation_opt']) - float(energies['neutral'])
  else:
    ip_adiabatic = "N/A"

  print("[ DONE ]")

  print("")
  print(''.center(50, '-'))
  print("{:>30} {:<20}".format("Koopmans' approx. (-HOMO): ", "{:.5f}".format(ip_koopmans)))
  print("{:>30} {:<20}".format("Vertical: ", "{:.5f}".format(ip_vertical) if ip_vertical != "N/A" else ip_vertical))
  print("{:>30} {:<20}".format("Adiabatic: ", "{:.5f}".format(ip_adiabatic) if ip_adiabatic != "N/A" else ip_adiabatic))
  print(''.center(50, '-'))

  # ========================================================= #
  # Add information to output CSV file                        #
  # ========================================================= #

  # Define line

  ip_line = {
    "Molecule" : source_name,
    "HOMO energy" : energies['homo'],
    "Neutral molecule energy" : energies['neutral'],
    "Cation energy" : energies.get('cation',"N/A"),
    "Cation energy (optimized)" : energies.get('cation_opt',"N/A"),
    "IP (Koopmans)" : ip_koopmans,
    "IP (vertical)" : ip_vertical,
    "IP (adiabatic)" : ip_adiabatic
  }

  # Define the CSV header

  csv_header = list(ip_line.keys())

  # If the output CSV file does not exist, create it and add the new line

  if not os.path.exists(output):

    print ("{:<40}".format('\nCreating the output CSV file ...'), end="")

    with open(output, 'w', newline='', encoding='utf-8') as csvfile:

      csv_writer = csv.DictWriter(csvfile, fieldnames=csv_header, delimiter=';', quoting=csv.QUOTE_MINIMAL)
      csv_writer.writeheader()
      csv_writer.writerow(ip_line) 

    print("[ DONE ]")

  # If the output CSV file already exists, just add the new line or update the file, using a temporary file as intermediary (see https://stackoverflow.com/questions/46126082/how-to-update-rows-in-a-csv-file for reference)

  else:

    tempfile = NamedTemporaryFile(mode='w', delete=False, encoding='utf-8')

    with open(output, 'r', encoding='utf-8') as csvfile, tempfile:

      csv_reader = csv.DictReader(csvfile, fieldnames=csv_header, delimiter=';')
      csv_writer = csv.DictWriter(tempfile, fieldnames=csv_header, delimiter=';', quoting=csv.QUOTE_MINIMAL)

      # Rewrite each line of our initial CSV file into the temp file, and update the values if our molecule was already present

      found = False

      for line in csv_reader:
      
        if line['Molecule'] == source_name:
          print("{:<40}".format("\nUpdating values for %s ..." % source_name), end ='')
          found = True
          line = ip_line
      
        csv_writer.writerow(line)
      
      # If our molecule was not already mentioned, add the new line

      if not found:
        print("{:<40}".format("\nAdding new line to output file ..."), end ='')
        csv_writer.writerow(ip_line)

    # Replace the initial CSV file by the temp file

    shutil.move(tempfile.name, output)

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

# If this script is executed through the command line, call the main function (see https://realpython.com/python-main-function/ for details)

if __name__ == "__main__":
  main()     
