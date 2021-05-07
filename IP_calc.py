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

def check_abspath(path:str,context:str,type="either",SkipError=False):
    """Checks if a path towards a file or directory exists and is of the correct type. If it is, the function returns its absolute version.

    Parameters
    ----------
    path : str
        The path towards the file or directory you want to test.
    context : str
        Message to show on screen to give more information in case of an exception (e.g. the role of the directory or file that was checked, where the checked path was given, etc.).
    type : str, optional
        The type of element for which you would like to test the path (file, directory or either).
        By default, checks if the path leads to either a file or a directory (type = either).
    SkipError : bool, optional
        By default, IP_Error exceptions will be caught and will cause the function to exit the script.
        Specify True to skip the error treatment and simply raise the exception.
    
    Returns
    -------
    abspath : str
        Normalized absolute version of the path.

    Raises
    ------
    ValueError
        If the specified type when calling the function is not "file", "directory" or "either".
    IP_Error
        If the type does not match what is given in the path, or if the path does not exist.
    """

    # For more information on try/except structures, see https://www.tutorialsteacher.com/python/exception-handling-in-python
    try:

      if type not in ["file","directory","either"]:
        raise ValueError ("The specified type for which the check_abspath function has been called is not one of 'file', 'directory' or 'either'")
      if not os.path.exists(path):
        raise IP_Error ("ERROR: %s does not seem to exist." % path)
      elif type == "file":
        if not os.path.isfile(path):
          raise IP_Error ("ERROR: %s is not a file" % path)
      elif type == "directory":
        if not os.path.isdir(path):
          raise IP_Error ("ERROR: %s is not a directory" % path)
      elif type == "either":
        if not os.path.isdir(path) and not os.path.isfile(path):
          raise IP_Error ("ERROR: %s is neither a file nor a directory" % path)

    except IP_Error as error:

        print("Something went wrong when checking the path ", path)
        print("Context: ",context)
        if not SkipError:
          print(error)
          exit(1)
        else:
          raise

    # If everything went well, get the normalized absolute version of the path
    
    abspath = os.path.abspath(path)

    return abspath

def gaussian_dft(source_content:list):
    """Parses the content of a Gaussian DFT calculation log file, looking to extract the HOMO energy and the SCF energy of the neutral molecule, as well as the SCF energy of the cation. The calculation must have performed a geometry optimization or a single point calculation of the neutral molecule first, followed by a a geometry optimization or a single point calculation of the cation. If a geometry optimization of the cation was performed, both the initial and final SCF energies of the cation will be extracted.
    Those values will be used by the main function to compute the ionization potentials of the molecule in three different ways: via Koopmans' theorem, in a vertical way (no changes in geometry) and in an adiabatic way (after having optimized the geometry of the cation)

    Parameters
    ----------
    source_content : list
        Content of the Gaussian log file. Each element of the list is a line of the file.
    
    Returns
    -------
    energies : dict
        The extracted information of the source file. It contains three keys and their associated values: ``homo``, ``neutral`` and ``cation`` where
        
        - ``homo`` is the energy of the HOMO (in Hartree)
        - ``neutral`` is the energy of the neutral molecule (in Hartree)
        - ``cation`` is the energy of the cation (in Hartree)

        If a geometry calculation of the cation was performed, an additionnal ``cation_opt`` key will be added, with a value corresponding to the energy of the cation in its optimized geometry (in Hartree)

    Raises
    ------
    IP_Error
        If some of the needed values are missing or unknown.

    """

    # ========================================================= #
    #                        HOMO energy                        #
    # ========================================================= #

    homo = False

    # Define the expression patterns for the lines containing information about the HOMO energy.

    homo_rx = {

      # Look for the "Alpha  occ. eigenvalues --   -0.35890  -0.35890  -0.35890  -0.28056  -0.28056" type of line (and store the last value of the line)
      'value': re.compile(r'^\s*Alpha  occ\. eigenvalues -- \s*(?:-?\d+\.\d+\s*)*\s*(?<value>-?\d+\.\d+)$'),

      # No need to go further than the "Normal termination of Gaussian" line (which signals the end of the neutral calculation)
      'end': re.compile(r'^\s*Normal termination of Gaussian')

    }

    # Parse the source file to get the information

    for line in source_content:
        
        # Compare each line to all the homo_rx expressions

        for key, rx in homo_rx.items():

          matching_line = rx.match(line)

          # If the line matches one of the rx patterns, act accordingly

          if matching_line is not None:

            if key == 'value':
              # Store the last numerical value of the line (this value will be overwritten everytime a new line of this type is encountered)
              homo = float(matching_line.group('value'))
    
            elif key == 'end':
              # Exit the loop once we reach the end of the neutral calculation, the last "homo" value stored in memory is the one corresponding to the HOMO.
              break         

    # Raise an exception if the HOMO energy has not been found

    if not homo:
      raise IP_Error ("ERROR: Unable to find the HOMO energy in the source file")

    # ========================================================= #
    #            SCF energy of the neutral molecule             #
    # ========================================================= #

    neutral = False

    # Define the expression patterns for the lines containing information about the SCF energy.

    neutral_rx = {

      # Look for the "SCF Done:  E(RB3LYP) =  -1454.81079575     A.U. after   11 cycles" type of line (and store the energy value)
      'value': re.compile(r'^\s*SCF Done:\s*E\(\w*\)\s*=\s*(?<value>-?\d+\.\d+)\s*A\.U\. after\s*\d+\s*cycles$'),

      # No need to go further than the "Normal termination of Gaussian" line (which signals the end of the neutral calculation)
      'end': re.compile(r'^\s*Normal termination of Gaussian')

    }

    # Parse the source file to get the information

    for line in source_content:
        
        # Compare each line to all the neutral_rx expressions

        for key, rx in neutral_rx.items():

          matching_line = rx.match(line)

          # If the line matches one of the rx patterns, act accordingly

          if matching_line is not None:

            if key == 'value':
              # Store the energy value of the line (this value will be overwritten everytime a new line of this type is encountered)
              neutral = float(matching_line.group('value'))
    
            elif key == 'end':
              # Exit the loop once we reach the end of the neutral calculation, the last "neutral" value stored in memory is the one corresponding to the SCF energy of the optimized geometry (for a single point calculation, only one value will be encountered).
              break         

    # Raise an exception if the SCF energy of the neutral molecule has not been found

    if not neutral:
      raise IP_Error ("ERROR: Unable to find the SCF energy of the neutral molecule in the source file")

    # ========================================================= #
    #               SCF energy(ies) of the cation               #
    # ========================================================= #

    cation_opt = False
    section_found = False
    first = True

    # Define the expression patterns for the lines containing information about the SCF energy.

    cation_rx = {

      # Look for the "Normal termination of Gaussian" line (which signals the end of the neutral calculation but also the end of the cation calculation)
      'limit': re.compile(r'^\s*Normal termination of Gaussian'),

      # Look for the "SCF Done:  E(RB3LYP) =  -1454.81079575     A.U. after   11 cycles" type of line (and store the energy value)
      'value': re.compile(r'^\s*SCF Done:\s*E\(\w*\)\s*=\s*(?<value>-?\d+\.\d+)\s*A\.U\. after\s*\d+\s*cycles$')

    }

    # Parse the source file to get the information

    for line in source_content:

      # Define when the section begins and ends

      if not section_found:
        if cation_rx['limit'].match(line): # The line matches the limit pattern
          section_found = True
          continue

      elif section_found and cation_rx['limit'].match(line):
        # If the section was already found and a second occurrence of the limit pattern is found, it is time to leave
        break

      # Process the section lines
  
      else:
  
        matching_line = cation_rx['value'].match(line)

        # If the line matches our pattern, store the energy value of the line 

        if matching_line is not None:

          if first:
            cation = float(matching_line.group('value')) # Only store the first occurrence of SCF energy in the cation variable
            first = False
        
          cation_opt = float(matching_line.group('value')) # This value will be overwritten everytime a new line of this type is encountered

    # Raise an exception if the section has not been found

    if not cation_opt:
      raise IP_Error ("ERROR: Unable to find the SCF energy of the cation in the source file")

    # ========================================================= #
    #                     Return the values                     #
    # ========================================================= #   

    energies = {
        "homo": homo,
        "neutral": neutral,
        "cation": cation
    }

    if cation != cation_opt:

        energies.update({
            "cation_opt": cation_opt
        })

    print("")
    print(''.center(50, '-'))
    print("{:<30} {:<20}".format("HOMO energy: ", "{:.5e}".format(homo)))
    print("{:<30} {:<20}".format("Neutral molecule energy: ", "{:.5e}".format(neutral)))
    print("{:<30} {:<20}".format("Cation energy: ", "{:.5e}".format(cation)))
    if cation != cation_opt:
        print("{:<30} {:<20}".format("Cation energy (optimized): ", "{:.5e}".format(cation_opt)))
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
required.add_argument("-o", "--output", type=str, help="Path towards the CSV file that will contain the resulting values.", required=True)
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

    print ("{:<40} {:<100}".format('\nOutput CSV file:',output))

    # ========================================================= #
    # Check arguments                                           #
    # ========================================================= #

    # Check the existence of the source file, then get its name and the name of the directory where it is located

    source = check_abspath(source,"Command line argument -s / --source","file")
    print ("{:<40} {:<100}".format('\nSource file:',source))

    source_filename = os.path.basename(source)
    source_name = str(source_filename.split('.')[0]) # Getting rid of the format extension to get the name of the source

    # Check the existence of the parsing function that will fetch the relevant information from the source file

    if (parsing_fct) not in dir() or not callable(parsing_fct):
      raise IP_Error ("ERROR: There is no parsing function named %s defined in this script." % parsing_fct)

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

    print ("{:<80}".format('\nLoading the source file ...'), end="")
    with open(source, 'r') as source_file:
      source_content = source_file.read().splitlines()
    print("[ DONE ]")

    # Cleaning up the source file from surrounding spaces and blank lines

    source_content = list(map(str.strip, source_content))   # Remove leading & trailing blank/spaces
    source_content = list(filter(None, source_content))     # Remove blank lines/no char

    # ========================================================= #
    # Call the parsing function                                 #
    # ========================================================= #

    print ("{:<50} {:<100}".format('\nParsing function:',parsing_fct))

    print ("\nParsing the source file ...")

    # Call the parsing function (fixed for now, but might change in the future)

    energies = eval(parsing_fct)(source_content)

    # Check the energies dictionary

    required_keys = frozenset({"homo", "neutral", "cation"})

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

  # Compute the IP value according to Koopmans' theorem (opposite of the HOMO energy)

  ip_koopmans = - float(energies['homo'])

  # Compute the vertical IP value

  ip_vertical =  float(energies['cation']) - float(energies['neutral'])

  # If cation_opt was available, compute the adiabatic IP value

  if energies.get('cation_opt'):
      ip_adiabatic = float(energies['cation_opt']) - float(energies['neutral'])
  else:
      ip_adiabatic = "N/A"

  # ========================================================= #
  # Add information to output CSV file                        #
  # ========================================================= #

  # Define line

  ip_line = {
      "Molecule" : source_name,
      "HOMO energy" : energies['homo'],
      "Neutral molecule energy" : energies['neutral'],
      "Cation energy" : energies['cation'],
      "Cation energy (optimized)" : energies.get('cation_opt'),
      "IP (Koopmans)" : ip_koopmans,
      "IP (vertical)" : ip_vertical,
      "IP (adiabatic)" : ip_adiabatic
  }

  # Define if we have to write the header or not (only write it if the file does not exist or is empty)

  csv_header = list(ip_line.keys())
  write_header = True

  if (os.path.exists(output) and os.path.isfile(output)):
    with open(output, 'r', encoding='utf-8') as f:
      write_header = (not f.readline()) # If the file is empty, write_header = True. Otherwise, write_header = False

  # Open the CSV file in 'Append' mode and add the new line (+ write header if required)

  with open(output, 'a', newline='', encoding='utf-8') as final_f:

    csv_writer = csv.DictWriter(final_f, fieldnames=csv_header, delimiter=';', quoting=csv.QUOTE_MINIMAL)

    if write_header:
      csv_writer.writeheader()

    csv_writer.writerow(ip_line)    

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
