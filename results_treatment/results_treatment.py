#!/usr/bin/env python3

################################################################################################################################################
##                                                             Results Treatment                                                              ##
##                                                                                                                                            ##
##                              This script scans one or more molecule directories containing all the information                             ##
##                    obtained through CHAINS and the various programs and generates the corresponding tables and graphs.                     ##
##                                                                                                                                            ##
##                    /!\ In order to run, this script requires Python 3.5+ as well as YAML, Jinja2 and GNUplot 5+. /!\                       ##
##                                          /!\ Ask your cluster(s) administrator(s) if needed. /!\                                           ##
################################################################################################################################################

import argparse
import contextlib
import csv
import os
import re
import shutil
import sys
from inspect import getsourcefile

import jinja2  # Only needed in the ind_results subscript, it is loaded here to check if your python installation does support jinja2
import numpy
import yaml

import results_errors

# =================================================================== #
# =================================================================== #
#                        FUNCTIONS DEFINITIONS                        #
# =================================================================== #
# =================================================================== #

def energy_unit_conversion(value:float,init:str,target:str) -> float:
    """|  Converts an energy value from an initial unit to a target unit by using atomic units of energy (Hartree) as an intermediary.
    |  Currently supported units: Hartree, cm\ :sup:`-1`\ , eV, nm, Hz and Joules

    Parameters
    ----------
    value : float
        The energy value we need to convert.
    init : str
        The unit of the value we need to convert.
    target : str
        The unit we must convert the value to.
    
    Returns
    -------
    conv_value : float
        The converted energy value.
    """

    # Define the dictionary of conversion factors, from atomic units (Hartree) to any unit you want. - Taken from the NIST website (https://physics.nist.gov/)

    conv_factors = {
      # 1 Hartree equals:
      "Ha" : 1,
      "cm-1" : 2.1947463136320e+05,
      "eV" : 27.211386245988,
      "nm" : 2.1947463136320e+05 / 1e+07,
      "Hz" : 6.579683920502e+15,
      "J" : 4.3597447222071e-18
      }

    # Put everything in lower cases, to make it case insensitive

    init_low = init.lower()
    target_low = target.lower()
    conv_factors_low = dict((key.lower(), value) for key, value in conv_factors.items())

    # Check if the desired units are supported

    if init_low not in conv_factors_low.keys():
      raise results_errors.ResultsError ("ERROR: The unit of the value you want to convert (%s) is currently not supported. Supported values include: %s" % (init, ', '.join(unit for unit in conv_factors.keys())))
    elif target_low not in conv_factors_low.keys():
      raise results_errors.ResultsError ("ERROR: The unit you want to convert the value to (%s) is currently not supported. Supported values include: %s" % (target, ', '.join(unit for unit in conv_factors.keys())))
    
    # Convert the value

    conv_value = (value / conv_factors_low[init_low]) * conv_factors_low[target_low]

    return conv_value

def import_path(fullpath:str):
    """ 
    Imports a file with full path specification. Allows one to import from anywhere, something __import__ does not do. Taken from https://stackoverflow.com/questions/72852/how-to-do-relative-imports-in-python

    Parameters
    ----------
    fullpath : str
        Full path towards the file you want to import

    Returns
    -------
    module
        The loaded file
    """

    # Split the path and filename (and remove extension of the filename)
    path, filename = os.path.split(fullpath)
    filename = os.path.splitext(filename)[0]

    # Add path to sys.path in order to be able to load the module, then remove it
    sys.path.insert(0, path)
    module = __import__(filename)
    del sys.path[0]

    return module

# =================================================================== #
# =================================================================== #
#                       COMMAND LINE ARGUMENTS                        #
# =================================================================== #
# =================================================================== #

# Define the arguments needed for the script (here they are defined as named arguments rather than positional arguments, check https://stackoverflow.com/questions/24180527/argparse-required-arguments-listed-under-optional-arguments for more info).

parser = argparse.ArgumentParser(add_help=False, description="For one or more molecule directories, this script reads the results files and generates the corresponding tables and graphs.")

required = parser.add_argument_group('Required arguments')
required.add_argument("-o","--out_dir", type=str, help="Path to the directory where you want to store the graphs and tables.", required=True)

mol_inp = parser.add_mutually_exclusive_group(required=True)
mol_inp.add_argument("-s","--single", type=str, help="Molecule directory containing the results files that need to be processed.")
mol_inp.add_argument("-m","--multiple", type=str, help="Directory containing multiple molecule directories.")

optional = parser.add_argument_group('Optional arguments')
optional.add_argument('-h','--help',action='help',default=argparse.SUPPRESS,help='Show this help message and exit')
optional.add_argument('-cf', '--config', type=str, help="Path to the YAML configuration file, default is this_script_directory/results_config.yml")

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
    print("EXECUTION OF THE RESULTS TREATMENT BEGINS NOW".center(columns))
    print("")
    print("".center(columns,"*"))

    # ========================================================= #
    # Read command line arguments                               #
    # ========================================================= #

    args = parser.parse_args()

    # Required arguments

    out_dir = args.out_dir                   # Directory where the graphs and tables will be created

    single_mol = args.single                 # Molecule directory containing the results files that need to be processed.
    multiple_mol = args.multiple             # Directory containing multiple molecule directories.

    # Optional arguments

    config_file = args.config                # YAML configuration file

    # ========================================================= #
    # Define codes directory                                    #
    # ========================================================= #

    # Determined by getting the path to the directory of this script

    code_dir = os.path.dirname(os.path.realpath(os.path.abspath(getsourcefile(lambda:0))))

    print ("{:<40} {:<100}".format('\nCodes directory:',code_dir))

    # ========================================================= #
    # Check and load the YAML configuration file                #
    # ========================================================= #

    if config_file: 
      config_file = results_errors.check_abspath(config_file,"Command line argument -cf / --config","file")
    else:
      # If no value has been provided through the command line, take the results_config.yml file in the same directory as this script 
      config_file = os.path.join(code_dir, "results_config.yml")

    print ("{:<40} {:<99}".format('\nLoading the configuration file',config_file + " ..."), end="")
    with open(config_file, 'r') as f_config:
      config = yaml.load(f_config, Loader=yaml.FullLoader)
    print('%12s' % "[ DONE ]")

    # ========================================================= #
    # Check molecule directories                                #
    # ========================================================= #

    if multiple_mol:

      multiple_mol = results_errors.check_abspath(multiple_mol,"Command line argument -m / --multiple","directory")
      mol_inp_path = multiple_mol

      print("{:<40} {:<99}".format("\nLooking for every molecule directory in", mol_inp_path + " ..."), end="")

      # We need to look for directories in the multiple_mol directory (see https://stackoverflow.com/questions/800197/how-to-get-all-of-the-immediate-subdirectories-in-python for reference).
      mol_inp_list = [dir.name for dir in os.scandir(mol_inp_path) if dir.is_dir()]

      if mol_inp_list == []:
        raise results_errors.ResultsError ("ERROR: Can't find any directory in %s" % mol_inp_path)
      
      print('%12s' % "[ DONE ]")

    else:

      single_mol = results_errors.check_abspath(single_mol,"Command line argument -s / --single","directory")
      print ("{:<40} {:<100}".format('\nMolecule directory:',single_mol))

      mol_inp_path = os.path.dirname(single_mol)
      mol_name = os.path.basename(single_mol)
      mol_inp_list = [mol_name]
 
    # ========================================================= #
    # Check other arguments                                     #
    # ========================================================= #

    out_dir = results_errors.check_abspath(out_dir,"Command line argument -o / --out_dir","directory")
    print ("{:<40} {:<100}".format('\nOutput directory:',out_dir))

    # ========================================================= #
    # Determine other important variables                       #
    # ========================================================= #

    quality_treshold = float(config["other"]["quality_treshold"])

    # Import the geom_scan.py file of ABIN LAUNCHER

    chains_path = os.path.dirname(code_dir)  
    geom_scan_path = os.path.join(chains_path,"abin_launcher","geom_scan.py")
    geom_scan = import_path(geom_scan_path)

  # ========================================================= #
  # Exception handling for the preparation step               #
  # ========================================================= #

  except results_errors.ResultsError as error:
    print("")
    print(error)
    exit(-1)

  # =================================================================== #
  # =================================================================== #
  #                   FILES MANIPULATION & GENERATION                   #
  # =================================================================== #
  # =================================================================== #

  for mol_name in mol_inp_list:

    mol_dir = os.path.join(mol_inp_path, mol_name)

    console_message = "Start procedure for the molecule " + mol_name
    print("")
    print(''.center(len(console_message)+11, '*'))
    print(console_message.center(len(console_message)+10))
    print(''.center(len(console_message)+11, '*'))

    # ========================================================= #
    # ========================================================= #
    #                 IDENTIFYING THE MOLECULE                  #
    # ========================================================= #
    # ========================================================= #

    # For more information on try/except structures, see https://www.tutorialsteacher.com/python/exception-handling-in-python
    try:

      print ("{:<140}".format('\nIdentifying molecule ...'), end="")

      # ========================================================= #
      # TAG of the molecule                                       #
      # ========================================================= #

      pattern = re.compile(r'^(?P<tag>[a-zA-Z]+\d+)-.*$')
      tag_finder = pattern.match(mol_name)

      tag = tag_finder.group('tag').upper()

      # ========================================================= #
      # Constitutive atoms                                        #
      # ========================================================= #

      # Check the optimized geometry file

      gaussian_dir = os.path.join(mol_dir, "GAUSSIAN")
      opt_geom_file = results_errors.check_abspath(os.path.join(gaussian_dir, mol_name + ".xyz"),"Optimized geometry file","file",True)

      # Load the XYZ geometry file

      with open(opt_geom_file, 'r') as mol_file:
        mol_content = mol_file.read().splitlines()

      # Call the scanning function from ABIN LAUNCHER but without its standard output (https://stackoverflow.com/questions/2828953/silence-the-stdout-of-a-function-in-python-without-trashing-sys-stdout-and-resto)

      with open(os.devnull, 'w') as devnull:
        with contextlib.redirect_stdout(devnull):
          file_data = geom_scan.xyz_scan(mol_content)

      chemical_formula = file_data["chemical_formula"]

      # Get the non-H atom types

      atoms = list(chemical_formula.keys())
      atoms.remove("H")
      atoms.sort()
      atoms.insert(0, atoms.pop(atoms.index("Si")))

      mol_type = ''.join(atoms)

      # ========================================================= #
      # "Pretty" name                                             #
      # ========================================================= #

      # If only 1 atom of that type, omit the number

      for atom,number in chemical_formula.items():
        if number == 1:
          chemical_formula[atom] = ""

      # Si on the front

      pretty_name_fr = "Si" + str(chemical_formula["Si"])

      # H on the end (if present)

      if "H" in chemical_formula:
        pretty_name_end = "H" + str(chemical_formula["H"])
      else:
        pretty_name_end = ""

      # Rest in the middle, in alphabetical order

      filtered_formula = {atom:number for atom, number in chemical_formula.items() if (atom != "Si" and atom != "H")}

      pretty_name_mid = ''.join('{}{}'.format(*atom) for atom in sorted(filtered_formula.items()))

      # Get the final name, correctly formatted

      pretty_name = pretty_name_fr + pretty_name_mid + pretty_name_end

      print('%12s' % "[ DONE ]")

    # ========================================================= #
    # Exception handling for the identification                 #
    # ========================================================= #

    except results_errors.ResultsError as error:
      print(error)
      print("Skipping %s molecule" % mol_name)
      continue

    # ========================================================= #
    # ========================================================= #
    #                 CHARACTERIZATION RESULTS                  #
    # ========================================================= #
    # ========================================================= #

    # For more information on try/except structures, see https://www.tutorialsteacher.com/python/exception-handling-in-python
    try:
    
      print ("{:<140}".format('\nComputing characterization values ...'), end="")

      # ========================================================= #
      # Check the data directory and its files                    #
      # ========================================================= #

      data_dir = results_errors.check_abspath(os.path.join(mol_dir,"CONTROL","data"),"Data directory created by control_launcher.py","directory",True)

      states_file = results_errors.check_abspath(os.path.join(data_dir, "states.csv"),"States file","file",True)
      mime_file = results_errors.check_abspath(os.path.join(data_dir, "mime"),"MIME file","file",True)
      momdip_file = results_errors.check_abspath(os.path.join(data_dir, "momdip_mtx"),"Transition dipole moments matrix file","file",True)
      transitions_file = results_errors.check_abspath(os.path.join(data_dir, "transitions.csv"),"Transitions file","file",True)

      eigenvalues_file = results_errors.check_abspath(os.path.join(data_dir, "eigenvalues"),"Eigenvalues file","file",True)
      eigenvectors_file = results_errors.check_abspath(os.path.join(data_dir, "eigenvectors"),"Eigenvectors file","file",True)
      transpose_file = results_errors.check_abspath(os.path.join(data_dir, "transpose"),"Transpose of the eigenvectors file","file",True)
      momdip_es_file = results_errors.check_abspath(os.path.join(data_dir, "momdip_es_mtx"),"Transition dipole moments matrix file, in the eigenstates basis set","file",True)

      # ========================================================= #
      # Load the files                                            #
      # ========================================================= #

      # CSV-type lists

      with open(states_file, 'r', newline='') as csv_file:

        states_content = csv.DictReader(csv_file, delimiter=';')
        states_list = list(states_content)
        states_header = states_content.fieldnames

      with open(transitions_file, 'r', newline='') as csv_file:

        transitions_content = csv.DictReader(csv_file, delimiter=';')
        transitions_list = list(transitions_content)
        transitions_header = transitions_content.fieldnames

      # Matrices (load with numpy then convert to list for easier handling)

      mime = numpy.loadtxt(mime_file).tolist()
      momdip = numpy.loadtxt(momdip_file).tolist()
      eigenvalues = numpy.loadtxt(eigenvalues_file).tolist()
      eigenvectors = numpy.loadtxt(eigenvectors_file).tolist()
      transpose = numpy.loadtxt(transpose_file).tolist()
      momdip_es = numpy.loadtxt(momdip_es_file).tolist()

      # Convert from strings to integer or floats when necessary

      eigenvalues = [float(value) for value in eigenvalues]
      momdip_es = [[float(moment) for moment in line] for line in momdip_es]

      for state in states_list:
        state["Energy (Ha)"] = float(state["Energy (Ha)"])
        state['Number'] = int(state['Number'])

      # ========================================================= #
      # Computing values                                          #
      # ========================================================= #

      # Relativistic shift
      # ==================

      shifts = [] # Initialize the list of relativistic shifts (one for each state)

      # Calculate each shift

      for eigenvalue in eigenvalues:
        
        shift = {}

        number = eigenvalues.index(eigenvalue)

        if number == 0:
          continue # Skip the ground state
        
        shift["Zero order state number"] = number
        shift["Zero order state energy"] = states_list[number]["Energy (Ha)"]

        shift["Eigenstate number"] = number
        shift["Eigenstate energy"] = eigenvalues[number]

        shift["Shift (Ha)"] = shift["Zero order state energy"] - shift["Eigenstate energy"]
        shift["Shift (%)"] = abs(shift["Shift (Ha)"]) / shift["Zero order state energy"]

        shifts.append(shift)

      # Calculate the mean and standard deviation

      shift_mean = numpy.mean([shift["Shift (%)"] for shift in shifts])
      shift_std = numpy.std([shift["Shift (%)"] for shift in shifts])

      # Brightest transition dipole moment
      # ==================================

      momdip_abs = [abs(moment) for moment in momdip_es[0]]
      brightest_number = momdip_abs.index(max(momdip_abs))
      brightest_mom = momdip_es[0][brightest_number]

      # Singlet-triplet transition dipole moment
      # ========================================

      # Identify singlet (bright) and triplet (dark) zero order states

      bright_list = []
      dark_list = []

      for state in states_list:
        if state['Number'] == 0:
          continue # exclude the ground state
        elif state['Type'].lower() == "dark":
          dark_list.append(state['Number'])
        elif state['Type'].lower() == "bright":
          bright_list.append(state['Number'])

      # Get each singlet-triplet transition dipole moment

      st_momdip_list = []

      for bright in bright_list:
        for dark in dark_list:
          st_momdip_list.append(momdip_es[bright][dark])

      # Calculate the mean and standard deviation

      st_momdip_mean = numpy.mean(st_momdip_list)
      st_momdip_std = numpy.std(st_momdip_list)

      # ========================================================= #
      # Writing values in the corresponding CSV file              #
      # ========================================================= #

      # Define line

      mol_line = {
        "Serie" : tag,
        "Type" : mol_type,
        "Molecule" : pretty_name,
        "Nb excited states" : len(eigenvalues) - 1,
        "Mean relativistic shift" : "{:.10e}".format(shift_mean),
        "Shift standard deviation" : "{:.10e}".format(shift_std),
        "Brightest dipole moment (a.u.)" : "{:.10e}".format(brightest_mom),
        "Mean singlet-triplet dipole moment (a.u.)" : "{:.10e}".format(st_momdip_mean),
        "S-T dipole moment std deviation (a.u.)" : "{:.10e}".format(st_momdip_std)
      }

      # Define the CSV file

      charac_csv_file = os.path.join(out_dir,"charac_results.csv")
      charac_csv_header = list(mol_line.keys())

      # Define if we have to write the header or not (only write it if the file does not exist or is empty)

      write_header = True

      if (os.path.exists(charac_csv_file) and os.path.isfile(charac_csv_file)):
        with open(charac_csv_file, 'r') as f:
          write_header = (not f.readline()) # If the file is empty, write_header = True. Otherwise, write_header = False

      # Open the CSV file in 'Append' mode and add the new line (+ write header if required)

      with open(charac_csv_file, 'a', newline='') as final_f:

        csv_writer = csv.DictWriter(final_f, fieldnames=charac_csv_header, delimiter=';', quoting=csv.QUOTE_MINIMAL)

        if write_header:
          csv_writer.writeheader()

        csv_writer.writerow(mol_line)      

      print('%12s' % "[ DONE ]")

    except results_errors.ResultsError as error:
      print(error)
      print("Skipping %s molecule" % mol_name)
      continue

    console_message = "End of procedure for the molecule " + mol_name
    print("")
    print(''.center(len(console_message)+10, '*'))
    print(console_message.center(len(console_message)+10))
    print(''.center(len(console_message)+10, '*'))

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
