#!/usr/bin/env python3

################################################################################################################################################
##                                                              Results Compiler                                                              ##
##                                                                                                                                            ##
##                              This script scans one or more molecule directories containing all the information                             ##
##                    obtained through CHAINS and the various programs and compiles the results into a single YAML file.                      ##
##                                                                                                                                            ##
##                                 /!\ In order to run, this script requires Python 3.5+ as well as YAML /!\                                  ##
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
      "nm" : 2.1947463136320e+05 * 1e+07,
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

parser = argparse.ArgumentParser(add_help=False, description="For one or more molecule directories, this script reads the results files and compiles the results into a single YAML file.")

required = parser.add_argument_group('Required arguments')
required.add_argument("-o","--out_yml", type=str, help="Path to the YAML file in which you want to compile the results. If the file already exists, it will be updated with the new content.", required=True)

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

    out_yml = args.out_yml                   # YAML file in which the results will be compiled

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
    # Load important files                                      #
    # ========================================================= #

    # YAML configuration file
    # =======================

    if config_file: 
      config_file = results_errors.check_abspath(config_file,"Command line argument -cf / --config","file")
    else:
      # If no value has been provided through the command line, take the results_config.yml file in the same directory as this script 
      config_file = os.path.join(code_dir, "results_config.yml")

    print ("{:<40} {:<99}".format('\nLoading the configuration file',config_file + " ..."), end="")
    with open(config_file, 'r') as f_config:
      config = yaml.load(f_config, Loader=yaml.FullLoader)
    print('%12s' % "[ DONE ]")

    # CHAINS YAML configuration file
    # ==============================

    chains_path = os.path.dirname(code_dir) 
    chains_config_file = results_errors.check_abspath(os.path.join(chains_path,"configs","chains_config.yml"),"CHAINS configuration YAML file","file")

    print ("{:<140}".format("\nLoading CHAINS configuration YAML file ..."), end="")
    with open(chains_config_file, 'r') as chains:
      chains_config = yaml.load(chains, Loader=yaml.FullLoader)
    print('%12s' % "[ DONE ]")

    # Ionization potentials CSV file
    # ==============================

    ip_file = results_errors.check_abspath(chains_config['ip_file'],"Ionization potentials CSV file","file")

    print ("{:<140}".format("\nLoading ionization potentials CSV file ..."), end="")
    with open(ip_file, 'r', newline='') as csv_file:
      ip_content = csv.DictReader(csv_file, delimiter=';')
      ip_list = list(ip_content)
    print('%12s' % "[ DONE ]")

    # ABIN LAUNCHER's geom_scan.py file
    # =================================

    geom_scan_path = os.path.join(chains_path,"abin_launcher","geom_scan.py")

    print ("{:<140}".format("\nImporting ABIN LAUNCHER's geom_scan.py file ..."), end="")
    geom_scan = import_path(geom_scan_path)
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
    # Determine other important variables                       #
    # ========================================================= #

    out_yml = os.path.abspath(out_yml)
    print ("{:<40} {:<100}".format('\nOutput YAML file:',out_yml))

    comp_results = {} # Dictionary consisting of multiples dictionaries containing data about each molecule (1 dict per molecule)

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
      # Group (constitutive atoms)                                #
      # ========================================================= #

      # Check the optimized geometry file

      gaussian_dir = os.path.join(mol_dir, "GAUSSIAN")
      opt_geom_file = results_errors.check_abspath(os.path.join(gaussian_dir, mol_name + ".xyz"),"Optimized geometry file","file")

      # Load the XYZ geometry file

      with open(opt_geom_file, 'r') as mol_file:
        mol_content = mol_file.read().splitlines()

      # Call the scanning function from ABIN LAUNCHER but without its standard output (https://stackoverflow.com/questions/2828953/silence-the-stdout-of-a-function-in-python-without-trashing-sys-stdout-and-resto)

      with open(os.devnull, 'w') as devnull:
        with contextlib.redirect_stdout(devnull):
          file_data = geom_scan.xyz_scan(mol_content)

      chemical_formula = file_data["chemical_formula"]

      # Get the non-H atom types, and put Si at the front

      atoms = list(chemical_formula.keys())
      atoms.remove("H")
      atoms.sort()
      atoms.insert(0, atoms.pop(atoms.index("Si")))

      mol_group = ''.join(atoms)

      # ========================================================= #
      # Total number of atoms                                     #
      # ========================================================= #

      nb_atoms = sum(file_data['chemical_formula'].values())

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
      # Store ID info                                             #
      # ========================================================= #

      comp_results[mol_name] = { "ID" : 
          { "TAG" : tag,
            "Group" : mol_group,
            "Name" : pretty_name,
            "Nb atoms" : nb_atoms
          }
        }

      # Add the number of Si atoms (for pure Si QDs only)

      if mol_group == "Si":
        comp_results[mol_name]['ID'].update({"Nb Si atoms" : chemical_formula["Si"] if chemical_formula["Si"] != "" else 1})

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
    
      # Get the path towards the data directory containing the relevant files

      data_dir = results_errors.check_abspath(os.path.join(mol_dir,"CONTROL","data"),"Data directory created by control_launcher.py","directory")

      # ========================================================= #
      # Kohn-Sham orbitals                                        #
      # ========================================================= #

      print ("{:<140}".format('\nFetching KS orbitals values ...'), end="")

      # Load the QCHEM output file

      qchem_file = results_errors.check_abspath(os.path.join(data_dir, mol_name + ".out"),"QCHEM output file","file")

      with open(qchem_file, 'r') as out_file:
        qchem_content = out_file.read().splitlines()

      qchem_content = list(map(str.strip, qchem_content))   # Remove leading & trailing blank/spaces
      qchem_content = list(filter(None, qchem_content))     # Remove blank lines/no char

      # Initialize some variables

      section_found = False
      occ_orb = []
      virt_orb = []

      # Define the expression patterns for the lines containing information about the orbitals
    
      orb_rx = {

        # Pattern for finding the "Orbital Energies (a.u.)" line (which marks the start of the section)
        'start': re.compile(r'^\s*Orbital Energies \(a\.u\.\)'),

        # Pattern for finding lines looking like '-- Occupied --'
        'occupied': re.compile(r'^\s*-- Occupied --\s*$'),

        # Pattern for finding lines looking like '-- Virtual --'
        'virtual': re.compile(r'^\s*-- Virtual --\s*$'),

        # Pattern for finding lines looking like ' -5.286  -5.276  -3.648  -3.648  -3.648  -3.648  -3.647  -3.647' or '********-101.5017-101.5016 -66.1592 -66.1591 -66.1591 -66.1251 -66.114'
        # The number of values is unknown but we do not need to capture the values
        'energies': re.compile(r'^\s*\**\s*(?:-?\d+\.\d+\s*)+$'),

        # Pattern for finding the "Beta MOs" line (which marks the end of the section, no need to get both sets of orbitals)
        'end': re.compile(r'^\s*Beta MOs')

      }

      # Parse the qchem output file to get the information

      for line in qchem_content:

        # Define when the section begins and ends (ignore beta orbitals)

        if not section_found:
          if orb_rx['start'].match(line):
            section_found = True
      
        elif section_found and orb_rx['end'].match(line):
          break

        # Store the type of the orbital (occupied / virtual)

        elif orb_rx['occupied'].match(line):
          orb_type = "Occupied"

        elif orb_rx['virtual'].match(line):
          orb_type = "Virtual"

        # If the line matches our energies pattern, iterate over the values and store them, according to their type (https://stackoverflow.com/questions/39668228/python-regex-catch-variable-number-of-groups)

        elif orb_rx['energies'].match(line):
          for match in re.finditer(r'(-?\d+\.\d+)', line):
            value = match.group(1)
            if orb_type == "Occupied":
              occ_orb.append(float(value))
            elif orb_type == "Virtual":
              virt_orb.append(float(value))

      # Raise an exception if the section has not been found

      if not section_found:
        raise results_errors.ResultsError ("ERROR: Unable to find the 'Orbital Energies (a.u.)' section in the QCHEM output file")

      # Raise an exception if the energies of the orbitals have not been found

      if occ_orb == []:
        raise results_errors.ResultsError ("ERROR: Unable to find the energies of the occupied orbitals in the QCHEM output file")

      if virt_orb == []:
        raise results_errors.ResultsError ("ERROR: Unable to find the energies of the virtual orbitals in the QCHEM output file")

      # Store the energies of the orbitals as a big string where values are delimited by ";"

      comp_results[mol_name].update({"QCHEM KS Orbitals" :
          { "Occupied" : ";".join(map(str,occ_orb)),
            "Virtual" : ";".join(map(str,virt_orb))
          }
        })

      print('%12s' % "[ DONE ]")

      # ========================================================= #
      # Compute energy gaps                                       #
      # ========================================================= #

      print ("{:<140}".format('\nComputing energy gaps ...'), end="")

      # HOMO-LUMO gap
      # =============

      homo_value = max(occ_orb)
      lumo_value = min(virt_orb)
      hl_gap = lumo_value - homo_value

      # Optical gap
      # ===========

      # Load the states list

      states_file = results_errors.check_abspath(os.path.join(data_dir, "states.csv"),"States file","file")

      with open(states_file, 'r', newline='') as csv_file:

        states_content = csv.DictReader(csv_file, delimiter=';')
        states_list = list(states_content)
        states_header = states_content.fieldnames

      # Get the energy of the first bright excited state

      bright_energies = [float(state['Energy (Ha)']) for state in states_list if state['Type'].lower() == "bright"]
      opt_gap = min(bright_energies)

      # Store data
      # ==========

      comp_results[mol_name].update({ "Energy gaps (Ha)" : 
          { "HOMO-LUMO" : hl_gap,
            "Optical" : opt_gap
          }
        })

      print('%12s' % "[ DONE ]")

      # ========================================================= #
      # Fetch ionization potentials                               #
      # ========================================================= #

      print ("{:<140}".format('\nFetching ionization potentials values ...'), end="")

      for line in ip_list:
        if line['Molecule'] == mol_name:
          comp_results[mol_name].update({ "IPs (Ha)" : 
            { "Koopmans" : float(line['IP (Koopmans)']),
              "Vertical" : float(line['IP (vertical)']) if line['IP (vertical)'] != "N/A" else line['IP (vertical)'],
              "Adiabatic" : float(line['IP (adiabatic)']) if line['IP (adiabatic)'] != "N/A" else line['IP (adiabatic)']
            }
          })
          break

      print('%12s' % "[ DONE ]")

    # ========================================================= #
    # Exception handling for the characterization results       #
    # ========================================================= #

    except results_errors.ResultsError as error:
      print(error)
      print("Skipping %s molecule" % mol_name)
      continue

    console_message = "End of procedure for the molecule " + mol_name
    print("")
    print(''.center(len(console_message)+10, '*'))
    print(console_message.center(len(console_message)+10))
    print(''.center(len(console_message)+10, '*'))

  # =================================================================== #
  # =================================================================== #
  #                      YAML FILE CREATION/UPDATE                      #
  # =================================================================== #
  # =================================================================== #

  # Create a custom class for dumping our data into the YAML format

  class CustomDumper(yaml.SafeDumper):

    # Insert blank lines between top-level objects (source: https://github.com/yaml/pyyaml/issues/127)
    def write_line_break(self, data=None):
      super().write_line_break(data)

      if len(self.indents) == 1:
        super().write_line_break()

    # Avoid writing anchors and aliases (source: https://ttl255.com/yaml-anchors-and-aliases-and-how-to-disable-them/)
    def ignore_aliases(self, data):
      return True

  # If the file already exists, update the data

  if os.path.exists(out_yml):

    with open(out_yml, 'r') as old_file:
      old_comp_results = yaml.load(old_file, Loader=yaml.FullLoader)

    old_comp_results.update(comp_results)
    comp_results = old_comp_results

  # Write the file (overwriting the old file if it already existed)

  with open(out_yml, 'w') as file:
    yaml.dump(comp_results, file, Dumper=CustomDumper, sort_keys=False)

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
