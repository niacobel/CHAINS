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
import math
import os
import re
import shutil
import sys
from inspect import getsourcefile

import numpy as np
import yaml
from scipy.spatial import ConvexHull, distance

import results_common

# =================================================================== #
# =================================================================== #
#                        FUNCTIONS DEFINITIONS                        #
# =================================================================== #
# =================================================================== #

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
      config_file = results_common.check_abspath(config_file,"Command line argument -cf / --config","file")
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
    chains_config_file = results_common.check_abspath(os.path.join(chains_path,"configs","chains_config.yml"),"CHAINS configuration YAML file","file")

    print ("{:<140}".format("\nLoading CHAINS configuration YAML file ..."), end="")
    with open(chains_config_file, 'r') as chains:
      chains_config = yaml.load(chains, Loader=yaml.FullLoader)
    print('%12s' % "[ DONE ]")

    # Ionization potentials CSV file
    # ==============================

    ip_file = results_common.check_abspath(chains_config['ip_file'],"Ionization potentials CSV file","file")

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

    # CONTROL LAUNCHER's source_parser.py file
    # =================================

    source_parser_path = os.path.join(chains_path,"control_launcher","source_parser.py")

    print ("{:<140}".format("\nImporting CONTROL LAUNCHER's source_parser.py file ..."), end="")
    source_parser = import_path(source_parser_path)
    print('%12s' % "[ DONE ]")

    # ========================================================= #
    # Check molecule directories                                #
    # ========================================================= #

    if multiple_mol:

      multiple_mol = results_common.check_abspath(multiple_mol,"Command line argument -m / --multiple","directory")
      mol_inp_path = multiple_mol

      print("{:<40} {:<99}".format("\nLooking for every molecule directory in", mol_inp_path + " ..."), end="")

      # We need to look for directories in the multiple_mol directory (see https://stackoverflow.com/questions/800197/how-to-get-all-of-the-immediate-subdirectories-in-python for reference).
      mol_inp_list = [dir.name for dir in os.scandir(mol_inp_path) if dir.is_dir()]

      if mol_inp_list == []:
        raise results_common.ResultsError ("ERROR: Can't find any directory in %s" % mol_inp_path)
      
      print('%12s' % "[ DONE ]")

    else:

      single_mol = results_common.check_abspath(single_mol,"Command line argument -s / --single","directory")
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

  except results_common.ResultsError as error:
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
    #             NAME & STRUCTURE OF THE MOLECULE              #
    # ========================================================= #
    # ========================================================= #

    # For more information on try/except structures, see https://www.tutorialsteacher.com/python/exception-handling-in-python
    try:

      print ("{:<140}".format('\nIdentifying name and structure of the molecule ...'), end="")

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
      opt_geom_file = results_common.check_abspath(os.path.join(gaussian_dir, mol_name + ".xyz"),"Optimized geometry file","file")

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

      if mol_group == 'SiC':
        if chemical_formula["H"] == 3*chemical_formula["C"]:
          mol_group = "SiC (100 %)"

      # ========================================================= #
      # Total number of atoms                                     #
      # ========================================================= #

      nb_atoms = sum(file_data['chemical_formula'].values())

      # ========================================================= #
      # "Pretty" name and "LaTeX" name                            #
      # ========================================================= #

      # If only 1 atom of that type, omit the number

      for atom,number in chemical_formula.items():
        if number == 1:
          chemical_formula[atom] = ""

      # Si on the front

      pretty_name_fr = "Si" + str(chemical_formula["Si"])
      latex_name_fr = "Si_{" + str(chemical_formula["Si"]) + "}"

      # H on the end (if present)

      if "H" in chemical_formula:
        pretty_name_end = "H" + str(chemical_formula["H"])
        latex_name_end = "H_{" + str(chemical_formula["H"]) + "}"
      else:
        pretty_name_end = ""
        latex_name_end = ""

      # Rest in the middle, in alphabetical order

      filtered_formula = {atom:number for atom, number in chemical_formula.items() if (atom != "Si" and atom != "H")}

      pretty_name_mid = ''.join([str(atom) + str(number) for atom, number in sorted(filtered_formula.items())])
      latex_name_mid = ''.join([str(atom) + "_{" + str(number) + "}" for atom, number in sorted(filtered_formula.items())])

      # Get the final name, correctly formatted

      pretty_name = pretty_name_fr + pretty_name_mid + pretty_name_end
      latex_name = latex_name_fr + latex_name_mid + latex_name_end

      print('%12s' % "[ DONE ]")

      # ========================================================= #
      # Size of the molecule                                      #
      # ========================================================= #

      # Create the list that will contain all coordinates

      coord_list = []

      # Define the pattern to extract coordinates from the XYZ file lines

      coord_pattern = re.compile(r'^\s*[a-zA-Z]{1,3}\s+(?P<coord_x>-?\d+\.\d+)\s+(?P<coord_y>-?\d+\.\d+)\s+(?P<coord_z>-?\d+\.\d+)\s*$')

      # Extract coordinates (and convert them from Angstrom to nanometers)

      for line in file_data['atomic_coordinates']:

        match = coord_pattern.match(line)

        x1 = float(match.group('coord_x'))/10
        y1 = float(match.group('coord_y'))/10
        z1 = float(match.group('coord_z'))/10

        coord_list.append([x1,y1,z1])

      # Compute the convex hull of the molecule and get the list of points consituting it

      points = np.array(coord_list)
      hull = ConvexHull(points)
      pts_hull = points[hull.vertices,:]

      # Compute the average distance between the points of the hull and their centroid, which will reflect the size of the molecule

      centroid = np.array(np.mean(pts_hull,axis=0),ndmin=2)
      dist = distance.cdist(pts_hull,centroid)
      size = 2*np.mean(dist)

      # Get the volume of the hull, as another way to reflect the size of the molecule

      volume = hull.volume

      # ========================================================= #
      # Store name and structure info                             #
      # ========================================================= #

      comp_results[mol_name] = { 
        "ID" : 
          { "TAG" : tag,
            "Group" : mol_group,
            "Name" : pretty_name,
            "LateX Name" : latex_name
          },
        "Structure":
          { "Nb atoms" : nb_atoms,
            "Size (nm)" : float(size),
            "Volume (nm^3)" : float(volume)
          }
        }

      # Add the number of Si atoms (for pure Si QDs only)

      if mol_group == "Si":
        comp_results[mol_name]['Structure'].update({"Nb Si atoms" : chemical_formula["Si"] if chemical_formula["Si"] != "" else 1})

    # ========================================================= #
    # Exception handling for the identification                 #
    # ========================================================= #

    except results_common.ResultsError as error:
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
    
      print ("{:<140}".format('\nCharacterization results'))

      # ========================================================= #
      # Load the QCHEM output file                                #
      # ========================================================= #

      print ("{:<133}".format('\n\tLoading QCHEM output file ...'), end="")

      qchem_dir = os.path.join(mol_dir, "QCHEM")
      qchem_file = results_common.check_abspath(os.path.join(qchem_dir, mol_name + ".out"),"QCHEM output file","file")

      with open(qchem_file, 'r') as out_file:
        qchem_content = out_file.read().splitlines()

      qchem_content = list(map(str.strip, qchem_content))   # Remove leading & trailing blank/spaces
      qchem_content = list(filter(None, qchem_content))     # Remove blank lines/no char

      print('%12s' % "[ DONE ]")

      # ========================================================= #
      # Get symmetry group                                        #
      # ========================================================= #

      print ("{:<133}".format('\n\tFetching molecular point group ...'), end="")

      # Initialize some variables

      sym_group = None
      sym_subgroup = None

      # Define the expression patterns for the lines containing information about the symmetry

      sym_rx = {
      
        # Pattern for finding the "Molecular Point Group                 Td    NOp = 24" type of line
        'group_line': re.compile(r'^\s*Molecular Point Group\s*(?P<group>[a-zA-Z0-9]+)\s*NOp\s*=\s*\d*\s*$'),

        # Pattern for finding the "Largest Abelian Subgroup              D2    NOp =  4" type of line
        'subgroup_line': re.compile(r'^\s*Largest Abelian Subgroup\s*(?P<subgroup>[a-zA-Z0-9]+)\s*NOp\s*=\s*\d*\s*$')

      }

      # Parse the qchem output file to get the information

      for line in qchem_content:

        # Store the data when encountered

        if sym_rx['group_line'].match(line):
          sym_group = sym_rx['group_line'].match(line).group('group')

        elif sym_rx['subgroup_line'].match(line):
          sym_subgroup = sym_rx['subgroup_line'].match(line).group('subgroup')
          break # This line will necessary follow the previous one so no need to go further

      # Raise an exception if the symmetry has not been found

      if not sym_group:
        raise results_common.ResultsError ("ERROR: Unable to find the molecule point group in the QCHEM output file")

      if not sym_subgroup:
        raise results_common.ResultsError ("ERROR: Unable to find the largest abelian subgroup in the QCHEM output file")

      # Add the data to the structure subdictionary for this molecule

      comp_results[mol_name]['Structure'].update({"Symmetry" : sym_group + " (" + sym_subgroup + ")"})

      print('%12s' % "[ DONE ]")

      # ========================================================= #
      # Fetch permanent dipole moment                             #
      # ========================================================= #

      print ("{:<133}".format('\n\tFetching permanent dipole moment value ...'), end="")

      # Initialize some variables

      section_found = False
      mu = "N/A"

      # Define the expression patterns for the lines containing information about the permanent dipole moment
    
      mu_rx = {

        # Pattern for finding the "Dipole Moment (Debye)" line (which marks the start of the section)
        'start': re.compile(r'^\s*Dipole Moment \(Debye\)\s*$'),

        # Pattern for finding lines looking like '       Tot       0.2618' (and capture the value)
        'value': re.compile(r'^\s*Tot\s+(?P<mu>\d+\.\d+)\s*$'),

        # Pattern for finding the "-----------------" line (which marks the end of the section)
        'end': re.compile(r'^\s*-+\s*$')

      }

      # Parse the qchem output file to get the information

      for line in qchem_content:

        # Define when the section begins and ends (ignore beta orbitals)

        if not section_found:
          if mu_rx['start'].match(line):
            section_found = True
      
        elif section_found and mu_rx['end'].match(line):
          break

        # Store the value

        elif mu_rx['value'].match(line):
          mu = float(mu_rx['value'].match(line).group('mu'))

      # Raise an exception if the section has not been found

      if not section_found:
        raise results_common.ResultsError ("ERROR: Unable to find the 'Dipole Moment (Debye)' line in the QCHEM output file")

      # Raise an exception if the energies of the orbitals have not been found

      if mu == "N/A":
        raise results_common.ResultsError ("ERROR: Unable to find the value of the permanent dipole moment in the QCHEM output file")

      # Store the value

      comp_results[mol_name]['Structure'].update({"Permanent dipole moment (Debye)" : mu})

      print('%12s' % "[ DONE ]")

      # ========================================================= #
      # Kohn-Sham orbitals                                        #
      # ========================================================= #

      print ("{:<133}".format('\n\tFetching KS orbitals values ...'), end="")

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
        raise results_common.ResultsError ("ERROR: Unable to find the 'Orbital Energies (a.u.)' section in the QCHEM output file")

      # Raise an exception if the energies of the orbitals have not been found

      if occ_orb == []:
        raise results_common.ResultsError ("ERROR: Unable to find the energies of the occupied orbitals in the QCHEM output file")

      if virt_orb == []:
        raise results_common.ResultsError ("ERROR: Unable to find the energies of the virtual orbitals in the QCHEM output file")

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

      print ("{:<133}".format('\n\tComputing energy gaps ...'), end="")

      # HOMO-LUMO gap
      # =============

      homo_value = max(occ_orb)
      lumo_value = min(virt_orb)
      hl_gap = lumo_value - homo_value

      # Optical gap
      # ===========

      # Call the parsing function from CONTROL LAUNCHER but without its standard output (https://stackoverflow.com/questions/2828953/silence-the-stdout-of-a-function-in-python-without-trashing-sys-stdout-and-resto)

      with open(os.devnull, 'w') as devnull:
        with contextlib.redirect_stdout(devnull):
          system = source_parser.qchem_tddft(qchem_content)

      states_list = system['states_list']

      # Get the list of bright states and sort them by ascending number

      bright_states = [state for state in states_list if state['type'].lower() == "bright"]
      bright_states.sort(key=lambda state: state['number'])

      # Get the energy of the first bright state (with a nonzero transition dipole moment)

      gs_number = [state['number'] for state in system['states_list'] if state['label'] == 'S0'][0]

      for bright_state in bright_states:

        momdip = 0
        for momdip_key in system['momdip_mtx']:
          momdip += system['momdip_mtx'][momdip_key][gs_number][bright_state['number']]**2
        momdip = math.sqrt(momdip)    

        if momdip != 0:
          opt_gap = bright_state['energy']
          break

      # Singlet-Triplet gap
      # ===================

      bright_energies = [float(state['energy']) for state in states_list if state['type'].lower() == "bright"]
      dark_energies = [float(state['energy']) for state in states_list if state['type'].lower() == "dark"]
      st_gap = abs( min(dark_energies) - min(bright_energies) )

      # Store data
      # ==========

      comp_results[mol_name].update({ "Energy gaps (Ha)" : 
          { "HOMO-LUMO" : hl_gap,
            "Optical" : opt_gap,
            "Singlet-Triplet" : st_gap
          }
        })

      print('%12s' % "[ DONE ]")

      # ========================================================= #
      # Fetch transition dipole moments                           #
      # ========================================================= #

      print ("{:<133}".format('\n\tFetching transition dipole moments with the ground state ...'), end="")

      # Initialize the dictionary that will contain the values

      comp_results[mol_name]["Transition dipole moments (au)"] = {}

      # Define the pairs of states that will be considered (here we consider the ground state with each of the first three non-degenerate singlets)

      pairs = []

      singlet_states = [state for state in states_list if state['label'].startswith('S') and state['label'] != 'S0']
      singlet_states.sort(key=lambda state: state['number'])

      prev_energy = 0

      for state in singlet_states:
        if state['energy'] == prev_energy:
          # Skip the singlet if it has the same energy as the previous one
          continue
        else:
          pairs.append(('S0',state['label']))
          prev_energy = state['energy']
          if len(pairs) == 3:
            # Stop after the third value
            break

      # Iterate over the pairs of states

      for pair in pairs:

        # Fetch the number of the involved states

        first_number = [state['number'] for state in system['states_list'] if state['label'] == pair[0]][0]
        second_number = [state['number'] for state in system['states_list'] if state['label'] == pair[1]][0]

        # Compute the module of the transition dipole moment

        momdip = 0
        for momdip_key in system['momdip_mtx']:
          momdip += system['momdip_mtx'][momdip_key][first_number][second_number]**2
        momdip = math.sqrt(momdip)          

        # Store the value

        comp_results[mol_name]["Transition dipole moments (au)"].update({
          pair[0] + '-' + pair[1] : momdip
          })

      print('%12s' % "[ DONE ]")      

      # ========================================================= #
      # Fetch ionization potentials                               #
      # ========================================================= #

      print ("{:<133}".format('\n\tFetching ionization potentials values ...'), end="")

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

    except results_common.ResultsError as error:
      print(error)
      print("Skipping %s molecule" % mol_name)
      continue

    # ========================================================= #
    # ========================================================= #
    #                      CONTROL RESULTS                      #
    # ========================================================= #
    # ========================================================= #

    # For more information on try/except structures, see https://www.tutorialsteacher.com/python/exception-handling-in-python
    try:
    
      print ("{:<140}".format('\nControl results'))

      # Check the data directory

      data_dir = results_common.check_abspath(os.path.join(mol_dir,"CONTROL","data"),"Data directory created by control_launcher.py","directory")

      # Load the eigenstates list

      eigenstates_file = results_common.check_abspath(os.path.join(data_dir, "eigenstates.csv"),"Eigenstates file","file")

      with open(eigenstates_file, 'r', newline='') as csv_file:

        eigenstates_content = csv.DictReader(csv_file, delimiter=';')
        eigenstates_list = list(eigenstates_content)

      # Check the eigenstates list

      required_keys = ['Number','Label','Energy (Ha)']
      results_common.check_keys(required_keys,eigenstates_list,"Eigenstates list file at %s" % eigenstates_file)

      # Store general information about the control system

      comp_results[mol_name].update({ 
          "Control" : 
            { "Nb excited states" : len(eigenstates_list) - 1,
              "Transitions" : []
            }
        })

      # Load the transitions list

      transitions_file = results_common.check_abspath(os.path.join(data_dir, "transitions.csv"),"Transitions file","file")

      with open(transitions_file, 'r', newline='') as csv_file:

        transitions_content = csv.DictReader(csv_file, delimiter=';')
        transitions_list = list(transitions_content)

      # Check the transitions list

      required_keys = ['Label','Transition dipole moments matrix','Initial file name','Target file name']
      results_common.check_keys(required_keys,transitions_list,"Transitions list file at %s" % transitions_file)

      # Load the transition dipole moment matrices

      momdip_es_mtx = {}
      for momdip_key in system['momdip_mtx']:
        file = results_common.check_abspath(os.path.join(data_dir, "momdip_es_mtx_" + momdip_key),"Transition dipole moment matrix associated with the %s key" % momdip_key,"file")
        momdip_es_mtx[momdip_key] = np.loadtxt(file)

      # Look for directories in the results CONTROL directory (see https://stackoverflow.com/questions/800197/how-to-get-all-of-the-immediate-subdirectories-in-python for reference).

      control_dir = os.path.join(mol_dir,"CONTROL")

      dir_list_all = [dir.name for dir in os.scandir(control_dir) if dir.is_dir()]

      # Only keep the 'transition_config'-type directories

      dir_list = []

      for transition in transitions_list:

        # Get the initial and target state number

        init_label = transition['Initial file name'][:-2] # [:-2] to remove the trailing "_1"
        init_number = eigenstates_list.index(next(eigenstate for eigenstate in eigenstates_list if eigenstate['Label'] == init_label))

        target_label = transition['Target file name'][:-2] # [:-2] to remove the trailing "_1"
        target_number = eigenstates_list.index(next(eigenstate for eigenstate in eigenstates_list if eigenstate['Label'] == target_label))

        # Get the orientation of the laser (momdip_key) and the transition dipole moment of this specific transition
        
        momdip_key = transition["Transition dipole moments matrix"]
        momdip = float(momdip_es_mtx[momdip_key][init_number][target_number])

        # Get the transition energy

        energy = abs(float(eigenstates_list[init_number]['Energy (Ha)']) - float(eigenstates_list[target_number]['Energy (Ha)']))

        # Iterate over the directories

        for dirname in dir_list_all:
        
          # For more information on try/except structures, see https://www.tutorialsteacher.com/python/exception-handling-in-python
          try:

            # Define the 'transition_config' regex and apply it to the dirname

            pattern = re.compile(r"^" + re.escape(transition["Label"]) + r"_(?P<config>.*)$")
            matching_dir = pattern.match(dirname)

            # If it is a '<transition>_<config>' directory, collect the data

            if matching_dir is not None:

              dir_list.append(dirname)
              config = pattern.match(dirname).group('config')

              print ("{:<133}".format("\n\tTreating the '%s' transition with the '%s' config ..." % (transition["Label"], config)), end="")

              # Open the PCP results file to get the fidelity of the pulse

              iter_file = results_common.check_abspath(os.path.join(control_dir, dirname, "PCP/obj.res"),"PCP Iterations QOCT-GRAD results file","file")

              with open(iter_file, 'r') as iter_f:
                iter_content = iter_f.read()

              # Define the expression patterns for the lines of the iterations file
              # For example "      0     1  1sec |Proba_moy  0.693654D-04 |Fidelity(U)  0.912611D-01 |Chp  0.531396D-04 -0.531399D-04 |Aire -0.202724D-03 |Fluence  0.119552D-03 |Recou(i)  0.693654D-04 |Tr_dist(i) -0.384547D-15 |Tr(rho)(i)  0.100000D+01 |Tr(rho^2)(i)  0.983481D+00 |Projector  0.100000D+01"
              rx_iter_line = re.compile(r"^\s+\d+\s+\d+\s+\d+sec\s\|Proba_moy\s+\d\.\d+D[+-]\d+\s\|Fidelity\(U\)\s+(?P<fidelity>\d\.\d+D[+-]\d+)\s\|Chp\s+\d\.\d+D[+-]\d+\s+-?\d\.\d+D[+-]\d+\s\|Aire\s+-?\d\.\d+D[+-]\d+\s\|Fluence\s+\d\.\d+D[+-]\d+\s\|Recou\(i\)\s+\d\.\d+D[+-]\d+\s\|Tr_dist\(i\)\s+-?\d\.\d+D[+-]\d+\s\|Tr\(rho\)\(i\)\s+\d\.\d+D[+-]\d+\s\|Tr\(rho\^2\)\(i\)\s+\d\.\d+D[+-]\d+\s\|Projector\s+\d\.\d+D[+-]\d+")

              # Get the fidelity

              iter_data = rx_iter_line.match(iter_content)
              if iter_data is not None:
                fidelity_raw = iter_data.group("fidelity")
                fidelity = float(re.compile(r'(\d*\.\d*)[dD]([-+]?\d+)').sub(r'\1E\2', fidelity_raw)) # Replace the possible d/D from Fortran double precision float format with an "E", understandable by Python)
              else:
                raise results_common.ResultsError ("ERROR: Unable to get the fidelity from the file %s" % iter_file) 

              # Store information specific to this transition

              comp_results[mol_name]['Control']['Transitions'].append({
                "Label" : transition["Label"],
                "Config" : config,
                "Initial state" : init_label, 
                "Target state" : target_label,
                "Energy (Ha)" : energy,
                "Orientation" : momdip_key,
                "Transition dipole moment (a.u.)" : momdip,
                "Fidelity" : fidelity
              })

              print('%12s' % "[ DONE ]")

          except results_common.ResultsError as error:
            print(error)
            print("Skipping %s directory" % dirname)
            continue

      if dir_list == []:
        raise results_common.ResultsError ("ERROR: Can't find any '<transition>_<config>' directory in %s" % os.path.join(mol_dir,"CONTROL"))
   
    # ========================================================= #
    #    Exception handling for the control results             #
    # ========================================================= #

    except results_common.ResultsError as error:
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

  # If the file already exists and is not empty, update the data

  if os.path.exists(out_yml) and os.stat(out_yml).st_size > 0:

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
