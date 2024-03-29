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
from itertools import combinations

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

    # CONTROL LAUNCHER's modelling_fcts.py file
    # =========================================

    modelling_fcts_path = os.path.join(chains_path,"control_launcher","modelling_fcts.py")

    print ("{:<140}".format("\nImporting CONTROL LAUNCHER's modelling_fcts.py file ..."), end="")
    modelling_fcts = import_path(modelling_fcts_path)
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

      if tag_finder is None:
        raise results_common.ResultsError ("ERROR: This directory does not seem to be a molecule directory.")

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

      if tag.startswith("B"):

        # Special bonus group is indicated by a tag that starts with B
        mol_group = 'Bonus'

      else:

        mol_group = ''.join(atoms)

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
      latex_name_fr = "Si$_{" + str(chemical_formula["Si"]) + "}$"

      # H on the end (if present)

      if "H" in chemical_formula:
        pretty_name_end = "H" + str(chemical_formula["H"])
        latex_name_end = "H$_{" + str(chemical_formula["H"]) + "}$"
      else:
        pretty_name_end = ""
        latex_name_end = ""

      # Rest in the middle, in alphabetical order

      filtered_formula = {atom:number for atom, number in chemical_formula.items() if (atom != "Si" and atom != "H")}

      pretty_name_mid = ''.join([str(atom) + str(number) for atom, number in sorted(filtered_formula.items())])
      latex_name_mid = ''.join([str(atom) + "$_{" + str(number) + "}$" for atom, number in sorted(filtered_formula.items())])

      # Get the final name, correctly formatted

      pretty_name = pretty_name_fr + pretty_name_mid + pretty_name_end
      latex_name = latex_name_fr + latex_name_mid + latex_name_end

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
      # Size of the molecule before geometry optimization         #
      # ========================================================= #

      # Check the original geometry file

      gaussian_dir = os.path.join(mol_dir, "GAUSSIAN")
      ori_geom_file = results_common.check_abspath(os.path.join(gaussian_dir, mol_name + "_ori.xyz"),"Original geometry file","file")

      # Load the XYZ geometry file

      with open(ori_geom_file, 'r') as mol_file:
        ori_mol_content = mol_file.read().splitlines()

      # Call the scanning function from ABIN LAUNCHER but without its standard output (https://stackoverflow.com/questions/2828953/silence-the-stdout-of-a-function-in-python-without-trashing-sys-stdout-and-resto)

      with open(os.devnull, 'w') as devnull:
        with contextlib.redirect_stdout(devnull):
          ori_file_data = geom_scan.xyz_scan(ori_mol_content)

      # Create the list that will contain all coordinates

      ori_coord_list = []

      # Define the pattern to extract coordinates from the XYZ file lines

      coord_pattern = re.compile(r'^\s*[a-zA-Z]{1,3}\s+(?P<coord_x>-?\d+\.\d+)\s+(?P<coord_y>-?\d+\.\d+)\s+(?P<coord_z>-?\d+\.\d+)\s*$')

      # Extract coordinates (and convert them from Angstrom to nanometers)

      for line in ori_file_data['atomic_coordinates']:

        match = coord_pattern.match(line)

        x1 = float(match.group('coord_x'))/10
        y1 = float(match.group('coord_y'))/10
        z1 = float(match.group('coord_z'))/10

        ori_coord_list.append([x1,y1,z1])

      # Compute the convex hull of the molecule and get the list of points consituting it

      points = np.array(ori_coord_list)
      hull = ConvexHull(points)
      pts_hull = points[hull.vertices,:]

      # Compute the average distance between the points of the hull and their centroid, which will reflect the size of the molecule

      centroid = np.array(np.mean(pts_hull,axis=0),ndmin=2)
      dist = distance.cdist(pts_hull,centroid)
      ori_size = 2*np.mean(dist)

      # Get the volume of the hull, as another way to reflect the size of the molecule

      ori_volume = hull.volume

      print('%12s' % "[ DONE ]")

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
            "Volume (nm^3)" : float(volume),
            "Original size (nm)" : float(ori_size),
            "Original volume (nm^3)" : float(ori_volume),
            "Coordinates" : file_data['atomic_coordinates']
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

      # Call the modelling function from CONTROL LAUNCHER but without its standard output (https://stackoverflow.com/questions/2828953/silence-the-stdout-of-a-function-in-python-without-trashing-sys-stdout-and-resto)

      with open(os.devnull, 'w') as devnull:
        with contextlib.redirect_stdout(devnull):
          system = modelling_fcts.qchem_tddft(qchem_content)

      zero_states_list = system['zero_states_list']

      # Get the list of singlet states and sort them by ascending number

      singlet_states = [state for state in zero_states_list if state['label'].startswith('S') and state['label'] != 'S0']
      singlet_states.sort(key=lambda state: state['number'])

      # Get the energy of the first bright state (with a nonzero transition dipole moment)

      gs_number = [state['number'] for state in zero_states_list if state['label'] == 'S0'][0]
      opt_gap = 0.0

      for singlet in singlet_states:

        momdip = 0
        for momdip_key in system['momdip_o_mtx']:
          momdip += system['momdip_o_mtx'][momdip_key][gs_number][singlet['number']]**2
        momdip = math.sqrt(momdip)    

        if momdip != 0:
          opt_gap = singlet['energy']
          break

      # Singlet-Triplet gap
      # ===================

      singlet_energies = [float(state['energy']) for state in zero_states_list if state['label'].startswith('S') and state['label'] != 'S0']
      triplet_energies = [float(state['energy']) for state in zero_states_list if state['label'].startswith('T')]
      st_gap = abs( min(triplet_energies) - min(singlet_energies) )

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
      # Fetch transition dipole moments                           #
      # ========================================================= #

      print ("{:<133}".format('\n\tFetching transition dipole moments ...'), end="")

      # Initialize the dictionary that will contain the values

      comp_results[mol_name]["Transition dipole moments (au)"] = {}

      # Define the pairs of states

      pairs = []

      singlet_states = [state for state in zero_states_list if state['label'].startswith('S')]
      singlet_states.sort(key=lambda state: state['number'])

      pairs.extend(list(combinations(singlet_states, 2))) # Look for all possible pairs of singlet states (see https://www.geeksforgeeks.org/python-all-possible-pairs-in-list/)

      triplet_states = []
      for state in zero_states_list:
        # Only add one substate of each triplet (no need to add ms=0 and ms=1 and ms=-1 as they all have the same transition dipole moments)
        if state['label'].startswith('T') and state['qchem_number'] not in [state['qchem_number'] for state in triplet_states]:
          triplet_states.append(state)
      triplet_states.sort(key=lambda state: state['number'])

      pairs.extend(list(combinations(triplet_states, 2)))

      # Iterate over the pairs of states

      for pair in pairs:

        # Compute the module of the transition dipole moment

        momdip = 0
        for momdip_key in system['momdip_o_mtx']:
          momdip += system['momdip_o_mtx'][momdip_key][pair[0]['number']][pair[1]['number']]**2
        momdip = math.sqrt(momdip)          

        # Store the value

        comp_results[mol_name]["Transition dipole moments (au)"].update({
          # Use partition to only keep the 'Tx' part of the 'Tx(ms=y)' labels
          pair[0]['label'].partition("(")[0] + '_' + pair[1]['label'].partition("(")[0] : momdip
          })

      print('%12s' % "[ DONE ]")

      # ========================================================= #
      # Fetch permanent dipole moments                            #
      # ========================================================= #

      print ("{:<133}".format('\n\tFetching permanent dipole moments ...'), end="")

      # Initialize the dictionary that will contain the values

      comp_results[mol_name]["Permanent dipole moments (au)"] = {}

      # Iterate over the states

      for state in zero_states_list:

        # Compute the module of the permanent dipole moment

        momdip = 0
        for momdip_key in system['momdip_o_mtx']:
          momdip += system['momdip_o_mtx'][momdip_key][state['number']][state['number']]**2
        momdip = math.sqrt(momdip)          

        # Store the value

        comp_results[mol_name]["Permanent dipole moments (au)"].update({
          # Use partition to only keep the 'Tx' part of the 'Tx(ms=y)' labels
          state['label'].partition("(")[0] : momdip
          })

      print('%12s' % "[ DONE ]")

      # ========================================================= #
      # Fetch spin-orbit couplings                                #
      # ========================================================= #

      print ("{:<133}".format('\n\tFetching SOC values ...'), end="")

      # Initialize the dictionary that will contain the values

      comp_results[mol_name]["SOC values (au)"] = {}

      # Iterate over the pairs of states

      pairs = [(zero_states_list[s1], zero_states_list[s2]) for s1 in range(len(zero_states_list)) for s2 in range(s1+1,len(zero_states_list))]
      pairs = [pair for pair in pairs if pair[0]['label'].startswith('T') or pair[1]['label'].startswith('T')]

      for pair in pairs:
        
        # Compute the module of the soc value

        soc = float(np.absolute(system['mime'][pair[0]['number']][pair[1]['number']]))
    
        # Store the value

        comp_results[mol_name]["SOC values (au)"].update({
          # Use partition to only keep the 'Tx' part of the 'Tx(ms=y)' labels
          pair[0]['label'] + '_' + pair[1]['label'] : soc
          })

      print('%12s' % "[ DONE ]")

      # ========================================================= #
      # Fetch zero states list                                    #
      # ========================================================= #

      print ("{:<133}".format('\n\tFetching zero states list ...'), end="")

      # Initialize the dictionary that will contain the values

      comp_results[mol_name]["Zero states list"] = {}

      # Filter the states (only one substate of each triplet)

      min_zero_states_list = []

      for state in zero_states_list:
        # Only add one substate of each triplet (no need to add ms=0 and ms=1 and ms=-1 as they all have the same transition dipole moments)
        if state['label'].startswith('T') and state['qchem_number'] not in [state['qchem_number'] for state in min_zero_states_list]:
          min_zero_states_list.append(state)
        elif state['label'].startswith('S'):
          min_zero_states_list.append(state)

      # Iterate over the states

      for state in min_zero_states_list:
      
        # Define the values
        # ~~~~~~~~~~~~~~~~~

        # Define important variables

        zero_state = {}
        sing_numbers = [state['number'] for state in zero_states_list if state['number'] != 0 and state['label'].startswith('S')]

        # Get the energies

        zero_state["Energy (Ha)"] = state['energy']

        # For singlet excited states, get the GS transition dipole moment

        if state['number'] != 0 and state['label'].startswith('S'):

          zero_state["GS transition dipole moment (au)"] = {
              "X" : float(system['momdip_o_mtx']["X"][0][state['number']]),
              "Y" : float(system['momdip_o_mtx']["Y"][0][state['number']]),
              "Z" : float(system['momdip_o_mtx']["Z"][0][state['number']])
              }

        # For triplet excited states, get the SOC values

        if state['label'].startswith('T'):

          # /!\ The order of those values is the same as the order defined in the modelling function (which created the zero_states_list)

          zero_state["GS SOC values (au)"] = {
              "ms=0" :  float(np.absolute(system['mime'][0][state['number']])),
              "ms=1" :  float(np.absolute(system['mime'][0][state['number']+1])),
              "ms=-1" : float(np.absolute(system['mime'][0][state['number']+2]))
              }

          zero_state["Average singlet SOC values (au)"] = {
              "ms=0" :  float(np.mean([np.absolute(system['mime'][sing_number][state['number']]) for sing_number in sing_numbers])),
              "ms=1" :  float(np.mean([np.absolute(system['mime'][sing_number][state['number']+1]) for sing_number in sing_numbers])),
              "ms=-1" : float(np.mean([np.absolute(system['mime'][sing_number][state['number']+2]) for sing_number in sing_numbers]))
              }

        # Store the values

        comp_results[mol_name]["Zero states list"].update({state['label'].partition("(")[0] : zero_state})

      print('%12s' % "[ DONE ]") 

      # ========================================================= #
      # Fetch relativistic states list                            #
      # ========================================================= #

      print ("{:<133}".format('\n\tFetching relativistic states list ...'), end="")

      # Initialize the dictionary that will contain the values

      comp_results[mol_name]["Relativistic states list"] = {}

      # Iterate over the states

      for state in system['states_list']:
      
        # Store the value

        comp_results[mol_name]["Relativistic states list"].update({
          state['label'] : {
            "Energy (Ha)": state['energy'],
            "GS percentage": float(state['gs_percent']),
            "Singlet percentage": float(state['sing_percent']),
            "Triplet percentage": float(state['trip_percent']),
            "GS transition dipole moment (au)" : {
              "X" : float(system['momdip_mtx']["X"][0][state['number']]) if state['number'] != 0 else None,
              "Y" : float(system['momdip_mtx']["Y"][0][state['number']]) if state['number'] != 0 else None,
              "Z" : float(system['momdip_mtx']["Z"][0][state['number']]) if state['number'] != 0 else None
              }
            }
          })

      print('%12s' % "[ DONE ]") 

      # ========================================================= #
      # Compute average and sum of non relavistic dipole moments  #
      # ========================================================= #

      print ("{:<133}".format('\n\tComputing average and sum of non relavistic dipole moments ...'), end="")

      # Initialization of the variables

      momdip_gs_s_list = []
      momdip_s_s_list = []
      momdip_t_t_list = []
      momdip_p_list = []
      
      # Iterate over the matrix to get the information

      it = np.nditer(system['momdip_o_mtx']['X'], flags=['multi_index'])
    
      for momdip in it:

        # Define the pair of states

        state_1 = it.multi_index[0]
        label_1 = [state['label'] for state in system['zero_states_list'] if state_1 == state['number']][0]
        energy_1 = [state['energy'] for state in system['zero_states_list'] if state_1 == state['number']][0]
        state_2 = it.multi_index[1]
        label_2 = [state['label'] for state in system['zero_states_list'] if state_2 == state['number']][0]
        energy_2 = [state['energy'] for state in system['zero_states_list'] if state_2 == state['number']][0]

        if state_1 > state_2:
          continue # No need to treat both triangles of the matrix

        # Compute the total length of the vector

        value = 0
        for momdip_key in system['momdip_o_mtx']:
          value += system['momdip_o_mtx'][momdip_key][state_1][state_2]**2
        value = math.sqrt(value)

        # Add the value to the corresponding list

        if state_1 == state_2:
          momdip_p_list.append(value)
        elif label_1.startswith('S0') and label_2.startswith('S'):
          momdip_gs_s_list.append(value)
        elif label_1.startswith('S') and label_2.startswith('S'):
          if energy_1 != energy_2:
            momdip_s_s_list.append(value)
        elif label_1.startswith('T') and label_2.startswith('T'):
          if energy_1 != energy_2:
            momdip_t_t_list.append(value)
  
      # Compute the total matrix

      total_mtx = np.zeros((len(system['zero_states_list']), len(system['zero_states_list'])))
      for momdip_key in system['momdip_o_mtx']:
        total_mtx = np.add(total_mtx,np.square(system['momdip_o_mtx'][momdip_key]))
      total_mtx = np.sqrt(total_mtx)

      # Store the values

      comp_results[mol_name]["Average NR dipole moments (au)"] = {
        "Permanent" : float(np.mean(momdip_p_list)),
        "GS-S" : float(np.mean(momdip_gs_s_list)),
        "S-S" : float(np.mean(momdip_s_s_list)),
        "T-T" : float(np.mean(momdip_t_t_list)),
        "Total" : float(total_mtx.mean())
        }

      comp_results[mol_name]["Sum NR dipole moments (au)"] = {
        "Permanent" : float(sum(momdip_p_list)),
        "GS-S" : float(sum(momdip_gs_s_list)),
        "S-S" : float(sum(momdip_s_s_list)),
        "T-T" : float(sum(momdip_t_t_list)),
        "Total" : float(total_mtx.sum())
        }

      print('%12s' % "[ DONE ]")

      # ========================================================= #
      # Compute average and sum of relavistic dipole moments      #
      # ========================================================= #

      print ("{:<133}".format('\n\tComputing average and sum of relavistic dipole moments ...'), end="")

      # Initialization of the variables

      momdip_gs_s_list = []
      momdip_gs_t_list = []
      momdip_s_s_list = []
      momdip_s_t_list = []
      momdip_t_t_list = []
      momdip_p_list = []
      
      # Iterate over the matrix to get the information

      it = np.nditer(system['momdip_mtx']['X'], flags=['multi_index'])
    
      for momdip in it:

        # Define the pair of states and their dominant multiplicity

        state_1 = it.multi_index[0]
        if state_1 == 0:
          multi_1 = 'GS'
        elif system['states_list'][state_1]['trip_percent'] < 0.5:
          multi_1 = 'S'
        else:
          multi_1 = 'T'
        
        state_2 = it.multi_index[1]
        if state_2 == 0:
          multi_2 = 'GS'
        elif system['states_list'][state_2]['trip_percent'] < 0.5:
          multi_2 = 'S'
        else:
          multi_2 = 'T'

        if state_1 > state_2:
          continue # No need to treat both triangles of the matrix

        # Compute the total length of the vector

        value = 0
        for momdip_key in system['momdip_mtx']:
          value += system['momdip_mtx'][momdip_key][state_1][state_2]**2
        value = math.sqrt(value) 

        # Add the value to the corresponding list

        if state_1 == state_2:
          momdip_p_list.append(value)
        elif multi_1 == 'GS' and multi_2 == 'S':
          momdip_gs_s_list.append(value)
        elif multi_1 == 'GS' and multi_2 == 'T':
          momdip_gs_t_list.append(value)
        elif multi_1 == 'S' and multi_2 == 'S':
          momdip_s_s_list.append(value)
        elif multi_1 == 'T' and multi_2 == 'T':
          momdip_t_t_list.append(value)
        elif (multi_1 == 'S' and multi_2 == 'T') or (multi_1 == 'T' and multi_2 == 'S'):
          momdip_s_t_list.append(value)

      # Compute the total matrix

      total_mtx = np.zeros((len(system['states_list']), len(system['states_list'])))
      for momdip_key in system['momdip_mtx']:
        total_mtx = np.add(total_mtx,np.square(system['momdip_mtx'][momdip_key]))
      total_mtx = np.sqrt(total_mtx)

      # Store the values

      comp_results[mol_name]["Average relativistic dipole moments (au)"] = {
        "Permanent" : float(np.mean(momdip_p_list)),
        "GS-S" : float(np.mean(momdip_gs_s_list)),
        "GS-T" : float(np.mean(momdip_gs_t_list)),
        "S-S" : float(np.mean(momdip_s_s_list)),
        "S-T" : float(np.mean(momdip_s_t_list)),
        "T-T" : float(np.mean(momdip_t_t_list)),
        "Total" : float(total_mtx.mean())
        }

      comp_results[mol_name]["Sum relativistic dipole moments (au)"] = {
        "Permanent" : float(sum(momdip_p_list)),
        "GS-S" : float(sum(momdip_gs_s_list)),
        "GS-T" : float(sum(momdip_gs_t_list)),
        "S-S" : float(sum(momdip_s_s_list)),
        "S-T" : float(sum(momdip_s_t_list)),
        "T-T" : float(sum(momdip_t_t_list)),
        "Total" : float(total_mtx.sum())
        }

      # Add the total average in each direction

      for momdip_key in system['momdip_mtx']:
        mtx_dir = system['momdip_mtx'][momdip_key]
        np.fill_diagonal(mtx_dir, 0)
        comp_results[mol_name]["Average relativistic dipole moments (au)"].update({
          momdip_key : float(mtx_dir.mean())
        })

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

      # ========================================================= #
      # Alpha / Duration parameters search                        #
      # ========================================================= #

      print ("{:<133}".format('\n\tFetching alpha and duration parameters search results ...'))      

      # Check the main directory

      aldu_dir = results_common.check_abspath(os.path.join(mol_dir,"CONTROL","aldu_param"),"Main directory for the alpha and duration parameters search","directory")

      # Check the data directory

      data_dir = results_common.check_abspath(os.path.join(aldu_dir,"data"),"Data directory created by control_launcher.py","directory")

      # Initialize the results dictionary

      comp_results[mol_name].update({ 
          "Control" : 
            { "Alpha/duration parameters search" : {}
            }
        })

      # Look for directories in the main directory (see https://stackoverflow.com/questions/800197/how-to-get-all-of-the-immediate-subdirectories-in-python for reference).

      dir_list_all = [dir.name for dir in os.scandir(aldu_dir) if dir.is_dir()]

      # Iterate over the directories

      for dirname in dir_list_all:
      
        # For more information on try/except structures, see https://www.tutorialsteacher.com/python/exception-handling-in-python
        try:

          # Define the '<transition>' regex and apply it to the dirname

          trans_pattern = re.compile(r"^(?P<key>[XYZ])_(?P<init>[pST0-9-]+)_(?P<target>[pST0-9-]+)_*([a-zA-Z0-9-_]*)$")
          matching_dir = trans_pattern.match(dirname)

          # If it is a directory named after a transition, collect the data

          if matching_dir is None:
            continue

          trans_dir = os.path.join(aldu_dir, dirname)
          momdip_key = matching_dir.group('key')

          print ("{:<126}".format("\n\t\tTreating the '%s' transition ..." % dirname), end="")

          # Create an entry for this transition

          comp_results[mol_name]['Control']['Alpha/duration parameters search'].update({ 
              dirname : 
                { "Polarisation" : momdip_key,
                  "Initial state" : matching_dir.group('init'),
                  "Target state" : matching_dir.group('target'),
                  "Values" : []
                }
            })

          # Load the CSV file

          aldu_file = results_common.check_abspath(os.path.join(trans_dir,"aldu_comp_results.csv"),"Alpha/duration parameters search compiled results CSV file","file")
          with open(aldu_file, 'r', newline='') as csv_file:
            aldu_content = csv.DictReader(csv_file, delimiter=';')
            aldu_list = list(aldu_content)
      
          # Store information specific to this calculation

          for line in aldu_list:
            for key in line.keys():
              if key == 'Best':
                line['Best'] = str(line['Best'])
              else:
                line[key] = float(line[key])

            # Add dict(line) in order to ensure line is not an OrderedDict, incompatible with YAML
            comp_results[mol_name]['Control']['Alpha/duration parameters search'][dirname]['Values'].append(dict(line))

          print('%12s' % "[ DONE ]")

        except results_common.ResultsError as error:
          print(error)
          print("Skipping %s directory" % dirname)
          continue

      # ========================================================= #
      # Constraints variation                                     #
      # ========================================================= #

      print ("{:<133}".format('\n\tFetching constraints variation results ...'))      

      try:

        # Check the main directory

        convar_dir = results_common.check_abspath(os.path.join(mol_dir,"CONTROL","const_var"),"Main directory for the constraints variation","directory")

        # Check the data directory

        data_dir = results_common.check_abspath(os.path.join(convar_dir,"data"),"Data directory created by control_launcher.py","directory")

        # Initialize the results dictionary

        comp_results[mol_name]['Control'].update({
          "Constraints variation" : {}
          })

        # Look for directories in the main directory (see https://stackoverflow.com/questions/800197/how-to-get-all-of-the-immediate-subdirectories-in-python for reference).

        dir_list_all = [dir.name for dir in os.scandir(convar_dir) if dir.is_dir()]

        # Iterate over the directories

        for dirname in dir_list_all:
        
          # For more information on try/except structures, see https://www.tutorialsteacher.com/python/exception-handling-in-python
          try:

            # Define the '<transition>' regex and apply it to the dirname

            trans_pattern = re.compile(r"^(?P<key>[XYZ])_(?P<init>[pST0-9-]+)_(?P<target>[pST0-9-]+)_*([a-zA-Z0-9-_]*)$")
            matching_dir = trans_pattern.match(dirname)

            # If it is a directory named after a transition, collect the data

            if matching_dir is None:
              continue

            trans_dir = os.path.join(convar_dir, dirname)
            momdip_key = matching_dir.group('key')

            print ("{:<126}".format("\n\t\tTreating the '%s' transition ..." % dirname), end="")

            # Create an entry for this transition

            comp_results[mol_name]['Control']['Constraints variation'].update({ 
                dirname : 
                  { "Polarisation" : momdip_key,
                    "Values" : []
                  }
              })

            # Load the CSV file

            convar_file = results_common.check_abspath(os.path.join(trans_dir,"convar_comp_results.csv"),"Constraints variation compiled results CSV file","file")
            with open(convar_file, 'r', newline='') as csv_file:
              convar_content = csv.DictReader(csv_file, delimiter=';')
              convar_list = list(convar_content)
        
            # Store information specific to this calculation

            for line in convar_list:
              for key in line.keys():
                if key == 'Best':
                  line['Best'] = str(line['Best'])
                else:
                  line[key] = float(line[key])

              # Add dict(line) in order to ensure line is not an OrderedDict, incompatible with YAML
              comp_results[mol_name]['Control']['Constraints variation'][dirname]['Values'].append(dict(line))

            print('%12s' % "[ DONE ]")

          except results_common.ResultsError as error:
            print(error)
            print("Skipping %s directory" % dirname)
            continue

      except results_common.ResultsError as error:
        print(error)
        print("Skipping constraints variation treatment")

      # ========================================================= #
      # Frequency filters variation                               #
      # ========================================================= #

      print ("{:<133}".format('\n\tFetching frequency filters variation results ...'))      

      # Check the main directory

      filt_freq_dir = results_common.check_abspath(os.path.join(mol_dir,"CONTROL","filt_freq"),"Main directory for the frequency filters variation","directory")

      # Check the data directory

      data_dir = results_common.check_abspath(os.path.join(filt_freq_dir,"data"),"Data directory created by control_launcher.py","directory")

      # Initialize the results dictionary

      comp_results[mol_name]['Control'].update({
        "Frequency filters variation" : {}
        })

      # Look for directories in the main directory (see https://stackoverflow.com/questions/800197/how-to-get-all-of-the-immediate-subdirectories-in-python for reference).

      dir_list_all = [dir.name for dir in os.scandir(filt_freq_dir) if dir.is_dir()]

      # Iterate over the directories

      for dirname in dir_list_all:
      
        # For more information on try/except structures, see https://www.tutorialsteacher.com/python/exception-handling-in-python
        try:

          # Define the '<transition>' regex and apply it to the dirname

          trans_pattern = re.compile(r"^(?P<key>[XYZ])_(?P<init>[pST0-9-]+)_(?P<target>[pST0-9-]+)_*([a-zA-Z0-9-_]*)$")
          matching_dir = trans_pattern.match(dirname)

          # If it is a directory named after a transition, collect the data

          if matching_dir is None:
            continue

          trans_dir = os.path.join(filt_freq_dir, dirname)
          momdip_key = matching_dir.group('key')

          print ("{:<126}".format("\n\t\tTreating the '%s' transition ..." % dirname), end="")

          # Create an entry for this transition

          comp_results[mol_name]['Control']['Frequency filters variation'].update({ 
              dirname : 
                { "Polarisation" : momdip_key,
                  "Values" : []
                }
            })

          # Load the CSV file

          filt_freq_file = results_common.check_abspath(os.path.join(trans_dir,"filt_freq_comp_results.csv"),"Frequency filters variation compiled results CSV file","file")
          with open(filt_freq_file, 'r', newline='') as csv_file:
            filt_freq_content = csv.DictReader(csv_file, delimiter=';')
            filt_freq_list = list(filt_freq_content)
      
          # Store information specific to this calculation

          for line in filt_freq_list:
            for key in line.keys():
              if key == 'Best':
                line['Best'] = str(line['Best'])
              else:
                line[key] = float(line[key])

            # Add dict(line) in order to ensure line is not an OrderedDict, incompatible with YAML
            comp_results[mol_name]['Control']['Frequency filters variation'][dirname]['Values'].append(dict(line))

          print('%12s' % "[ DONE ]")

        except results_common.ResultsError as error:
          print(error)
          print("Skipping %s directory" % dirname)
          continue

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

    print ("{:<140}".format('\nUpdating the output YAML file ...'), end="")
    with open(out_yml, 'r') as old_file:
      old_comp_results = yaml.load(old_file, Loader=yaml.FullLoader)

    old_comp_results.update(comp_results)
    comp_results = old_comp_results
  
  else:

    print ("{:<140}".format('\nCreating the output YAML file ...'), end="")

  # Write the file (overwriting the old file if it already existed)

  with open(out_yml, 'w') as file:
    yaml.dump(comp_results, file, Dumper=CustomDumper, sort_keys=False)

  print('%12s' % "[ DONE ]")

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
