#!/usr/bin/env python3

################################################################################################################################################
##                                                          Results Plotting Script                                                           ##
##                                                                                                                                            ##
##                                   This script extracts the relevant information from a single YAML file                                    ##
##                                             and produces the various needed tables and graphs.                                             ##
##                                                                                                                                            ##
##                       /!\ In order to run, this script requires Python 3.5+ as well as YAML, matplotlib, Jinja2 /!\                        ##
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
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

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

# =================================================================== #
# =================================================================== #
#                       COMMAND LINE ARGUMENTS                        #
# =================================================================== #
# =================================================================== #

# Define the arguments needed for the script (here they are defined as named arguments rather than positional arguments, check https://stackoverflow.com/questions/24180527/argparse-required-arguments-listed-under-optional-arguments for more info).

parser = argparse.ArgumentParser(add_help=False, description="This script extracts the relevant information from a single YAML file and produces the various needed tables and graphs.")

required = parser.add_argument_group('Required arguments')
required.add_argument("-i","--inp_yml", type=str, help="Path to the YAML file containing the results.", required=True)
required.add_argument("-o","--out_dir", type=str, help="Path to the directory where you want to store the graphs and tables.", required=True)

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
    print("EXECUTION OF THE RESULTS PLOTTING SCRIPT BEGINS NOW".center(columns))
    print("")
    print("".center(columns,"*"))

    # ========================================================= #
    # Read command line arguments                               #
    # ========================================================= #

    args = parser.parse_args()

    # Required arguments

    inp_yml = args.inp_yml                   # YAML file containing the results
    out_dir = args.out_dir                   # Directory where the graphs and tables will be stored

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
    # Check and load the YAML input file                        #
    # ========================================================= #

    inp_yml = results_errors.check_abspath(inp_yml,"Command line argument -i / --inp_yml","file")

    print ("{:<40} {:<99}".format('\nLoading the YAML input file',inp_yml + " ..."), end="")
    with open(inp_yml, 'r') as f_yml:
      yml_content = yaml.load(f_yml, Loader=yaml.FullLoader)
    print('%12s' % "[ DONE ]")

    # ========================================================= #
    # Check other arguments                                     #
    # ========================================================= #

    out_dir = results_errors.check_abspath(out_dir,"Command line argument -o / --out_dir","directory")
    print ("{:<40} {:<100}".format('\nOutput directory:',out_dir))

  # ========================================================= #
  # Exception handling for the preparation step               #
  # ========================================================= #

  except results_errors.ResultsError as error:
    print("")
    print(error)
    exit(-1)

  # =================================================================== #
  # =================================================================== #
  #                          SORTING THE DATA                           #
  # =================================================================== #
  # =================================================================== #

  # For more information on try/except structures, see https://www.tutorialsteacher.com/python/exception-handling-in-python
  try:

    print ("{:<140}".format('\nSorting the data ...'), end="")

    # Get all the different types of molecules

    mol_types = [yml_content[mol]['ID']['Type'] for mol in yml_content]
    mol_types = list(dict.fromkeys(mol_types)) # Remove duplicates (https://www.w3schools.com/python/python_howto_remove_duplicates.asp)

    # Regroup the data by type (constitutive atoms) (#!only include spheroids for now)

    yml_sorted = {}

    for mol_type in mol_types:
      yml_sorted[mol_type] = {mol:value for (mol,value) in yml_content.items() if yml_content[mol]['ID']['Type'] == mol_type and yml_content[mol]['ID']['TAG'].startswith("S")}

    # Separate the reference (Si) from the other types

    si_data = yml_sorted.pop("Si")
    mol_types.remove("Si")

    print('%12s' % "[ DONE ]")

  # ========================================================= #
  # Exception handling for sorting the data                   #
  # ========================================================= #

  except results_errors.ResultsError as error:
    print("")
    print(error)
    exit(-1)

  # =================================================================== #
  # =================================================================== #
  #                        PLOTTING ENERGY GAPS                         #
  # =================================================================== #
  # =================================================================== #

  # For more information on try/except structures, see https://www.tutorialsteacher.com/python/exception-handling-in-python
  try:

    section_title = "1. Energy gaps"

    print("")
    print(''.center(len(section_title)+10, '*'))
    print(section_title.center(len(section_title)+10))
    print(''.center(len(section_title)+10, '*'))

    # Get the values for Si
    # =====================

    print ("{:<140}".format('\nFetching the values for Si QDs ...'), end="")

    si_gaps = []

    for mol in si_data:
      if si_data[mol].get('Energy_gaps'):
        nb_si_atoms = si_data[mol]['ID']['Nb_Si_atoms']
        hl_gap = energy_unit_conversion(si_data[mol]['Energy_gaps']['HOMO-LUMO'],"ha","ev")
        opt_gap = energy_unit_conversion(si_data[mol]['Energy_gaps']['Optical'],"ha","ev")
        si_gaps.append((nb_si_atoms,hl_gap,opt_gap))

    si_gaps.sort(key=lambda tup: tup[0]) # Sort the values by number of Si atoms

    print('%12s' % "[ DONE ]")

    # Iterate over each molecule type and compare it to Si
    # ====================================================

    for mol_type in mol_types:

      print ("{:<140}".format('\nTreating the values for %s QDs ...' % mol_type), end="")

      # Get the values
      # ==============

      gaps = []

      for mol in yml_sorted[mol_type]:
        if yml_sorted[mol_type][mol].get('Energy_gaps'):

          # Get the gaps values

          hl_gap = energy_unit_conversion(yml_sorted[mol_type][mol]['Energy_gaps']['HOMO-LUMO'],"ha","ev")
          opt_gap = energy_unit_conversion(yml_sorted[mol_type][mol]['Energy_gaps']['Optical'],"ha","ev")

          # Use the TAG key to find the corresponding Si QD and fetch its number of Si atoms

          for si_mol in si_data:
            if si_data[si_mol]['ID']['TAG'] == yml_sorted[mol_type][mol]['ID']['TAG']:
              nb_si_atoms = si_data[si_mol]['ID']['Nb_Si_atoms']

          # Store the data for this molecule

          gaps.append((nb_si_atoms,hl_gap,opt_gap))

      gaps.sort(key=lambda tup: tup[0]) # Sort the values by number of Si atoms

      # Plot the graphs
      # ===============

      plt.style.use('seaborn-colorblind')

      fig, ax = plt.subplots()

      # Plot the Si values

      ax.plot([mol[0] for mol in si_gaps],[mol[1] for mol in si_gaps],marker='.',linestyle='--',label='H-L gaps - Si (DFT)')
      ax.plot([mol[0] for mol in si_gaps],[mol[2] for mol in si_gaps],marker='^',markersize=4,linestyle='--',label='Optical gaps - Si (TD-DFT)')

      # Plot the specific type value

      ax.plot([mol[0] for mol in gaps],[mol[1] for mol in gaps],marker='.',linestyle='-',label='H-L gaps - %s (DFT)' % mol_type)
      ax.plot([mol[0] for mol in gaps],[mol[2] for mol in gaps],marker='^',markersize=4,linestyle='-',label='Optical gaps - %s (TD-DFT)' % mol_type)

      # Add the legend and titles

      ax.set_title('Energy gaps: Si vs %s' % mol_type)
      ax.set_xlabel("Number of Si atoms in the Si QD")
      ax.set_ylabel('Energy (eV)')
      ax.legend()

      # Set other parameters

      ax.tick_params(top=False, right=False)
      ax.xaxis.set_minor_locator(AutoMinorLocator(2))
      ax.yaxis.set_minor_locator(AutoMinorLocator(2))

      plt.tight_layout()
      plt.grid(True,which='both',linestyle='--')

      # Save the figure and clear the axes to prepare for the next plot

      plt.savefig(os.path.join(out_dir,'gaps','gaps_Si_vs_%s.png' % mol_type),dpi=200)
      plt.cla()

      print('%12s' % "[ DONE ]")

      print(mol_type,gaps)
   
  # ========================================================= #
  # Exception handling for plotting the energy gaps           #
  # ========================================================= #

  except results_errors.ResultsError as error:
    print("")
    print(error)
    exit(-1)
