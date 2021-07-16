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
import os
import shutil
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

  print ("{:<140}".format('\nSorting the data ...'), end="")

  # Get all the different groups of molecules

  mol_groups = [yml_content[mol]['ID']['Group'] for mol in yml_content]
  mol_groups = list(dict.fromkeys(mol_groups)) # Remove duplicates (https://www.w3schools.com/python/python_howto_remove_duplicates.asp)

  # Regroup the data by group (constitutive atoms) (#!only include spheroids for now)

  yml_sorted = {}

  for mol_group in mol_groups:
    yml_sorted[mol_group] = {mol:value for (mol,value) in yml_content.items() if yml_content[mol]['ID']['Group'] == mol_group and yml_content[mol]['ID']['TAG'].startswith("S")}

  print('%12s' % "[ DONE ]")

  # =================================================================== #
  # =================================================================== #
  #                        PLOTTING ENERGY GAPS                         #
  # =================================================================== #
  # =================================================================== #

  section_title = "1. Energy gaps"

  print("")
  print(''.center(len(section_title)+10, '*'))
  print(section_title.center(len(section_title)+10))
  print(''.center(len(section_title)+10, '*'))

  # Create the directory that will contain the graphs
  # =================================================

  gaps_dir = os.path.join(out_dir,'gaps')
  os.makedirs(gaps_dir, exist_ok = True)

  # Get the values for Si
  # =====================

  print ("{:<140}".format('\nFetching the values for Si QDs ...'), end="")

  si_gaps = []

  for mol in yml_sorted['Si']:
    if yml_sorted['Si'][mol].get('Energy gaps (Ha)'):
      nb_si_atoms = yml_sorted['Si'][mol]['ID']['Nb Si atoms']
      hl_gap = energy_unit_conversion(yml_sorted['Si'][mol]['Energy gaps (Ha)']['HOMO-LUMO'],"ha","ev")
      opt_gap = energy_unit_conversion(yml_sorted['Si'][mol]['Energy gaps (Ha)']['Optical'],"ha","ev")
      si_gaps.append((nb_si_atoms,hl_gap,opt_gap))

  si_gaps.sort(key=lambda tup: tup[0]) # Sort the values by number of Si atoms

  print('%12s' % "[ DONE ]")

  # Iterate over each molecule group and compare it to Si
  # ====================================================

  for mol_group in mol_groups:

    if mol_group != 'Si':

      print ("{:<140}".format('\nTreating the values for %s QDs ...' % mol_group), end="")

      # Get the values
      # ==============

      gaps = []

      for mol in yml_sorted[mol_group]:
        if yml_sorted[mol_group][mol].get('Energy gaps (Ha)'):

          # Get the gaps values

          hl_gap = energy_unit_conversion(yml_sorted[mol_group][mol]['Energy gaps (Ha)']['HOMO-LUMO'],"ha","ev")
          opt_gap = energy_unit_conversion(yml_sorted[mol_group][mol]['Energy gaps (Ha)']['Optical'],"ha","ev")

          # Use the TAG key to find the corresponding Si QD and fetch its number of Si atoms

          for si_mol in yml_sorted['Si']:
            if yml_sorted['Si'][si_mol]['ID']['TAG'] == yml_sorted[mol_group][mol]['ID']['TAG']:
              nb_si_atoms = yml_sorted['Si'][si_mol]['ID']['Nb Si atoms']

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

      # Plot the specific group value

      ax.plot([mol[0] for mol in gaps],[mol[1] for mol in gaps],marker='.',linestyle='-',label='H-L gaps - %s (DFT)' % mol_group)
      ax.plot([mol[0] for mol in gaps],[mol[2] for mol in gaps],marker='^',markersize=4,linestyle='-',label='Optical gaps - %s (TD-DFT)' % mol_group)

      # Add the legend and titles

      ax.set_title('Energy gaps: Si vs %s' % mol_group)
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

      plt.savefig(os.path.join(gaps_dir,'gaps_Si_vs_%s.png' % mol_group),dpi=200)
      plt.cla()

      print('%12s' % "[ DONE ]")
   
  # =================================================================== #
  # =================================================================== #
  #                      PLOTTING ORBITAL ENERGIES                      #
  # =================================================================== #
  # =================================================================== #

  section_title = "2. Orbital energies"

  print("")
  print(''.center(len(section_title)+10, '*'))
  print(section_title.center(len(section_title)+10))
  print(''.center(len(section_title)+10, '*'))

  # Create the directory that will contain the graphs
  # =================================================

  orb_dir = os.path.join(out_dir,'orbitals')
  os.makedirs(orb_dir, exist_ok = True)

  # Iterate over each group
  # =======================

  for mol_group in mol_groups:

    print ("{:<140}".format('\nTreating the values for %s QDs ...' % mol_group), end="")

    # Initialize the list of tuples that will contain all the values

    occ_orbitals = []
    virt_orbitals = []

    # Smallest molecule (reference)
    # =============================

    # Identify the smallest molecule of the group that will act as reference

    min_atoms =  min([yml_sorted[mol_group][mol]['ID']['Nb atoms'] for mol in yml_sorted[mol_group] if yml_sorted[mol_group][mol].get('QCHEM KS Orbitals')])
    ref_mol = list(filter(lambda mol: yml_sorted[mol_group][mol]['ID']['Nb atoms'] == min_atoms, yml_sorted[mol_group]))[0]

    # Convert the orbital values from string to list and take the top 5 levels

    occ_values = sorted(map(float,yml_sorted[mol_group][ref_mol]['QCHEM KS Orbitals']['Occupied'].split(';')),reverse=True)[0:5]
    virt_values = sorted(map(float,yml_sorted[mol_group][ref_mol]['QCHEM KS Orbitals']['Virtual'].split(';')))[0:5]

    # Filter the extreme values (values higher than the double of the mean)

    occ_values = [val for val in occ_values if abs(val) < (2 * abs(sum(occ_values)/len(occ_values)))]
    virt_values = [val for val in virt_values if abs(val) < (2 * abs(sum(virt_values)/len(virt_values)))]

    # Shift the values by centering them around the HOMO-LUMO gap

    hl_gap = min(virt_values) - max(occ_values)
    shift = max(occ_values) + (hl_gap/2)
    occ_values = [val - shift for val in sorted(occ_values)]
    virt_values = [val - shift for val in virt_values]

    # Store the values in the list as tuples (x = 1 since it is the smallest molecule of the group)

    for value in occ_values:
      occ_orbitals.append((1,value))
    for value in virt_values:
      virt_orbitals.append((1,value))

    # Define the upper and lower limits that will be used to filter the orbital energies of the bigger molecules of the group

    upper_limit = max(virt_values)
    lower_limit = min(occ_values)

    # Rest of the group
    # =================

    # Sort the molecules by size

    sorted_mol = sorted(yml_sorted[mol_group].keys(), key=lambda mol: yml_sorted[mol_group][mol]['ID']['Nb atoms'])

    # Iterate over each molecule

    for mol in yml_sorted[mol_group]:
      if mol != ref_mol and yml_sorted[mol_group][mol].get('QCHEM KS Orbitals'):

        # Determine the x value

        x_value = sorted_mol.index(mol) + 1

        # Convert the orbital values from string to list

        occ_values = sorted(map(float,yml_sorted[mol_group][mol]['QCHEM KS Orbitals']['Occupied'].split(';')),reverse=True)
        virt_values = sorted(map(float,yml_sorted[mol_group][mol]['QCHEM KS Orbitals']['Virtual'].split(';')))

        # Shift the values by centering them around the HOMO-LUMO gap

        hl_gap = min(virt_values) - max(occ_values)
        shift = max(occ_values) + (hl_gap/2)
        occ_values = [val - shift for val in sorted(occ_values)]
        virt_values = [val - shift for val in virt_values]

        # Filter the values based on the group limits

        occ_values = [val for val in occ_values if val >= lower_limit]
        virt_values = [val for val in virt_values if val <= upper_limit]

        # Store the values in the list as tuples

        for value in occ_values:
          occ_orbitals.append((x_value,value))
        for value in virt_values:
          virt_orbitals.append((x_value,value))
    
    # Plot the values
    # ===============

    plt.style.use('seaborn-colorblind')

    fig, ax = plt.subplots()

    # Define the X labels and ticks

    xlabels = [yml_sorted[mol_group][mol]['ID']['Name'] for mol in sorted_mol]
    xticks = list(range(1,len(xlabels)+1))

    # Plot the values

    ax.scatter([orb[0] for orb in occ_orbitals],[orb[1] for orb in occ_orbitals],label='Occupied',color='red',marker='_',s=900)
    ax.scatter([orb[0] for orb in virt_orbitals],[orb[1] for orb in virt_orbitals],label='Virtual',color='blue',marker='_',s=900)

    # Add the legend and titles

    ax.set_title('Orbital levels for the %s QDs' % mol_group)
    ax.set_ylabel('Energy (Ha)')

    # Set other parameters

    ax.tick_params(top=False, right=False)
    plt.xticks(ticks=xticks,labels=xlabels)
    plt.tight_layout()

    # Save the figure and clear the axes to prepare for the next plot

    plt.savefig(os.path.join(orb_dir,'orb_%s.png' % mol_group),dpi=200)
    plt.cla()

    print('%12s' % "[ DONE ]")

# =================================================================== #
# =================================================================== #
#                          CALL MAIN FUNCTION                         #
# =================================================================== #
# =================================================================== #

# If this script is executed through the command line, call the main function (see https://realpython.com/python-main-function/ for details)

if __name__ == "__main__":
    main()   