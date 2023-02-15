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
import csv
import os
import re
import shutil
from inspect import getsourcefile

import matplotlib.mathtext
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from matplotlib.cm import ScalarMappable
from matplotlib.ticker import AutoMinorLocator
from scipy import constants

import results_common

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
#                        FUNCTIONS DEFINITIONS                        #
# =================================================================== #
# =================================================================== #

def diff_percent(value_1:float, value_2:float):
    """ 
    Computes the difference percentage of value_2 compared to value_1

    Parameters
    ----------
    value_1 : float
        Value of reference
    value_2 : float
        Value for which the difference needs to be computed

    Returns
    -------
    diff_percent : float
        The difference percentage of value_2 compared to value_1
    """

    diff_per = ((value_2 - value_1) / value_1) * 100

    return diff_per

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
      config_file = results_common.check_abspath(config_file,"Command line argument -cf / --config","file")
    else:
      # If no value has been provided through the command line, take the results_config.yml file in the same directory as this script 
      config_file = os.path.join(code_dir, "results_config.yml")

    print ("{:<40} {:<99}".format('\nLoading the configuration file',config_file + " ..."), end="")
    with open(config_file, 'r') as f_config:
      config = yaml.load(f_config, Loader=yaml.FullLoader)
    print('%12s' % "[ DONE ]")

    res_dpi = 400

    # ========================================================= #
    # Check and load the YAML input file                        #
    # ========================================================= #

    inp_yml = results_common.check_abspath(inp_yml,"Command line argument -i / --inp_yml","file")

    print ("{:<40} {:<99}".format('\nLoading the YAML input file',inp_yml + " ..."), end="")
    with open(inp_yml, 'r') as f_yml:
      yml_content = yaml.load(f_yml, Loader=yaml.FullLoader)
    print('%12s' % "[ DONE ]")

    # ========================================================= #
    # Check other arguments                                     #
    # ========================================================= #

    out_dir = results_common.check_abspath(out_dir,"Command line argument -o / --out_dir","directory")
    print ("{:<40} {:<100}".format('\nOutput directory:',out_dir))

    section_count = 0

  # ========================================================= #
  # Exception handling for the preparation step               #
  # ========================================================= #

  except results_common.ResultsError as error:
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

  # Put the Si group at the front for easier handling

  mol_groups.insert(0, mol_groups.pop(mol_groups.index('Si')))

  # Regroup the data by group (constitutive atoms) (#!only include spheroids for now)

  yml_sorted = {}

  for mol_group in mol_groups:
    yml_sorted[mol_group] = {mol:value for (mol,value) in yml_content.items() if yml_content[mol]['ID']['Group'] == mol_group and ((yml_content[mol]['ID']['TAG'].startswith("S") and not yml_content[mol]['ID']['TAG'] == 'S00') or yml_content[mol]['ID']['TAG'].startswith("B"))}
    # yml_sorted[mol_group] = {mol:value for (mol,value) in yml_content.items() if yml_content[mol]['ID']['Group'] == mol_group and yml_content[mol]['ID']['TAG'].startswith("S") and not yml_content[mol]['ID']['TAG'] == 'S00'}

  # mol_groups.pop(mol_groups.index("Bonus"))

  print('%12s' % "[ DONE ]")

  # Create a subdirectory for each group

  for mol_group in mol_groups:
    mol_group_dir = os.path.join(out_dir,mol_group)
    os.makedirs(mol_group_dir, exist_ok=True)

    # Create a subsubdirectory for each molecule in that group

    for mol in yml_sorted[mol_group]:
      mol_dir = os.path.join(mol_group_dir,mol)
      os.makedirs(mol_dir, exist_ok=True)

  # =================================================================== #
  # =================================================================== #
  #                          OPTIMIZED GEOMETRIES                       #
  # =================================================================== #
  # =================================================================== #

  section_count += 1
  section_title = str(section_count) + ". Optimized geometries"

  print("")
  print(''.center(len(section_title)+10, '*'))
  print(section_title.center(len(section_title)+10))
  print(''.center(len(section_title)+10, '*'))

  # Iterate over each molecule group
  # ================================

  for mol_group in mol_groups:

    print ("{:<140}".format('\nTreating the geometries for %s clusters ...' % mol_group), end="")

    # Iterate over each molecule
    # ==========================

    for mol in yml_sorted[mol_group]:

      geom_file = os.path.join(out_dir,mol_group,mol,"geometry.xyz")

      with open(geom_file, 'w+', encoding='utf-8') as f:
        f.write("%s\n\n" % len(yml_sorted[mol_group][mol]['Structure']['Coordinates']))
        for line in yml_sorted[mol_group][mol]['Structure']['Coordinates']:
          f.write("%s\n" % line)
      
    print('%12s' % "[ DONE ]")

  # =================================================================== #
  # =================================================================== #
  #                     SIZE AND ATOMS RELATIONSHIP                     #
  # =================================================================== #
  # =================================================================== #

  section_count += 1
  section_title = str(section_count) + ". Sizes"

  print("")
  print(''.center(len(section_title)+10, '*'))
  print(section_title.center(len(section_title)+10))
  print(''.center(len(section_title)+10, '*'))

  # Get the values
  # ==============

  print ("{:<140}".format('\nTreating the values for Si clusters ...'), end="")

  si_sizes = []

  for mol in yml_sorted['Si']:
    if yml_sorted['Si'][mol].get('Structure'):

      si_sizes.append({
        "Molécule": yml_sorted['Si'][mol]['ID']['Name'],
        "Nombre d'atomes de Si" : yml_sorted['Si'][mol]['Structure']['Nb Si atoms'],
        "Diamètre (nm)": yml_sorted['Si'][mol]['Structure']['Size (nm)']
        })

  si_sizes.sort(key=lambda mol: mol['Diamètre (nm)']) # Sort the values by size

  # Store the values in a CSV file
  # ==============================

  csv_header = list(si_sizes[0].keys())

  with open(os.path.join(out_dir,'Si','sizes.csv'), 'w', newline='', encoding='utf-8') as csvfile:

    csv_writer = csv.DictWriter(csvfile, fieldnames=csv_header, delimiter=';', quoting=csv.QUOTE_MINIMAL)
    csv_writer.writeheader()

    for mol in si_sizes:
      csv_writer.writerow(mol) 

  print('%12s' % "[ DONE ]")

  # =================================================================== #
  # =================================================================== #
  #                        PLOTTING ENERGY GAPS                         #
  # =================================================================== #
  # =================================================================== #

  section_count += 1
  section_title = str(section_count) + ". Energy gaps"

  print("")
  print(''.center(len(section_title)+10, '*'))
  print(section_title.center(len(section_title)+10))
  print(''.center(len(section_title)+10, '*'))

  # ========================================================= #
  # Values for Si group                                       #
  # ========================================================= #

  # Get the values for Si
  # =====================

  print ("{:<140}".format('\nTreating the values for Si clusters ...'), end="")

  si_gaps = []

  for mol in yml_sorted['Si']:
    if yml_sorted['Si'][mol].get('Energy gaps (Ha)'):

      sym_label = str(yml_sorted['Si'][mol]['Structure'].get('Symmetry')).partition(" (")[0]
      sym_latex = "$" + sym_label[0] + "_{" + sym_label[1:] + "}$"

      si_gaps.append({
        "Molécule": yml_sorted['Si'][mol]['ID']['LateX Name'],
        "Diamètre (nm)": yml_sorted['Si'][mol]['Structure']['Size (nm)'],
        "Symétrie" : sym_latex,
        "Gap H-L (eV)": results_common.energy_unit_conversion(yml_sorted['Si'][mol]['Energy gaps (Ha)']['HOMO-LUMO'],"ha","ev"),
        "Gap optique (eV)": results_common.energy_unit_conversion(yml_sorted['Si'][mol]['Energy gaps (Ha)']['Optical'],"ha","ev"),
        "Gap S-T (eV)": results_common.energy_unit_conversion(yml_sorted['Si'][mol]['Energy gaps (Ha)']['Singlet-Triplet'],"ha","ev")
        })

  si_gaps.sort(key=lambda mol: mol['Diamètre (nm)']) # Sort the values by size

  # Store the values in a CSV file
  # ==============================

  csv_header = list(si_gaps[0].keys())

  with open(os.path.join(out_dir,'Si','gaps.csv'), 'w', newline='', encoding='utf-8') as csvfile:

    csv_writer = csv.DictWriter(csvfile, fieldnames=csv_header, delimiter=';', quoting=csv.QUOTE_MINIMAL)
    csv_writer.writeheader()

    for mol in si_gaps:
      csv_writer.writerow(mol) 

  # Create the LaTeX table
  # ======================

  gap_data = [{key:data[key] for key in data if key != 'Gap S-T (eV)'} for data in si_gaps]

  df = pd.DataFrame(gap_data)

  file_content = df.to_latex(
      index=False,
      float_format="%.2f",
      column_format=r"*{5}{>{\centering\arraybackslash}m{0.17\textwidth}}",
      escape=False,
      na_rep="-")

  file_content = file_content.splitlines()

  file_content[2] = r"\multirow{2}{*}{Molécule} & \multirow{2}{*}{Diamètre (nm)} & \multirow{2}{*}{Symétrie} & \multicolumn{2}{c}{Gap (eV)}\\"
  file_content.insert(3, r" &  &  & H-L & optique\\")

  with open(os.path.join(out_dir,'Si',"gaps.tex"), 'w+', encoding='utf-8') as f:
    f.write("\n".join(file_content))

  # Plot the optical and HOMO-LUMO gaps graphs
  # ==========================================

  plt.style.use('seaborn-colorblind')

  fig, ax = plt.subplots()

  # Plot the values

  ax.plot([mol['Diamètre (nm)'] for mol in si_gaps],[mol['Gap H-L (eV)'] for mol in si_gaps],marker='.',linestyle='--',label='Gap H-L (DFT/B3LYP/def2-SVP)')
  ax.plot([mol['Diamètre (nm)'] for mol in si_gaps],[mol['Gap optique (eV)'] for mol in si_gaps],marker='^',markersize=4,linestyle='--',label='Gap optique (TD-DFT/B3LYP/def2-SVP)')

  # Add the legend and titles

  #ax.set_title('Energy gaps for Si QDs')
  ax.set_xlabel("Diamètre (nm)")
  ax.set_ylabel("Énergie (eV)")
  ax.legend()

  # Set other parameters

  ax.tick_params(top=False, right=False)
  ax.xaxis.set_minor_locator(AutoMinorLocator(2))
  ax.yaxis.set_minor_locator(AutoMinorLocator(2))

  plt.tight_layout()
  plt.grid(True,which='both',linestyle='--')

  # Save the file and close the figure

  plt.savefig(os.path.join(out_dir,'Si','gaps.png'),dpi=res_dpi)
  plt.close()

  # Plot the singlet-triplet gaps graphs
  # ====================================

  plt.style.use('seaborn-colorblind')

  fig, ax = plt.subplots()

  # Plot the values

  ax.plot([mol['Diamètre (nm)'] for mol in si_gaps],[mol['Gap S-T (eV)'] for mol in si_gaps],marker='.',linestyle='--')

  # Add the legend and titles

  #ax.set_title('Singlet-Triplet gaps for Si QDs')
  ax.set_xlabel("Diamètre (nm)")
  ax.set_ylabel(r'$\Delta E_{ST}~(eV)$')

  # Set other parameters

  ax.tick_params(top=False, right=False)
  ax.xaxis.set_minor_locator(AutoMinorLocator(2))
  ax.yaxis.set_minor_locator(AutoMinorLocator(2))

  plt.tight_layout()
  plt.grid(True,which='both',linestyle='--')

  params = {'mathtext.default': 'regular' }          
  plt.rcParams.update(params)

  # Save the file and close the figure

  plt.savefig(os.path.join(out_dir,'Si','st_gaps.png'),dpi=res_dpi)
  plt.close()

  print('%12s' % "[ DONE ]")

  # ========================================================= #
  # Values for other groups compared to Si                    #
  # ========================================================= #

  # Iterate over each molecule group
  # ================================

  for mol_group in mol_groups:

    if mol_group != 'Si':

      print ("{:<140}".format('\nTreating the values for %s clusters ...' % mol_group), end="")

      # Get the values
      # ==============

      gaps = []

      for mol in yml_sorted[mol_group]:
        if yml_sorted[mol_group][mol].get('Energy gaps (Ha)'):

          sym_label = str(yml_sorted[mol_group][mol]['Structure'].get('Symmetry')).partition(" (")[0]
          sym_latex = "$" + sym_label[0] + "_{" + sym_label[1:] + "}$"

          gaps.append({
              "Molécule": yml_sorted[mol_group][mol]['ID']['LateX Name'],
              "Diamètre (nm)": yml_sorted[mol_group][mol]['Structure']['Size (nm)'],
              "Symétrie" : sym_latex,
              "Gap H-L (eV)": results_common.energy_unit_conversion(yml_sorted[mol_group][mol]['Energy gaps (Ha)']['HOMO-LUMO'],"ha","ev"),
              "Gap optique (eV)": results_common.energy_unit_conversion(yml_sorted[mol_group][mol]['Energy gaps (Ha)']['Optical'],"ha","ev"),
              "Gap S-T (eV)": results_common.energy_unit_conversion(yml_sorted[mol_group][mol]['Energy gaps (Ha)']['Singlet-Triplet'],"ha","ev")
            })

      gaps.sort(key=lambda mol: mol['Diamètre (nm)']) # Sort the values by size

      # Store the values in a CSV file
      # ==============================

      csv_header = list(gaps[0].keys())

      with open(os.path.join(out_dir,mol_group,'gaps.csv'), 'w', newline='', encoding='utf-8') as csvfile:

        csv_writer = csv.DictWriter(csvfile, fieldnames=csv_header, delimiter=';', quoting=csv.QUOTE_MINIMAL)
        csv_writer.writeheader()

        for mol in gaps:
          csv_writer.writerow(mol) 

      # Create the LaTeX table
      # ======================

      gap_data = [{key:data[key] for key in data if key != 'Gap S-T (eV)'} for data in gaps]

      df = pd.DataFrame(gap_data)

      file_content = df.to_latex(
          index=False,
          float_format="%.2f",
          column_format=r"*{5}{>{\centering\arraybackslash}m{0.17\textwidth}}",
          escape=False,
          na_rep="-")

      file_content = file_content.splitlines()

      file_content[2] = r"\multirow{2}{*}{Molécule} & \multirow{2}{*}{Diamètre (nm)} & \multirow{2}{*}{Symétrie} & \multicolumn{2}{c}{Gap (eV)}\\"
      file_content.insert(3, r" &  &  & H-L & optique\\")

      with open(os.path.join(out_dir,mol_group,"gaps.tex"), 'w+', encoding='utf-8') as f:
        f.write("\n".join(file_content))

      # Plot the optical and HOMO-LUMO gaps graphs
      # ==========================================

      plt.style.use('seaborn-colorblind')

      fig, ax = plt.subplots()

      # Plot the Si values

      ax.plot([mol['Diamètre (nm)'] for mol in si_gaps],[mol['Gap H-L (eV)'] for mol in si_gaps],marker='.',linestyle='--',label='Gap H-L - Si (DFT)')
      ax.plot([mol['Diamètre (nm)'] for mol in si_gaps],[mol['Gap optique (eV)'] for mol in si_gaps],marker='^',markersize=4,linestyle='--',label='Gap optique - Si (TD-DFT)')

      # Plot the specific group value

      ax.plot([mol['Diamètre (nm)'] for mol in gaps],[mol['Gap H-L (eV)'] for mol in gaps],marker='.',linestyle='-',label='Gap H-L - %s (DFT)' % mol_group)
      ax.plot([mol['Diamètre (nm)'] for mol in gaps],[mol['Gap optique (eV)'] for mol in gaps],marker='^',markersize=4,linestyle='-',label='Gap optique - %s (TD-DFT)' % mol_group)

      # Add the legend and titles

      #ax.set_title('Energy gaps for Si QDs vs %s QDs' % mol_group)
      ax.set_xlabel("Diamètre (nm)")
      ax.set_ylabel("Énergie (eV)")
      ax.legend()

      # Set other parameters

      ax.tick_params(top=False, right=False)
      ax.xaxis.set_minor_locator(AutoMinorLocator(2))
      ax.yaxis.set_minor_locator(AutoMinorLocator(2))

      plt.tight_layout()
      plt.grid(True,which='both',linestyle='--')

      # Save the file and close the figure

      plt.savefig(os.path.join(out_dir,mol_group,'gaps.png'),dpi=res_dpi)
      plt.close()

      # Plot the singlet-triplet gaps graphs
      # ====================================

      plt.style.use('seaborn-colorblind')

      fig, ax = plt.subplots()

      # Plot the Si values

      ax.plot([mol['Diamètre (nm)'] for mol in si_gaps],[mol['Gap S-T (eV)'] for mol in si_gaps],marker='.',linestyle='--',label='Gap S-T - Si (TD-DFT)')

      # Plot the specific group value

      ax.plot([mol['Diamètre (nm)'] for mol in gaps],[mol['Gap S-T (eV)'] for mol in gaps],marker='.',linestyle='-',label='Gap S-T - %s (TD-DFT)' % mol_group)

      # Add the legend and titles

      #ax.set_title('Singlet-Triplet gaps for Si QDs vs %s QDs' % mol_group)
      ax.set_xlabel("Diamètre (nm)")
      ax.set_ylabel(r'$\Delta E_{ST}~(eV)$')
      ax.legend()

      # Set other parameters

      ax.tick_params(top=False, right=False)
      ax.xaxis.set_minor_locator(AutoMinorLocator(2))
      ax.yaxis.set_minor_locator(AutoMinorLocator(2))

      plt.tight_layout()
      plt.grid(True,which='both',linestyle='--')

      params = {'mathtext.default': 'regular' }          
      plt.rcParams.update(params)

      # Save the file and close the figure

      plt.savefig(os.path.join(out_dir,mol_group,'st_gaps.png'),dpi=res_dpi)
      plt.close()

      print('%12s' % "[ DONE ]")

  # =================================================================== #
  # =================================================================== #
  #                      PLOTTING ORBITAL ENERGIES                      #
  # =================================================================== #
  # =================================================================== #

  section_count += 1
  section_title = str(section_count) + ". Orbital energies"

  print("")
  print(''.center(len(section_title)+10, '*'))
  print(section_title.center(len(section_title)+10))
  print(''.center(len(section_title)+10, '*'))

  # Iterate over each group
  # =======================

  for mol_group in mol_groups:

    print ("{:<140}".format('\nTreating the values for %s clusters ...' % mol_group), end="")

    # Initialize the list of tuples that will contain all the values

    occ_orbitals = []
    virt_orbitals = []

    # Define the orbitals cutoff (based on the smallest molecule of the reference group)
    # ==========================

    if mol_group == "Si":

      # Identify the smallest molecule of the reference group

      min_atoms =  min([yml_sorted[mol_group][mol]['Structure']['Nb atoms'] for mol in yml_sorted[mol_group] if yml_sorted[mol_group][mol].get('QCHEM KS Orbitals')])
      ref_mol = list(filter(lambda mol: yml_sorted[mol_group][mol]['Structure']['Nb atoms'] == min_atoms, yml_sorted[mol_group]))[0]

      # Convert the orbital values from string to list and take the top 10 levels

      occ_values = sorted(map(float,yml_sorted[mol_group][ref_mol]['QCHEM KS Orbitals']['Occupied'].split(';')),reverse=True)[0:10]
      virt_values = sorted(map(float,yml_sorted[mol_group][ref_mol]['QCHEM KS Orbitals']['Virtual'].split(';')))[0:10]

      # Filter the extreme values (values higher than the triple of the mean)

      occ_values = [val for val in occ_values if abs(val) < (3 * abs(sum(occ_values)/len(occ_values)))]
      virt_values = [val for val in virt_values if abs(val) < (3 * abs(sum(virt_values)/len(virt_values)))]

      # Shift the values by centering them around the HOMO-LUMO gap

      hl_gap = min(virt_values) - max(occ_values)
      shift = max(occ_values) + (hl_gap/2)
      occ_values = [val - shift for val in sorted(occ_values)]
      virt_values = [val - shift for val in virt_values]

      # Define the upper and lower limits that will be used to filter the orbital energies of the bigger molecules of the group

      upper_limit = max(virt_values)
      lower_limit = min(occ_values)

    # Get the values
    # ==============

    # Filter the molecules for which the data is available

    mols = [mol for mol in yml_sorted[mol_group] if yml_sorted[mol_group][mol].get('QCHEM KS Orbitals')]

    # Sort the molecules by size

    sorted_mol = sorted(mols, key=lambda mol: yml_sorted[mol_group][mol]['Structure']['Nb atoms'])

    # Iterate over each molecule

    for mol in sorted_mol:

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

    if mol_group == "Si":
      xlabels = [yml_sorted[mol_group][mol]['ID']['LateX Name'] for mol in sorted_mol]
      xticks = list(range(1,len(xlabels)+1))
    else:
      # Build a list with the Si version and this group version of each molecule side by side, identified by their tag
      xlabels = []
      for group_mol in sorted_mol:
        tag = group_mol.partition("-")[0]
        for si_mol in si_mols:
          if si_mol.partition("-")[0] == tag:
            xlabels.append(yml_sorted['Si'][si_mol]['ID']['LateX Name'])
            xlabels.append(yml_sorted[mol_group][group_mol]['ID']['LateX Name'])
            shift = si_mols.index(si_mol) - sorted_mol.index(group_mol)
      xticks = list(range(1,len(xlabels)+1))

    # Plot the values

    if mol_group == 'Si':

      ax.scatter([orb[0] for orb in occ_orbitals],[results_common.energy_unit_conversion(orb[1],"ha","ev") for orb in occ_orbitals],label='Occupées',color='red',marker='_',s=900)
      ax.scatter([orb[0] for orb in virt_orbitals],[results_common.energy_unit_conversion(orb[1],"ha","ev") for orb in virt_orbitals],label='Virtuelles',color='blue',marker='_',s=900)

    else:

      for x in range(0, len(sorted_mol)):
        ax.scatter([orb[0] - shift + x for orb in si_occ_orbitals if orb[0] == x + 1 + shift],[results_common.energy_unit_conversion(orb[1],"ha","ev") for orb in si_occ_orbitals if orb[0] == x + 1 + shift],color='orange',marker='_',s=450)
        ax.scatter([orb[0] + x + 1 for orb in occ_orbitals if orb[0] == x + 1],[results_common.energy_unit_conversion(orb[1],"ha","ev") for orb in occ_orbitals if orb[0] == x + 1],color='red',marker='_',s=450)
        ax.scatter([orb[0] - shift + x for orb in si_virt_orbitals if orb[0] == x + 1 + shift],[results_common.energy_unit_conversion(orb[1],"ha","ev") for orb in si_virt_orbitals if orb[0] == x + 1 + shift],color='green',marker='_',s=450)
        ax.scatter([orb[0] + x + 1 for orb in virt_orbitals if orb[0] == x + 1],[results_common.energy_unit_conversion(orb[1],"ha","ev") for orb in virt_orbitals if orb[0] == x + 1],color='blue',marker='_',s=450)

    # Add the legend and titles

    ax.set_ylabel("Énergie (eV)")

    # Set other parameters

    ax.tick_params(top=False, right=False, bottom=False)
    plt.xticks(ticks=xticks, labels=xlabels, rotation=-45)
    plt.axhline(y=0, color='grey', linestyle='--')
    plt.tight_layout()

    params = {'mathtext.default': 'regular' }          
    plt.rcParams.update(params)

    # Save the file and close the figure

    plt.savefig(os.path.join(out_dir,mol_group,'orb.png'),dpi=res_dpi)
    plt.close()

    # Save the results for the Si group, for comparison with the other groups

    if mol_group == "Si":
      si_occ_orbitals = occ_orbitals
      si_virt_orbitals = virt_orbitals
      si_mols = sorted_mol

    print('%12s' % "[ DONE ]")

  # =================================================================== #
  # =================================================================== #
  #                   PLOTTING IONIZATION POTENTIALS                    #
  # =================================================================== #
  # =================================================================== #

  section_count += 1
  section_title = str(section_count) + ". Ionization potentials"

  print("")
  print(''.center(len(section_title)+10, '*'))
  print(section_title.center(len(section_title)+10))
  print(''.center(len(section_title)+10, '*'))

  # Iterate over each molecule group
  # ================================

  for mol_group in mol_groups:

      print ("{:<140}".format('\nTreating the values for %s clusters ...' % mol_group), end="")

      # Get the values
      # ==============

      ips_all = []
      ips_vert = []
      ips_adiab = []

      for mol in yml_sorted[mol_group]:
        if yml_sorted[mol_group][mol].get('IPs (Ha)'):

          # Get the size values

          size = yml_sorted[mol_group][mol]['Structure']['Size (nm)']

          # Get the IPs values and add them to their corresponding lists (separate lists because all sizes are not necessarily represented for each IP)

          ip_vert = results_common.energy_unit_conversion(yml_sorted[mol_group][mol]['IPs (Ha)']['Vertical'],"ha","ev")
          ips_vert.append((size,ip_vert))

          ip_adiab = yml_sorted[mol_group][mol]['IPs (Ha)']['Adiabatic']
          if ip_adiab != 'N/A':
            ip_adiab = results_common.energy_unit_conversion(ip_adiab,"ha","ev")
            ips_adiab.append((size,ip_adiab))
          else:
            ip_adiab = None
  
          # Compute the differences between the two IPs

          if ip_adiab:
            diff_va = ip_vert - ip_adiab
          else:
            diff_va = None

          # Store the data for this molecule

          ips_all.append({
            "Molécule": yml_sorted[mol_group][mol]['ID']['LateX Name'],
            "Diamètre (nm)": size,
            "P.I. vertical (eV)": ip_vert,
            "P.I. adiabatique (eV)": ip_adiab,
            r"$\Delta_{VA}$ (eV)" : diff_va
            })

      # If they were no values for this group, skip it

      if ips_all == []:
        continue

      # Sort the values by size

      ips_vert.sort(key=lambda tup: tup[0])
      ips_adiab.sort(key=lambda tup: tup[0])

      ips_all.sort(key=lambda mol: mol['Diamètre (nm)'])

      # Store the values in a CSV file
      # ==============================

      csv_header = list(ips_all[0].keys())

      with open(os.path.join(out_dir,mol_group,'ips.csv'), 'w', newline='', encoding='utf-8') as csvfile:

        csv_writer = csv.DictWriter(csvfile, fieldnames=csv_header, delimiter=';', quoting=csv.QUOTE_MINIMAL)
        csv_writer.writeheader()

        for mol in ips_all:
          csv_writer.writerow(mol) 

      # Create the LaTeX table
      # ======================

      df = pd.DataFrame(ips_all)

      file_content = df.to_latex(
          index=False,
          float_format="%.2f",
          column_format=r"*{" + str(len(ips_all[0])) + r"}{>{\centering\arraybackslash}m{0.17\textwidth}}",
          escape=False,
          na_rep="-")

      file_content = file_content.splitlines()

      file_content[2] = r"\multirow{2}{*}{Molécule} & \multirow{2}{*}{Diamètre (nm)} & \multicolumn{2}{c}{Potentiels d'ionisation (eV)} & \multirow{2}{*}{$\Delta_{VA}$ (eV)} \\"
      file_content.insert(3, r"  &  & vertical & adiabatique &  \\")

      with open(os.path.join(out_dir,mol_group,"ips.tex"), 'w+', encoding='utf-8') as f:
        f.write("\n".join(file_content))

      # Plot the graphs
      # ===============

      plt.style.use('seaborn-colorblind')

      fig, ax = plt.subplots()

      # Plot the values

      ax.plot([mol[0] for mol in ips_vert],[mol[1] for mol in ips_vert],marker='.',linestyle='-',label='Vertical')
      ax.plot([mol[0] for mol in ips_adiab],[mol[1] for mol in ips_adiab],marker='.',linestyle='-',label='Adiabatique')

      # Add the legend and titles

      #ax.set_title('Ionization potentials for %s clusters' % mol_group)
      ax.set_xlabel("Diamètre (nm)")
      ax.set_ylabel("Énergie (eV)")
      ax.legend()

      # Set other parameters

      ax.tick_params(top=False, right=False)
      ax.xaxis.set_minor_locator(AutoMinorLocator(2))
      ax.yaxis.set_minor_locator(AutoMinorLocator(2))

      plt.tight_layout()
      plt.grid(True,which='both',linestyle='--')

      # Save the file and close the figure

      plt.savefig(os.path.join(out_dir,mol_group,'ips.png'),dpi=res_dpi)
      plt.close()

      print('%12s' % "[ DONE ]")

  # =================================================================== #
  # =================================================================== #
  #                  PLOTTING CHARACTERIZATION RESULTS                  #
  # =================================================================== #
  # =================================================================== #

  section_count += 1
  section_title = str(section_count) + ". Characterization results"

  print("")
  print(''.center(len(section_title)+10, '*'))
  print(section_title.center(len(section_title)+10))
  print(''.center(len(section_title)+10, '*'))

  # Iterate over each molecule group
  # ================================

  for mol_group in mol_groups:

      print ("{:<140}".format('\nTreating the values for %s clusters ...' % mol_group), end="")

      # Get the values
      # ==============

      gs_momdip_ind = []
      gs_soc_ind = []

      nr_momdip = []
      soc_results = []
      mixing_results = []
      rel_momdip = []
      avg_dir_momdip = []

      for mol in yml_sorted[mol_group]:
        if yml_sorted[mol_group][mol].get('Transition dipole moments (au)') and yml_sorted[mol_group][mol].get('Relativistic states list'):

          # States list
          # ===========

          # Non relativistic states list

          zero_states = []

          for state in yml_sorted[mol_group][mol]["Zero states list"]:
            if state != "S0":
              zero_states.append({
                "Label" : state,
                "Energie (cm" + r'$^{-1}$' + ")" : results_common.energy_unit_conversion(yml_sorted[mol_group][mol]["Zero states list"][state]['Energy (Ha)'],"ha","cm-1")
              })

          # Relativistic states list

          rel_states = []

          for state in yml_sorted[mol_group][mol]["Relativistic states list"]:
            if state != "GS":
              rel_states.append({
                "Label" : state,
                "Energie (cm" + r'$^{-1}$' + ")" : results_common.energy_unit_conversion(yml_sorted[mol_group][mol]["Relativistic states list"][state]['Energy (Ha)'],"ha","cm-1")
              })

          # Mixing singlet-triplet percentages

          mix_perc = []

          for state in yml_sorted[mol_group][mol]["Relativistic states list"]:
            if state != "GS":
              mix_perc.append({
                "Label" : state,
                r'$\%_{GS}$' : yml_sorted[mol_group][mol]["Relativistic states list"][state]['GS percentage'] * 100,
                r'$\%_{S}$' : yml_sorted[mol_group][mol]["Relativistic states list"][state]['Singlet percentage'] * 100,
                r'$\%_{T}$' : yml_sorted[mol_group][mol]["Relativistic states list"][state]['Triplet percentage'] * 100
              })

          # Store the states list in CSV files
          # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          # Non relativistic states list

          csv_header = list(zero_states[0].keys())

          with open(os.path.join(out_dir,mol_group,mol,'zero_states.csv'), 'w', newline='', encoding='utf-8') as csvfile:

            csv_writer = csv.DictWriter(csvfile, fieldnames=csv_header, delimiter=';', quoting=csv.QUOTE_MINIMAL)
            csv_writer.writeheader()

            for state in zero_states:
              csv_writer.writerow(state) 

          # Relativistic states list

          csv_header = list(rel_states[0].keys())

          with open(os.path.join(out_dir,mol_group,mol,'rel_states.csv'), 'w', newline='', encoding='utf-8') as csvfile:

            csv_writer = csv.DictWriter(csvfile, fieldnames=csv_header, delimiter=';', quoting=csv.QUOTE_MINIMAL)
            csv_writer.writeheader()

            for state in rel_states:
              csv_writer.writerow(state) 

          # Mixing singlet-triplet percentages

          csv_header = list(mix_perc[0].keys())

          with open(os.path.join(out_dir,mol_group,mol,'mix_perc.csv'), 'w', newline='', encoding='utf-8') as csvfile:

            csv_writer = csv.DictWriter(csvfile, fieldnames=csv_header, delimiter=';', quoting=csv.QUOTE_MINIMAL)
            csv_writer.writeheader()

            for state in mix_perc:
              csv_writer.writerow(state) 

          # Create the states LaTeX tables
          # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          # The zero (non-relativistic) states data will be split into two tables, that will be presented side-by-side in the thesis

          zero_states_1 = zero_states[:len(zero_states)//2] # First half
          zero_states_2 = zero_states[len(zero_states)//2:] # Second half 

          for idx, frame in enumerate([zero_states_1, zero_states_2]):
            df = pd.DataFrame(frame)
            with open(os.path.join(out_dir,mol_group,mol,"zero_states_%s.tex" % (idx+1)), 'w+', encoding='utf-8') as f:
              f.write(df.to_latex(
                index=False,
                formatters=[None, results_common.format_num(2, "f")],
                column_format="cc",
                escape=False,
                na_rep="-"))

          # The relativistic states data will be split into two tables, that will be presented side-by-side in the thesis

          rel_states_1 = rel_states[:len(rel_states)//2] # First half
          rel_states_2 = rel_states[len(rel_states)//2:] # Second half 

          for idx, frame in enumerate([rel_states_1, rel_states_2]):
            df = pd.DataFrame(frame)
            with open(os.path.join(out_dir,mol_group,mol,"rel_states_%s.tex" % (idx+1)), 'w+', encoding='utf-8') as f:
              f.write(df.to_latex(
                index=False,
                formatters=[None, results_common.format_num(2, "f")],
                column_format="cc",
                escape=False,
                na_rep="-"))

          # Ground sate values
          # ==================

         # Non relativistic transition dipole moments with the ground state

          gs_momdip_line = {
            "Molécule": yml_sorted[mol_group][mol]['ID']['LateX Name'],
            "Diamètre (nm)": yml_sorted[mol_group][mol]['Structure']['Size (nm)']
          }

          for state in yml_sorted[mol_group][mol]["Zero states list"]:
            if state != "S0" and state.startswith("S"):
              gs_momdip = 0
              for momdip_key in yml_sorted[mol_group][mol]["Zero states list"][state]['GS transition dipole moment (au)']:
                gs_momdip += yml_sorted[mol_group][mol]["Zero states list"][state]['GS transition dipole moment (au)'][momdip_key] ** 2
              gs_momdip = np.sqrt(gs_momdip)

              gs_momdip_line.update({state : gs_momdip})

          gs_momdip_ind.append(gs_momdip_line)

          # Total SOC with the ground state

          gs_soc_line = {
            "Molécule": yml_sorted[mol_group][mol]['ID']['LateX Name'],
            "Diamètre (nm)": yml_sorted[mol_group][mol]['Structure']['Size (nm)']
          }

          for state in yml_sorted[mol_group][mol]["Zero states list"]:
            if state.startswith("T"):
              soc = 0
              for ms_value in yml_sorted[mol_group][mol]["Zero states list"][state]['GS SOC values (au)']:
                soc += yml_sorted[mol_group][mol]["Zero states list"][state]['GS SOC values (au)'][ms_value] ** 2
              soc = results_common.energy_unit_conversion(np.sqrt(soc),"ha","cm-1")

              gs_soc_line.update({state : soc})

          gs_soc_ind.append(gs_soc_line)

          # Average dipole moments
          # ======================

          # Define the labels for each type of dipole moment

          gs_momdip_label = r'$\bar{\mu}_{S0 \to S}$'
          ss_momdip_label = r'$\bar{\mu}_{S \to S}$'
          tt_momdip_label = r'$\bar{\mu}_{T \to T}$'
          p_momdip_label = r'$\bar{\mu}_p$'
          all_momdip_label = r'$\bar{\mu}_g$'

          # Store the values for this molecule

          nr_momdip.append({
              "Molécule": yml_sorted[mol_group][mol]['ID']['LateX Name'],
              "TAG": yml_sorted[mol_group][mol]['ID']['TAG'],
              "Diamètre (nm)": yml_sorted[mol_group][mol]['Structure']['Size (nm)'],
              gs_momdip_label : yml_sorted[mol_group][mol]["Average NR dipole moments (au)"]["GS-S"],
              ss_momdip_label : yml_sorted[mol_group][mol]["Average NR dipole moments (au)"]["S-S"],
              tt_momdip_label : yml_sorted[mol_group][mol]["Average NR dipole moments (au)"]["T-T"],
              p_momdip_label : yml_sorted[mol_group][mol]["Average NR dipole moments (au)"]["Permanent"],
              all_momdip_label : yml_sorted[mol_group][mol]["Average NR dipole moments (au)"]["Total"]
              })

          # Average SOC values
          # ==================

          # Average soc values with the ground state

          gs_soc_label = r"$\overline{soc}_{S0 \to T'}$"
          gs_soc = [float(yml_sorted[mol_group][mol]['SOC values (au)'][key]) for key in yml_sorted[mol_group][mol]['SOC values (au)'] if key.startswith('S0')]

          mean_gs_soc = np.mean(gs_soc)

          # Average soc values between singlets and triplet substates

          st_soc_label = r"$\overline{soc}_{S \to T'}$"
          st_soc = [float(yml_sorted[mol_group][mol]['SOC values (au)'][key]) for key in yml_sorted[mol_group][mol]['SOC values (au)'] if 'S' in key and not key.startswith('S0')]

          mean_st_soc = np.mean(st_soc)

          # Average soc values between triplet substates

          tt_soc_label = r"$\overline{soc}_{T' \to T'}$"
          tt_soc = [float(yml_sorted[mol_group][mol]['SOC values (au)'][key]) for key in yml_sorted[mol_group][mol]['SOC values (au)'] if not 'S' in key]

          mean_tt_soc = np.mean(tt_soc)

          # Global average of SOC values

          all_soc_label = r'$\overline{soc}_g$'
          all_soc = gs_soc + st_soc + tt_soc

          mean_all_soc = np.mean(all_soc)

          # Store the values for this molecule

          soc_results.append({
              "Molécule": yml_sorted[mol_group][mol]['ID']['LateX Name'],
              "TAG": yml_sorted[mol_group][mol]['ID']['TAG'],
              "Diamètre (nm)": yml_sorted[mol_group][mol]['Structure']['Size (nm)'],
              gs_soc_label : results_common.energy_unit_conversion(mean_gs_soc,"ha","cm-1"),
              st_soc_label : results_common.energy_unit_conversion(mean_st_soc,"ha","cm-1"),
              tt_soc_label : results_common.energy_unit_conversion(mean_tt_soc,"ha","cm-1"),
              all_soc_label : results_common.energy_unit_conversion(mean_all_soc,"ha","cm-1")
              })

          # Average mixing percentages
          # ==========================

          # Average singlet percentages in pseudo-triplet relativistic states

          sing_percents_label = r'$\bar{\%}_{S}$ (\%)'
          sing_percents = [yml_sorted[mol_group][mol]["Relativistic states list"][state]["Singlet percentage"] for state in yml_sorted[mol_group][mol]["Relativistic states list"] if yml_sorted[mol_group][mol]["Relativistic states list"][state]["Triplet percentage"] >= 0.5]
          mean_percent = np.mean(sing_percents)

          # Average excitation energy of non relativistic states

          energies_label = r'$\bar{E}_{NR}$ (cm$^{-1}$)'
          energies = [yml_sorted[mol_group][mol]["Zero states list"][state]["Energy (Ha)"] for state in yml_sorted[mol_group][mol]["Zero states list"]]

          mean_energies = np.mean(energies)

          # Compute the average ratio between S-T SOC and triplet energies

          ratio_label = r'$\overline{C/\Delta E}$'
          ratio_list = []

          for soc_label in yml_sorted[mol_group][mol]['SOC values (au)']:
            if 'S' in soc_label and not soc_label.startswith('S0'):

              state_1_label = soc_label.partition("_")[0]
              state_2_label = soc_label.partition("_")[2]
              soc_value = yml_sorted[mol_group][mol]['SOC values (au)'][soc_label]

              energy_1 = yml_sorted[mol_group][mol]["Zero states list"][state_1_label.partition("(")[0]]['Energy (Ha)']
              energy_2 = yml_sorted[mol_group][mol]["Zero states list"][state_2_label.partition("(")[0]]['Energy (Ha)']
              delta_e = abs(energy_1 - energy_2)

              ratio_list.append(soc_value / delta_e)

          mean_ratio = np.mean(ratio_list)

          # Store the values for this molecule

          mixing_results.append({
              "Molécule": yml_sorted[mol_group][mol]['ID']['LateX Name'],
              "Diamètre (nm)": yml_sorted[mol_group][mol]['Structure']['Size (nm)'],
              sing_percents_label : mean_percent*100,
              energies_label : results_common.energy_unit_conversion(mean_energies,"ha","cm-1"),
              st_soc_label + r' (cm$^{-1}$)' : results_common.energy_unit_conversion(mean_st_soc,"ha","cm-1"),
              ratio_label : mean_ratio
              })

          # Average relativistic dipole moments
          # ===================================

          gs_ps_momdip_label = r'$\bar{\mu}_{GS \to pS}$'
          gs_pt_momdip_label = r'$\bar{\mu}_{GS \to pT}$'
          ps_ps_momdip_label = r'$\bar{\mu}_{pS \to pS}$'
          pt_pt_momdip_label = r'$\bar{\mu}_{pT \to pT}$'
          ps_pt_momdip_label = r'$\bar{\mu}_{pS \to pT}$'
          rp_momdip_label = r'$\bar{\mu}_{p}$'

          # Store the values for this molecule

          rel_momdip.append({
              "Molécule": yml_sorted[mol_group][mol]['ID']['LateX Name'],
              "TAG": yml_sorted[mol_group][mol]['ID']['TAG'],
              "Diamètre (nm)": yml_sorted[mol_group][mol]['Structure']['Size (nm)'],
              gs_ps_momdip_label : yml_sorted[mol_group][mol]["Average relativistic dipole moments (au)"]["GS-S"],
              ps_ps_momdip_label : yml_sorted[mol_group][mol]["Average relativistic dipole moments (au)"]["S-S"],
              pt_pt_momdip_label : yml_sorted[mol_group][mol]["Average relativistic dipole moments (au)"]["T-T"],
              #rp_momdip_label : yml_sorted[mol_group][mol]["Average relativistic dipole moments (au)"]["Permanent"],
              gs_pt_momdip_label : yml_sorted[mol_group][mol]["Average relativistic dipole moments (au)"]["GS-T"],
              ps_pt_momdip_label : yml_sorted[mol_group][mol]["Average relativistic dipole moments (au)"]["S-T"]
              })

          # Average relativistic dipole moments in each direction
          # =====================================================

          x_momdip_label = r'$\bar{\mu}_{X}$'
          y_momdip_label = r'$\bar{\mu}_{Y}$'
          z_momdip_label = r'$\bar{\mu}_{Z}$'

          # Store the values for this molecule

          avg_dir_momdip.append({
              "Molécule": yml_sorted[mol_group][mol]['ID']['LateX Name'],
              "Diamètre (nm)": yml_sorted[mol_group][mol]['Structure']['Size (nm)'],
              x_momdip_label : yml_sorted[mol_group][mol]["Average relativistic dipole moments (au)"]["X"],
              y_momdip_label : yml_sorted[mol_group][mol]["Average relativistic dipole moments (au)"]["Y"],
              z_momdip_label : yml_sorted[mol_group][mol]["Average relativistic dipole moments (au)"]["Z"]
              })

          # Plot the average non relativistic vs relativistic dipole moments pie charts
          # ===========================================================================

          plt.style.use('seaborn-colorblind')

          fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(nrows=2, ncols=2, gridspec_kw={'height_ratios': [3, 1]})

          # Define the number of roots (useful for restoring sum of values based on average values)

          nb_roots = len([state for state in yml_sorted[mol_group][mol]["Zero states list"] if state != "S0" and state.startswith("S")])

          # Non relativistic values
          # ~~~~~~~~~~~~~~~~~~~~~~~

          # Define the non relativistic values (multiply the averages by the number of times the values appear)

          gs_s_sum = yml_sorted[mol_group][mol]["Sum NR dipole moments (au)"]["GS-S"]
          s_s_sum = yml_sorted[mol_group][mol]["Sum NR dipole moments (au)"]["S-S"]
          t_t_sum = yml_sorted[mol_group][mol]["Sum NR dipole moments (au)"]["T-T"]

          total_sum = gs_s_sum + s_s_sum + t_t_sum # used to compute the percentages for the label

          y_values = np.array([
            gs_s_sum,
            s_s_sum,
            t_t_sum
          ])

          # Define the labels

          y_labels = [
            r'$\sum \mu_{S0 \to S}$' + " - " + np.format_float_positional((gs_s_sum / total_sum) * 100, precision=1, pad_left=3, unique=False) + " %",
            r'$\sum \mu_{S \to S}$' + "  - " + np.format_float_positional((s_s_sum / total_sum) * 100, precision=1, pad_left=3, unique=False) + " %",
            r'$\sum \mu_{T \to T}$' + "  - " + np.format_float_positional((t_t_sum / total_sum) * 100, precision=1, pad_left=3, unique=False) + " %"
          ]

          # Plot the relativistic values

          pie = ax1.pie(y_values, wedgeprops={'linewidth': 1.0, 'edgecolor': 'white'})

          # Add the legend and titles

          ax3.axis("off")
          ax3.legend(pie[0], y_labels, loc="center")

          # Relativistic values
          # ~~~~~~~~~~~~~~~~~~~

          # Define the relativistic values (multiply the averages by the number of times the values appear)

          gs_ps_sum = yml_sorted[mol_group][mol]["Sum relativistic dipole moments (au)"]["GS-S"]
          gs_pt_sum = yml_sorted[mol_group][mol]["Sum relativistic dipole moments (au)"]["GS-T"]
          ps_ps_sum = yml_sorted[mol_group][mol]["Sum relativistic dipole moments (au)"]["S-S"]
          pt_pt_sum = yml_sorted[mol_group][mol]["Sum relativistic dipole moments (au)"]["T-T"]
          ps_pt_sum = yml_sorted[mol_group][mol]["Sum relativistic dipole moments (au)"]["S-T"]

          total_sum = gs_ps_sum + gs_pt_sum + ps_ps_sum + pt_pt_sum + ps_pt_sum # used to compute the percentages for the label

          y_values = np.array([
            gs_ps_sum,
            ps_ps_sum,
            pt_pt_sum,
            gs_pt_sum,
            ps_pt_sum
          ])

          # Define the labels

          y_labels = [
            r'$\sum \mu_{GS \to pS}$' + " - " + np.format_float_positional((gs_ps_sum / total_sum) * 100, precision=1, pad_left=3, unique=False) + " %",
            r'$\sum \mu_{pS \to pS}$' + " - " + np.format_float_positional((ps_ps_sum / total_sum) * 100, precision=1, pad_left=3, unique=False) + " %",
            r'$\sum \mu_{pT \to pT}$' + " - " + np.format_float_positional((pt_pt_sum / total_sum) * 100, precision=1, pad_left=3, unique=False) + " %",
            r'$\sum \mu_{GS \to pT}$' + " - " + np.format_float_positional((gs_pt_sum / total_sum) * 100, precision=1, pad_left=3, unique=False) + " %",
            r'$\sum \mu_{pS \to pT}$' + " - " + np.format_float_positional((ps_pt_sum / total_sum) * 100, precision=1, pad_left=3, unique=False) + " %"
          ]

          # Plot the values

          explosion_values = [0, 0, 0, 0.2, 0.2]
          pie2 = ax2.pie(y_values, explode = explosion_values, wedgeprops={'linewidth': 1.0, 'edgecolor': 'white'})

          # Add the legend and titles

          ax4.axis("off")
          ax4.legend(pie2[0], y_labels, loc="center")

          # Set other parameters

          plt.tight_layout()

          params = {'mathtext.default': 'regular' }          
          plt.rcParams.update(params)

          # Save the file and close the figure

          plt.savefig(os.path.join(out_dir,mol_group,mol,'momdip_pies.png'),dpi=res_dpi)
          plt.close()  

      # Sort the values by size

      gs_momdip_ind.sort(key=lambda mol: mol['Diamètre (nm)'])
      gs_soc_ind.sort(key=lambda mol: mol['Diamètre (nm)'])         

      nr_momdip.sort(key=lambda mol: mol['Diamètre (nm)'])
      soc_results.sort(key=lambda mol: mol['Diamètre (nm)'])
      mixing_results.sort(key=lambda mol: mol['Diamètre (nm)'])
      rel_momdip.sort(key=lambda mol: mol['Diamètre (nm)'])
      avg_dir_momdip.sort(key=lambda mol: mol['Diamètre (nm)'])

      # Compare the values with Si group
      # ================================

      # Store the values for the Si group

      if mol_group == 'Si':
        si_nr_momdip = nr_momdip
        si_rel_momdip = rel_momdip
        si_soc_results = soc_results

      # Compare the values for the other groups

      else:

        # Difference in transition dipole moments

        diff_nr_momdip = []
        for group_mol in nr_momdip:
          for si_mol in si_nr_momdip:
            if si_mol['TAG'] == group_mol['TAG']:
              diff_nr_momdip.append({
                  "Molécule": group_mol["Molécule"],
                  r"$\Delta$ " + gs_momdip_label : diff_percent(si_mol[gs_momdip_label], group_mol[gs_momdip_label]),
                  r"$\Delta$ " + ss_momdip_label : diff_percent(si_mol[ss_momdip_label], group_mol[ss_momdip_label]),
                  r"$\Delta$ " + tt_momdip_label : diff_percent(si_mol[tt_momdip_label], group_mol[tt_momdip_label]),
                # r"$\Delta$ " + p_momdip_label : diff_percent(si_mol[p_momdip_label], group_mol[p_momdip_label]),
                  r"$\Delta$ " + all_momdip_label : diff_percent(si_mol[all_momdip_label], group_mol[all_momdip_label])
                  })
              break

        # Difference in SOC values

        diff_soc_results = []
        for group_mol in soc_results:
          for si_mol in si_soc_results:
            if si_mol['TAG'] == group_mol['TAG']:
              diff_soc_results.append({
                  "Molécule": group_mol["Molécule"],
                  r"$\Delta$ " + gs_soc_label : diff_percent(si_mol[gs_soc_label], group_mol[gs_soc_label]),
                  r"$\Delta$ " + st_soc_label : diff_percent(si_mol[st_soc_label], group_mol[st_soc_label]),
                  r"$\Delta$ " + tt_soc_label : diff_percent(si_mol[tt_soc_label], group_mol[tt_soc_label]),
                  r"$\Delta$ " + all_soc_label : diff_percent(si_mol[all_soc_label], group_mol[all_soc_label])
                  })
              break

        # Difference in relativistic transition dipole moments

        diff_rel_momdip = []
        for group_mol in rel_momdip:
          for si_mol in si_rel_momdip:
            if si_mol['TAG'] == group_mol['TAG']:
              diff_rel_momdip.append({
                  "Molécule": group_mol["Molécule"],
                  r"$\Delta$ " + gs_ps_momdip_label : diff_percent(si_mol[gs_ps_momdip_label], group_mol[gs_ps_momdip_label]),
                  r"$\Delta$ " + ps_ps_momdip_label : diff_percent(si_mol[ps_ps_momdip_label], group_mol[ps_ps_momdip_label]),
                  r"$\Delta$ " + pt_pt_momdip_label : diff_percent(si_mol[pt_pt_momdip_label], group_mol[pt_pt_momdip_label]),
                # r"$\Delta$ " + rp_momdip_label : diff_percent(si_mol[rp_momdip_label], group_mol[rp_momdip_label]),
                  r"$\Delta$ " + gs_pt_momdip_label : diff_percent(si_mol[gs_pt_momdip_label], group_mol[gs_pt_momdip_label]),
                  r"$\Delta$ " + ps_pt_momdip_label : diff_percent(si_mol[ps_pt_momdip_label], group_mol[ps_pt_momdip_label])
                  })
              break

      # Store the values in CSV files
      # =============================

      # Ground state values for transition dipole moments

      csv_header = list(gs_momdip_ind[0].keys())

      with open(os.path.join(out_dir,mol_group,'gs_momdip_ind.csv'), 'w', newline='', encoding='utf-8') as csvfile:

        csv_writer = csv.DictWriter(csvfile, fieldnames=csv_header, delimiter=';', quoting=csv.QUOTE_MINIMAL)
        csv_writer.writeheader()

        for mol in gs_momdip_ind:
          csv_writer.writerow(mol) 

      # Ground state values for total SOC

      csv_header = list(gs_soc_ind[0].keys())

      with open(os.path.join(out_dir,mol_group,'gs_soc_ind.csv'), 'w', newline='', encoding='utf-8') as csvfile:

        csv_writer = csv.DictWriter(csvfile, fieldnames=csv_header, delimiter=';', quoting=csv.QUOTE_MINIMAL)
        csv_writer.writeheader()

        for mol in gs_soc_ind:
          csv_writer.writerow(mol) 

      # Average transition dipole moments

      csv_header = list(nr_momdip[0].keys())

      with open(os.path.join(out_dir,mol_group,'nr_momdip.csv'), 'w', newline='', encoding='utf-8') as csvfile:

        csv_writer = csv.DictWriter(csvfile, fieldnames=csv_header, delimiter=';', quoting=csv.QUOTE_MINIMAL)
        csv_writer.writeheader()

        for mol in nr_momdip:
          csv_writer.writerow(mol)

      # Average SOC values

      csv_header = list(soc_results[0].keys())

      with open(os.path.join(out_dir,mol_group,'soc_results.csv'), 'w', newline='', encoding='utf-8') as csvfile:

        csv_writer = csv.DictWriter(csvfile, fieldnames=csv_header, delimiter=';', quoting=csv.QUOTE_MINIMAL)
        csv_writer.writeheader()

        for mol in soc_results:
          csv_writer.writerow(mol)

      # Average mixing percentages

      csv_header = list(mixing_results[0].keys())

      with open(os.path.join(out_dir,mol_group,'mixing_results.csv'), 'w', newline='', encoding='utf-8') as csvfile:

        csv_writer = csv.DictWriter(csvfile, fieldnames=csv_header, delimiter=';', quoting=csv.QUOTE_MINIMAL)
        csv_writer.writeheader()

        for mol in mixing_results:
          csv_writer.writerow(mol) 

      # Average relativistic dipole moments

      all_rel_momdip = []
      for idx, line in enumerate(rel_momdip):
        all_rel_momdip.append({**line, **avg_dir_momdip[idx]}) # Merge rel_momdip and avg_dir_momdip

      csv_header = list(all_rel_momdip[0].keys())

      with open(os.path.join(out_dir,mol_group,'all_rel_momdip.csv'), 'w', newline='', encoding='utf-8') as csvfile:

        csv_writer = csv.DictWriter(csvfile, fieldnames=csv_header, delimiter=';', quoting=csv.QUOTE_MINIMAL)
        csv_writer.writeheader()

        for mol in all_rel_momdip:
          csv_writer.writerow(mol) 

      if mol_group not in ["Si","Bonus"]:

        # Difference in transition dipole moments

        csv_header = list(diff_nr_momdip[0].keys())

        with open(os.path.join(out_dir,mol_group,'diff_nr_momdip.csv'), 'w', newline='', encoding='utf-8') as csvfile:

          csv_writer = csv.DictWriter(csvfile, fieldnames=csv_header, delimiter=';', quoting=csv.QUOTE_MINIMAL)
          csv_writer.writeheader()

          for mol in diff_nr_momdip:
            csv_writer.writerow(mol)

        # Difference in SOC values

        csv_header = list(diff_soc_results[0].keys())

        with open(os.path.join(out_dir,mol_group,'diff_soc_results.csv'), 'w', newline='', encoding='utf-8') as csvfile:

          csv_writer = csv.DictWriter(csvfile, fieldnames=csv_header, delimiter=';', quoting=csv.QUOTE_MINIMAL)
          csv_writer.writeheader()

          for mol in diff_soc_results:
            csv_writer.writerow(mol)

        # Difference in relativistic transition dipole moments

        csv_header = list(diff_rel_momdip[0].keys())

        with open(os.path.join(out_dir,mol_group,'diff_rel_momdip.csv'), 'w', newline='', encoding='utf-8') as csvfile:

          csv_writer = csv.DictWriter(csvfile, fieldnames=csv_header, delimiter=';', quoting=csv.QUOTE_MINIMAL)
          csv_writer.writeheader()

          for mol in diff_rel_momdip:
            csv_writer.writerow(mol)

      # Create the LaTeX tables
      # =======================

      # Ground state values for transition dipole moments

      gs_momdip_ind = [{key:data[key] for key in data if key != 'Diamètre (nm)'} for data in gs_momdip_ind]

      df = pd.DataFrame(gs_momdip_ind)

      file_content = df.to_latex(
          index=False,
          float_format="%.3f",
          column_format= r">{\centering\arraybackslash}m{0.15\textwidth}" + r"*{" + str(len(gs_momdip_ind[0])-1) + r"}{>{\centering\arraybackslash}m{0.10\textwidth}}",
          escape=False,
          na_rep="-")

      file_content = file_content.splitlines()

      file_content[2] = r"\multirow{2}{*}{Molécule} & \multicolumn{" + str(len(gs_momdip_ind[0])-1) + r"}{c}{Moment dipolaire de transition $|\vec{\mu}_{ij}|$ depuis S0 vers ... (u.a.)}\\"
      file_content.insert(3, "".join(" & %s" % state for state in list(gs_momdip_ind[0].keys()) if state != "Molécule") + r"\\")
      file_content.insert(3, r"\cline{2-" + str(len(gs_momdip_ind[0])) + "}")

      with open(os.path.join(out_dir,mol_group,"gs_momdip_ind.tex"), 'w+', encoding='utf-8') as f:
        f.write("\n".join(file_content))

      # Ground state values for total SOC

      gs_soc_ind = [{key:data[key] for key in data if key != 'Diamètre (nm)'} for data in gs_soc_ind]

      df = pd.DataFrame(gs_soc_ind)

      file_content = df.to_latex(
          index=False,
          float_format="%.3f",
          column_format=r">{\centering\arraybackslash}m{0.15\textwidth}" + r"*{" + str(len(gs_soc_ind[0])-1) + r"}{>{\centering\arraybackslash}m{0.10\textwidth}}",
          escape=False,
          na_rep="-")

      file_content = file_content.splitlines()

      file_content[2] = r"\multirow{2}{*}{Molécule} & \multicolumn{" + str(len(gs_soc_ind[0])-1) + r"}{c}{Couplage spin-orbite total de S0 avec ... (cm$^{-1}$)}\\"
      file_content.insert(3, "".join(" & %s" % state for state in list(gs_soc_ind[0].keys()) if state != "Molécule") + r"\\")
      file_content.insert(3, r"\cline{2-" + str(len(gs_soc_ind[0])) + "}")

      with open(os.path.join(out_dir,mol_group,"gs_soc_ind.tex"), 'w+', encoding='utf-8') as f:
        f.write("\n".join(file_content))

      # Average transition dipole moments

      df = pd.DataFrame([{key:data[key] for key in data if key not in ['Diamètre (nm)','TAG']} for data in nr_momdip])

      file_content = df.to_latex(
          index=False,
          float_format="%.3f",
          column_format=r">{\centering\arraybackslash}m{0.15\textwidth}" + r"*{5}{>{\centering\arraybackslash}m{0.12\textwidth}}",
          escape=False,
          na_rep="-")

      file_content = file_content.splitlines()

      file_content.insert(2, r"\multirow{2}{*}{Molécule} & \multicolumn{5}{c}{Moyennes des moments dipolaires (u.a.)}\\")
      file_content[3] = (file_content[3].lstrip()).lstrip("Molécule")
      file_content.insert(3, r"\cline{2-6}")

      with open(os.path.join(out_dir,mol_group,"nr_momdip.tex"), 'w+', encoding='utf-8') as f:
        f.write("\n".join(file_content))

      # Average SOC values

      df = pd.DataFrame([{key:data[key] for key in data if key not in ['Diamètre (nm)','TAG']} for data in soc_results])

      file_content = df.to_latex(
          index=False,
          float_format="%.3f",
          column_format=r">{\centering\arraybackslash}m{0.15\textwidth}" + r"*{5}{>{\centering\arraybackslash}m{0.15\textwidth}}",
          escape=False,
          na_rep="-")

      file_content = file_content.splitlines()

      file_content.insert(2, r"\multirow{2}{*}{Molécule} & \multicolumn{4}{c}{Moyennes des couplages spin-orbite (cm$^{-1}$)}\\")
      file_content[3] = (file_content[3].lstrip()).lstrip("Molécule")
      file_content.insert(3, r"\cline{2-5}")

      with open(os.path.join(out_dir,mol_group,"soc_results.tex"), 'w+', encoding='utf-8') as f:
        f.write("\n".join(file_content))

      # Average mixing percentages

      mixing_results = [{key:data[key] for key in data if key != 'Diamètre (nm)'} for data in mixing_results]

      df = pd.DataFrame(mixing_results)
      with open(os.path.join(out_dir,mol_group,"mixing_results.tex"), 'w+', encoding='utf-8') as f:
        f.write(df.to_latex(
          index=False,
          formatters=[None, results_common.format_num(3, "f"), results_common.format_num(2, "f"), results_common.format_num(3, "f"), results_common.format_num(2, "e")],
          column_format="lcccc",
          escape=False,
          na_rep="-"))

      # Average relativistic dipole moments

      df = pd.DataFrame([{key:data[key] for key in data if key not in ['Diamètre (nm)','TAG']} for data in rel_momdip])

      file_content = df.to_latex(
          index=False,
          float_format="%.2e",
          column_format=r">{\centering\arraybackslash}m{0.15\textwidth}" + r"*{5}{>{\centering\arraybackslash}m{0.12\textwidth}}",
          escape=False,
          na_rep="-")

      file_content = file_content.splitlines()

      file_content.insert(2, r"\multirow{2}{*}{Molécule} & \multicolumn{5}{c}{Moyennes des moments dipolaires (u.a.)}\\")
      file_content[3] = (file_content[3].lstrip()).lstrip("Molécule")
      file_content.insert(3, r"\cline{2-6}")

      with open(os.path.join(out_dir,mol_group,"rel_momdip.tex"), 'w+', encoding='utf-8') as f:
        f.write("\n".join(file_content))

      if mol_group not in ["Si","Bonus"]:

        # Difference in transition dipole moments

        df = pd.DataFrame(diff_nr_momdip)

        file_content = df.to_latex(
            index=False,
            float_format="%.1f",
            column_format=r">{\centering\arraybackslash}m{0.2\textwidth}" + r"*{4}{>{\centering\arraybackslash}m{0.15\textwidth}}",
            escape=False,
            na_rep="-")

        file_content = file_content.splitlines()

        file_content.insert(2, r"\multirow{2}{*}{Molécule} & \multicolumn{4}{c}{Pourcentages d'écart des moyennes (\%)}\\")
        file_content[3] = (file_content[3].lstrip()).lstrip("Molécule")
        file_content.insert(3, r"\cline{2-5}")

        with open(os.path.join(out_dir,mol_group,"diff_nr_momdip.tex"), 'w+', encoding='utf-8') as f:
          f.write("\n".join(file_content))        

        # Difference in SOC values

        df = pd.DataFrame(diff_soc_results)

        file_content = df.to_latex(
            index=False,
            float_format="%.1f",
            column_format=r">{\centering\arraybackslash}m{0.14\textwidth}" + r"*{5}{>{\centering\arraybackslash}m{0.16\textwidth}}",
            escape=False,
            na_rep="-")

        file_content = file_content.splitlines()

        file_content.insert(2, r"\multirow{2}{*}{Molécule} & \multicolumn{4}{c}{Pourcentages d'écart des moyennes (\%)}\\")
        file_content[3] = (file_content[3].lstrip()).lstrip("Molécule")
        file_content.insert(3, r"\cline{2-5}")

        with open(os.path.join(out_dir,mol_group,"diff_soc_results.tex"), 'w+', encoding='utf-8') as f:
          f.write("\n".join(file_content))

        # Difference in relativistic dipole moments

        df = pd.DataFrame(diff_rel_momdip)

        file_content = df.to_latex(
            index=False,
            float_format="%.1f",
            column_format=r">{\centering\arraybackslash}m{0.14\textwidth}" + r"*{5}{>{\centering\arraybackslash}m{0.16\textwidth}}",
            escape=False,
            na_rep="-")

        file_content = file_content.splitlines()

        file_content.insert(2, r"\multirow{2}{*}{Molécule} & \multicolumn{5}{c}{Pourcentages d'écart des moyennes (\%)}\\")
        file_content[3] = (file_content[3].lstrip()).lstrip("Molécule")
        file_content.insert(3, r"\cline{2-6}")

        with open(os.path.join(out_dir,mol_group,"diff_rel_momdip.tex"), 'w+', encoding='utf-8') as f:
          f.write("\n".join(file_content))

      # Plot the average transition dipole moments graph
      # ================================================

      plt.style.use('seaborn-colorblind')

      fig, ax = plt.subplots()

      # Define the values

      barwidth = 0.1

      br1 = np.arange(len(nr_momdip))
      br2 = [bar + barwidth for bar in br1]
      br3 = [bar + barwidth for bar in br2]
      br4 = [bar + barwidth for bar in br3]
      br5 = [bar + barwidth for bar in br4]

      # Define the X labels and ticks

      xlabels = [mol['Molécule'] for mol in nr_momdip]
      xticks = [r + barwidth for r in range(len(xlabels))]

      # Plot the values

      ax.bar(br1, [mol[gs_momdip_label] for mol in nr_momdip], width = barwidth, edgecolor ='grey', label = gs_momdip_label)
      ax.bar(br2, [mol[ss_momdip_label] for mol in nr_momdip], width = barwidth, edgecolor ='grey', label = ss_momdip_label)
      ax.bar(br3, [mol[tt_momdip_label] for mol in nr_momdip], width = barwidth, edgecolor ='grey', label = tt_momdip_label)
      ax.bar(br4, [mol[p_momdip_label] for mol in nr_momdip], width = barwidth, edgecolor ='grey', label = p_momdip_label)
      ax.bar(br5, [mol[all_momdip_label] for mol in nr_momdip], width = barwidth, edgecolor ='grey', label = all_momdip_label)

      # Add the legend and titles

      ax.set_ylabel('Moyenne des moments dipolaires (u.a.)')
      ax.legend()

      # Set other parameters

      ax.tick_params(top=False, right=False, bottom=False)
      ax.set_axisbelow(True) # axes and grid beneath the plots
      plt.xticks(ticks=xticks, labels=xlabels, rotation=-45)
      plt.tight_layout()
      plt.grid(True, which='both', axis='y', linestyle='--')

      params = {'mathtext.default': 'regular' }          
      plt.rcParams.update(params)

      # Save the file and close the figure

      plt.savefig(os.path.join(out_dir,mol_group,'nr_momdip.png'),dpi=res_dpi)
      plt.close()  

      # Plot the average SOC values graph
      # =================================

      plt.style.use('seaborn-colorblind')

      fig, ax = plt.subplots()

      # Define the values

      barwidth = 0.15

      br1 = np.arange(len(soc_results))
      br2 = [bar + barwidth for bar in br1]
      br3 = [bar + barwidth for bar in br2]
      br4 = [bar + barwidth for bar in br3]

      # Define the X labels and ticks

      xlabels = [mol['Molécule'] for mol in soc_results]
      xticks = [r + barwidth for r in range(len(xlabels))]

      # Plot the values

      ax.bar(br1, [mol[gs_soc_label] for mol in soc_results], width = barwidth, edgecolor ='grey', label = gs_soc_label)
      ax.bar(br2, [mol[st_soc_label] for mol in soc_results], width = barwidth, edgecolor ='grey', label = st_soc_label)
      ax.bar(br3, [mol[tt_soc_label] for mol in soc_results], width = barwidth, edgecolor ='grey', label = tt_soc_label)
      ax.bar(br4, [mol[all_soc_label] for mol in soc_results], width = barwidth, edgecolor ='grey', label = all_soc_label)

      # Add the legend and titles

      ax.set_ylabel('Moyenne des couplages spin-orbite ' + r'(cm$^{-1}$)')
      ax.legend()

      # Set other parameters

      ax.tick_params(top=False, right=False, bottom=False)
      ax.set_axisbelow(True) # axes and grid beneath the plots
      plt.xticks(ticks=xticks, labels=xlabels, rotation=-45)
      plt.tight_layout()
      plt.grid(True, which='both', axis='y', linestyle='--')

      params = {'mathtext.default': 'regular' }          
      plt.rcParams.update(params)

      # Save the file and close the figure

      plt.savefig(os.path.join(out_dir,mol_group,'soc_results.png'),dpi=res_dpi)
      plt.close() 

      """
      # Plot the mixing percentage graphs
      # =================================

      plt.style.use('seaborn-colorblind')

      fig, ax = plt.subplots()

      # Plot the values

      ax.plot([mol['Diamètre (nm)'] for mol in mixing_results],[mol[sing_percents_label] for mol in mixing_results],marker='.',linestyle='-')

      # Add the legend and titles

      #ax.set_title('Ionization potentials for %s clusters' % mol_group)
      ax.set_xlabel("Diamètre (nm)")
      ax.set_ylabel("Pourcentage moyen d'états singulets (%)")

      # Set other parameters

      ax.tick_params(top=False, right=False)
      ax.xaxis.set_minor_locator(AutoMinorLocator(2))
      ax.yaxis.set_minor_locator(AutoMinorLocator(2))

      plt.tight_layout()
      plt.grid(True,which='both',linestyle='--')

      # Save the file and close the figure

      plt.savefig(os.path.join(out_dir,mol_group,'mixing_percent.png'),dpi=res_dpi)
      plt.close()
      """

      # Plot the average relativistic dipole moments graph
      # ==================================================

      plt.style.use('seaborn-colorblind')

      fig, ax = plt.subplots()

      # Define the values

      barwidth = 0.15

      br1 = np.arange(len(rel_momdip))
      br2 = [bar + barwidth for bar in br1]

      # Define the X labels and ticks

      xlabels = [mol['Molécule'] for mol in rel_momdip]
      xticks = [r + barwidth for r in range(len(xlabels))]

      # Plot the values

      ax.bar(br1, [mol[gs_pt_momdip_label] for mol in rel_momdip], width = barwidth, edgecolor ='grey', label = gs_pt_momdip_label)
      ax.bar(br2, [mol[ps_pt_momdip_label] for mol in rel_momdip], width = barwidth, edgecolor ='grey', label = ps_pt_momdip_label)

      # Add the legend and titles

      ax.set_ylabel('Moyenne des moments dipolaires (u.a.)')
      ax.legend()

      # Set other parameters

      # ax.set_yscale('log')

      ax.tick_params(top=False, right=False, bottom=False)
      ax.set_axisbelow(True) # axes and grid beneath the plots
      plt.xticks(ticks=xticks, labels=xlabels, rotation=-45)
      plt.tight_layout()
      plt.grid(True, which='both', axis='y', linestyle='--')

      params = {'mathtext.default': 'regular' }          
      plt.rcParams.update(params)

      # Save the file and close the figure

      plt.savefig(os.path.join(out_dir,mol_group,'rel_momdip.png'),dpi=res_dpi)
      plt.close()  

      # Plot the average relativistic dipole moments in each direction
      # ==============================================================

      plt.style.use('seaborn-colorblind')

      fig, ax = plt.subplots()

      # Define the values

      barwidth = 0.15

      br1 = np.arange(len(avg_dir_momdip))
      br2 = [bar + barwidth for bar in br1]
      br3 = [bar + barwidth for bar in br2]

      # Define the X labels and ticks

      xlabels = [mol['Molécule'] for mol in avg_dir_momdip]
      xticks = [r + barwidth for r in range(len(xlabels))]

      # Plot the values

      ax.bar(br1, [mol[x_momdip_label] for mol in avg_dir_momdip], color='red', width = barwidth, edgecolor ='grey', label = x_momdip_label)
      ax.bar(br2, [mol[y_momdip_label] for mol in avg_dir_momdip], color='green', width = barwidth, edgecolor ='grey', label = y_momdip_label)
      ax.bar(br3, [mol[z_momdip_label] for mol in avg_dir_momdip], color='blue', width = barwidth, edgecolor ='grey', label = z_momdip_label)

      # Add the legend and titles

      ax.set_ylabel('Moyenne des moments dipolaires (u.a.)')
      ax.legend()

      # Set other parameters

      ax.tick_params(top=False, right=False, bottom=False)
      ax.set_axisbelow(True) # axes and grid beneath the plots
      plt.xticks(ticks=xticks, labels=xlabels, rotation=-45)
      plt.tight_layout()
      plt.grid(True, which='both', axis='y', linestyle='--')

      params = {'mathtext.default': 'regular' }          
      plt.rcParams.update(params)

      # Save the file and close the figure

      plt.savefig(os.path.join(out_dir,mol_group,'avg_dir_momdip.png'),dpi=res_dpi)
      plt.close()  

      print('%12s' % "[ DONE ]")

  # =================================================================== #
  # =================================================================== #
  #                   PLOTTING ALDU PARAMETERS SEARCH                   #
  # =================================================================== #
  # =================================================================== #

  section_count += 1
  section_title = str(section_count) + ". Alpha / Duration parameters search"

  print("")
  print(''.center(len(section_title)+10, '*'))
  print(section_title.center(len(section_title)+10))
  print(''.center(len(section_title)+10, '*'))

  # Iterate over each molecule
  # ==========================

  for mol_group in mol_groups:

    transitions = []
    best_pulses = [] # List of dictionaries containing the data of the best pulses
    dir_sens = []    # List of dictionaries containing the data about the directional sensibility of the best pulses

    for mol in yml_sorted[mol_group]:
      if yml_sorted[mol_group][mol].get('Control'):
        if yml_sorted[mol_group][mol]['Control'].get('Alpha/duration parameters search'):

          print ("{:<140}".format('\nTreating the %s molecule ...' % mol))

          # Create a subdirectory for the graphs

          aldu_param_dir = os.path.join(out_dir,mol_group,mol,"aldu_param")
          os.makedirs(aldu_param_dir, exist_ok=True)

          # Initialize a variable for the directional sensibility

          best_perf = 0

          # Iterate over each transition

          for transition in yml_sorted[mol_group][mol]['Control']['Alpha/duration parameters search']:

            print ("{:<133}".format('\n\tTreating the %s transition ...' % transition), end="")

            trans_dir = os.path.join(aldu_param_dir, transition)
            os.makedirs(trans_dir, exist_ok=True)

            # Identify the transition
            # =======================

            init_state = transition.split("_")[1]
            target_state = transition.split("_")[2]
            main_polarisation = yml_sorted[mol_group][mol]['Control']['Alpha/duration parameters search'][transition]['Polarisation']

            transitions.append({
              "Molécule": yml_sorted[mol_group][mol]['ID']['LateX Name'],
              "Diamètre (nm)": yml_sorted[mol_group][mol]['Structure']['Size (nm)'],
              "Polarisation": main_polarisation,
              "Etat initial" : init_state,
              "Etat cible" : target_state
            })

            # Get the values
            # ==============

            values_list = yml_sorted[mol_group][mol]['Control']['Alpha/duration parameters search'][transition]['Values']

            alpha_list = []
            duration_list = []
            fluence_list = []
            efficiency_list = []

            for data in values_list:

              data['Fluence (J/m^2)'] = results_common.energy_unit_conversion(data.pop('Fluence'),'ha','j')/(constants.value('atomic unit of length')**2)

              alpha_list.append(data['Alpha'])
              duration_list.append(data['Duration (ps)'])
              fluence_list.append(data['Fluence (J/m^2)'])
              efficiency_list.append(data['Efficiency'])

              if data['Best'] == 'True':

                best_pulses.append({
                  "Molécule": yml_sorted[mol_group][mol]['ID']['LateX Name'],
                  "Diamètre (nm)": yml_sorted[mol_group][mol]['Structure']['Size (nm)'],
                  "Polarisation" : main_polarisation,
                  r'$\alpha_{0}$': data['Alpha'],
                  "Durée (ps)": data['Duration (ps)'],
                  "Performance": data[main_polarisation + '_Projector'],
                  "Fluence (J/m"+r'$^2$' + ")": data['Fluence (J/m^2)']
                })
            
                max_perf = max(data['X_Projector'],data['Y_Projector'],data['Z_Projector'])
                
                if max_perf > best_perf:

                  best_perf = max_perf
                  old_line = [line for line in dir_sens if line['Molécule'] == yml_sorted[mol_group][mol]['ID']['LateX Name']]

                  if old_line != []:
                    dir_sens.remove(old_line[0])
                  
                  dir_sens.append({
                    "Molécule": yml_sorted[mol_group][mol]['ID']['LateX Name'],
                    "Diamètre (nm)": yml_sorted[mol_group][mol]['Structure']['Size (nm)'],
                    "Performance (X)": data['X_Projector'],
                    "Performance rel. (X)": (data['X_Projector'] / max_perf) * 100,
                    "Performance (Y)": data['Y_Projector'],
                    "Performance rel. (Y)": (data['Y_Projector'] / max_perf) * 100,
                    "Performance (Z)": data['Z_Projector'],
                    "Performance rel. (Z)": (data['Z_Projector'] / max_perf) * 100
                  })

            # Prepare the data (as shown on https://stackoverflow.com/questions/54437559/numpy-meshgrid-from-unordered-x-y-z-data)
            # ~~~~~~~~~~~~~~~~

            alpha_array = np.asarray(alpha_list)
            duration_array = np.asarray(duration_list)
            fluence_array = np.asarray(fluence_list)
            efficiency_array = np.asarray(efficiency_list)

            # Sort coordinates and reshape in grid
            
            len_alpha = len(list(dict.fromkeys(alpha_list))) # Remove duplicates (https://www.w3schools.com/python/python_howto_remove_duplicates.asp)
            len_dur = len(list(dict.fromkeys(duration_list)))
            idx = np.lexsort((duration_array, alpha_array)).reshape(len_alpha, len_dur)

            # Plot the fluence data
            # ~~~~~~~~~~~~~~~~~~~~~

            plt.style.use('ggplot')

            fig, ax = plt.subplots()

            aldu_plot = ax.scatter(alpha_array[idx], duration_array[idx], c=fluence_array[idx], cmap=matplotlib.cm.get_cmap('viridis_r'))

            # Add the titles

            ax.set_xlabel('Facteur de pénalité ' + r'$\alpha_{0}$')
            ax.set_ylabel('Durée du champ (ps)')

            # Set the ticks for alpha (only major ticks corresponding to the data)

            # ax.set_xscale("log")
            # ax.set_xticks(list(dict.fromkeys(alpha_list)))
            # ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
            # ax.minorticks_off()

            ax.set_xscale("log")
            xlabels = ["{:.0f}".format(value) for value in list(dict.fromkeys(alpha_list))]
            ax.set_xticks(ticks = list(dict.fromkeys(alpha_list)), labels = xlabels, rotation=-30)
            ax.minorticks_off()

            # Set other parameters

            ax.set_yticks(list(dict.fromkeys(duration_list)))

            cb = fig.colorbar(ScalarMappable(norm=aldu_plot.norm, cmap=aldu_plot.cmap))
            cb.ax.invert_yaxis()
            cb.set_label("Fluence (J/m" + r'$^{2}$' + ")")

            plt.tight_layout()

            params = {'mathtext.default': 'regular' }          
            plt.rcParams.update(params)

            plt.setp(ax.get_xticklabels(), rotation=-30)

            # Save the file and close the figure

            plot_filename = "fluence.png"

            plt.savefig(os.path.join(trans_dir,plot_filename),dpi=res_dpi,bbox_inches='tight')
            plt.close()

            # Plot the efficiency score
            # ~~~~~~~~~~~~~~~~~~~~~~~~~

            plt.style.use('ggplot')

            fig, ax = plt.subplots()

            aldu_plot = ax.scatter(alpha_array[idx], duration_array[idx], c=efficiency_array[idx], cmap='viridis')

            # Add the titles

            ax.set_xlabel('Facteur de pénalité ' + r'$\alpha_{0}$')
            ax.set_ylabel('Durée du champ (ps)')

            # Set the ticks for alpha (only major ticks corresponding to the data)

            # ax.set_xscale("log")
            # ax.set_xticks(list(dict.fromkeys(alpha_list)))
            # ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
            # ax.minorticks_off()

            ax.set_xscale("log")
            xlabels = ["{:.0f}".format(value) for value in list(dict.fromkeys(alpha_list))]
            ax.set_xticks(ticks = list(dict.fromkeys(alpha_list)), labels = xlabels, rotation=-30)
            ax.minorticks_off()

            # Set other parameters

            ax.set_yticks(list(dict.fromkeys(duration_list)))

            cb = fig.colorbar(ScalarMappable(norm=aldu_plot.norm, cmap=aldu_plot.cmap))
            cb.set_label("Score d'efficacité")

            plt.tight_layout()

            params = {'mathtext.default': 'regular' }          
            plt.rcParams.update(params)

            plt.setp(ax.get_xticklabels(), rotation=-30)

            # Save the file and close the figure

            plot_filename = "efficiency.png"

            plt.savefig(os.path.join(trans_dir,plot_filename),dpi=res_dpi,bbox_inches='tight')
            plt.close()

            # Plot the performance data
            # ~~~~~~~~~~~~~~~~~~~~~~~~~

            # Iterate over each polarisation

            for polarisation in ['X','Y','Z']:

              overlap_list = []

              for data in values_list:
                overlap_list.append(data[polarisation + '_Projector'])
              
              overlap_array = np.asarray(overlap_list)

              plt.style.use('ggplot')

              fig, ax = plt.subplots()

              vmax=min(max(overlap_list),1)

              aldu_plot = ax.scatter(alpha_array[idx], duration_array[idx], c=overlap_array[idx], cmap='viridis', vmax=vmax)

              # Add the titles

              """
              if polarisation == main_polarisation:
                ax.set_title('Parameters search for the %s molecule,\n with a polarization along the %s axis\n' % (yml_sorted[mol_group][mol]['ID']['LateX Name'],polarisation))
              else:
                ax.set_title('Parameters search for the %s molecule,\n with an alternate polarization along the %s axis\n' % (yml_sorted[mol_group][mol]['ID']['LateX Name'],polarisation)) 
              """

              ax.set_xlabel('Facteur de pénalité ' + r'$\alpha_{0}$')
              ax.set_ylabel('Durée du champ (ps)')

              # Set the ticks for alpha

              # ax.set_xscale("log")
              # ax.set_xticks(list(dict.fromkeys(alpha_list)))
              # ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
              # ax.minorticks_off()

              ax.set_xscale("log")
              xlabels = ["{:.0f}".format(value) for value in list(dict.fromkeys(alpha_list))]
              ax.set_xticks(ticks = list(dict.fromkeys(alpha_list)), labels = xlabels, rotation=-30)
              ax.minorticks_off()

              # Set other parameters

              ax.set_yticks(list(dict.fromkeys(duration_list)))

              #cb = fig.colorbar(ScalarMappable(norm=aldu_plot.norm, cmap=aldu_plot.cmap),ticks=np.linspace(vmin, vmax, num=2))
              cb = fig.colorbar(ScalarMappable(norm=aldu_plot.norm, cmap=aldu_plot.cmap))
             
              cb.set_label("Indice de Performance")

              plt.tight_layout()

              params = {'mathtext.default': 'regular' }          
              plt.rcParams.update(params)

              plt.setp(ax.get_xticklabels(), rotation=-30)

              # Save the file and close the figure

              if polarisation == main_polarisation:
                plot_filename = "performance.png"
              else:
                plot_filename = "performance_alt" + polarisation + ".png"

              plt.savefig(os.path.join(trans_dir,plot_filename),dpi=res_dpi,bbox_inches='tight')
              plt.close()

            # Store the values in a CSV file
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            csv_header = list(values_list[0].keys())

            with open(os.path.join(trans_dir,'aldu_param_results.csv'), 'w', newline='', encoding='utf-8') as csvfile:

              csv_writer = csv.DictWriter(csvfile, fieldnames=csv_header, delimiter=';', quoting=csv.QUOTE_MINIMAL)
              csv_writer.writeheader()

              for data in values_list:
                csv_writer.writerow(data) 

            print('%12s' % "[ DONE ]")

    # Create the LaTeX tables
    # =======================

    transitions.sort(key=lambda mol: mol['Diamètre (nm)']) # Sort the values by size
    transitions = [{key:data[key] for key in data if key != 'Diamètre (nm)'} for data in transitions]
    df = pd.DataFrame(transitions)

    with open(os.path.join(out_dir,mol_group,"transitions.tex"), 'w+', encoding='utf-8') as f:
      f.write(df.to_latex(
        index=False,
        column_format="lccc",
        escape=False,
        na_rep="-"))

    if best_pulses != []:

      best_pulses.sort(key=lambda mol: mol['Diamètre (nm)']) # Sort the values by size
      best_pulses = [{key:data[key] for key in data if key != 'Diamètre (nm)'} for data in best_pulses]
      df = pd.DataFrame(best_pulses)

      with open(os.path.join(out_dir,mol_group,"best_pulses.tex"), 'w+', encoding='utf-8') as f:
        f.write(df.to_latex(
          index=False,
          formatters=[None, None, results_common.format_num(0, "f"), results_common.format_num(1, "f"), results_common.format_num(4, "f"), results_common.format_num(2, "f")],
          column_format="lccccc",
          escape=False,
          na_rep="-"))

    if dir_sens != []:
   
      dir_sens.sort(key=lambda mol: mol['Diamètre (nm)']) # Sort the values by size
      dir_sens = [{key:data[key] for key in data if key != 'Diamètre (nm)'} for data in dir_sens]
      df = pd.DataFrame(dir_sens)

      with open(os.path.join(out_dir,mol_group,"dir_sens.tex"), 'w+', encoding='utf-8') as f:
        f.write(df.to_latex(
          index=False,
          formatters=[None, results_common.format_num(2, "e"), results_common.format_num(2, "f"), results_common.format_num(2, "e"), results_common.format_num(2, "f"), results_common.format_num(2, "e"), results_common.format_num(2, "f")],
          column_format="lcccccc",
          escape=False,
          na_rep="-"))

    # Plot the data for directional sensibility of the best pulses
    # ============================================================

    plt.rcParams.update(plt.rcParamsDefault)

    plt.style.use('seaborn-colorblind')

    fig, ax = plt.subplots()

    # Set bars

    barwidth = 0.15

    br1 = np.arange(len(dir_sens))
    br2 = [bar + barwidth for bar in br1]
    br3 = [bar + barwidth for bar in br2]

    # Define the X labels and ticks

    xlabels = [mol['Molécule'] for mol in dir_sens]
    xticks = [r + barwidth for r in range(len(xlabels))]

    # Plot the values

    ax.bar(br1, [mol['Performance (X)'] for mol in dir_sens], color='red', width = barwidth, edgecolor ='grey', label ='X')
    ax.bar(br2, [mol['Performance (Y)'] for mol in dir_sens], color='green', width = barwidth, edgecolor ='grey', label ='Y')
    ax.bar(br3, [mol['Performance (Z)'] for mol in dir_sens], color='blue', width = barwidth, edgecolor ='grey', label ='Z')

    # Add the legend and titles

    ax.set_ylabel('Indice de performance')
    ax.legend()

    # Set other parameters

    ax.tick_params(top=False, right=False, bottom=False)
    ax.set_axisbelow(True) # axes and grid beneath the plots
    plt.xticks(ticks=xticks, labels=xlabels, rotation=-45)
    plt.tight_layout()
    plt.grid(True,which='both', axis='y', linestyle='--')

    params = {'mathtext.default': 'regular' }          
    plt.rcParams.update(params)

    # Save the file and close the figure

    plt.savefig(os.path.join(out_dir,mol_group,'dir_sens.png'),dpi=res_dpi)
    plt.close()    

  # =================================================================== #
  # =================================================================== #
  #                   PLOTTING CONSTRAINTS VARIATION                    #
  # =================================================================== #
  # =================================================================== #

  section_count += 1
  section_title = str(section_count) + ". Constraints variation"

  print("")
  print(''.center(len(section_title)+10, '*'))
  print(section_title.center(len(section_title)+10))
  print(''.center(len(section_title)+10, '*'))

  # Iterate over each molecule
  # ==========================

  for mol_group in mol_groups:

    ratio_const = [] # List of dictionaries containing the max and min values for each constraint

    for mol in yml_sorted[mol_group]:
      if yml_sorted[mol_group][mol].get('Control'):
        if yml_sorted[mol_group][mol]['Control'].get('Constraints variation'):

          print ("{:<140}".format('\nTreating the %s molecule ...' % mol))

          # Create a subdirectory for the graphs

          const_var_dir = os.path.join(out_dir,mol_group,mol,"const_var")
          os.makedirs(const_var_dir, exist_ok=True)

          # Iterate over each transition

          for transition in yml_sorted[mol_group][mol]['Control']['Constraints variation']:

            print ("{:<133}".format('\n\tTreating the %s transition ...' % transition), end="")

            trans_dir = os.path.join(const_var_dir, transition)
            os.makedirs(trans_dir, exist_ok=True)

            # Get the values
            # ==============

            values_list = yml_sorted[mol_group][mol]['Control']['Constraints variation'][transition]['Values']

            if values_list == []:
              print('empty')
              continue

            main_polarisation = yml_sorted[mol_group][mol]['Control']['Constraints variation'][transition]['Polarisation']
            cent_freq = yml_sorted[mol_group][mol]['Control']['Constraints variation'][transition]['Values'][0]['Central frequency (cm-1)']

            fluence_list = []
            window_list = []
            overlap_list = []

            for data in values_list:

              window_list.append(data['Window (cm-1)'])
              fluence_list.append(data['Fluence (J/m^2)'])
              overlap_list.append(data[main_polarisation + '_Projector'])

            ratio_const.append({
              "Molécule": yml_sorted[mol_group][mol]['ID']['LateX Name'],
              "Diamètre (nm)": yml_sorted[mol_group][mol]['Structure']['Size (nm)'],
              "Polarisation" : main_polarisation,
              "Fréquence centrale (cm" + r"$^{-1}$" + ")" : cent_freq,
              "Fenêtre optimale (cm" + r"$^{-1}$" + ")" : max(list(dict.fromkeys(window_list))),
              "Fenêtre idéale (cm" + r"$^{-1}$" + ")" : min(list(dict.fromkeys(window_list))),
              "Fenêtre rapport" : max(list(dict.fromkeys(window_list))) / min(list(dict.fromkeys(window_list))),
              "Fluence optimale (J/m"+r'$^2$' + ")":  max(list(dict.fromkeys(fluence_list))),
              "Fluence idéale (J/m"+r'$^2$' + ")":  min(list(dict.fromkeys(fluence_list))),
              "Fluence rapport" : max(list(dict.fromkeys(fluence_list))) / min(list(dict.fromkeys(fluence_list)))
            })

            # Prepare the data (as shown on https://stackoverflow.com/questions/54437559/numpy-meshgrid-from-unordered-x-y-z-data)
            # ~~~~~~~~~~~~~~~~

            fluence_array = np.asarray(fluence_list)
            window_array = np.asarray(window_list)
            overlap_array = np.asarray(overlap_list)

            # Sort coordinates and reshape in grid
            
            len_flu = len(list(dict.fromkeys(fluence_list))) # Remove duplicates (https://www.w3schools.com/python/python_howto_remove_duplicates.asp)
            len_win = len(list(dict.fromkeys(window_list)))
            idx = np.lexsort((fluence_array, window_array)).reshape(len_flu, len_win)

            # Plot the performance data
            # ~~~~~~~~~~~~~~~~~~~~~~~~~

            plt.style.use('ggplot')

            fig, ax = plt.subplots()

            vmax=min(max(overlap_list),1)

            convar_plot = ax.scatter(fluence_array[idx], window_array[idx], c=overlap_array[idx], cmap='viridis', vmax=vmax)

            # Add the titles

            """
            ax.set_title('Parameters search for the %s molecule,\n with a polarization along the %s axis\n' % (yml_sorted[mol_group][mol]['ID']['LateX Name'],polarisation))
            """

            ax.set_xlabel("Fluence (J/m" + r'$^{2}$' + ")")
            ax.set_ylabel("Largeur spectrale (cm" + r'$^{-1}$' + ")")

            # Set the ticks for alpha

            ax.set_xscale("log")
            xlabels = ["{:.2e}".format(value) for value in list(dict.fromkeys(fluence_list))]
            ax.set_xticks(ticks = list(dict.fromkeys(fluence_list)), labels = xlabels, rotation=-30)
            ax.minorticks_off()
            #ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
            #ax.ticklabel_format(axis='x',style='sci',useMathText=True)

            # Set other parameters

            ax.set_yticks(list(dict.fromkeys(window_list)))

            #cb = fig.colorbar(ScalarMappable(norm=convar_plot.norm, cmap=convar_plot.cmap),ticks=np.linspace(vmin, vmax, num=2))
            cb = fig.colorbar(ScalarMappable(norm=convar_plot.norm, cmap=convar_plot.cmap))
            
            cb.set_label("Indice de Performance")

            plt.tight_layout()

            params = {'mathtext.default': 'regular' }          
            plt.rcParams.update(params)

            #♣plt.setp(ax.get_xticklabels(), rotation=-30)

            # Save the file and close the figure

            plt.savefig(os.path.join(trans_dir,"performance.png"),dpi=res_dpi)
            plt.close()

            # Store the values in a CSV file
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            csv_header = list(values_list[0].keys())

            with open(os.path.join(trans_dir,'convar_comp_results.csv'), 'w', newline='', encoding='utf-8') as csvfile:

              csv_writer = csv.DictWriter(csvfile, fieldnames=csv_header, delimiter=';', quoting=csv.QUOTE_MINIMAL)
              csv_writer.writeheader()

              for data in values_list:
                csv_writer.writerow(data) 

            print('%12s' % "[ DONE ]")

    # Create the LaTeX tables
    # =======================

    if ratio_const != []:

      ratio_const.sort(key=lambda mol: mol['Diamètre (nm)']) # Sort the values by size
      ratio_const = [{key:data[key] for key in data if key != 'Diamètre (nm)'} for data in ratio_const]
      df = pd.DataFrame(ratio_const)

      file_content = df.to_latex(
          index=False,
          formatters=[None, None, results_common.format_num(0, "f"), results_common.format_num(0, "f"), results_common.format_num(0, "f"), results_common.format_num(0, "f"), results_common.format_num(2, "e"), results_common.format_num(2, "e"), results_common.format_num(2, "e")],
          column_format=r">{\centering\arraybackslash}m{0.12\textwidth}" + r">{\centering\arraybackslash}m{0.05\textwidth}" + r">{\centering\arraybackslash}m{0.10\textwidth}" + r"*{6}{>{\centering\arraybackslash}m{0.10\textwidth}}",
          escape=False,
          na_rep="-")

      file_content = file_content.splitlines()

      file_content[2] = r"\multirow{2}{*}{Molécule} & \multirow{2}{*}{Pol.} & $\omega_0$ & \multicolumn{3}{c}{Largeur spectrale (cm$^{-1}$)} & \multicolumn{3}{c}{Fluence (J/m$^2$)}\\"
      file_content.insert(3, r"  &  & (cm$^{-1}$) & optim. & idéale & rapport & optim. & idéale & rapport \\")
      file_content.insert(3, r"\cline{4-9}")

      with open(os.path.join(out_dir,mol_group,"ratio_const.tex"), 'w+', encoding='utf-8') as f:
        f.write("\n".join(file_content))

  # =================================================================== #
  # =================================================================== #
  #                PLOTTING FREQUENCY FILTERS VARIATION                 #
  # =================================================================== #
  # =================================================================== #

  section_count += 1
  section_title = str(section_count) + ". Frequency filters variation"

  print("")
  print(''.center(len(section_title)+10, '*'))
  print(section_title.center(len(section_title)+10))
  print(''.center(len(section_title)+10, '*'))

  # Iterate over each molecule
  # ==========================

  for mol_group in mol_groups:

    filt_freq_results = []

    for mol in yml_sorted[mol_group]:
      if yml_sorted[mol_group][mol].get('Control'):
        if yml_sorted[mol_group][mol]['Control'].get('Frequency filters variation'):

          print ("{:<140}".format('\nTreating the %s molecule ...' % mol))

          # Create a subdirectory for the graphs

          filt_freq_dir = os.path.join(out_dir,mol_group,mol,"filt_freq")
          os.makedirs(filt_freq_dir, exist_ok=True)

          # Iterate over each transition

          for transition in yml_sorted[mol_group][mol]['Control']['Frequency filters variation']:

            print ("{:<133}".format('\n\tTreating the %s transition ...' % transition), end="")

            trans_dir = os.path.join(filt_freq_dir, transition)
            os.makedirs(trans_dir, exist_ok=True)

            # Get the values
            # ==============

            filt_freq_ind = [] # List of dictionaries containing the results for each frequency filter
            values_list = yml_sorted[mol_group][mol]['Control']['Frequency filters variation'][transition]['Values']

            if values_list == []:
              print('empty')
              continue

            main_polarisation = yml_sorted[mol_group][mol]['Control']['Frequency filters variation'][transition]['Polarisation']
            cent_freq = yml_sorted[mol_group][mol]['Control']['Frequency filters variation'][transition]['Values'][0]['Central frequency (cm-1)']

            window_label = "Largeur spectrale (cm" + r"$^{-1}$" + ")"
            perf_label = "Indice de performance"
            fluence_label = "Fluence (J/m"+r'$^2$' + ")"

            opt_wind = max([data['Window (cm-1)'] for data in values_list])
            id_wind = min([data['Window (cm-1)'] for data in values_list])

            for data in values_list:

              # Get the data for the individual tables

              data['Fluence (J/m^2)'] = results_common.energy_unit_conversion(data.pop('Fluence'),'ha','j')/(constants.value('atomic unit of length')**2)

              filt_freq_ind.append({
                window_label : data['Window (cm-1)'],
                perf_label : data[main_polarisation + '_Projector'],
                fluence_label :  data['Fluence (J/m^2)']
              })

              # Get the representative data for the summary tables

              if data['Window (cm-1)'] == opt_wind:
                opt_flu = data['Fluence (J/m^2)']
                opt_perf = data[main_polarisation + '_Projector']

              elif data['Window (cm-1)'] == id_wind:
                id_flu = data['Fluence (J/m^2)']
                id_perf = data[main_polarisation + '_Projector']

            # Sort the individual values by window size

            filt_freq_ind.sort(key=lambda mol: mol[window_label], reverse=True)

            # Plot the data
            # =============

            plt.rcParams.update(plt.rcParamsDefault)

            plt.style.use('seaborn-colorblind')

            fig = plt.figure()
            ax1 = fig.add_subplot(111)

            # Define the X labels and ticks

            xlabels = [str(round(data[window_label])) for data in filt_freq_ind]
            xticks = range(1,len(xlabels)+1)

            # Plot the values

            plot_1 = ax1.bar(xticks, [data[fluence_label] for data in filt_freq_ind], width = 0.2, color="red", edgecolor ='grey', label = "Fluence")
            ax1.set_xticks(ticks=xticks, labels=xlabels)
            ax2 = ax1.twinx() # Two Y axes on the same graph
            plot_2 = ax2.plot(xticks, [data[perf_label] for data in filt_freq_ind], label = "Performance")

            # Add the legend (https://stackoverflow.com/questions/5484922/secondary-axis-with-twinx-how-to-add-to-legend)

            fig.legend(loc="upper center", ncol=2)

            # Add the titles

            ax1.set_xlabel(window_label)
            ax1.set_ylabel(fluence_label)
            ax2.set_ylabel(perf_label)

            # Set other parameters

            #ax1.tick_params(top=False, right=False, bottom=False)
            ax1.set_axisbelow(True) # axes and grid beneath the plots
            plt.tight_layout()
            plt.grid(False)

            plt.tight_layout()

            params = {'mathtext.default': 'regular' }          
            plt.rcParams.update(params)

            # Save the file and close the figure

            plt.savefig(os.path.join(trans_dir,"perf_vs_flu.png"),dpi=res_dpi)
            plt.close()

            # Store the values in a CSV file
            # ==============================

            csv_header = list(values_list[0].keys())

            with open(os.path.join(trans_dir,'filt_freq_comp_results.csv'), 'w', newline='', encoding='utf-8') as csvfile:

              csv_writer = csv.DictWriter(csvfile, fieldnames=csv_header, delimiter=';', quoting=csv.QUOTE_MINIMAL)
              csv_writer.writeheader()

              for data in values_list:
                csv_writer.writerow(data) 

            # Create the LaTeX table
            # ======================

            df = pd.DataFrame(filt_freq_ind)
            with open(os.path.join(trans_dir,"filt_freq_ind.tex"), 'w+', encoding='utf-8') as f:
              f.write(df.to_latex(
                index=False,
                formatters=[results_common.format_num(0, "f"), results_common.format_num(4, "f"), results_common.format_num(2, "f")],
                column_format="ccc",
                escape=False,
                na_rep="-"))

            # Prepare the data for the summary table
            # ======================================

            filt_freq_results.append({
              "Molécule": yml_sorted[mol_group][mol]['ID']['LateX Name'],
              "Diamètre (nm)": yml_sorted[mol_group][mol]['Structure']['Size (nm)'],
              "Polarisation" : main_polarisation,
              "Fréquence centrale (cm^-1)" : cent_freq,
              "Fenêtre optimale (cm^-1)" : opt_wind,
              "Fluence optimale (J/m^2)":  opt_flu,
              "Performance optimale" : opt_perf,
              "Fenêtre idéale (cm^-1)" : id_wind,
              "Fluence idéale (J/m^2)":  id_flu,
              "Performance idéale" : id_perf
            })

            print('%12s' % "[ DONE ]")

    # Create the csv file for the summary table
    # =========================================

    if filt_freq_results != []:

      filt_freq_results.sort(key=lambda mol: mol['Diamètre (nm)'])

      # Average mixing percentages

      csv_header = list(filt_freq_results[0].keys())

      with open(os.path.join(out_dir,mol_group,'filt_freq_results.csv'), 'w', newline='', encoding='utf-8') as csvfile:

        csv_writer = csv.DictWriter(csvfile, fieldnames=csv_header, delimiter=';', quoting=csv.QUOTE_MINIMAL)
        csv_writer.writeheader()

        for mol in filt_freq_results:
          csv_writer.writerow(mol) 
   
    # Create the LaTeX summary table
    # ==============================

    if filt_freq_results != []:

      filt_freq_results = [{key:data[key] for key in data if key != 'Diamètre (nm)'} for data in filt_freq_results]

      df = pd.DataFrame(filt_freq_results)

      file_content = df.to_latex(
          index=False,
          formatters=[None, None, results_common.format_num(0, "f"), results_common.format_num(0, "f"), results_common.format_num(2, "e"), results_common.format_num(3, "f"), results_common.format_num(0, "f"), results_common.format_num(2, "e"), results_common.format_num(3, "f")],
          column_format=r">{\centering\arraybackslash}m{0.12\textwidth}" + r">{\centering\arraybackslash}m{0.05\textwidth}" + r">{\centering\arraybackslash}m{0.10\textwidth}" + r"*{6}{>{\centering\arraybackslash}m{0.10\textwidth}}",
          escape=False,
          na_rep="-")

      file_content = file_content.splitlines()

      file_content[2] = r"\multirow{3}{*}{Molécule} & \multirow{3}{*}{Pol.} & \multirow{2}{*}{$\omega_0$} & \multicolumn{3}{c}{Largeur optimale} & \multicolumn{3}{c}{Largeur idéale} \\"
      file_content.insert(3, r"  &  &             & Valeur      & Fluence   & \multirow{2}{*}{Perf.} & Valeur      & Fluence   & \multirow{2}{*}{Perf.} \\")
      file_content.insert(4, r"  &  & (cm$^{-1}$) & (cm$^{-1}$) & (J/m$^2$) &                        & (cm$^{-1}$) & (J/m$^2$) &                        \\")
      file_content.insert(3, r"\cline{4-9}")

      with open(os.path.join(out_dir,mol_group,"filt_freq_results.tex"), 'w+', encoding='utf-8') as f:
        f.write("\n".join(file_content))

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
