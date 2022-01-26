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
import shutil
from inspect import getsourcefile

import matplotlib.mathtext
import matplotlib.pyplot as plt
import yaml
from matplotlib.ticker import AutoMinorLocator

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

  # Regroup the data by group (constitutive atoms) (#!only include spheroids for now)

  yml_sorted = {}

  for mol_group in mol_groups:
    yml_sorted[mol_group] = {mol:value for (mol,value) in yml_content.items() if yml_content[mol]['ID']['Group'] == mol_group and yml_content[mol]['ID']['TAG'].startswith("S")}

  print('%12s' % "[ DONE ]")

  # =================================================================== #
  # =================================================================== #
  #                     SIZE AND ATOMS RELATIONSHIP                     #
  # =================================================================== #
  # =================================================================== #

  section_title = "0. Sizes"

  print("")
  print(''.center(len(section_title)+10, '*'))
  print(section_title.center(len(section_title)+10))
  print(''.center(len(section_title)+10, '*'))

  # Get the values
  # ==============

  print ("{:<140}".format('\nTreating the values for Si QDs ...'), end="")

  si_sizes = []

  for mol in yml_sorted['Si']:
    if yml_sorted['Si'][mol].get('Structure'):

      si_sizes.append({
        "Molecule": yml_sorted['Si'][mol]['ID']['Name'],
        "Number of Si atoms" : yml_sorted['Si'][mol]['Structure']['Nb Si atoms'],
        "Original diameter (nm)": yml_sorted['Si'][mol]['Structure']['Original size (nm)'],
        "Diameter (nm)": yml_sorted['Si'][mol]['Structure']['Size (nm)']
        })

  si_sizes.sort(key=lambda mol: mol['Diameter (nm)']) # Sort the values by size

  # Store the values in a CSV file
  # ==============================

  csv_header = list(si_sizes[0].keys())

  with open(os.path.join(out_dir,'Si_sizes.csv'), 'w', newline='', encoding='utf-8') as csvfile:

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

  section_title = "1. Energy gaps"

  print("")
  print(''.center(len(section_title)+10, '*'))
  print(section_title.center(len(section_title)+10))
  print(''.center(len(section_title)+10, '*'))

  # ========================================================= #
  # Values for Si group                                       #
  # ========================================================= #

  # Get the values for Si
  # =====================

  print ("{:<140}".format('\nTreating the values for Si QDs ...'), end="")

  si_gaps = []

  for mol in yml_sorted['Si']:
    if yml_sorted['Si'][mol].get('Energy gaps (Ha)'):

      si_gaps.append({
        "Molecule": yml_sorted['Si'][mol]['ID']['Name'],
        "Diameter (nm)": yml_sorted['Si'][mol]['Structure']['Size (nm)'],
        "H-L gap (eV)": results_common.energy_unit_conversion(yml_sorted['Si'][mol]['Energy gaps (Ha)']['HOMO-LUMO'],"ha","ev"),
        "Optical gap (eV)": results_common.energy_unit_conversion(yml_sorted['Si'][mol]['Energy gaps (Ha)']['Optical'],"ha","ev"),
        "S-T gap (eV)": results_common.energy_unit_conversion(yml_sorted['Si'][mol]['Energy gaps (Ha)']['Singlet-Triplet'],"ha","ev")
        })

  si_gaps.sort(key=lambda mol: mol['Diameter (nm)']) # Sort the values by size

  # Store the values in a CSV file
  # ==============================

  csv_header = list(si_gaps[0].keys())

  with open(os.path.join(out_dir,'Si_gaps.csv'), 'w', newline='', encoding='utf-8') as csvfile:

    csv_writer = csv.DictWriter(csvfile, fieldnames=csv_header, delimiter=';', quoting=csv.QUOTE_MINIMAL)
    csv_writer.writeheader()

    for mol in si_gaps:
      csv_writer.writerow(mol) 

  # Plot the optical and HOMO-LUMO gaps graphs
  # ==========================================

  plt.style.use('seaborn-colorblind')

  fig, ax = plt.subplots()

  # Plot the values

  ax.plot([mol['Diameter (nm)'] for mol in si_gaps],[mol['H-L gap (eV)'] for mol in si_gaps],marker='.',linestyle='--',label='H-L gaps (DFT/B3LYP/def2-SVP)')
  ax.plot([mol['Diameter (nm)'] for mol in si_gaps],[mol['Optical gap (eV)'] for mol in si_gaps],marker='^',markersize=4,linestyle='--',label='Optical gaps (TD-DFT/B3LYP/def2-SVP)')

  # Add the legend and titles

  ax.set_title('Energy gaps for Si QDs')
  ax.set_xlabel("Diameter of the QD (nm)")
  ax.set_ylabel('Energy (eV)')
  ax.legend()

  # Set other parameters

  ax.tick_params(top=False, right=False)
  ax.xaxis.set_minor_locator(AutoMinorLocator(2))
  ax.yaxis.set_minor_locator(AutoMinorLocator(2))

  plt.tight_layout()
  plt.grid(True,which='both',linestyle='--')

  # Save the file and close the figure

  plt.savefig(os.path.join(out_dir,'Si_gaps.png'),dpi=200)
  plt.close()

  # Plot the singlet-triplet gaps graphs
  # ====================================

  plt.style.use('seaborn-colorblind')

  fig, ax = plt.subplots()

  # Plot the values

  ax.plot([mol['Diameter (nm)'] for mol in si_gaps],[mol['S-T gap (eV)'] for mol in si_gaps],marker='.',linestyle='--')

  # Add the legend and titles

  ax.set_title('Singlet-Triplet gaps for Si QDs')
  ax.set_xlabel("Diameter of the QD (nm)")
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

  plt.savefig(os.path.join(out_dir,'Si_st_gaps.png'),dpi=200)
  plt.close()

  print('%12s' % "[ DONE ]")

  # ========================================================= #
  # Values for other groups compared to Si                    #
  # ========================================================= #

  # Iterate over each molecule group
  # ================================

  for mol_group in mol_groups:

    if mol_group != 'Si':

      print ("{:<140}".format('\nTreating the values for %s QDs ...' % mol_group), end="")

      # Get the values
      # ==============

      gaps = []

      for mol in yml_sorted[mol_group]:
        if yml_sorted[mol_group][mol].get('Energy gaps (Ha)'):

          gaps.append({
            "Molecule": yml_sorted[mol_group][mol]['ID']['Name'],
            "Diameter (nm)": yml_sorted[mol_group][mol]['Structure']['Size (nm)'],
            "H-L gap (eV)": results_common.energy_unit_conversion(yml_sorted[mol_group][mol]['Energy gaps (Ha)']['HOMO-LUMO'],"ha","ev"),
            "Optical gap (eV)": results_common.energy_unit_conversion(yml_sorted[mol_group][mol]['Energy gaps (Ha)']['Optical'],"ha","ev"),
            "S-T gap (eV)": results_common.energy_unit_conversion(yml_sorted[mol_group][mol]['Energy gaps (Ha)']['Singlet-Triplet'],"ha","ev")
            })

      gaps.sort(key=lambda mol: mol['Diameter (nm)']) # Sort the values by size

      # Store the values in a CSV file
      # ==============================

      csv_header = list(gaps[0].keys())

      with open(os.path.join(out_dir,'%s_gaps.csv' % mol_group), 'w', newline='', encoding='utf-8') as csvfile:

        csv_writer = csv.DictWriter(csvfile, fieldnames=csv_header, delimiter=';', quoting=csv.QUOTE_MINIMAL)
        csv_writer.writeheader()

        for mol in gaps:
          csv_writer.writerow(mol) 

      # Plot the optical and HOMO-LUMO gaps graphs
      # ==========================================

      plt.style.use('seaborn-colorblind')

      fig, ax = plt.subplots()

      # Plot the Si values

      ax.plot([mol['Diameter (nm)'] for mol in si_gaps],[mol['H-L gap (eV)'] for mol in si_gaps],marker='.',linestyle='--',label='H-L gaps - Si (DFT)')
      ax.plot([mol['Diameter (nm)'] for mol in si_gaps],[mol['Optical gap (eV)'] for mol in si_gaps],marker='^',markersize=4,linestyle='--',label='Optical gaps - Si (TD-DFT)')

      # Plot the specific group value

      ax.plot([mol['Diameter (nm)'] for mol in gaps],[mol['H-L gap (eV)'] for mol in gaps],marker='.',linestyle='-',label='H-L gaps - %s (DFT)' % mol_group)
      ax.plot([mol['Diameter (nm)'] for mol in gaps],[mol['Optical gap (eV)'] for mol in gaps],marker='^',markersize=4,linestyle='-',label='Optical gaps - %s (TD-DFT)' % mol_group)

      # Add the legend and titles

      ax.set_title('Energy gaps for Si QDs vs %s QDs' % mol_group)
      ax.set_xlabel("Diameter of the QD (nm)")
      ax.set_ylabel('Energy (eV)')
      ax.legend()

      # Set other parameters

      ax.tick_params(top=False, right=False)
      ax.xaxis.set_minor_locator(AutoMinorLocator(2))
      ax.yaxis.set_minor_locator(AutoMinorLocator(2))

      plt.tight_layout()
      plt.grid(True,which='both',linestyle='--')

      # Save the file and close the figure

      plt.savefig(os.path.join(out_dir,'%s_gaps.png' % mol_group),dpi=200)
      plt.close()

      # Plot the singlet-triplet gaps graphs
      # ====================================

      plt.style.use('seaborn-colorblind')

      fig, ax = plt.subplots()

      # Plot the Si values

      ax.plot([mol['Diameter (nm)'] for mol in si_gaps],[mol['S-T gap (eV)'] for mol in si_gaps],marker='.',linestyle='--',label='S-T gaps - Si (TD-DFT)')

      # Plot the specific group value

      ax.plot([mol['Diameter (nm)'] for mol in gaps],[mol['S-T gap (eV)'] for mol in gaps],marker='.',linestyle='-',label='S-T gaps - %s (TD-DFT)' % mol_group)

      # Add the legend and titles

      ax.set_title('Singlet-Triplet gaps for Si QDs vs %s QDs' % mol_group)
      ax.set_xlabel("Diameter of the QD (nm)")
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

      plt.savefig(os.path.join(out_dir,'%s_st_gaps.png' % mol_group),dpi=200)
      plt.close()

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

    min_atoms =  min([yml_sorted[mol_group][mol]['Structure']['Nb atoms'] for mol in yml_sorted[mol_group] if yml_sorted[mol_group][mol].get('QCHEM KS Orbitals')])
    ref_mol = list(filter(lambda mol: yml_sorted[mol_group][mol]['Structure']['Nb atoms'] == min_atoms, yml_sorted[mol_group]))[0]

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

    # Filter the molecules for which the data is available

    mols = [mol for mol in yml_sorted[mol_group] if yml_sorted[mol_group][mol].get('QCHEM KS Orbitals')]

    # Sort the molecules by size

    sorted_mol = sorted(mols, key=lambda mol: yml_sorted[mol_group][mol]['Structure']['Nb atoms'])

    # Iterate over each molecule

    for mol in sorted_mol:
      if mol != ref_mol:

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

    xlabels = ["$" + yml_sorted[mol_group][mol]['ID']['LateX Name'] + "$" for mol in sorted_mol]
    xticks = list(range(1,len(xlabels)+1))

    # Plot the values

    ax.scatter([orb[0] for orb in occ_orbitals],[orb[1] for orb in occ_orbitals],label='Occupied',color='red',marker='_',s=900)
    ax.scatter([orb[0] for orb in virt_orbitals],[orb[1] for orb in virt_orbitals],label='Virtual',color='blue',marker='_',s=900)

    # Add the legend and titles

    ax.set_title('Orbital energies for the %s QDs' % mol_group)
    ax.set_ylabel('Energy (Ha)')

    # Set other parameters

    ax.tick_params(top=False, right=False, bottom=False)
    plt.xticks(ticks=xticks, labels=xlabels, rotation=-45)
    plt.axhline(y=0, color='grey', linestyle='--')
    plt.tight_layout()

    params = {'mathtext.default': 'regular' }          
    plt.rcParams.update(params)

    # Save the file and close the figure

    plt.savefig(os.path.join(out_dir,'%s_orb.png' % mol_group),dpi=200)
    plt.close()

    print('%12s' % "[ DONE ]")

  # =================================================================== #
  # =================================================================== #
  #                   PLOTTING IONIZATION POTENTIALS                    #
  # =================================================================== #
  # =================================================================== #

  section_title = "3. Ionization potentials"

  print("")
  print(''.center(len(section_title)+10, '*'))
  print(section_title.center(len(section_title)+10))
  print(''.center(len(section_title)+10, '*'))

  # Iterate over each molecule group
  # ================================

  for mol_group in mol_groups:

      print ("{:<140}".format('\nTreating the values for %s QDs ...' % mol_group), end="")

      # Get the values
      # ==============

      ips_all = []
      ips_koop = []
      ips_vert = []
      ips_adiab = []

      for mol in yml_sorted[mol_group]:
        if yml_sorted[mol_group][mol].get('IPs (Ha)'):

          # Get the size values

          size = yml_sorted[mol_group][mol]['Structure']['Size (nm)']

          # Get the IPs values and add them to their corresponding lists (separate lists because all sizes are not necessarily represented for each IP)

          ip_koop = results_common.energy_unit_conversion(yml_sorted[mol_group][mol]['IPs (Ha)']['Koopmans'],"ha","ev")
          ips_koop.append((size,ip_koop))

          ip_vert = yml_sorted[mol_group][mol]['IPs (Ha)']['Vertical']
          if ip_vert != 'N/A':
             ip_vert = results_common.energy_unit_conversion(ip_vert,"ha","ev")
             ips_vert.append((size,ip_vert))

          ip_adiab = yml_sorted[mol_group][mol]['IPs (Ha)']['Adiabatic']
          if ip_adiab != 'N/A':
             ip_adiab = results_common.energy_unit_conversion(ip_adiab,"ha","ev")
             ips_adiab.append((size,ip_adiab))       

          # Store the data for this molecule

          ips_all.append({
            "Molecule": yml_sorted[mol_group][mol]['ID']['Name'],
            "Diameter (nm)": size,
            "Koopmans (eV)": ip_koop,
            "Vertical (eV)": ip_vert,
            "Adiabatic (eV)": ip_adiab
            })

      # Sort the values by size

      ips_koop.sort(key=lambda tup: tup[0])
      ips_vert.sort(key=lambda tup: tup[0])
      ips_adiab.sort(key=lambda tup: tup[0])

      ips_all.sort(key=lambda mol: mol['Diameter (nm)'])

      # Store the values in a CSV file
      # ==============================

      csv_header = list(ips_all[0].keys())

      with open(os.path.join(out_dir,'%s_ips.csv' % mol_group), 'w', newline='', encoding='utf-8') as csvfile:

        csv_writer = csv.DictWriter(csvfile, fieldnames=csv_header, delimiter=';', quoting=csv.QUOTE_MINIMAL)
        csv_writer.writeheader()

        for mol in ips_all:
          csv_writer.writerow(mol) 

      # Plot the graphs
      # ===============

      plt.style.use('seaborn-colorblind')

      fig, ax = plt.subplots()

      # Plot the values

      ax.plot([mol[0] for mol in ips_koop],[mol[1] for mol in ips_koop],marker='.',linestyle='-',label='Koopmans')
      ax.plot([mol[0] for mol in ips_vert],[mol[1] for mol in ips_vert],marker='.',linestyle='-',label='Vertical')
      ax.plot([mol[0] for mol in ips_adiab],[mol[1] for mol in ips_adiab],marker='.',linestyle='-',label='Adiabatic')

      # Add the legend and titles

      ax.set_title('Ionization potentials for %s QDs' % mol_group)
      ax.set_xlabel("Diameter of the QD (nm)")
      ax.set_ylabel('Energy (eV)')
      ax.legend()

      # Set other parameters

      ax.tick_params(top=False, right=False)
      ax.xaxis.set_minor_locator(AutoMinorLocator(2))
      ax.yaxis.set_minor_locator(AutoMinorLocator(2))

      plt.tight_layout()
      plt.grid(True,which='both',linestyle='--')

      # Save the file and close the figure

      plt.savefig(os.path.join(out_dir,'%s_ips.png' % mol_group),dpi=200)
      plt.close()

      print('%12s' % "[ DONE ]")

  # =================================================================== #
  # =================================================================== #
  #                     PLOTTING HIGHEST FIDELITIES                     #
  # =================================================================== #
  # =================================================================== #

  section_title = "4. Fidelities"

  print("")
  print(''.center(len(section_title)+10, '*'))
  print(section_title.center(len(section_title)+10))
  print(''.center(len(section_title)+10, '*'))

  # Treat the values for Si
  # =======================

  print ("{:<140}".format('\nTreating the values for Si QDs ...'), end="")

  si_fids_all = []
  si_fids_spgw = []
  si_fids_mgw = []

  for mol in yml_sorted['Si']:
    if yml_sorted['Si'][mol].get('Control'):

      # Get the size values

      size = yml_sorted['Si'][mol]['Structure']['Size (nm)']

      # Get the fidelities values for each config and add them to their corresponding lists (separate lists because all sizes are not necessarily represented for each type of config)

      values_spgw = [transition['Fidelity'] for transition in yml_sorted['Si'][mol]['Control']['Transitions'] if not transition['Label'].startswith('R_') and transition['Config'] == 'opc_filters']

      if values_spgw != []:
        max_spgw = max(values_spgw)
        si_fids_spgw.append((size,max_spgw))
      else:
        max_spgw = "N/A"

      values_mgw = [transition['Fidelity'] for transition in yml_sorted['Si'][mol]['Control']['Transitions'] if not transition['Label'].startswith('R_') and transition['Config'] == 'mgw']

      if values_mgw != []:
        max_mgw = max(values_mgw)
        si_fids_mgw.append((size,max_mgw))
      else:
        max_mgw = "N/A"

      # Store the data for this molecule

      si_fids_all.append({
        "Molecule": yml_sorted['Si'][mol]['ID']['Name'],
        "Diameter (nm)": size,
        "Super Gaussian filter": max_spgw,
        "Multi Gaussian filter": max_mgw
        })

      # Sort the values by size

      si_fids_spgw.sort(key=lambda tup: tup[0])
      si_fids_mgw.sort(key=lambda tup: tup[0])

      si_fids_all.sort(key=lambda mol: mol['Diameter (nm)'])

  # Store the values in a CSV file
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  csv_header = list(si_fids_all[0].keys())

  with open(os.path.join(out_dir,'Si_fids.csv'), 'w', newline='', encoding='utf-8') as csvfile:

    csv_writer = csv.DictWriter(csvfile, fieldnames=csv_header, delimiter=';', quoting=csv.QUOTE_MINIMAL)
    csv_writer.writeheader()

    for mol in si_fids_all:
      csv_writer.writerow(mol) 

  # Plot the graphs
  # ~~~~~~~~~~~~~~~

  plt.style.use('seaborn-colorblind')

  fig, ax = plt.subplots()

  ax.plot([mol[0] for mol in si_fids_spgw],[mol[1] for mol in si_fids_spgw],marker='.',linestyle='--',label='Si (Super Gaussian)')
  ax.plot([mol[0] for mol in si_fids_mgw],[mol[1] for mol in si_fids_mgw],marker='.',linestyle='--',label='Si (Multi Gaussian)')

  # Add the legend and titles

  ax.set_title('Highest fidelities: Si')
  ax.set_xlabel("Diameter of the QD (nm)")
  ax.set_ylabel('Fidelity')
  ax.legend()

  # Set other parameters

  ax.tick_params(top=False, right=False)
  ax.xaxis.set_minor_locator(AutoMinorLocator(2))
  ax.yaxis.set_minor_locator(AutoMinorLocator(2))

  plt.tight_layout()
  plt.grid(True,which='both',linestyle='--')

  # Save the file and close the figure

  plt.savefig(os.path.join(out_dir,'Si_fids.png'),dpi=200)
  plt.close()

  print('%12s' % "[ DONE ]")

  # Iterate over each molecule group and compare it to Si
  # ====================================================

  for mol_group in mol_groups:

    if mol_group != 'Si':

      print ("{:<140}".format('\nTreating the values for %s QDs ...' % mol_group), end="")

      # Get the values
      # ==============

      fids_all = []
      fids_spgw = []
      fids_mgw = []

      for mol in yml_sorted[mol_group]:
        if yml_sorted[mol_group][mol].get('Control'):

          # Get the size values

          size = yml_sorted[mol_group][mol]['Structure']['Size (nm)']

          # Get the fidelities values for each config and add them to their corresponding lists (separate lists because all sizes are not necessarily represented for each type of config)

          values_spgw = [transition['Fidelity'] for transition in yml_sorted[mol_group][mol]['Control']['Transitions'] if not transition['Label'].startswith('R_') and transition['Config'] == 'opc_filters']

          if values_spgw != []:
            max_spgw = max(values_spgw)
            fids_spgw.append((size,max_spgw))
          else:
            max_spgw = "N/A"

          values_mgw = [transition['Fidelity'] for transition in yml_sorted[mol_group][mol]['Control']['Transitions'] if not transition['Label'].startswith('R_') and transition['Config'] == 'mgw']

          if values_mgw != []:
            max_mgw = max(values_mgw)
            fids_mgw.append((size,max_mgw))
          else:
            max_mgw = "N/A"

          # Store the data for this molecule

          fids_all.append({
            "Molecule": yml_sorted[mol_group][mol]['ID']['Name'],
            "Diameter (nm)": size,
            "Super Gaussian filter": max_spgw,
            "Multi Gaussian filter": max_mgw
            })

          # Sort the values by size

          fids_spgw.sort(key=lambda tup: tup[0])
          fids_mgw.sort(key=lambda tup: tup[0])

          fids_all.sort(key=lambda mol: mol['Diameter (nm)'])

      # Store the values in a CSV file

      csv_header = list(fids_all[0].keys())

      with open(os.path.join(out_dir,'%s_fids.csv' % mol_group), 'w', newline='', encoding='utf-8') as csvfile:

        csv_writer = csv.DictWriter(csvfile, fieldnames=csv_header, delimiter=';', quoting=csv.QUOTE_MINIMAL)
        csv_writer.writeheader()

        for mol in fids_all:
          csv_writer.writerow(mol) 

      # Plot the graphs
      # ===============

      plt.style.use('seaborn-colorblind')

      fig, ax = plt.subplots()

      # Plot the Si values

      ax.plot([mol[0] for mol in si_fids_spgw],[mol[1] for mol in si_fids_spgw],marker='.',linestyle='--',label='Si (Super Gaussian)')
      ax.plot([mol[0] for mol in si_fids_mgw],[mol[1] for mol in si_fids_mgw],marker='.',linestyle='--',label='Si (Multi Gaussian)')

      # Plot the specific group value

      ax.plot([mol[0] for mol in fids_spgw],[mol[1] for mol in fids_spgw],marker='.',linestyle='-',label= mol_group + ' (Super Gaussian)')
      ax.plot([mol[0] for mol in fids_mgw],[mol[1] for mol in fids_mgw],marker='.',linestyle='-',label= mol_group + ' (Multi Gaussian)')

      # Add the legend and titles

      ax.set_title('Highest fidelities: Si vs %s' % mol_group)
      ax.set_xlabel("Diameter of the QD (nm)")
      ax.set_ylabel('Fidelity')
      ax.legend()

      # Set other parameters

      ax.tick_params(top=False, right=False)
      ax.xaxis.set_minor_locator(AutoMinorLocator(2))
      ax.yaxis.set_minor_locator(AutoMinorLocator(2))

      plt.tight_layout()
      plt.grid(True,which='both',linestyle='--')

      # Save the file and close the figure

      plt.savefig(os.path.join(out_dir,'%s_fids.png' % mol_group),dpi=200)
      plt.close()

      print('%12s' % "[ DONE ]") 
  
  # =================================================================== #
  # =================================================================== #
  #                  PLOTTING TRANSITION DIPOLE MOMENT                  #
  # =================================================================== #
  # =================================================================== #

  section_title = "5. Transition dipole moments"

  print("")
  print(''.center(len(section_title)+10, '*'))
  print(section_title.center(len(section_title)+10))
  print(''.center(len(section_title)+10, '*'))

  # Iterate over each molecule group
  # ================================

  for mol_group in mol_groups:

      print ("{:<140}".format('\nTreating the values for %s QDs ...' % mol_group), end="")

      # Get the values
      # ==============

      momdips = []
      first_momdips = []
      second_momdips = []
      third_momdips = []

      for mol in yml_sorted[mol_group]:
        if yml_sorted[mol_group][mol].get('Transition dipole moments (au)'):

          # Get the size values

          size = yml_sorted[mol_group][mol]['Structure']['Size (nm)']

          # Get the values and convert them from atomic units to Debye (conversion factor from https://link.springer.com/content/pdf/bbm%3A978-3-319-89972-5%2F1.pdf)
          # Then add them to their corresponding lists (separate lists because some molecules do not have three different values)

          conv_factor = 2.541746
          data = yml_sorted[mol_group][mol]['Transition dipole moments (au)']
          sortedkeys = sorted(data, key=str.lower) # Sort the keys alphabetically to ensure increasing order of states (S1 then S2, then S3, ...)

          first_value = data[sortedkeys[0]] * conv_factor
          first_momdips.append((size,first_value))
          second_value = "N/A"
          third_value = "N/A"

          if len(sortedkeys) >= 2:
            second_value = data[sortedkeys[1]] * conv_factor
            second_momdips.append((size,second_value))
          if len(sortedkeys) >= 3:
            third_value = data[sortedkeys[2]] * conv_factor
            third_momdips.append((size,third_value))
    
          # Store the data for this molecule

          momdips.append({
            "Molecule": yml_sorted[mol_group][mol]['ID']['Name'],
            "Diameter (nm)": size,
            "First value (D)": first_value,
            "Second value (D)": second_value,
            "Third value (D)": third_value,
            "Comment" : "  ".join(sortedkeys)
            })

      # Sort the values by size

      first_momdips.sort(key=lambda tup: tup[0])
      second_momdips.sort(key=lambda tup: tup[0])
      third_momdips.sort(key=lambda tup: tup[0])

      momdips.sort(key=lambda mol: mol['Diameter (nm)'])

      # Store the values in a CSV file
      # ==============================

      csv_header = list(momdips[0].keys())

      with open(os.path.join(out_dir,'%s_momdips.csv' % mol_group), 'w', newline='', encoding='utf-8') as csvfile:

        csv_writer = csv.DictWriter(csvfile, fieldnames=csv_header, delimiter=';', quoting=csv.QUOTE_MINIMAL)
        csv_writer.writeheader()

        for mol in momdips:
          csv_writer.writerow(mol) 

      # Plot the graphs
      # ===============

      plt.style.use('seaborn-colorblind')

      fig, ax = plt.subplots()

      # Plot the values

      ax.plot([mol[0] for mol in first_momdips],[mol[1] for mol in first_momdips],marker='.',linestyle='-',label="First singlet")
      ax.plot([mol[0] for mol in second_momdips],[mol[1] for mol in second_momdips],marker='.',linestyle='-',label="Second singlet")
      ax.plot([mol[0] for mol in third_momdips],[mol[1] for mol in third_momdips],marker='.',linestyle='-',label="Third singlet")

      # Add the legend and titles

      ax.set_title('Transition dipole moments for %s QDs' % mol_group)
      ax.set_xlabel("Diameter of the QD (nm)")
      ax.set_ylabel("Dipole moment (D)")
      ax.legend()

      # Set other parameters

      ax.tick_params(top=False, right=False)
      ax.xaxis.set_minor_locator(AutoMinorLocator(2))
      ax.yaxis.set_minor_locator(AutoMinorLocator(2))

      plt.tight_layout()
      plt.grid(True,which='both',linestyle='--')

      # Save the file and close the figure

      plt.savefig(os.path.join(out_dir,'%s_momdips.png' % mol_group),dpi=200)
      plt.close()

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
