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
from matplotlib.ticker import AutoMinorLocator
import matplotlib.mathtext

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
  #                        PLOTTING ENERGY GAPS                         #
  # =================================================================== #
  # =================================================================== #

  section_title = "1. Energy gaps"

  print("")
  print(''.center(len(section_title)+10, '*'))
  print(section_title.center(len(section_title)+10))
  print(''.center(len(section_title)+10, '*'))

  # Get the values for Si
  # =====================

  print ("{:<140}".format('\nFetching the values for Si QDs ...'), end="")

  si_gaps = []

  for mol in yml_sorted['Si']:
    if yml_sorted['Si'][mol].get('Energy gaps (Ha)'):
      size = yml_sorted['Si'][mol]['Structure']['Size (nm)']
      hl_gap = results_common.energy_unit_conversion(yml_sorted['Si'][mol]['Energy gaps (Ha)']['HOMO-LUMO'],"ha","ev")
      opt_gap = results_common.energy_unit_conversion(yml_sorted['Si'][mol]['Energy gaps (Ha)']['Optical'],"ha","ev")
      st_gap = results_common.energy_unit_conversion(yml_sorted['Si'][mol]['Energy gaps (Ha)']['Singlet-Triplet'],"ha","ev")
      si_gaps.append((size,hl_gap,opt_gap,st_gap))

  si_gaps.sort(key=lambda tup: tup[0]) # Sort the values by size

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

          size = yml_sorted[mol_group][mol]['Structure']['Size (nm)']
          hl_gap = results_common.energy_unit_conversion(yml_sorted[mol_group][mol]['Energy gaps (Ha)']['HOMO-LUMO'],"ha","ev")
          opt_gap = results_common.energy_unit_conversion(yml_sorted[mol_group][mol]['Energy gaps (Ha)']['Optical'],"ha","ev")
          st_gap = results_common.energy_unit_conversion(yml_sorted[mol_group][mol]['Energy gaps (Ha)']['Singlet-Triplet'],"ha","ev")

          # Store the data for this molecule

          gaps.append((size,hl_gap,opt_gap,st_gap))

      gaps.sort(key=lambda tup: tup[0]) # Sort the values by size

      # Plot the optical and HOMO-LUMO gaps graphs
      # ==========================================

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

      ax.plot([mol[0] for mol in si_gaps],[mol[3] for mol in si_gaps],marker='.',linestyle='--',label='S-T gaps - Si (TD-DFT)')

      # Plot the specific group value

      ax.plot([mol[0] for mol in gaps],[mol[3] for mol in gaps],marker='.',linestyle='-',label='S-T gaps - %s (TD-DFT)' % mol_group)

      # Add the legend and titles

      ax.set_title('Singlet-Triplet gaps: Si vs %s' % mol_group)
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
    plt.xticks(ticks=xticks, labels=xlabels, rotation=45)
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

      ips_koop = []
      ips_vert = []
      ips_adiab = []

      for mol in yml_sorted[mol_group]:
        if yml_sorted[mol_group][mol].get('IPs (Ha)'):

          # Get the size values

          size = yml_sorted[mol_group][mol]['Structure']['Size (nm)']

          # Get the IPs values and add them to their corresponding lists

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

      # Sort the values by number of Si atoms

      ips_koop.sort(key=lambda tup: tup[0])
      ips_vert.sort(key=lambda tup: tup[0])
      ips_adiab.sort(key=lambda tup: tup[0])

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
  #                  PLOTTING TRANSITION DIPOLE MOMENT                  #
  # =================================================================== #
  # =================================================================== #

  section_title = "4. Transition dipole moment"

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

      first_momdip_list = []

      for mol in yml_sorted[mol_group]:
        if yml_sorted[mol_group][mol]['Structure'].get('First non-zero transition dipole moment (au)'):

          # Get the value

          size = yml_sorted[mol_group][mol]['Structure']['Size (nm)']
          s0_s1_momdip = yml_sorted[mol_group][mol]['Structure']['S0-S1 transition dipole moment (au)']
          first_momdip = yml_sorted[mol_group][mol]['Structure']['First non-zero transition dipole moment (au)']
          
          # Store the data for this molecule

          first_momdip_list.append((size,s0_s1_momdip,first_momdip))

      first_momdip_list.sort(key=lambda tup: tup[0]) # Sort the values by size

      # Plot the singlet-triplet gaps graphs
      # ====================================

      plt.style.use('seaborn-colorblind')

      fig, ax = plt.subplots()

      # Plot the values

      ax.plot([mol[0] for mol in first_momdip_list],[mol[1] for mol in first_momdip_list],marker='.',linestyle='-',label="S0-S1 µ")
      ax.plot([mol[0] for mol in first_momdip_list],[mol[2] for mol in first_momdip_list],marker='.',linestyle='-',label="First non-zero µ")

      # Add the legend and titles

      ax.set_title('First transition dipole moments for %s group' % mol_group)
      ax.set_xlabel("Diameter of the QD (nm)")
      ax.set_ylabel("Dipole moment (au)")
      ax.legend()

      # Set other parameters

      ax.tick_params(top=False, right=False)
      ax.xaxis.set_minor_locator(AutoMinorLocator(2))
      ax.yaxis.set_minor_locator(AutoMinorLocator(2))

      plt.tight_layout()
      plt.grid(True,which='both',linestyle='--')

      # Save the file and close the figure

      plt.savefig(os.path.join(out_dir,'%s_first_momdip.png' % mol_group),dpi=200)
      plt.close()

      print('%12s' % "[ DONE ]")

# =================================================================== #
# =================================================================== #
#                          CALL MAIN FUNCTION                         #
# =================================================================== #
# =================================================================== #

# If this script is executed through the command line, call the main function (see https://realpython.com/python-main-function/ for details)

if __name__ == "__main__":
    main()   