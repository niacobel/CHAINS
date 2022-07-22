#!/usr/bin/env python3

################################################################################################################################################
##                                                     Additional Results Plotting Script                                                     ##
##                                                                                                                                            ##
##                                 This script produces specific tables and graphs based on specific results.                                 ##
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
import yaml
from matplotlib.cm import ScalarMappable
from matplotlib.ticker import AutoMinorLocator
from scipy import constants
from tabulate import tabulate
from cycler import cycler
import results_common

# =================================================================== #
# =================================================================== #
#                       COMMAND LINE ARGUMENTS                        #
# =================================================================== #
# =================================================================== #

# Define the arguments needed for the script (here they are defined as named arguments rather than positional arguments, check https://stackoverflow.com/questions/24180527/argparse-required-arguments-listed-under-optional-arguments for more info).

parser = argparse.ArgumentParser(add_help=False, description="This script extracts the relevant information from a single YAML file and produces the various needed tables and graphs.")

required = parser.add_argument_group('Required arguments')
required.add_argument("-i","--inp_dir", type=str, help="Path to the directory containing the results.", required=True)
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
    print("EXECUTION OF THE ADDITIONAL RESULTS PLOTTING SCRIPT BEGINS NOW".center(columns))
    print("")
    print("".center(columns,"*"))

    # ========================================================= #
    # Read command line arguments                               #
    # ========================================================= #

    args = parser.parse_args()

    # Required arguments

    inp_dir = args.inp_dir                   # YAML file containing the results
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
    # Check other arguments                                     #
    # ========================================================= #

    inp_dir = results_common.check_abspath(inp_dir,"Command line argument -i / --inp_dir","directory")
    print ("{:<40} {:<100}".format('\nInput directory:',inp_dir))

    out_dir = results_common.check_abspath(out_dir,"Command line argument -o / --out_dir","directory")
    print ("{:<40} {:<100}".format('\nOutput directory:',out_dir))

    section_count = 0

    # ========================================================= #
    # Define markers and colours cycle for Matplotlib's Pyplot  #
    # ========================================================= #

    # See https://stackoverflow.com/questions/13091649/unique-plot-marker-for-each-plot-in-matplotlib

    dots_cycler = (
      cycler(marker=['o', 'd', '^', 'p']) +
      cycler(color=['black', 'orange', 'purple','darkturquoise']))

    line_cycler = (
      cycler(marker=['.', '+', 'x', '*']) +
      cycler(color=['blue', 'red', 'green', 'magenta']) +
      cycler(linestyle=['-','--',':','-.']))

  # ========================================================= #
  # Exception handling for the preparation step               #
  # ========================================================= #

  except results_common.ResultsError as error:
    print("")
    print(error)
    exit(-1)

  # =================================================================== #
  # =================================================================== #
  #                     SIZE AND ATOMS RELATIONSHIP                     #
  # =================================================================== #
  # =================================================================== #

  section_count += 1
  section_title = str(section_count) + ". Sizes vs number of atoms"

  print("")
  print(''.center(len(section_title)+10, '*'))
  print(section_title.center(len(section_title)+10))
  print(''.center(len(section_title)+10, '*'))

  # Get the values
  # ==============

  # Our values

  sizes_file = results_common.check_abspath(os.path.join(inp_dir,"sizes.csv"),"Si clusters sizes CSV file","file")

  print ("{:<140}".format("\nLoading Si clusters sizes CSV file ..."), end="")
  with open(sizes_file, 'r', newline='') as csv_file:
    sizes_content = csv.DictReader(csv_file, delimiter=';')
    sizes_list = list(sizes_content)
  print('%12s' % "[ DONE ]")

  # Puzder et al's function (https://aip.scitation.org/doi/pdf/10.1063/1.1504707)
  
  puzder_x = np.linspace(0, 150, 100)
  puzder_y = 0.54310*(( (3*puzder_x) / (4*np.pi) )**(1/3)) # D = a * [(3/4pi)N]^1/3

  # Litterature values

  sizes_litt_file = results_common.check_abspath(os.path.join(inp_dir,"sizes_litt.yml"),"Sizes litterature values YAML file","file")

  print ("{:<50} {:<89}".format('\nLoading the sizes litterature values YAML file',sizes_litt_file + " ..."), end="")
  with open(sizes_litt_file, 'r') as f_yml:
    sizes_litt = yaml.load(f_yml, Loader=yaml.FullLoader)
  print('%12s' % "[ DONE ]")

  # Plot the graph
  # ==============

  print ("{:<140}".format('\nPlotting the graph ...'), end="")

  plt.style.use('seaborn-colorblind')

  fig, ax = plt.subplots()

  # Plot the values

  ax.plot(puzder_x,puzder_y,linestyle='-',color='green',label="Puzder et al.'s function")

  ax.set_prop_cycle(line_cycler)

  ax.plot([float(mol["Nombre d'atomes de Si"]) for mol in sizes_list],[float(mol['Diamètre (nm)']) for mol in sizes_list],marker='x',label='Values of this work (optimized)')
  ax.plot([float(mol["Nombre d'atomes de Si"]) for mol in sizes_list],[float(mol['Diamètre original (nm)']) for mol in sizes_list],marker='^',label='Values of this work (not optimized)')

  ax.set_prop_cycle(dots_cycler)

  for article in sizes_litt:
    values_x = list(sizes_litt[article]['Values'].keys())
    valuess_y = list(sizes_litt[article]['Values'].values())
    ax.plot(values_x,valuess_y,label=sizes_litt[article]['Label'],linestyle=' ')

  # Add the legend and titles

  ax.set_xlabel('Number of Si atoms')
  ax.set_ylabel("Diameter of the cluster (nm)")
  ax.legend()

  # Set other parameters

  ax.tick_params(top=False, right=False)
  ax.xaxis.set_minor_locator(AutoMinorLocator(2))
  ax.yaxis.set_minor_locator(AutoMinorLocator(2))

  plt.tight_layout()
  plt.grid(True,which='both',linestyle='--')

  # Save the file and close the figure

  plt.savefig(os.path.join(out_dir,'sizes.png'),dpi=300)
  plt.close()

  print('%12s' % "[ DONE ]")

  # =================================================================== #
  # =================================================================== #
  #                             ENERGY GAPS                             #
  # =================================================================== #
  # =================================================================== #

  section_count += 1
  section_title = str(section_count) + ". Energy gaps"

  print("")
  print(''.center(len(section_title)+10, '*'))
  print(section_title.center(len(section_title)+10))
  print(''.center(len(section_title)+10, '*'))

  # Get the values
  # ==============

  # Our values

  funct_data = {}
  functs = ["B3LYP","B3PW91","PBE0"]

  for funct in functs:

    # Load each functional file

    funct_file = results_common.check_abspath(os.path.join(inp_dir,"gaps_%s.csv" % funct.lower()),"%s gaps CSV file" % funct,"file")

    print ("{:<140}".format("\nLoading %s gaps CSV file ..." % funct), end="")
    with open(funct_file, 'r', newline='', encoding='utf-8') as csv_file:
      funct_content = csv.DictReader(csv_file, delimiter=';')
      funct_list = list(funct_content)
    print('%12s' % "[ DONE ]")

    # Replace commas by point (French decimal to US decimal)

    for line in funct_list:
      for key, value in line.items():
        line[key] = value.replace(",",".")

    # Exclude empty lines

    funct_data[funct] = [line for line in funct_list if line['Diametre (nm)'] != '']

  # Niatz and Zdetsis's function (https://pubs.acs.org/doi/10.1021/acs.jpcc.6b02955)
  
  zdetsis_x = np.linspace(0.6, 2.5, 100)
  zdetsis_y = 1.33 + (41.8 / (zdetsis_x*10))

  # Litterature values

  gaps_litt_file = results_common.check_abspath(os.path.join(inp_dir,"gaps_litt.yml"),"Gaps litterature values YAML file","file")

  print ("{:<50} {:<89}".format('\nLoading the gaps litterature values YAML file',gaps_litt_file + " ..."), end="")
  with open(gaps_litt_file, 'r') as f_yml:
    gaps_litt = yaml.load(f_yml, Loader=yaml.FullLoader)
  print('%12s' % "[ DONE ]")

  # Plot the graph
  # ==============

  print ("{:<140}".format('\nPlotting the graph ...'), end="")

  plt.style.use('seaborn-colorblind')

  fig, ax = plt.subplots()

  # Plot the values

  #ax.plot(zdetsis_x,zdetsis_y,linestyle='-',color='r',label="Niatz and Zdetsis's function")

  ax.set_prop_cycle(line_cycler)

  for funct in functs:
    ax.plot([float(mol['Diametre (nm)']) for mol in funct_data[funct]],[float(mol['Gap optique (eV)']) for mol in funct_data[funct]],linestyle='-',label=funct)

  ax.set_prop_cycle(dots_cycler)

  for article in gaps_litt:
    values_x = list(gaps_litt[article]['Values'].keys())
    valuess_y = list(gaps_litt[article]['Values'].values())
    ax.plot(values_x,valuess_y,label=gaps_litt[article]['Label'],linestyle=' ')

  # Annotate a set of points to indicate which point corresponds to which molecule
  # See https://queirozf.com/entries/add-labels-and-text-to-matplotlib-plots-annotation-examples

  for diam, gap in zip([float(mol['Diametre (nm)']) for mol in funct_data['B3LYP']],[float(mol['Gap optique (eV)']) for mol in funct_data['B3LYP']]):

    molecule = [mol['Molecule'] for mol in funct_data['B3LYP'] if float(mol['Diametre (nm)']) == diam][0]

    plt.annotate(molecule,                    # this is the text
                 (diam, gap),                 # these are the coordinates to position the label
                 textcoords="offset points",  # how to position the text
                 xytext=(0,5),                # distance from text to points (x,y)
                 ha='left')                  # horizontal alignment can be left, right or center

  # Add the legend and titles

  ax.set_xlabel('Diameter of the cluster (nm)')
  ax.set_ylabel("Optical gap (eV)")
  ax.legend()

  # Set other parameters

  ax.tick_params(top=False, right=False)
  ax.xaxis.set_minor_locator(AutoMinorLocator(2))
  ax.yaxis.set_minor_locator(AutoMinorLocator(2))

  plt.tight_layout()
  plt.grid(True,which='both',linestyle='--')

  # Save the file and close the figure

  plt.savefig(os.path.join(out_dir,'gaps.png'),dpi=300)
  plt.close()

  print('%12s' % "[ DONE ]")

  # =================================================================== #
  # =================================================================== #
  #                        IONIZATION POTENTIALS                        #
  # =================================================================== #
  # =================================================================== #

  section_count += 1
  section_title = str(section_count) + ". Ionization potentials"

  print("")
  print(''.center(len(section_title)+10, '*'))
  print(section_title.center(len(section_title)+10))
  print(''.center(len(section_title)+10, '*'))

  # Get the values
  # ==============

  # Our values

  ips_file = results_common.check_abspath(os.path.join(inp_dir,"ips.csv"),"Si IPs CSV file","file")

  print ("{:<140}".format("\nLoading Si IPs CSV file ..."), end="")
  with open(ips_file, 'r', newline='') as csv_file:
    ips_content = csv.DictReader(csv_file, delimiter=';')
    ips_list = list(ips_content)
  print('%12s' % "[ DONE ]")

  # Melnikov et al's function (https://journals.aps.org/prb/abstract/10.1103/PhysRevB.69.113305)

  #melnikov_x = np.linspace(0.5, 2, 100)
  #melnikov_y = 4.8 + (44.4/((melnikov_x*(1e-9)/constants.value('Bohr radius'))**1.2))

  # Litterature values

  ips_litt_file = results_common.check_abspath(os.path.join(inp_dir,"ips_litt.yml"),"IPs litterature values YAML file","file")

  print ("{:<50} {:<89}".format('\nLoading the IPs litterature values YAML file',ips_litt_file + " ..."), end="")
  with open(ips_litt_file, 'r') as f_yml:
    ips_litt = yaml.load(f_yml, Loader=yaml.FullLoader)
  print('%12s' % "[ DONE ]")

  # Plot the graph
  # ==============

  print ("{:<140}".format('\nPlotting the graph ...'), end="")

  plt.style.use('seaborn-colorblind')

  fig, ax = plt.subplots()

  # Plot the values

  ax.plot([float(mol['Diamètre (nm)']) for mol in ips_list],[float(mol[r'$-\text{E}_{\text{HOMO}}$~(eV)']) for mol in ips_list],marker='.',linestyle='-',label=r'$-E_{HOMO}$')
  ax.plot([float(mol['Diamètre (nm)']) for mol in ips_list],[float(mol['P.I. vertical (eV)']) for mol in ips_list],marker='.',linestyle='-',label='Vertical')
  ax.plot([float(mol['Diamètre (nm)']) for mol in ips_list if mol['P.I. adiabatique (eV)'] != ""],[float(mol['P.I. adiabatique (eV)']) for mol in ips_list if mol['P.I. adiabatique (eV)'] != ""],marker='.',linestyle='-',label='Adiabatic')
  #ax.plot(melnikov_x,melnikov_y,linestyle='-',color='r',label="Melnikov et al.'s function")

  # Add the legend and titles

  #ax.set_title('Ionization potentials for %s clusters' % mol_group)
  ax.set_xlabel("Diameter of the cluster (nm)")
  ax.set_ylabel('Energy (eV)')
  ax.legend()

  # Set other parameters

  ax.tick_params(top=False, right=False)
  ax.xaxis.set_minor_locator(AutoMinorLocator(2))
  ax.yaxis.set_minor_locator(AutoMinorLocator(2))

  plt.tight_layout()
  plt.grid(True,which='both',linestyle='--')

  # Save the file and close the figure

  plt.savefig(os.path.join(out_dir,'ips.png'),dpi=300)
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
