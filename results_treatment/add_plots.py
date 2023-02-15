#!/usr/bin/env python3

################################################################################################################################################
##                                                     Additional Results Plotting Script                                                     ##
##                                                                                                                                            ##
##                                 This script produces specific tables and graphs based on specific results.                                 ##
##                                                                                                                                            ##
##                              /!\ In order to run, this script requires Python 3.5+ as well as matplotlib /!\                               ##
##                                          /!\ Ask your cluster(s) administrator(s) if needed. /!\                                           ##
################################################################################################################################################

import argparse
import csv
import os
import shutil
from inspect import getsourcefile

import matplotlib.pyplot as plt
import numpy as np
import yaml
from cycler import cycler
from matplotlib.ticker import AutoMinorLocator
from scipy import constants

import results_common

# =================================================================== #
# =================================================================== #
#                        FUNCTIONS DEFINITIONS                        #
# =================================================================== #
# =================================================================== #

def compute_size(nb_atoms:int):
    """ 
    Computes the diameter of a silicon cluster by using the number of silicon atoms and our own fitting equation:
    $d (\text{nm}) = 0.3738 (N_a)^0.3276$ 
    where $N_a$ is the number of silicon atoms in the cluster

    Parameters
    ----------
    nb_atoms : int
        Number of silicon atoms in the cluster

    Returns
    -------
    size : float
        The computed diameter, in nanometers
    """

    size = 0.3738 * (nb_atoms ** 0.3276)

    return size

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

    res_dpi = 400

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
      cycler(marker=['o', 'd', '^', 'p', 'v', 's']) +
      cycler(color=['black', 'orange', 'purple', 'darkturquoise', 'magenta', 'lawngreen']))

    line_cycler = (
      cycler(marker=['.', '+', 'x', '*']) +
      cycler(color=['blue', 'red', 'green', 'gold']) +
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

  ax.plot(puzder_x,puzder_y,linestyle='-',color='green',label="Fonction de Puzder et al. (théor.)")

  ax.set_prop_cycle(line_cycler)

  ax.plot([float(mol["Nombre d'atomes de Si"]) for mol in sizes_list],[float(mol['Diamètre (nm)']) for mol in sizes_list],marker='x',label='Valeurs de ce travail (optimisées)')
  ax.plot([float(mol["Nombre d'atomes de Si"]) for mol in sizes_list],[float(mol['Diamètre original (nm)']) for mol in sizes_list],marker='^',label='Valeurs de ce travail (originales)')

  ax.set_prop_cycle(dots_cycler)

  for article in sizes_litt:
    values_x = list(sizes_litt[article]['Values'].keys())
    valuess_y = list(sizes_litt[article]['Values'].values())
    ax.plot(values_x,valuess_y,label=sizes_litt[article]['Label'],linestyle=' ')

  # Add the legend and titles

  ax.set_xlabel("Nombre d'atomes de Si")
  ax.set_ylabel("Diamètre (nm)")
  ax.legend()

  # Set other parameters

  ax.tick_params(top=False, right=False)
  ax.xaxis.set_minor_locator(AutoMinorLocator(2))
  ax.yaxis.set_minor_locator(AutoMinorLocator(2))

  plt.tight_layout()
  plt.grid(True,which='both',linestyle='--')

  # Save the file and close the figure

  plt.savefig(os.path.join(out_dir,'sizes.png'),dpi=res_dpi)
  plt.close()

  print('%12s' % "[ DONE ]")

  # =================================================================== #
  # =================================================================== #
  #                      ENERGY GAPS vs FUNCTIONALS                     #
  # =================================================================== #
  # =================================================================== #

  section_count += 1
  section_title = str(section_count) + ". Energy gaps vs functionals"

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
    x_axis = gaps_litt[article].get('X_axis')
    if x_axis == 'atoms':
      values_x = [compute_size(nb_atoms) for nb_atoms in list(gaps_litt[article]['Values'].keys())]
    else:
      values_x = list(gaps_litt[article]['Values'].keys())
    valuess_y = list(gaps_litt[article]['Values'].values())
    ax.plot(values_x,valuess_y,label=gaps_litt[article]['Label'],linestyle=' ')

  # Annotate a set of points to indicate which point corresponds to which molecule
  # See https://queirozf.com/entries/add-labels-and-text-to-matplotlib-plots-annotation-examples

  for molecule in [mol['Molecule'] for mol in funct_data['B3LYP']]:

    max_gap = -float('inf')

    for funct in functs:
      gap = [float(mol['Gap optique (eV)']) for mol in funct_data[funct] if mol['Molecule'] == molecule][0]
      if gap > max_gap:
        max_gap = gap
        chosen_funct = funct

    diam = [float(mol['Diametre (nm)']) for mol in funct_data[chosen_funct] if mol['Molecule'] == molecule][0]

    if molecule == "Si$_{5}$H$_{12}$":
      txt_offset = (50,10)
    elif molecule == "Si$_{71}$H$_{84}$":
      txt_offset = (25,50)
    elif molecule == "Si$_{99}$H$_{100}$":
      txt_offset = (60,10)
    else:
      txt_offset = (25,40)

    ax.annotate(molecule,
                xy=(diam, gap), xycoords='data',
                xytext=txt_offset, textcoords='offset points',
                arrowprops=dict(facecolor='black', arrowstyle="->"),
                horizontalalignment='right', verticalalignment='top')

  # Add the legend and titles

  ax.set_xlabel("Diamètre (nm)")
  ax.set_ylabel("Gap optique (eV)")
  ax.legend()

  # Set other parameters

  ax.tick_params(top=False, right=False)
  ax.xaxis.set_minor_locator(AutoMinorLocator(2))
  ax.yaxis.set_minor_locator(AutoMinorLocator(2))

  plt.tight_layout()
  plt.grid(True,which='both',linestyle='--')

  # Save the file and close the figure

  plt.savefig(os.path.join(out_dir,'gaps.png'),dpi=res_dpi)
  plt.close()

  print('%12s' % "[ DONE ]")

  # =================================================================== #
  # =================================================================== #
  #                 ENERGY GAPS Si vs SiBP vs SiP vs SiB                #
  # =================================================================== #
  # =================================================================== #

  section_count += 1
  section_title = str(section_count) + ". Energy gaps Si/SiBP/SiB/SiP"

  print("")
  print(''.center(len(section_title)+10, '*'))
  print(section_title.center(len(section_title)+10))
  print(''.center(len(section_title)+10, '*'))

  # Get the values
  # ==============

  # Our values

  groups_data = {}
  groups = ["Si","SiBP","SiB","SiP"]

  for group in groups:

    # Load each functional file

    group_file = results_common.check_abspath(os.path.join(inp_dir,"gaps_%s.csv" % group.lower()),"%s gaps CSV file" % group,"file")

    print ("{:<140}".format("\nLoading %s gaps CSV file ..." % group), end="")
    with open(group_file, 'r', newline='', encoding='utf-8') as csv_file:
      group_content = csv.DictReader(csv_file, delimiter=';')
      group_list = list(group_content)
    print('%12s' % "[ DONE ]")

    # Replace commas by point (French decimal to US decimal)

    for line in group_list:
      for key, value in line.items():
        line[key] = value.replace(",",".")

    # Exclude empty lines

    groups_data[group] = [line for line in group_list if line['Diamètre (nm)'] != '']

  # Plot the graph
  # ==============

  print ("{:<140}".format('\nPlotting the graph ...'), end="")

  plt.style.use('seaborn-colorblind')

  fig, ax = plt.subplots()

  # Plot the values

  ax.set_prop_cycle(line_cycler)

  for group in groups:
    ax.plot([float(mol['Diamètre (nm)']) for mol in groups_data[group]],[float(mol['Gap optique (eV)']) for mol in groups_data[group]],linestyle='-',label=group)

  # Add the legend and titles

  ax.set_xlabel("Diamètre (nm)")
  ax.set_ylabel("Gap optique (eV)")
  ax.legend()

  # Set other parameters

  ax.tick_params(top=False, right=False)
  ax.xaxis.set_minor_locator(AutoMinorLocator(2))
  ax.yaxis.set_minor_locator(AutoMinorLocator(2))

  plt.tight_layout()
  plt.grid(True,which='both',linestyle='--')

  # Save the file and close the figure

  plt.savefig(os.path.join(out_dir,'gaps_sibp.png'),dpi=res_dpi)
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

  melnikov_x = np.linspace(0.5, 2, 100)
  melnikov_y = 4.8 + (44.4/((melnikov_x*(1e-9)/constants.value('Bohr radius'))**1.2))

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

  ax.set_prop_cycle(line_cycler)

  ax.plot([float(mol['Diamètre (nm)']) for mol in ips_list],[float(mol['P.I. vertical (eV)']) for mol in ips_list],label='Valeurs de ce travail (vertical)')
  ax.plot([float(mol['Diamètre (nm)']) for mol in ips_list if mol['P.I. adiabatique (eV)'] != ""],[float(mol['P.I. adiabatique (eV)']) for mol in ips_list if mol['P.I. adiabatique (eV)'] != ""],label='Valeurs de ce travail (adiabatique)')

  #ax.plot(melnikov_x,melnikov_y,label="Fonction de Melnikov et al. (PP-LDA")

  ax.set_prop_cycle(dots_cycler)

  for article in ips_litt:
    x_axis = ips_litt[article].get('X_axis')
    if x_axis == 'atoms':
      values_x = [compute_size(nb_atoms) for nb_atoms in list(ips_litt[article]['Values'].keys())]
    else:
      values_x = list(ips_litt[article]['Values'].keys())
    valuess_y = list(ips_litt[article]['Values'].values())
    ax.plot(values_x,valuess_y,label=ips_litt[article]['Label'],linestyle=' ')

  # Add the legend and titles

  ax.set_xlabel("Diamètre (nm)")
  ax.set_ylabel("Potentiel d'ionisation (eV)")
  ax.legend()

  # Set other parameters

  ax.tick_params(top=False, right=False)
  ax.xaxis.set_minor_locator(AutoMinorLocator(2))
  ax.yaxis.set_minor_locator(AutoMinorLocator(2))

  plt.tight_layout()
  plt.grid(True,which='both',linestyle='--')

  # Save the file and close the figure

  plt.savefig(os.path.join(out_dir,'ips.png'),dpi=res_dpi)
  plt.close()

  print('%12s' % "[ DONE ]")

  # =================================================================== #
  # =================================================================== #
  #                         GENERAL COMPARISONS                         #
  # =================================================================== #
  # =================================================================== #

  section_count += 1
  section_title = str(section_count) + ". General comparisons"

  print("")
  print(''.center(len(section_title)+10, '*'))
  print(section_title.center(len(section_title)+10))
  print(''.center(len(section_title)+10, '*'))

  # Plot the average pS -> pT transition dipole moments graph
  # =========================================================

  print ("{:<140}".format('\nPlotting the average pS -> pT transition dipole moments graph ...'), end="")

  plt.style.use('seaborn-colorblind')

  fig, ax = plt.subplots()

  # Define the values

  ps_pt_momdip_si =   [0.001494296, 0.007773517, 0.000407183, 0.009551477, 0.000596850, 0.000488429, 0.003129356]
  ps_pt_momdip_sige = [0.010676326, 0.067770619, 0.031168208, 0.111790101, 0.046776813, 0.026332329, 0]
  ps_pt_momdip_sibp = [0,           0.280754039, 0.230401132, 0.059531313, 0.493414296, 0.538615300, 0]

  barwidth = 0.15

  br1 = np.arange(len(ps_pt_momdip_si))
  br2 = [bar + barwidth for bar in br1]
  br3 = [bar + barwidth for bar in br2]

  # Define the X labels and ticks

  xlabels = [r'Si$_{5}$H$_{12}$',r'Si$_{17}$H$_{36}$',r'Si$_{29}$H$_{36}$',r'Si$_{47}$H$_{60}$',r'Si$_{71}$H$_{84}$',r'Si$_{87}$H$_{76}$',r'Si$_{99}$H$_{100}$']
  xticks = [r + barwidth for r in range(len(xlabels))]

  # Plot the values

  ax.bar(br1, ps_pt_momdip_si, width = barwidth, edgecolor ='grey', label = "Si")
  ax.bar(br2, ps_pt_momdip_sige, width = barwidth, edgecolor ='grey', label = "SiGe")
  ax.bar(br3, ps_pt_momdip_sibp, width = barwidth, edgecolor ='grey', label = "SiBP")

  # Add the legend and titles

  ax.set_ylabel('Moyenne des moments dipolaires ' + r'$\bar{\mu}_{pS \to pT}$' + ' (u.a.)')
  ax.set_yscale('log')
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

  plt.savefig(os.path.join(out_dir,'all_ps_pt_momdip.png'),dpi=res_dpi)
  plt.close()  

  print('%12s' % "[ DONE ]")

  # Plot the highest performance graph
  # ==================================

  print ("{:<140}".format('\nPlotting the highest performance graph ...'), end="")

  plt.style.use('seaborn-colorblind')

  fig, ax = plt.subplots()

  # Define the values

  high_perf_si =   [7.65E-05, 0.981357, 8.78E-08,    0.000103052, 0.000194045, 0.490496, 8.43E-08]
  high_perf_sige = [0.662969, 0.850243, 0.980076,    0.998173,    1,           0.998252, 0]
  high_perf_sibp = [0,        0.401346, 0.000100033, 0.396336,    0.0273536,   0.552744, 0]

  barwidth = 0.15

  br1 = np.arange(len(high_perf_si))
  br2 = [bar + barwidth for bar in br1]
  br3 = [bar + barwidth for bar in br2]

  # Define the X labels and ticks

  xlabels = [r'Si$_{5}$H$_{12}$',r'Si$_{17}$H$_{36}$',r'Si$_{29}$H$_{36}$',r'Si$_{47}$H$_{60}$',r'Si$_{71}$H$_{84}$',r'Si$_{87}$H$_{76}$',r'Si$_{99}$H$_{100}$']
  xticks = [r + barwidth for r in range(len(xlabels))]

  # Plot the values

  ax.bar(br1, high_perf_si, width = barwidth, edgecolor ='grey', label = "Si")
  ax.bar(br2, high_perf_sige, width = barwidth, edgecolor ='grey', label = "SiGe")
  ax.bar(br3, high_perf_sibp, width = barwidth, edgecolor ='grey', label = "SiBP")

  # Add the legend and titles

  ax.set_ylabel('Indice de performance')
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

  plt.savefig(os.path.join(out_dir,'highest_perf.png'),dpi=res_dpi)
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
