#!/usr/bin/env python3

################################################################################################################################################
##                                                           Pulses Plotting Script                                                           ##
##                                                                                                                                            ##
##                                     This script produces the graphs corresponding to specific pulses.                                      ##
##                                                                                                                                            ##
##                           /!\ In order to run, this script requires Python 3.5+ as well as YAML, matplotlib /!\                            ##
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

import matplotlib.mathtext
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from cycler import cycler
from matplotlib.ticker import AutoMinorLocator
from scipy import constants
from scipy.fft import rfft, rfftfreq

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

parser = argparse.ArgumentParser(add_help=False, description="This script produces the graphs corresponding to specific pulses.")

required = parser.add_argument_group('Required arguments')
required.add_argument("-o","--out_dir", type=str, help="Path to the directory where you want to store the graphs.", required=True)

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
    print("EXECUTION OF THE PULSES PLOTTING SCRIPT BEGINS NOW".center(columns))
    print("")
    print("".center(columns,"*"))

    # ========================================================= #
    # Read command line arguments                               #
    # ========================================================= #

    args = parser.parse_args()

    # Required arguments

    out_dir = args.out_dir                   # Directory where the graphs will be stored

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

    res_dpi = 400

    # ABIN LAUNCHER's geom_scan.py file
    # =================================

    chains_path = os.path.dirname(code_dir)
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

    # CONTROL LAUNCHER's transition_fcts.py file
    # ==========================================
    
    transition_fcts_path = os.path.join(chains_path,"control_launcher","transition_fcts.py")

    print ("{:<140}".format("\nImporting CONTROL LAUNCHER's transition_fcts.py file ..."), end="")
    transition_fcts = import_path(transition_fcts_path)
    print('%12s' % "[ DONE ]")

    # CONTROL LAUNCHER's control_common.py file
    # ==========================================
    
    control_common_path = os.path.join(chains_path,"control_launcher","control_common.py")

    print ("{:<140}".format("\nImporting CONTROL LAUNCHER's control_common.py file ..."), end="")
    control_common = import_path(control_common_path)
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
    # Check arguments                                           #
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
  #                          GRAPHS GENERATION                          #
  # =================================================================== #
  # =================================================================== #

  for mol_name in mol_inp_list:

    mol_dir = os.path.join(mol_inp_path, mol_name)

    console_message = "Start procedure for the molecule " + mol_name
    print("")
    print(''.center(len(console_message)+11, '*'))
    print(console_message.center(len(console_message)+10))
    print(''.center(len(console_message)+11, '*'))

    # =================================================================== #
    # =================================================================== #
    #                         MOLECULE INFORMATION                        #
    # =================================================================== #
    # =================================================================== #

    # For more information on try/except structures, see https://www.tutorialsteacher.com/python/exception-handling-in-python
    try:

      print ("{:<140}".format('\nFetching molecular information ...'))

      # Check the data directory

      data_dir = results_common.check_abspath(os.path.join(mol_dir,"CONTROL","aldu_param","data"),"Data directory created by control_launcher.py","directory")

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

      # ========================================================= #
      # Load the QCHEM output file                                #
      # ========================================================= #

      print ("{:<133}".format('\n\tLoading QCHEM output file ...'), end="")

      qchem_file = results_common.check_abspath(os.path.join(data_dir, mol_name + ".out"),"QCHEM output file","file")

      with open(qchem_file, 'r') as out_file:
        qchem_content = out_file.read().splitlines()

      qchem_content = list(map(str.strip, qchem_content))   # Remove leading & trailing blank/spaces
      qchem_content = list(filter(None, qchem_content))     # Remove blank lines/no char

      print('%12s' % "[ DONE ]")

      # ========================================================= #
      # Get the transitions list                                  #
      # ========================================================= #

      print ("{:<133}".format('\n\tLoading transitions list ...'), end="")

      # Call the modelling function from CONTROL LAUNCHER but without its standard output (https://stackoverflow.com/questions/2828953/silence-the-stdout-of-a-function-in-python-without-trashing-sys-stdout-and-resto)

      with open(os.devnull, 'w') as devnull:
        with contextlib.redirect_stdout(devnull):
          system = modelling_fcts.qchem_tddft(qchem_content)

      # Call the transition function from CONTROL LAUNCHER but without its standard output (https://stackoverflow.com/questions/2828953/silence-the-stdout-of-a-function-in-python-without-trashing-sys-stdout-and-resto)

      with open(os.devnull, 'w') as devnull:
        with contextlib.redirect_stdout(devnull):
          transitions_list = transition_fcts.brightest_to_darkest(system)

      print('%12s' % "[ DONE ]")

      # ========================================================= #
      # Define the frequency range for the FFT                    #
      # ========================================================= #

      print ("{:<133}".format('\n\tDefining frequency range for FFT ...'), end="")

      # Consider each pair of excited states

      upper_limit = 0

      for state_1 in range(1,len(system['states_list'])):
        for state_2 in range(state_1 + 1, len(system['states_list'])):

          # Compute the energy difference

          energy_diff = results_common.energy_unit_conversion(abs(system['states_list'][state_1]['energy'] - system['states_list'][state_2]['energy']),"ha","cm-1")

          # Store the transition energy if it is the minimum or maximum value so far

          if energy_diff > upper_limit:
            upper_limit = energy_diff

      lower_limit = 0
      upper_limit = 1.2 * upper_limit

      print('%12s' % "[ DONE ]")

    # ========================================================= #
    # Exception handling for the molecular information          #
    # ========================================================= #

    except results_common.ResultsError as error:
      print(error)
      print("Skipping %s molecule" % mol_name)
      continue

    # =================================================================== #
    # =================================================================== #
    #                            DATA TREATMENT                           #
    # =================================================================== #
    # =================================================================== #

    # Create a subdirectory for this molecule in the output directory

    out_mol_dir = os.path.join(out_dir, mol_group, mol_name)
    os.makedirs(out_mol_dir, exist_ok=True)

    # Iterate over each wanted profile

    for profile in ['aldu_param','const_var','filt_freq']:

      # For more information on try/except structures, see https://www.tutorialsteacher.com/python/exception-handling-in-python
      try:

        profile_dir = results_common.check_abspath(os.path.join(mol_dir,"CONTROL",profile),"Main directory for the %s profile" % profile,"directory")

        # Look for directories in the main directory (see https://stackoverflow.com/questions/800197/how-to-get-all-of-the-immediate-subdirectories-in-python for reference).

        dir_list_all = [dir.name for dir in os.scandir(profile_dir) if dir.is_dir()]

        # Only keep the 'transition_config'-type directories

        # Iterate over the directories

        for dirname in dir_list_all:
        
          # Define the 'transition_config' regex and apply it to the dirname

          trans_pattern = re.compile(r"^(?P<key>[XYZ])_(?P<init>[pST0-9-]+)_(?P<target>[pST0-9-]+)_*([a-zA-Z0-9-_]*)$")
          matching_dir = trans_pattern.match(dirname)

          # If it is a '<transition>_<config>' directory, collect the data and start generating the graphs

          if matching_dir is None:
            continue

          trans_dir = os.path.join(profile_dir, dirname)
          momdip_key = matching_dir.group('key')

          print ("{:<140}".format("\nTreating the %s transition with the '%s' profile ..." % (dirname, profile)))

          # Create a subsubdirectory for this profile in the output molecule subdirectory

          out_trans_dir = os.path.join(out_mol_dir,profile,dirname)
          os.makedirs(out_trans_dir, exist_ok=True)

          # Iterate over each pulse directory

          for sub_dir in os.listdir(trans_dir):

            # For more information on try/except structures, see https://www.tutorialsteacher.com/python/exception-handling-in-python
            try:

              if not sub_dir.startswith("pulse_") and sub_dir != 'best_pulse':
                continue

              print ("{:<133}".format("\n\tTreating the %s pulse directory ..." % (sub_dir)))

              pulse_dir = os.path.join(trans_dir, sub_dir)

              # Check key files

              iter_file = results_common.check_abspath(os.path.join(pulse_dir, "obj.res"),"Iterations QOCT-GRAD results file","file")
              guess_pulse_file = results_common.check_abspath(os.path.join(pulse_dir, "Pulse", "Pulse_init"),"Guess pulse file","file")
              pulse_file = results_common.check_abspath(os.path.join(pulse_dir, "Pulse", "Pulse_best"),"Best pulse file","file")
              pop_file = results_common.check_abspath(os.path.join(pulse_dir, "PCP_%s" % momdip_key, "pop1"),"PCP relativistic populations file","file")

              # =================================================================== #
              # =================================================================== #
              #                          PULSES TREATMENT                           #
              # =================================================================== #
              # =================================================================== #

              # Create a directory for this pulse in the profile subsubdirectory

              out_pulse_dir = os.path.join(out_trans_dir,sub_dir)
              os.makedirs(out_pulse_dir, exist_ok=True)

              # ========================================================= #
              # Temporal profiles                                         #
              # ========================================================= #

              print ("{:<126}".format('\n\t\tTreating the temporal profiles ...'), end="")

              plt.style.use('seaborn-colorblind')

              fig, (ax1,ax2) = plt.subplots(ncols=2, sharey=True, figsize=(12.8,4.8))

              # Guess pulse
              # ~~~~~~~~~~~

              # Import the pulse values

              guess_pulse = np.loadtxt(guess_pulse_file)
              time_au = guess_pulse[:,0]
              field_au = guess_pulse[:,1]
              time_step_au = time_au[1] - time_au[0]

              # Convert them from atomic units to SI units

              time = time_au * constants.value('atomic unit of time')
              guess_field = field_au * constants.value('atomic unit of electric field')
              time_step = time_step_au * constants.value('atomic unit of time')

              # Plot the temporal profile

              ax1.plot(time * 1e12,guess_field)

              # Add the legend and titles

              ax1.set_xlabel("Temps (ps)")
              ax1.set_ylabel("Champ électrique (V/m)")
              #ax_gpulse_time.set_title("Guess Pulse")

              # Best pulse
              # ~~~~~~~~~~

              # Import the pulse values (time is the same as for the guess pulse)

              pulse = np.loadtxt(pulse_file)
              field_au = pulse[:,1]

              # Convert them from atomic units to SI units

              best_field = field_au * constants.value('atomic unit of electric field')

              # Plot the temporal profile

              ax2.plot(time * 1e12,best_field)

              # Add the legend and titles

              ax2.set_xlabel("Temps (ps)")
              #ax_pulse_time.set_title("Final Pulse")

              # Save the file and close the figure

              plt.tight_layout()
              plt.savefig(os.path.join(out_pulse_dir,'temporal_pro.png'),dpi=res_dpi)
              plt.close()  

              print('%12s' % "[ DONE ]")

              # ========================================================= #
              # Spectral profiles                                         #
              # ========================================================= #

              print ("{:<126}".format('\n\t\tTreating the spectral profiles ...'), end="")

              fig, (ax1,ax2) = plt.subplots(ncols=2, sharey=True, figsize=(12.8,4.8))

              # Define the frequency range for the FFT
              # ======================================

              # if profile == 'aldu_param': # Uncomment this line and comment the next one if you want to zoom in on each specific frequency range
              if True:

                # Consider each pair of excited states

                upper_limit = 0

                for state_1 in range(1,len(system['states_list'])):
                  for state_2 in range(state_1 + 1, len(system['states_list'])):

                    # Compute the energy difference

                    energy_diff = results_common.energy_unit_conversion(abs(system['states_list'][state_1]['energy'] - system['states_list'][state_2]['energy']),"ha","cm-1")

                    # Store the transition energy if it is the maximum value so far

                    if energy_diff > upper_limit:
                      upper_limit = energy_diff

                lower_limit = 0
                upper_limit = 1.2 * upper_limit

              else:

                param_file_path = 0

                for nml_file in os.listdir(pulse_dir):
                  if nml_file.startswith("param") and nml_file.endswith('nml'):
                    param_file_path = os.path.join(pulse_dir, nml_file)
                    break
                
                if param_file_path == 0:
                  raise results_common.ResultsError ("ERROR: Can't find any parameters file in %s" % pulse_dir)

                with open(param_file_path, 'r') as param_file:
                  param_content = param_file.read().splitlines()

                window_pattern = re.compile(r'^\s+spectral_filter_fwhm\s+=\s+(?P<window>\d+\.\d+d[+-]\d+)')
                cent_freq_pattern = re.compile(r"^\s*spectral_filter_center\s+=\s+(?P<cent_freq>\d+\.\d+[dD][+-]?\d+)")

                for line in param_content:
                
                  if window_pattern.match(line): 
                    window = float((window_pattern.match(line).group('window')).replace('d','e'))

                  elif cent_freq_pattern.match(line): 
                    cent_freq = float((cent_freq_pattern.match(line).group('cent_freq')).replace('d','e'))

                upper_limit = cent_freq + window/2
                lower_limit = cent_freq - window/2

              # Guess pulse
              # ~~~~~~~~~~~

              # Compute the FFT (see tutorial at https://realpython.com/python-scipy-fft/)

              intensity = rfft(guess_field)
              freq = rfftfreq(len(time),time_step)
              freq = results_common.energy_unit_conversion(freq,'Hz','cm-1')

              # Plot the spectral profile

              ax1.plot(freq,np.abs(intensity))

              # Add the legend and titles

              ax1.set_xlim([lower_limit,upper_limit])
              ax1.set_xlabel("Nombre d'onde (cm-1)")
              ax1.set_ylabel("Intensité")

              # Set other parameters

              ax1.xaxis.set_minor_locator(AutoMinorLocator(2))
              ax1.grid(True, which='both', linestyle='--')

              # Best pulse
              # ~~~~~~~~~~

              # Compute the FFT (see tutorial at https://realpython.com/python-scipy-fft/)

              intensity = rfft(best_field)
              freq = rfftfreq(len(time),time_step)
              freq = results_common.energy_unit_conversion(freq,'Hz','cm-1')

              # Plot the spectral profile

              ax2.plot(freq,np.abs(intensity))

              # Add the legend and titles

              ax2.set_xlim([lower_limit,upper_limit])
              ax2.set_xlabel("Nombre d'onde (cm-1)")

              # Set other parameters

              ax2.xaxis.set_minor_locator(AutoMinorLocator(2))
              ax2.grid(True, which='both', linestyle='--')

              # Save the file and close the figure

              plt.tight_layout()
              plt.savefig(os.path.join(out_pulse_dir,'spectral_pro.png'),dpi=res_dpi)
              plt.close()

              # Max intensities
              # ~~~~~~~~~~~~~~~

              # Rank the maximum intensities

              list_int = list(np.abs(intensity))
              list_freq = list(freq)

              int_ranks = sorted( [(intens, idx) for (idx, intens) in enumerate(list_int)], reverse=True )

              # Fetch the 50 maximum intensities and their associated frequencies

              max_frequencies = []

              for intens, idx in int_ranks:

                max_frequencies.append({
                  "Frequency (cm$^{-1}$)" : list_freq[idx],
                  "Intensity" : intens
                })

                if len(max_frequencies) == 50:
                  break

              # Store the values in a CSV file

              csv_header = list(max_frequencies[0].keys())

              with open(os.path.join(out_pulse_dir,'max_frequencies.csv'), 'w', newline='', encoding='utf-8') as csvfile:

                csv_writer = csv.DictWriter(csvfile, fieldnames=csv_header, delimiter=';', quoting=csv.QUOTE_MINIMAL)
                csv_writer.writeheader()

                for line in max_frequencies:
                  csv_writer.writerow(line) 
      
              print('%12s' % "[ DONE ]")

              # ========================================================= #
              # Isolated best pulse (no comparison)                       #
              # ========================================================= #

              # Prepare an additionnal plot with just data about the optimal pulse

              fig_final, (ax_time,ax_freq) = plt.subplots(ncols=2, figsize=(12.8,4.8))

              # Plot the temporal profile

              ax_time.plot(time * 1e12,best_field)

              # Add the legend and titles

              ax_time.set_xlabel("Temps (ps)")
              ax_time.set_ylabel("Champ électrique (V/m)")

              # Plot the spectral profile

              ax_freq.plot(freq,np.abs(intensity))

              # Add the legend and titles

              ax_freq.set_xlim([lower_limit,upper_limit])
              ax_freq.set_xlabel("Nombre d'onde (cm-1)")
              ax_freq.set_ylabel("Intensité")

              # Set other parameters

              ax_freq.xaxis.set_minor_locator(AutoMinorLocator(2))
              ax_freq.grid(True, which='both', linestyle='--')

              # Save the file and close the figure

              plt.tight_layout()
              plt.savefig(os.path.join(out_pulse_dir,'opt_pulse.png'),dpi=res_dpi)
              plt.close()

              # =================================================================== #
              # =================================================================== #
              #                        POPULATIONS TREATMENT                        #
              # =================================================================== #
              # =================================================================== #

              print ("{:<126}".format('\n\t\tTreating the population values ...'), end="")

              # Populations graph (major)
              # =========================

              fig, (ax,axl) = plt.subplots(nrows=2, ncols=1, gridspec_kw={'height_ratios': [3, 1]})

              # Initialize the list of main states (will be used later for describing the transitions between those sates)

              main_states = []

              # Define markers and colours cycle for Matplotlib's Pyplot  
              # See https://stackoverflow.com/questions/13091649/unique-plot-marker-for-each-plot-in-matplotlib

              line_cycler = (cycler(linestyle=[(0, ()),(0, (5, 3)), (0, (1, 1))])
              * cycler(color=['blue','red','green','purple','black','gold','darkturquoise','peru','magenta','lightgreen']))

              # Load the populations file

              pop_file_content = np.loadtxt(pop_file)

              # Import the time values and convert them to seconds
              
              time_au = pop_file_content[:,0]
              time = time_au * constants.value('atomic unit of time')

              # Iterate over each state and plot the populations

              ax.set_prop_cycle(line_cycler)
              axl.set_prop_cycle(line_cycler)
              current_column = 0

              for state in system['states_list']:
                current_column += 1
                if not math.isclose(max(pop_file_content[:,current_column]),0, abs_tol=5e-2):
                  multip = "pT" if state['trip_percent'] >= 0.5 else "pS"
                  main_states.append((state['number'],state['label'],multip))
                  ax.plot(time * 1e12, pop_file_content[:,current_column], label=state['label'], linewidth=0.75)
                  axl.plot(0,0,label=state['label'])

              # Add the legend and titles

              ax.set_xlabel("Temps (ps)")
              ax.set_ylabel("Population")

              axl.axis("off")
              axl.legend(ncol=5, loc="center")

              # Set other parameters

              #ax.xaxis.set_minor_locator(AutoMinorLocator(2))
              ax.yaxis.set_minor_locator(AutoMinorLocator(2))

              ax.grid(True, which='both', linestyle='--')
              plt.tight_layout()

              # Save the file and close the figure

              plt.savefig(os.path.join(out_pulse_dir,'pcp_pop.png'), dpi=res_dpi, bbox_inches="tight")
              plt.close()

              # Populations graph (minor)
              # =========================

              fig, (ax,axl) = plt.subplots(nrows=2, ncols=1, gridspec_kw={'height_ratios': [3, 1]})

              # Iterate over each state and plot the populations

              ax.set_prop_cycle(line_cycler)
              axl.set_prop_cycle(line_cycler)
              current_column = 0

              for state in system['states_list']:
                current_column += 1
                if math.isclose(max(pop_file_content[:,current_column]),0, abs_tol=5e-2):
                  ax.plot(time * 1e12, pop_file_content[:,current_column], label=state['label'], linewidth=0.75)
                  axl.plot(0,0,label=state['label'])

              # Add the legend and titles

              ax.set_xlabel("Temps (ps)")
              ax.set_ylabel("Population")

              axl.axis("off")
              axl.legend(ncol=5, loc="center")

              # Set other parameters

              #ax.xaxis.set_minor_locator(AutoMinorLocator(2))
              ax.yaxis.set_minor_locator(AutoMinorLocator(2))

              ax.grid(True, which='both', linestyle='--')
              plt.tight_layout()

              # Save the file and close the figure

              plt.savefig(os.path.join(out_pulse_dir,'pcp_pop_minor.png'), dpi=res_dpi, bbox_inches="tight")
              plt.close()

              # Transition between the main involved states
              # ===========================================

              # Initialize the list of dictionaries

              main_transitions = []

              # Fetch the transition frequency and dipole moments between the main states
              
              for state_1 in main_states:
                for state_2 in main_states:

                  if state_1[0] > state_2[0]:
                    continue

                  energy_1 = results_common.energy_unit_conversion(system['states_list'][state_1[0]]['energy'],"ha","cm-1")
                  energy_2 = results_common.energy_unit_conversion(system['states_list'][state_2[0]]['energy'],"ha","cm-1")

                  if control_common.is_indistinguishable(energy_1, energy_2):
                    continue
                  
                  main_transitions.append({
                    "Transition" : state_1[1] + "-" + state_2[1],
                    "Fréquence (cm$^{-1}$)" : abs(energy_1 - energy_2),
                    "$\mu_%s$ (u.a.)" % momdip_key : system['momdip_mtx'][momdip_key][state_1[0]][state_2[0]],
                    "Type" : "$%s \to %s$" % (state_1[2],state_2[2]) if state_1[2] == state_2[2] else "$pS \to pT$"
                  })

              # Store the values in a CSV file
              # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

              if main_transitions != []:

                csv_header = list(main_transitions[0].keys())

                with open(os.path.join(out_pulse_dir,'main_transitions.csv'), 'w', newline='', encoding='utf-8') as csvfile:

                  csv_writer = csv.DictWriter(csvfile, fieldnames=csv_header, delimiter=';', quoting=csv.QUOTE_MINIMAL)
                  csv_writer.writeheader()

                  for trans in main_transitions:
                    csv_writer.writerow(trans) 

              # Create the LaTeX table
              # ======================

              if main_transitions != []:

                df = pd.DataFrame(main_transitions)
                with open(os.path.join(out_pulse_dir,"main_transitions.tex"), 'w+', encoding='utf-8') as f:
                  f.write(df.to_latex(
                    index=False,
                    formatters=[None, results_common.format_num(2, "f"), results_common.format_num(2, "e"), None],
                    column_format="ccccc",
                    escape=False,
                    na_rep="-"))

              print('%12s' % "[ DONE ]")

              # =================================================================== #
              # =================================================================== #
              #                        ITERATIONS TREATMENT                         #
              # =================================================================== #
              # =================================================================== #

              print ("{:<126}".format('\n\t\tTreating the iterations values ...'), end="")

              # Get the values
              # ==============

              # Load the iterations file

              with open(iter_file, 'r') as file:
                iter_content = file.read().splitlines()

              # Define the expression patterns for the lines of the iterations file
              # For example "      0     1  1sec |Proba_moy  0.693654D-04 |Fidelity(U)  0.912611D-01 |Chp  0.531396D-04 -0.531399D-04 |Aire -0.202724D-03 |Fluence  0.119552D-03 |Recou(i)  0.693654D-04 |Tr_dist(i) -0.384547D-15 |Tr(rho)(i)  0.100000D+01 |Tr(rho^2)(i)  0.983481D+00 |Projector  0.100000D+01"          

              rx_iter_line = re.compile(r"^\s+(?P<niter>\d+)\s+\d+\s+\d+(?:sec|min)\s\|Proba_moy\s+\d\.\d+D[+-]\d+\s\|Fidelity\(U\)\s+(?P<fidelity>\d\.\d+D[+-]\d+)\s\|Chp\s+\d\.\d+D[+-]\d+\s+-?\d\.\d+D[+-]\d+\s\|Aire\s+-?\d\.\d+D[+-]\d+\s\|Fluence\s+(?P<fluence>\d\.\d+D[+-]\d+)\s\|Recou\(i\)\s+(?P<overlap>\d\.\d+D[+-]\d+)\s\|Tr_dist\(i\)\s+-?\d\.\d+D[+-]\d+\s\|Tr\(rho\)\(i\)\s+\d\.\d+D[+-]\d+\s\|Tr\(rho\^2\)\(i\)\s+\d\.\d+D[+-]\d+\s\|Projector\s+(?P<projector>\d\.\d+D[+-]\d+)")
              
              # Extract the number of iterations and the fidelity for each line of the iterations file

              iterations = []
              projectors = []
              fluences = []

              for line in iter_content:

                iter_data = rx_iter_line.match(line)

                if iter_data is not None:

                  niter = int(iter_data.group("niter"))
                  iterations.append(niter)

                  projector_raw = iter_data.group("projector")
                  projector = float(re.compile(r'(\d*\.\d*)[dD]([-+]?\d+)').sub(r'\1E\2', projector_raw)) # Replace the possible d/D from Fortran double precision float format with an "E", understandable by Python)
                  projectors.append(projector) 

                  fluence_raw = iter_data.group("fluence")
                  fluence = float(re.compile(r'(\d*\.\d*)[dD]([-+]?\d+)').sub(r'\1E\2', fluence_raw)) # Replace the possible d/D from Fortran double precision float format with an "E", understandable by Python)
                  fluence = results_common.energy_unit_conversion(fluence,'ha','j')/(constants.value('atomic unit of length')**2)
                  fluences.append(fluence) 

              # Plot performance
              # ================

              # Plot the evolution of performance over iterations

              fig, ax = plt.subplots()
              ax.plot(iterations,projectors)

              # Add the legend and titles

              ax.set_xlabel("Nombre d'itérations")
              ax.set_ylabel("Indice de performance")

              # Save the file and close the figure

              plt.tight_layout()
              plt.savefig(os.path.join(out_pulse_dir,'perf_iter.png'),dpi=res_dpi)
              plt.close()

              # Plot fluence
              # ============

              # Plot the evolution of fluence over iterations

              fig, ax = plt.subplots()
              ax.plot(iterations,fluences)

              # Add the legend and titles

              ax.set_xlabel("Nombre d'itérations")
              ax.set_ylabel(r"Fluence (J/m$^2$)")

              # Save the file and close the figure

              params = {'mathtext.default': 'regular' }          
              plt.rcParams.update(params)

              plt.tight_layout()
              plt.savefig(os.path.join(out_pulse_dir,'flu_iter.png'),dpi=res_dpi)
              plt.close()

              print('%12s' % "[ DONE ]")

            # ========================================================= #
            # Exception handling for the pulse treatment                #
            # ========================================================= #

            except results_common.ResultsError as error:
              print(error)
              print("Skipping %s directory" % sub_dir)
              continue

      # ========================================================= #
      # Exception handling for the data treatment                 #
      # ========================================================= #

      except results_common.ResultsError as error:
        print(error)
        print("Skipping %s profile" % profile)
        continue

    console_message = "End of procedure for the molecule " + mol_name
    print("")
    print(''.center(len(console_message)+10, '*'))
    print(console_message.center(len(console_message)+10))
    print(''.center(len(console_message)+10, '*'))

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
