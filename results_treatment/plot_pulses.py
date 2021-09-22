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
import csv
import os
import re
import shutil
from inspect import getsourcefile

import matplotlib.mathtext
import matplotlib.pyplot as plt
import numpy as np
import yaml
from matplotlib.ticker import AutoMinorLocator
from scipy.fft import rfft, rfftfreq
from scipy import constants

import results_common

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
optional.add_argument('-qt', '--threshold', type=float, help="Quality threshold for the fidelity below which the graphs will not be created.")

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
    threshold = args.threshold               # Quality threshold below which graphs won't be created

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

    if threshold and (threshold < 0 or threshold > 1):
      raise results_common.ResultsError ("The value for the quality threshold (%s) must be comprised between 0 and 1." % threshold)

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

      print ("{:<140}".format('\nFetching molecular information ...'), end="")

      # Check the data directory

      data_dir = results_common.check_abspath(os.path.join(mol_dir,"CONTROL","data"),"Data directory created by control_launcher.py","directory")

      # Load the eigenstates list

      eigenstates_file = results_common.check_abspath(os.path.join(data_dir, "eigenstates.csv"),"Eigenstates file","file")

      with open(eigenstates_file, 'r', newline='') as csv_file:

        eigenstates_content = csv.DictReader(csv_file, delimiter=';')
        eigenstates_list = list(eigenstates_content)
        eigenstates_header = eigenstates_content.fieldnames

      # Check the eigenstates list

      required_keys = ['Number','Label']
      results_common.check_keys(required_keys,eigenstates_list,"Eigenstates list file at %s" % eigenstates_file)

      # Load the transitions list

      transitions_file = results_common.check_abspath(os.path.join(data_dir, "transitions.csv"),"Transitions file","file")

      with open(transitions_file, 'r', newline='') as csv_file:

        transitions_content = csv.DictReader(csv_file, delimiter=';')
        transitions_list = list(transitions_content)
        transitions_header = transitions_content.fieldnames

      # Check the transitions list

      required_keys = ['Label', 'Energy (Ha)', 'Initial state number', 'Target state number', 'Transition dipole moments matrix']
      results_common.check_keys(required_keys,transitions_list,"Transitions list file at %s" % transitions_file)

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

    # Look for directories in the results CONTROL directory (see https://stackoverflow.com/questions/800197/how-to-get-all-of-the-immediate-subdirectories-in-python for reference).

    control_dir = os.path.join(mol_dir,"CONTROL")

    dir_list_all = [dir.name for dir in os.scandir(control_dir) if dir.is_dir()]

    # Only keep the 'transition_config'-type directories

    dir_list = []

    for dirname in dir_list_all:
      
      # For more information on try/except structures, see https://www.tutorialsteacher.com/python/exception-handling-in-python
      try:

        for transition in transitions_list:

          # Define the 'transition_config' regex and apply it to the dirname

          pattern = re.compile(r"^" + re.escape(transition["Label"]) + r"_(?P<config>.*)$")
          matching_dir = pattern.match(dirname)

          # If it is a '<transition>_<config>' directory, collect the data and start generating the graphs

          if matching_dir is not None:

            dir_list.append(dirname)
            config = pattern.match(dirname).group('config')

            print ("{:<140}".format("\nTreating the %s transition with the '%s' config ..." % (transition["Label"], config)))

            # Check key files

            iter_file = results_common.check_abspath(os.path.join(control_dir, dirname, "obj.res"),"Iterations QOCT-GRAD results file","file")
            guess_pulse_file = results_common.check_abspath(os.path.join(control_dir, dirname, "Pulse", "Pulse_init"),"Guess pulse file","file")
            pulse_file = results_common.check_abspath(os.path.join(control_dir, dirname, "Pulse", "Pulse"),"Final pulse file","file")
            pop_file = results_common.check_abspath(os.path.join(control_dir, dirname, "PCP", "pop1"),"PCP eigenstates populations file","file")

            # =================================================================== #
            # =================================================================== #
            #                       QUALITY THRESHOLD CHECK                       #
            # =================================================================== #
            # =================================================================== #

            if threshold:

              print ("{:<133}".format('\n\tChecking quality ...'), end="")

              # Go straight to the last line of the iterations file (see https://stackoverflow.com/questions/46258499/read-the-last-line-of-a-file-in-python for reference)

              with open(iter_file, 'rb') as f:
                f.seek(-2, os.SEEK_END)
                while f.read(1) != b'\n':
                    f.seek(-2, os.SEEK_CUR)
                last_line = f.readline().decode()

              # Define the expression patterns for the lines of the iterations file
              # For example "    300     2  2sec |Proba_moy  0.000000E+00 |Fidelity(U)  0.000000E+00 |Chp  0.123802E+00 -0.119953E+00 |Aire  0.140871E-03 |Fluence  0.530022E+01 |Recou(i)  0.000000E+00 |Tr_dist(i) -0.500000E+00 |Tr(rho)(i)  0.100000E+01 |Tr(rho^2)(i)  0.100000E+01 |Projector  0.479527E-13"

              rx_iter_line = re.compile(r"^\s+(?P<niter>\d+)\s+\d+\s+\d+sec\s\|Proba_moy\s+\d\.\d+D[+-]\d+\s\|Fidelity\(U\)\s+(?P<fidelity>\d\.\d+D[+-]\d+)\s\|Chp\s+\d\.\d+D[+-]\d+\s+-?\d\.\d+D[+-]\d+\s\|Aire\s+-?\d\.\d+D[+-]\d+\s\|Fluence\s+\d\.\d+D[+-]\d+\s\|Recou\(i\)\s+\d\.\d+D[+-]\d+\s\|Tr_dist\(i\)\s+-?\d\.\d+D[+-]\d+\s\|Tr\(rho\)\(i\)\s+\d\.\d+D[+-]\d+\s\|Tr\(rho\^2\)\(i\)\s+\d\.\d+D[+-]\d+\s\|Projector\s+\d\.\d+D[+-]\d+")

              # Get the last fidelity from the last line

              iter_data = rx_iter_line.match(last_line)
              if iter_data is not None:
                fidelity_raw = iter_data.group("fidelity")
                fidelity = float(re.compile(r'(\d*\.\d*)[dD]([-+]?\d+)').sub(r'\1E\2', fidelity_raw)) # Replace the possible d/D from Fortran double precision float format with an "E", understandable by Python)
              else:
                raise results_common.ResultsError ("ERROR: Unable to get fidelity from the last line of %s" % iter_file) 

              # Compare the fidelity to the quality threshold

              if fidelity < threshold:
                print("The fidelity for this pulse (%s) is inferior to the specified quality threshold (%s)" % (fidelity,threshold))
                print("Skipping %s directory" % dirname)
                continue

              print('%12s' % "[ DONE ]")

            # =================================================================== #
            # =================================================================== #
            #                        GUESS PULSE TREATMENT                        #
            # =================================================================== #
            # =================================================================== #

            # Define the figure that will host the six graphs for this specific '<transition>_<config>' directory

            fig, ((ax_gpulse_time,ax_gpulse_freq),(ax_pulse_time,ax_pulse_freq),(ax_pop,ax_fidel)) = plt.subplots(nrows=3,ncols=2)

            print ("{:<133}".format('\n\tTreating the guess pulse values ...'), end="")

            # Import the pulse values

            guess_pulse = np.loadtxt(guess_pulse_file)
            time_au = guess_pulse[:,0]
            amplitude_au = guess_pulse[:,1]
            time_step_au = time_au[1] - time_au[0]

            # Convert them from atomic units to SI units

            time = time_au * constants.value('atomic unit of time')
            amplitude = amplitude_au * constants.value('atomic unit of electric field')
            time_step = time_step_au * constants.value('atomic unit of time')

            # Plot the temporal profile

            ax_gpulse_time.plot(time * 1e12,amplitude)
            ax_gpulse_freq.set_yscale('log')
            ax_gpulse_time.set_xlabel("Time (ps)")
            ax_gpulse_time.set_ylabel("Amplitude (V/m)")

            # Compute the FFT (see tutorial at https://realpython.com/python-scipy-fft/)

            intensity = rfft(amplitude)
            freq = rfftfreq(len(time),time_step)

            # Plot the spectral profile

            ax_gpulse_freq.plot(freq,np.abs(intensity))
            ax_gpulse_freq.set_xlabel("Frequencies (Hz)")
            ax_gpulse_freq.set_ylabel("Intensity")

            print('%12s' % "[ DONE ]")

            # =================================================================== #
            # =================================================================== #
            #                        FINAL PULSE TREATMENT                        #
            # =================================================================== #
            # =================================================================== #

            print ("{:<133}".format('\n\tTreating the final pulse values ...'), end="")

            # Import the pulse values

            pulse = np.loadtxt(pulse_file)
            time_au = pulse[:,0]
            amplitude_au = pulse[:,1]
            time_step_au = time_au[1] - time_au[0]

            # Convert them from atomic units to SI units

            time = time_au * constants.value('atomic unit of time')
            amplitude = amplitude_au * constants.value('atomic unit of electric field')
            time_step = time_step_au * constants.value('atomic unit of time')

            # Plot the temporal profile

            ax_pulse_time.plot(time * 1e12,amplitude)
            ax_pulse_time.set_xlabel("Time (ps)")
            ax_pulse_time.set_ylabel("Amplitude (V/m)")

            # Compute the FFT (see tutorial at https://realpython.com/python-scipy-fft/)

            intensity = rfft(amplitude)
            freq = rfftfreq(len(time),time_step)

            # Plot the spectral profile

            ax_pulse_freq.plot(freq,np.abs(intensity))
            ax_pulse_freq.set_yscale('log')
            ax_pulse_freq.set_xlabel("Frequencies (Hz)")
            ax_pulse_freq.set_ylabel("Intensity")

            print('%12s' % "[ DONE ]")

            # =================================================================== #
            # =================================================================== #
            #                        POPULATIONS TREATMENT                        #
            # =================================================================== #
            # =================================================================== #

            print ("{:<133}".format('\n\tTreating the post-controle values ...'), end="")

            # Load the populations file

            pop_file_content = np.loadtxt(pop_file)

            # Import the time values and convert them to seconds
            
            time_au = pop_file_content[:,0]
            time = time_au * constants.value('atomic unit of time')

            # Iterate over each state and plot the populations

            current_column = 0

            for eigenstate in eigenstates_list:
              current_column += 1
              ax_pop.plot(time * 1e12,pop_file_content[:,current_column],label=eigenstate['Label'])

            # Legend the graph

            ax_pop.set_xlabel("Time (ps)")
            ax_pop.set_ylabel("Population")
            ax_pop.legend()

            print('%12s' % "[ DONE ]")

            # =================================================================== #
            # =================================================================== #
            #                        FIDELITIES TREATMENT                         #
            # =================================================================== #
            # =================================================================== #

            print ("{:<133}".format('\n\tTreating the fidelities values ...'), end="")

            # Load the iterations file

            with open(iter_file, 'r') as file:
              iter_content = file.read().splitlines()

            # Define the expression patterns for the lines of the iterations file
            # For example "    300     2  2sec |Proba_moy  0.000000E+00 |Fidelity(U)  0.000000E+00 |Chp  0.123802E+00 -0.119953E+00 |Aire  0.140871E-03 |Fluence  0.530022E+01 |Recou(i)  0.000000E+00 |Tr_dist(i) -0.500000E+00 |Tr(rho)(i)  0.100000E+01 |Tr(rho^2)(i)  0.100000E+01 |Projector  0.479527E-13            

            rx_iter_line = re.compile(r"^\s+(?P<niter>\d+)\s+\d+\s+\d+sec\s\|Proba_moy\s+\d\.\d+D[+-]\d+\s\|Fidelity\(U\)\s+(?P<fidelity>\d\.\d+D[+-]\d+)\s\|Chp\s+\d\.\d+D[+-]\d+\s+-?\d\.\d+D[+-]\d+\s\|Aire\s+-?\d\.\d+D[+-]\d+\s\|Fluence\s+\d\.\d+D[+-]\d+\s\|Recou\(i\)\s+\d\.\d+D[+-]\d+\s\|Tr_dist\(i\)\s+-?\d\.\d+D[+-]\d+\s\|Tr\(rho\)\(i\)\s+\d\.\d+D[+-]\d+\s\|Tr\(rho\^2\)\(i\)\s+\d\.\d+D[+-]\d+\s\|Projector\s+\d\.\d+D[+-]\d+")

            # Extract the number of iterations and the fidelity for each line of the iterations file

            iterations = []
            fidelities = []

            for line in iter_content:

              iter_data = rx_iter_line.match(line)

              if iter_data is not None:
                niter = int(iter_data.group("niter"))
                iterations.append(niter)
                fidelity_raw = iter_data.group("fidelity")
                fidelity = float(re.compile(r'(\d*\.\d*)[dD]([-+]?\d+)').sub(r'\1E\2', fidelity_raw)) # Replace the possible d/D from Fortran double precision float format with an "E", understandable by Python)
                fidelities.append(fidelity) 

            # Plot the evolution of fidelity over iterations

            ax_fidel.plot(iterations,fidelities)
            ax_fidel.set_xlabel("Number of iterations")
            ax_fidel.set_ylabel("Fidelity")

            print('%12s' % "[ DONE ]")

            # =================================================================== #
            # =================================================================== #
            #                       CREATE THE FIGURE FILE                        #
            # =================================================================== #
            # =================================================================== #

            fig.set_size_inches(8, 10)
            plt.tight_layout()
            plt.savefig(os.path.join(out_dir,'%s_%s.png' % (mol_name,dirname)),dpi=200)
            plt.close()

      # ========================================================= #
      # Exception handling for the pulse treatment                #
      # ========================================================= #

      except results_common.ResultsError as error:
        print(error)
        print("Skipping %s directory" % dirname)
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