#!/usr/bin/env python3

########################################################################################################################################################
##                                                        QOCT-RA Input Builder & Job Launcher                                                        ##
##                                                                                                                                                    ##
##                                    This script prepares the input files needed to run the QOCT-RA program by extracting                            ##
##                                  information from a given source file and launches the corresponding jobs on the cluster.                          ##
##                                      Extended documentation is available at https://chains-ulb.readthedocs.io/                                     ##
##                                                                                                                                                    ##
## /!\ In order to run, this script requires Python 3.5+ as well as NumPy 1.14+, YAML and Jinja2. Ask your cluster(s) administrator(s) if needed. /!\ ##
########################################################################################################################################################

import argparse
import copy
import importlib
import os
import shutil
import sys
from collections import OrderedDict
from inspect import getsourcefile

import jinja2  # Only needed in the renderer subscript, it is loaded here to check if your python installation does support jinja2
import numpy as np
import yaml

# Subscripts (files that end with .py and must be placed in the same directory as this script)

import control_errors
import control_renderer
import source_parser
import transition_fcts

# =================================================================== #
# =================================================================== #
#                       COMMAND LINE ARGUMENTS                        #
# =================================================================== #
# =================================================================== #

# Define the arguments needed for the script (here they are defined as named arguments rather than positional arguments, check https://stackoverflow.com/questions/24180527/argparse-required-arguments-listed-under-optional-arguments for more info).

parser = argparse.ArgumentParser(add_help=False, description="This script prepares the input files needed to run the QOCT-RA program by extracting information from a given source file and launches the corresponding jobs on the cluster.")

required = parser.add_argument_group('Required arguments')
required.add_argument("-s","--source", type=str, help="Path to the source file containing all the necessary information that needs to be processed.", required=True)
required.add_argument("-p","--parsing_fct", type=str, help="Name of the parsing function that will be used to parse the source file, as defined in source_parser.py.", required=True)
required.add_argument('-cf', '--config', type=str, help="Path to the YAML configuration file.", required=True)
required.add_argument("-o","--out_dir", type=str, help="Path to the directory where you want to create the subdirectories for each job.", required=True)
required.add_argument('-cl', '--cluster_name', type=str, help="Name of the cluster where this script is running, as defined in the YAML clusters configuration file.", required=True)

optional = parser.add_argument_group('Optional arguments')
optional.add_argument('-h','--help',action='help',default=argparse.SUPPRESS,help='Show this help message and exit')
optional.add_argument("-ow","--overwrite",action="store_true",help="If a job subdirectory for the same source file already exists, remove it before creating a new one.")
optional.add_argument("-d","--dry_run",action="store_true",help="Do not launch the jobs, just create the files and directories.")

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
    print("EXECUTION OF THE QOCT-RA INPUT BUILDER & JOB LAUNCHER BEGINS NOW".center(columns))
    print("")
    print("".center(columns,"*"))

    # ========================================================= #
    # Read command line arguments                               #
    # ========================================================= #

    args = parser.parse_args()

    # Required arguments

    source = args.source                     # Source file containing all the necessary information
    parsing_fct = args.parsing_fct           # Name of the parsing function, as defined in source_parser.py
    out_dir = args.out_dir                   # Directory where all jobs subdirectories will be created
    config_file = args.config                # YAML configuration file
    cluster_name = args.cluster_name         # Name of the cluster where this script is running, as defined in the clusters configuration YAML file

    # Optional arguments

    overwrite = args.overwrite               # Flag for removing job subdirectories before creating a new one, if they have the same name
    dry_run = args.dry_run                   # Flag to not launch the jobs and just create the files

    # Other important variable

    profile = "qoctra"                         # Name of the qoctra key that appears in the clusters configuration YAML files

    # ========================================================= #
    # Define codes directory                                    #
    # ========================================================= #

    # Determined by getting the path to the directory of this script

    code_dir = os.path.dirname(os.path.realpath(os.path.abspath(getsourcefile(lambda:0))))

    print ("{:<40} {:<100}".format('\nCodes directory:',code_dir))

    # ========================================================= #
    # Check and load the YAML clusters configuration file       #
    # ========================================================= #

    # Loading the clusters_file 

    clusters_file = control_errors.check_abspath(os.path.join(code_dir,"clusters.yml"),"YAML clusters configuration file","file")
    print ("{:<80}".format('\nLoading the clusters configuration file "clusters.yml" ...'), end="")
    with open(clusters_file, 'r') as f_clusters:
      clusters_cfg = yaml.load(f_clusters, Loader=yaml.FullLoader)
    print('%12s' % "[ DONE ]")

    # Check the name of the cluster

    if cluster_name not in clusters_cfg:
      raise control_errors.ControlError ("ERROR: There is no information about the %s cluster in the clusters configuration file. Please add relevant information or change the cluster before proceeding further." % cluster_name.upper())

    print("\nThis script is running on the %s cluster" % cluster_name.upper())

    # Check if the submit_command key has been defined

    submit_command = clusters_cfg[cluster_name].get("submit_command")

    if submit_command is None:
      raise control_errors.ControlError ("ERROR: There is no defined submit_command for the %s cluster in the clusters configuration file." % cluster_name.upper()) 

    # Check if the profile exists 

    if "profiles" not in clusters_cfg[cluster_name]:
      raise control_errors.ControlError ('ERROR: There is no "profiles" key defined for the %s cluster in the clusters configuration file. Consult official documentation for details.' % cluster_name.upper())    

    if profile not in clusters_cfg[cluster_name]["profiles"]:
      raise control_errors.ControlError ('ERROR: There is no "%s" key defined for the %s cluster in the clusters configuration file. \nPlease add information for this profile to the clusters configuration file.' % (profile, cluster_name.upper()))

    # ========================================================= #
    # Establishing the different job scales                     #
    # ========================================================= #

    # Gather all the different job scales from the clusters configuration file in a temporary dictionary

    job_scales_tmp = clusters_cfg[cluster_name]['profiles'][profile].get('job_scales')

    if job_scales_tmp is None:
      raise control_errors.ControlError ("ERROR: There is no defined job_scales for the %s profile in the %s cluster in the clusters configuration file." % (profile, cluster_name.upper())) 

    # Defined the required keys in our job scales

    required_keys = frozenset({"label", "scale_limit", "time", "mem_per_cpu" })

    # Define an array for correct English spelling during printing

    special_numbers = {1:"st", 2:"nd", 3:"rd"}

    # Initialize the final dictionary where the job scales will be sorted by their upper limit

    job_scales = {}

    # Check the job scales

    for scale in job_scales_tmp:

      # Check if all the required keys are present

      for key in required_keys:
        if key not in scale:
          raise control_errors.ControlError ('ERROR: There is no defined "%s" key for the %s%s job scale of the %s profile in the %s cluster in the clusters configuration file.' % (key, job_scales_tmp.index(scale), ("th" if not job_scales_tmp.index(scale) in special_numbers else special_numbers[job_scales_tmp.index(scale)]), profile, cluster_name.upper()))           

      # Extract the scale upper limit from the job scales

      scale_limit = scale.pop('scale_limit')
      job_scales[scale_limit] = scale

    # Sort the different job scales by their upper limit and store them in the job_scales dictionary

    job_scales = OrderedDict(sorted(job_scales.items()))

    print("\nJob scales for %s on %s:" % (profile,cluster_name.upper()))
    print("")
    print(''.center(146, '-'))
    print ("{:<15} {:<20} {:<20} {:<20} {:<10} {:<20} {:<40}".format('Scale Limit','Label','Partition Name','Time','Cores','Mem per CPU (MB)','Delay Command'))
    print(''.center(146, '-'))
    for scale_limit, scale in job_scales.items():
      print ("{:<15} {:<20} {:<20} {:<20} {:<10} {:<20} {:<40}".format(scale_limit, scale['label'], scale.get('partition_name', "not specified"), scale['time'], scale.get('cores', "default"), scale['mem_per_cpu'], scale.get('delay_command', "not specified")))
    print(''.center(146, '-'))

    # ========================================================= #
    # Check other arguments                                     #
    # ========================================================= #

    out_dir = control_errors.check_abspath(out_dir,"Command line argument -o / --out_dir","directory")
    print ("{:<40} {:<100}".format('\nJobs main directory:',out_dir))

    source = control_errors.check_abspath(source,"Command line argument -s / --source","file")
    print ("{:<40} {:<100}".format('\nSource file:',source))

    if (parsing_fct) not in dir(source_parser) or not callable(getattr(source_parser, parsing_fct)):
      raise control_errors.ControlError ("ERROR: There is no parsing function named %s defined in source_parser.py." % parsing_fct)
    print ("{:<40} {:<100}".format('\nParsing function:',parsing_fct))

    # Check and load the configuration file for the information about the molecule

    config_file = control_errors.check_abspath(config_file,"Command line argument -cf / --config","file")
    print ("{:<40} {:<100}".format('\nConfiguration file:',config_file))
  
    print ("{:<80}".format('\nLoading the configuration file ...'), end='')
    with open(config_file, 'r') as f_config:
      config = yaml.load(f_config, Loader=yaml.FullLoader)
    print('%12s' % "[ DONE ]")

    # ========================================================= #
    # Important files and directories                           #
    # ========================================================= #

    # Get the name of the source file and the name of the directory where the source file is

    source_path = os.path.dirname(source)
    source_filename = os.path.basename(source)

    # Check if a directory already exists for that molecule

    mol_name = str(source_filename.split('.')[0]) # Getting rid of the format extension to get the name of the molecule

    if os.path.exists(os.path.join(out_dir,mol_name)) and not overwrite:
      raise control_errors.ControlError ("ERROR: A directory for the %s molecule (or source file) already exists in %s !" % (mol_name, out_dir))

    # Define the rendering function that will render the job script and the parameters file(s) (defined in control_renderer.py)

    render_fct = profile + "_render"

    if (render_fct) not in dir(control_renderer) or not callable(getattr(control_renderer, render_fct)):
      raise control_errors.ControlError ("ERROR: There is no function defined for the %s profile in renderer.py." % profile)

  # ========================================================= #
  # Exception handling for the preparation step               #
  # ========================================================= #

  except control_errors.ControlError as error:
    print("")
    print(error)
    exit(-1)

  # =================================================================== #
  # =================================================================== #
  #                       PARSING THE SOURCE FILE                       #
  # =================================================================== #
  # =================================================================== #

  # For more information on try/except structures, see https://www.tutorialsteacher.com/python/exception-handling-in-python
  try:

    section_title = "1. Parsing the source file"

    print("")
    print(''.center(len(section_title)+10, '*'))
    print(section_title.center(len(section_title)+10))
    print(''.center(len(section_title)+10, '*'))

    # ========================================================= #
    # Read the file                                             #
    # ========================================================= #

    print ("{:<80}".format('\nLoading %s file ... ' % source_filename), end="")
    with open(source, 'r') as source_file:
      source_content = source_file.read().splitlines()
    print('%12s' % "[ DONE ]")

    # Cleaning up the source file from surrounding spaces and blank lines

    source_content = list(map(str.strip, source_content))   # Remove leading & trailing blank/spaces
    source_content = list(filter(None, source_content))     # Remove blank lines/no char

    # Call the parsing function (defined in source_parser.py, see the documentation for more information)

    system = eval("source_parser." + parsing_fct)(source_content)

    print("\nThe source file has been succesfully parsed.")

    # ========================================================= #
    # MIME diagonalization                                      #
    # ========================================================= #

    print("{:<80}".format("\nDiagonalizing the MIME ..."), end="")
    # Using NumPy to diagonalize the matrix (see https://numpy.org/doc/stable/reference/generated/numpy.linalg.eig.html for reference)   
    system['eigenvalues'], system['eigenvectors'] = np.linalg.eig(system['mime'])
    print('%12s' % "[ DONE ]")

    # Sort the eigenvalues and associated eigenvectors (see https://stackoverflow.com/questions/8092920/sort-eigenvalues-and-associated-eigenvectors-after-using-numpy-linalg-eig-in-pyt for reference)

    idx = system['eigenvalues'].argsort()   
    system['eigenvalues'] = system['eigenvalues'][idx]
    system['eigenvectors'] = system['eigenvectors'][:,idx]

    # ========================================================= #
    # Eigenvalues                                               #
    # ========================================================= #

    # Converting the eigenvalues from cm-1 to ua, nm and eV
    eigenvalues_ua = system['eigenvalues'] / 219474.6313705
    eigenvalues_nm = 10000000 / system['eigenvalues']
    eigenvalues_ev = system['eigenvalues'] / 8065.6

    print("")
    print(''.center(40, '-'))
    print('Eigenvalues'.center(40, ' '))
    print(''.center(40, '-'))
    print("{:<10} {:<10} {:<10} {:<10}".format('cm-1','ua','eV','nm'))
    print(''.center(40, '-'))
    for val in range(len(system['eigenvalues'])):
      print("{:<9.2f} {:<1.4e} {:<8.4f} {:<8.4f}".format(system['eigenvalues'][val],eigenvalues_ua[val],eigenvalues_ev[val],eigenvalues_nm[val]))
    print(''.center(40, '-'))

    # ========================================================= #
    # Eigenvectors matrix                                       #
    # ========================================================= #

    print("\nEigenvectors matrix")
    print('')
    for vector in system['eigenvectors']:
      for val in vector:
        print(np.format_float_scientific(val,precision=5,unique=False,pad_left=2), end = " ")
      print('')

    # ========================================================= #
    # Eigenvectors transpose matrix                             #
    # ========================================================= #

    # Using NumPy to transpose the eigenvectors matrix (see https://numpy.org/doc/stable/reference/generated/numpy.transpose.html for reference)
    system['transpose'] = np.transpose(system['eigenvectors'])

    print("\nEigenvectors transpose matrix")
    print('')
    for vector in system['transpose']:
      for val in vector:
        print(np.format_float_scientific(val,precision=5,unique=False,pad_left=2), end = " ")
      print('')

    # ========================================================= #
    # Diagonalized MIME                                         #
    # ========================================================= #

    # Using NumPy to convert from the zero order basis set to the eigenstates basis set through a matrix product (see https://numpy.org/doc/stable/reference/generated/numpy.matmul.html#numpy.matmul for reference)
    system['mime_diag'] = np.matmul(np.matmul(system['transpose'],system['mime']),system['eigenvectors'])
        
    print("\nMIME in the eigenstates basis set (cm-1)")
    print('')
    for row in system['mime_diag']:
      for val in row:
        print(np.format_float_scientific(val,precision=5,unique=False,pad_left=2), end = " ")
      print('')

    # ========================================================= #
    # Dipole moment matrix in the eigenstates basis set         #
    # ========================================================= #

    # Using NumPy to convert from the zero order basis set to the eigenstates basis set through a matrix product (see https://numpy.org/doc/stable/reference/generated/numpy.matmul.html#numpy.matmul for reference)
    system['momdip_es_mtx'] = np.matmul(np.matmul(system['transpose'],system['momdip_mtx']),system['eigenvectors'])
        
    print("\nDipole moments matrix in the eigenstates basis set (ua)")
    print('')
    for row in system['momdip_es_mtx']:
      for val in row:
        print(np.format_float_scientific(val,precision=5,unique=False,pad_left=2), end = " ")
      print('')

  # ========================================================= #
  # Exception handling for the parsing process                #
  # ========================================================= #

  except control_errors.ControlError as error:
    print("")
    print(error)
    exit(-1)
  
  # ===================================================================
  # ===================================================================
  #                       FILES CREATION
  # ===================================================================
  # ===================================================================

  section_title = "2. Data files creation"

  print("")
  print(''.center(len(section_title)+10, '*'))
  print(section_title.center(len(section_title)+10))
  print(''.center(len(section_title)+10, '*'))

  # =========================================================
  # Creating the molecule directory and the data subdirectory
  # =========================================================

  mol_dir = os.path.join(out_dir,mol_name)

  if os.path.exists(mol_dir): # Overwrite was already checked previously, no need to check it again
    shutil.rmtree(mol_dir)

  os.makedirs(mol_dir)

  # Creating the data subdirectory where all the files describing the molecule will be stored
  data_dir = os.path.join(mol_dir,"data")
  os.makedirs(data_dir)

  print("\nThe data subdirectory has been created at %s" % data_dir)

  # Copying the config file and the source file into the data subdirectory
  
  shutil.copy(os.path.join(source_path,source_filename), data_dir)
  shutil.copy(config_file, data_dir)

  print("\nThe files %s and %s have been successfully copied into the data subdirectory." % (source_filename, os.path.basename(config_file)))
  print("")

  # =========================================================
  # Writing already calculated values to files
  # =========================================================

  # Values obtained through the parser script (in CSV format)

  """   print("{:<60}".format('\nCreating states.csv file ... '), end="")
  with open(os.path.join(data_dir, "states.csv"), "w") as f:
    print("Number;Multiplicity;Energy (cm-1);Label", file = f)
    for state in states_list:
      # Print every item in state, separated by ";"
      print(";".join(map(str,state)), file = f)
  print('%12s' % "[ DONE ]")

  print("{:<60}".format('\nCreating coupling_list.csv file ... '), end="")
  with open(os.path.join(data_dir, "coupling_list.csv"), "w") as f:
    print("State 1;State 2;Energy (cm-1)", file = f)
    for line in ori_coupling_list:
      # Print every item in line, separated by ";"
      print(";".join(map(str,line)), file = f)
  print('%12s' % "[ DONE ]")

  print("{:<60}".format('\nCreating momdip_list.csv file ... '), end="")
  with open(os.path.join(data_dir, "momdip_list.csv"), "w") as f:
    print("State 1;State 2;Dipole (a.u.)", file = f)
    for line in momdip_list:
      # Print every item in line, separated by ";"
      print(";".join(map(str,line)), file = f)
  print('%12s' % "[ DONE ]")
  """
  # MIME

  mime_file = config[profile]['created_files']['mime_file']
  print("{:<60}".format('\nCreating %s file ... ' % mime_file), end="")
  np.savetxt(os.path.join(data_dir,mime_file),system['mime'],fmt='% 18.10e')
  print('%12s' % "[ DONE ]")

  # Energies

  energies_file = config[profile]['created_files']['energies_file']
  print("{:<60}".format('\nCreating %s file ... ' % energies_file), end="")
  np.savetxt(os.path.join(data_dir,energies_file + '_cm-1'),system['eigenvalues'],fmt='%1.10e')
  np.savetxt(os.path.join(data_dir,energies_file + '_ua'),eigenvalues_ua,fmt='%1.10e')
  np.savetxt(os.path.join(data_dir,energies_file + '_nm'),eigenvalues_nm,fmt='%1.10e')
  np.savetxt(os.path.join(data_dir,energies_file + '_ev'),eigenvalues_ev,fmt='%1.10e')
  print('%12s' % "[ DONE ]")

  # Eigenvectors matrix and eigenvectors transpose matrix

  mat_et0 = config[profile]['created_files']['mat_et0']
  print("{:<60}".format('\nCreating %s file ... ' % mat_et0), end="")
  np.savetxt(os.path.join(data_dir,mat_et0),system['eigenvectors'],fmt='% 18.10e')
  print('%12s' % "[ DONE ]")

  mat_0te = config[profile]['created_files']['mat_0te']
  print ("{:<60}".format('\nCreating %s file ... ' % mat_0te), end="")
  np.savetxt(os.path.join(data_dir,mat_0te),system['transpose'],fmt='% 18.10e')
  print('%12s' % "[ DONE ]")

  # Dipole moments matrix

  momdip_0 = config[profile]['created_files']['momdip_zero']
  print("{:<60}".format('\nCreating %s file ... ' % momdip_0), end="")
  np.savetxt(os.path.join(data_dir,momdip_0),system['momdip_mtx'],fmt='% 18.10e')
  print('%12s' % "[ DONE ]")

  momdip_e = config[profile]['created_files']['momdip_eigen']
  print("{:<60}".format('\nCreating %s file ... ' % momdip_e), end="")
  np.savetxt(os.path.join(data_dir,momdip_e),system['momdip_es_mtx'],fmt='% 18.10e')	
  print('%12s' % "[ DONE ]")

  # =========================================================
  # Density matrices
  # =========================================================

  # Initial population

  init_file = config[profile]['created_files']['init_pop']
  print("{:<60}".format("\nCreating %s file ..." % init_file), end="") 

  init_pop = np.zeros((len(system['states_list']), len(system['states_list'])),dtype=complex)  # Quick init of a zero-filled matrix
  init_pop[0,0] = 1+0j # All the population is in the ground state at the beginning

  with open(os.path.join(data_dir, init_file + "_1"), "w") as f:
    for line in init_pop:
      for val in line:
        print('( {0.real:.2f} , {0.imag:.2f} )'.format(val), end = " ", file = f)
      print('', file = f)

  print('%12s' % "[ DONE ]")

  # Final population (dummy file but still needed by QOCT-RA)

  final_file = config[profile]['created_files']['final_pop']
  print("{:<60}".format("\nCreating %s file ..." % final_file), end="") 

  final_pop = np.zeros((len(system['states_list']), len(system['states_list'])),dtype=complex)  # Quick init of a zero-filled matrix

  with open(os.path.join(data_dir, final_file + "_1"), "w") as f:
    for line in final_pop:
      for val in line:
        print('( {0.real:.2f} , {0.imag:.2f} )'.format(val), end = " ", file = f)
      print('', file = f)

  print('%12s' % "[ DONE ]")

  # Projectors

  print('')
  proj_file = config[profile]['created_files']['projectors']
  target_state = config[profile]['target_state'] # The type of states that will be targeted by the control
  targets_list = [] # List of target states

  for state in system['states_list']:
    if state['Multiplicity'] == target_state:
      targets_list.append(state['Label'])
      print("{:<59}".format("Creating %s file ..." % (proj_file + state['Label'])), end="")
      proj = np.zeros((len(system['states_list']),len(system['states_list'])),dtype=complex)
      proj[state['Number'],state['Number']] = 1+0j
      with open(os.path.join(data_dir, proj_file + state['Label'] + "_1"), "w") as f:
        for line in proj:
          for val in line:
            print('( {0.real:.2f} , {0.imag:.2f} )'.format(val), end = " ", file = f)
          print('', file = f)
      print('%12s' % "[ DONE ]")

  # ===================================================================
  # ===================================================================
  #                     CALCULATION REQUIREMENTS
  # ===================================================================
  # ===================================================================

  section_title = "3. Calculation requirements"

  print("")
  print(''.center(len(section_title)+10, '*'))
  print(section_title.center(len(section_title)+10))
  print(''.center(len(section_title)+10, '*'))

  # Use the number of states to determine the job scale

  scale_index = len(system['eigenvalues'])

  print("")
  print(''.center(50, '-'))
  print("{:<20} {:<30}".format("Number of states: ", scale_index))

  # Job scale category definition

  jobscale = None

  for scale_limit in job_scales:
    if scale_index > scale_limit:
      continue
    else:
      jobscale = job_scales[scale_limit]
      jobscale_limit = scale_limit
      break

  if not jobscale:
    raise control_errors.ControlError ("ERROR: The number of states is too big for this cluster (%s). Please change cluster." % cluster_name.upper())

  # Obtaining the information associated to our job scale

  job_partition = jobscale.get('partition_name')
  job_walltime = jobscale['time']
  job_cores = jobscale.get('cores')
  job_mem_per_cpu = jobscale['mem_per_cpu']
  delay_command = jobscale.get("delay_command", '')

  print(''.center(50, '-'))
  print("{:<20} {:<30}".format("Cluster: ", cluster_name))
  print("{:<20} {:<30}".format("Job scale: ", jobscale["label"]))
  print("{:<20} {:<30}".format("Job scale limit: ", jobscale_limit))
  print(''.center(50, '-'))
  print("{:<20} {:<30}".format("Job partition: ", (job_partition or "not specified")))
  print("{:<20} {:<30}".format("Job walltime: ", job_walltime))
  print("{:<20} {:<30}".format("Number of cores: ", (job_cores or "default")))
  print("{:<20} {:<30}".format("Mem per CPU (MB): ", job_mem_per_cpu))
  print("{:<20} {:<30}".format("Delay command: ", ("not specified" if delay_command == '' else delay_command)))
  print(''.center(50, '-'))

  # ===================================================================
  # ===================================================================
  #         RENDERING OF THE PARAMETERS FILE AND JOBS LAUNCHING
  # ===================================================================
  # ===================================================================

  section_title = "4. Rendering of the parameters file and jobs launching"

  print("")
  print(''.center(len(section_title)+10, '*'))
  print(section_title.center(len(section_title)+10))
  print(''.center(len(section_title)+10, '*'))

  if not dry_run:
    job_count = 0   # Launched jobs counter, this number will be showed on the console screen at the end of the execution

  # For each projector, render the parameters file and run the corresponding job

  for target in targets_list:

    console_message = "Start procedure for the " + target + " target"
    print("")
    print(''.center(len(console_message)+11, '*'))
    print(console_message.center(len(console_message)+10))
    print(''.center(len(console_message)+11, '*'))

    # Define the job directory for that specific target

    job_dirname = proj_file + target
    job_dir = os.path.join(mol_dir,job_dirname)

    # Get the path to the jinja templates directory (a directory named "templates" in the same directory as this script)
    
    templates_dir = os.path.join(code_dir,"templates")

    # Build a dictionary that will contain all information related to the data directory

    data = {
      "path" : data_dir,
      "mime_file" : mime_file,
      "energies_file" : energies_file,
      "momdip_0" : momdip_0,
      "momdip_e" : momdip_e,
      "mat_et0" : mat_et0,
      "mat_0te" : mat_0te,
      "init_file" : init_file,
      "final_file" : final_file,
      "proj_file" : proj_file
    }

    # Build a dictionary that will contain all information related to the job

    job_specs = {
      "profile" : profile,
      "scale_index" : scale_index,
      "cluster_name" : cluster_name,
      "scale_label" : jobscale["label"],
      "scale_limit" : jobscale_limit,
      "partition" : job_partition,
      "walltime" : job_walltime,
      "mem_per_cpu" : job_mem_per_cpu
    }

    # Build a dictionary containing the additional variables that might be needed by the rendering function
    # If your rendering function needs anything else, you can add it in this dictionary

    misc = {  
        "code_dir" : code_dir,
        "templates_dir" : templates_dir,
        "mol_name" : mol_name,
        "config_name" : os.path.basename(config_file),
        "parsing_fct" : parsing_fct,
        "mol_dir" : mol_dir,
        "job_dirname" : job_dirname,
        "nb_targets" : len(targets_list),
        "target" : target
        }

    # Call the rendering function (defined in control_renderer.py, see the documentation for more information)

    rendered_content, rendered_script = eval("control_renderer." + render_fct)(clusters_cfg, config, system, data, job_specs, misc)

    # Create the job directory

    os.makedirs(job_dir)
    print("\nThe %s job directory has been created in %s" % (job_dirname,mol_dir))

    # Write the content of each rendered file into its own file with the corresponding filename

    for filename, file_content in rendered_content.items():
      rendered_file_path = os.path.join(job_dir, filename)
      with open(rendered_file_path, "w", encoding='utf-8') as result_file:
        result_file.write(file_content)
      print("\nThe %s file has been created into the job directory" % filename)
    
    # Launch the job
    
    if not dry_run:

      print("{:<80}".format("\nLaunching the job ..."), end="")
      os.chdir(job_dir)

      # Define the launch command

      launch_command = submit_command +  " " + delay_command + " " + rendered_script

      # Execute the command and get the command status

      retcode = os.system(launch_command)
      
      # If something went wrong when submitting the job, do not raise an exception and just quit the execution. It is likely a problem linked to the cluster.

      if retcode != 0 :
        print("Job submit encountered an issue")
        print("Aborting ...")
        exit(5)
    
      job_count += 1

      print('%12s' % "[ DONE ]")

    console_message = "End of procedure for the " + target + " target"
    print("")
    print(''.center(len(console_message)+10, '*'))
    print(console_message.center(len(console_message)+10))
    print(''.center(len(console_message)+10, '*'))

  if not dry_run:
    print("")
    if job_count == 1:
      print("One job has been succesfully launched.")
    elif job_count > 1:
      print("%s jobs have been succesfully launched." % job_count)
    else:
      print("WARNING: No job could be launched.")
      
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
