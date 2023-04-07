#!/usr/bin/env python3

########################################################################################################################################################
##                                                        QOCT-GRAD Input Builder & Job Launcher                                                      ##
##                                                                                                                                                    ##
##                                    This script prepares the input files needed to run the QOCT-GRAD program by extracting                          ##
##                                  information from a given source file and launches the corresponding jobs on the cluster.                          ##
##                                      Extended documentation is available at https://chains-ulb.readthedocs.io/                                     ##
##                                                                                                                                                    ##
## /!\ In order to run, this script requires Python 3.5+ as well as NumPy 1.14+, YAML and Jinja2. Ask your cluster(s) administrator(s) if needed. /!\ ##
########################################################################################################################################################

import argparse
import fnmatch
import os
import re
import shutil
import sys
from collections import OrderedDict
from inspect import getsourcefile

import jinja2  # Only needed in the renderer subscript, it is loaded here to check if your python installation does support jinja2
import numpy as np
import yaml
from scipy import constants

# Subscripts (files that end with .py and must be placed in the same directory as this script)

import control_common
import control_renderer
import modelling_fcts
import transition_fcts

# =================================================================== #
# =================================================================== #
#                       COMMAND LINE ARGUMENTS                        #
# =================================================================== #
# =================================================================== #

# Define the arguments needed for the script (here they are defined as named arguments rather than positional arguments, check https://stackoverflow.com/questions/24180527/argparse-required-arguments-listed-under-optional-arguments for more info).

parser = argparse.ArgumentParser(add_help=False, description="This script prepares the input files needed to run the QOCT-GRAD program by extracting information from a given source file and launches the corresponding jobs on the cluster.")

required = parser.add_argument_group('Required arguments')
required.add_argument("-s","--source", type=str, help="Path to the source file containing all the necessary information that needs to be processed.", required=True)
required.add_argument('-cf', '--config', type=str, help="Path to either a YAML configuration file or a directory containing multiple YAML configuration files, extension must be .yml or .yaml.", required=True)
required.add_argument('-cl', '--cluster_name', type=str, help="Name of the cluster where this script is running, as defined in the YAML clusters configuration file.", required=True)
required.add_argument("-p","--profile", type=str, help="Name of the profile you wish to run jobs with, as defined in the YAML clusters configuration file.", required=True)
required.add_argument("-o","--out_dir", type=str, help="Path to the directory where you want to create the directory for your source file containing the subdirectories for each job.", required=True)

optional = parser.add_argument_group('Optional arguments')
optional.add_argument('-h','--help',action='help',default=argparse.SUPPRESS,help='Show this help message and exit')
optional.add_argument("-ow","--overwrite",action="store_true",help="If a data directory for the same source file or a job directory for the same transition-configuration combination already exists, remove it before creating a new one.")
optional.add_argument("-d","--dry_run",action="store_true",help="Do not launch the jobs, just create the files and directories.")
optional.add_argument("-as","--arch_src",action="store_true",help="Archive the source file after it has been processed.")
optional.add_argument("-ac","--arch_cf",action="store_true",help="Archive the configuration files after they have been processed.")

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

    # Define an array for correct English spelling during printing

    special_numbers = {1:"st", 2:"nd", 3:"rd"}

    # Save a reference to the original standard output as it will be modified later on (see https://stackabuse.com/writing-to-a-file-with-pythons-print-function/ for reference)

    original_stdout = sys.stdout

    # Get the size of the terminal in order to have a prettier output, if you need something more robust, go check http://granitosaurus.rocks/getting-terminal-size.html

    columns, rows = shutil.get_terminal_size()

    # Output Header

    print("".center(columns,"*"))
    print("")
    print("EXECUTION OF THE QOCT-GRAD INPUT BUILDER & JOB LAUNCHER BEGINS NOW".center(columns))
    print("")
    print("".center(columns,"*"))

    # ========================================================= #
    # Read command line arguments                               #
    # ========================================================= #

    args = parser.parse_args()

    # Required arguments

    source = args.source                     # Source file containing all the necessary information
    config_inp = args.config                 # YAML configuration file or directory containing the YAML configuration files

    cluster_name = args.cluster_name         # Name of the cluster where this script is running, as defined in the clusters configuration YAML file
    profile = args.profile                   # Name of the profile for which files need to be created
    out_dir = args.out_dir                   # Directory where all jobs subdirectories will be created

    # Optional arguments

    overwrite = args.overwrite               # Flag for removing job subdirectories before creating a new one, if they have the same name
    dry_run = args.dry_run                   # Flag to not launch the jobs and just create the files

    arch_src = args.arch_src                 # Flag for archiving the source file after it has been processed
    arch_cf = args.arch_cf                   # Flag for archiving the configuration files after they have been processed

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

    clusters_file = control_common.check_abspath(os.path.join(code_dir,"clusters.yml"),"YAML clusters configuration file","file")
    print ("{:<80}".format('\nLoading the clusters configuration file "clusters.yml" ...'), end="")
    with open(clusters_file, 'r', encoding='utf-8') as f_clusters:
      clusters_cfg = yaml.load(f_clusters, Loader=yaml.FullLoader)
    print('%12s' % "[ DONE ]")

    # Check the name of the cluster

    if cluster_name not in clusters_cfg:
      raise control_common.ControlError ("ERROR: There is no information about the %s cluster in the clusters configuration file. Please add relevant information or change the cluster before proceeding further." % cluster_name)

    print ("{:<40} {:<100}".format('\nCluster name:',cluster_name))

    # Check if the submit_command key has been defined

    submit_command = clusters_cfg[cluster_name].get("submit_command")

    if submit_command is None:
      raise control_common.ControlError ("ERROR: There is no defined submit_command for the %s cluster in the clusters configuration file." % cluster_name) 

    # Check if the profile exists 

    if "profiles" not in clusters_cfg[cluster_name]:
      raise control_common.ControlError ('ERROR: There is no "profiles" key defined for the %s cluster in the clusters configuration file. Consult official documentation for details.' % cluster_name)    

    if profile not in clusters_cfg[cluster_name]["profiles"]:
      raise control_common.ControlError ("ERROR: The specified profile (%s) is unknown on this cluster. Possible profiles include: %s \nPlease use one of those, change cluster or add information for this profile to the clusters configuration file." % (profile, ', '.join(profile for profile in clusters_cfg[cluster_name]["profiles"].keys())))
    
    print ("{:<40} {:<100}".format('\nProfile:',profile))

    # Get the modelling function that will fetch the relevant information from the source file - defined in modelling_fcts.py and specified in the clusters configuration file

    modelling_fct = clusters_cfg[cluster_name]["profiles"][profile].get("modelling_function")

    if modelling_fct is None:
      raise control_common.ControlError ("ERROR: There is no defined modelling function for the %s profile in the %s cluster in the clusters configuration file." % (profile, cluster_name))
    if (modelling_fct) not in dir(modelling_fcts) or not callable(getattr(modelling_fcts, modelling_fct)):
      raise control_common.ControlError ("ERROR: There is no modelling function named %s defined in modelling_fcts.py." % modelling_fct)

    print ("{:<40} {:<100}".format('\nModelling function for that profile:',modelling_fct))

    # Get the transition function that will determine the initial and target states of the control procedure - defined in transition_fcts.py and specified in the clusters configuration file

    transition_fct = clusters_cfg[cluster_name]["profiles"][profile].get("transition_function")

    if transition_fct is None:
      raise control_common.ControlError ("ERROR: There is no defined transition function for the %s profile in the %s cluster in the clusters configuration file." % (profile, cluster_name))
    if (transition_fct) not in dir(transition_fcts) or not callable(getattr(transition_fcts, transition_fct)):
      raise control_common.ControlError ("ERROR: There is no transition function named %s defined in transition_fcts.py." % transition_fct)

    print ("{:<40} {:<100}".format('\nTransition function for that profile:',transition_fct))

    # Define the rendering function that will render the job script and the input file (depends on the profile)  - defined in control_renderer.py

    render_fct = clusters_cfg[cluster_name]["profiles"][profile].get("rendering_function")

    if render_fct is None:
      raise control_common.ControlError ("ERROR: There is no defined rendering function for the %s profile in the %s cluster in the clusters configuration file." % (profile, cluster_name))
    if (render_fct) not in dir(control_renderer) or not callable(getattr(control_renderer, render_fct)):
      raise control_common.ControlError ("ERROR: There is no rendering function named %s defined in renderer.py." % render_fct)

    print ("{:<40} {:<100}".format('\nRendering function for that profile:',render_fct))

    # ========================================================= #
    # Establishing the different job scales                     #
    # ========================================================= #

    # Gather all the different job scales from the clusters configuration file in a temporary dictionary

    job_scales_tmp = clusters_cfg[cluster_name]['profiles'][profile].get('job_scales')

    if job_scales_tmp is None:
      raise control_common.ControlError ("ERROR: There is no defined job scales for the %s profile in the %s cluster in the clusters configuration file." % (profile, cluster_name)) 

    # Initialize the final dictionary where the job scales will be sorted by their upper limit

    job_scales = {}

    # Check the job scales

    required_keys = ["label", "scale_limit", "time", "memory"]
    control_common.check_keys(required_keys,job_scales_tmp,"Job scales of the %s profile in the %s cluster in the clusters configuration file." % (profile, cluster_name))

    # Extract the scale upper limit from the job scales

    for scale in job_scales_tmp:
  
      scale_limit = scale.pop('scale_limit')
      job_scales[scale_limit] = scale

    # Sort the different job scales by their upper limit and store them in the job_scales dictionary

    job_scales = OrderedDict(sorted(job_scales.items()))

    print("\nJob scales for %s on %s:" % (profile,cluster_name))
    print("")
    print(''.center(136, '-'))
    print ("{:<15} {:<20} {:<20} {:<20} {:<20} {:<40}".format('Scale Limit','Label','Partition Name','Time','Memory (MB)','Delay Command'))
    print(''.center(136, '-'))
    for scale_limit, scale in job_scales.items():
      print ("{:<15} {:<20} {:<20} {:<20} {:<20} {:<40}".format(scale_limit, scale['label'], scale.get('partition_name', "not specified"), scale['time'], scale['memory'], scale.get('delay_command', "not specified")))
    print(''.center(136, '-'))

    # ========================================================= #
    # Check other arguments                                     #
    # ========================================================= #

    out_dir = control_common.check_abspath(out_dir,"Command line argument -o / --out_dir","directory")
    print ("{:<40} {:<100}".format('\nJobs main directory:',out_dir))

    # Check source file
    # =================

    # Check the existence of the source file, then get its name and the name of the directory where it is located

    source = control_common.check_abspath(source,"Command line argument -s / --source","file")
    print ("{:<40} {:<100}".format('\nSource file:',source))

    source_path = os.path.dirname(source)
    source_filename = os.path.basename(source)
    source_name = str(source_filename.split('.')[0]) # Getting rid of the format extension to get the name of the source

    # Check config file(s)
    # ====================

    config_inp = control_common.check_abspath(config_inp,"Command line argument -cf / --config")

    # If the argument config_inp is a directory, we need to look for every YAML configuration file in that directory.

    if os.path.isdir(config_inp):

      print("{:<40} {:<100}".format("\nLooking for .yml or .yaml files in", config_inp + " ..."), end="")

      config_inp_path = config_inp

      # Define which type of file we are looking for in a case-insensitive way (see https://gist.github.com/techtonik/5694830)

      rule = re.compile(fnmatch.translate("*.yml"), re.IGNORECASE)
      rule2 = re.compile(fnmatch.translate("*.yaml"), re.IGNORECASE)

      # Find all matching files in config_inp directory

      config_inp_list = [config for config in os.listdir(config_inp) if (rule.match(config) or rule2.match(config))]

      if config_inp_list == []:
        raise control_common.ControlError ("ERROR: Can't find any .yml or .yaml file in %s" % config_inp_path)

      print('%12s' % "[ DONE ]")

    # If given a single config file as argument, check its extension.

    else:

      print ("{:<40} {:<100}".format('\nConfiguration file:',config_inp))

      if os.path.isfile(config_inp) and os.path.splitext(config_inp)[-1].lower() != (".yml") and os.path.splitext(config_inp)[-1].lower() != (".yaml"):
        raise control_common.ControlError ("  ^ ERROR: This is not a YAML file (YAML file extension is either .yml or .yaml).")

      config_inp_path = os.path.dirname(config_inp)
      config_inp_file = os.path.basename(config_inp)
      config_inp_list = [config_inp_file]

    # Check dry run option
    # ====================

    if dry_run:
      print("\nThe dry run option has been enabled: the job files and directories will be created but the jobs will not be submitted to the job scheduler.")
    else:
      job_count = 0   # Launched jobs counter, this number will be shown on the console screen at the end of the execution

  # ========================================================= #
  # Exception handling for the preparation step               #
  # ========================================================= #

  except control_common.ControlError as error:
    print("")
    print(error)
    exit(-1)

  # =================================================================== #
  # =================================================================== #
  #                        SOURCE FILE TREATMENT                        #
  # =================================================================== #
  # =================================================================== #

  # For more information on try/except structures, see https://www.tutorialsteacher.com/python/exception-handling-in-python
  try:

    console_message = "Start of the source file treatment"
    print("")
    print(''.center(len(console_message)+11, '*'))
    print(console_message.center(len(console_message)+10))
    print(''.center(len(console_message)+11, '*'))

    # Create an output log file that will contain all the information about the source file treatment (system modelling and determining transitions)

    data_log_name = source_name + ".log"
    data_log = open(os.path.join(out_dir,data_log_name), 'w', encoding='utf-8')

    # Redirect standard output to the data_log file (see https://stackabuse.com/writing-to-a-file-with-pythons-print-function/ for reference)

    sys.stdout = data_log
    
    print ("{:<49} {:<100}".format('Cluster name:',cluster_name))      
    print ("{:<50} {:<100}".format('\nProfile:',profile))

    # =================================================================== #
    # =================================================================== #
    #                           SYSTEM MODELLING                          #
    # =================================================================== #
    # =================================================================== #

    section_title = "1. System modelling"

    print("")
    print(''.center(len(section_title)+10, '*'))
    print(section_title.center(len(section_title)+10))
    print(''.center(len(section_title)+10, '*'))

    # Console screen notification (we need to temporarily switch the standard outputs to show this message on screen and not in the log file)

    sys.stdout = original_stdout                                 
    print ("{:<80}".format('\nModelling the system ...'), end="")
    sys.stdout = data_log  

    # ========================================================= #
    # Load the source file                                      #
    # ========================================================= #

    print ("{:<50}".format('\nLoading %s file ... ' % source_filename), end="")
    with open(source, 'r') as source_file:
      source_content = source_file.read().splitlines()
    print("[ DONE ]")

    # Cleaning up the source file from surrounding spaces and blank lines

    source_content = list(map(str.strip, source_content))   # Remove leading & trailing blank/spaces
    source_content = list(filter(None, source_content))     # Remove blank lines/no char

    # ========================================================= #
    # Modelling function                                        #
    # ========================================================= #

    print ("{:<50} {:<100}".format('\nModelling function:',modelling_fct))

    # Call the modelling function (defined in modelling_fcts.py, see the documentation for more information)

    system = eval("modelling_fcts." + modelling_fct)(source_content)

    # Check the system dictionary

    if not isinstance(system, dict):
      raise control_common.ControlError ('ERROR: The "system" variable returned by the %s modelling function is not a dictionary.' % modelling_fct) 

    required_keys = ["states_list", "momdip_mtx"]
    control_common.check_keys(required_keys,system,"The 'system' dictionary returned by the %s modelling function." % modelling_fct)

    # Check the states list

    if not isinstance(system["states_list"], list):
      raise control_common.ControlError ('ERROR: The "states_list" value in the system dictionary returned by the %s modelling function is not a list.' % modelling_fct)

    required_keys = ["number", "label", "energy"]
    control_common.check_keys(required_keys,system["states_list"],"The 'states_list' list of the 'system' dictionary returned by the %s modelling function." % modelling_fct)

    control_common.is_consecutive([state['number'] for state in system['states_list']],"State numbers from the modelling function")

    # Check the dipole moments matrices

    if not isinstance(system["momdip_mtx"], dict):
      raise control_common.ControlError ('ERROR: The "momdip_mtx" value in the system dictionary returned by the %s modelling function is not a dictionary.' % modelling_fct)

    for momdip_key in system["momdip_mtx"]:
      if not isinstance(system["momdip_mtx"][momdip_key], (list, np.ndarray)):
        raise control_common.ControlError ('ERROR: The "%s" value in the "momdip_mtx" dictionary returned by the %s modelling function is neither a list nor a NumPy array.' % (momdip_key, modelling_fct))

    print("\nThe system has been succesfully modelled.")

    # Console screen notification (we need to temporarily switch the standard outputs to show this message on screen and not in the log file)
    
    sys.stdout = original_stdout                                 
    print('%12s' % "[ DONE ]")
    sys.stdout = data_log 

    # =================================================================== #
    # =================================================================== #
    #                             TRANSITIONS                             #
    # =================================================================== #
    # =================================================================== #

    section_title = "2. Determining the transitions"

    print("")
    print(''.center(len(section_title)+10, '*'))
    print(section_title.center(len(section_title)+10))
    print(''.center(len(section_title)+10, '*'))

    # Console screen notification (we need to temporarily switch the standard outputs to show this message on screen and not in the log file)

    sys.stdout = original_stdout                                 
    print ("{:<80}".format('\nDetermining the transitions ...'), end="")
    sys.stdout = data_log  

    # ========================================================= #
    # Transition function                                       #
    # ========================================================= #

    print ("{:<50} {:<100}".format('\nTransition function:',transition_fct))

    # Call the transition function (defined in transition_fcts.py, see the documentation for more information)

    transitions_list = eval("transition_fcts." + transition_fct)(system)

    # Check the transitions list

    if not isinstance(transitions_list, list):
      raise control_common.ControlError ('ERROR: The transitions_list returned by the %s transition function is not a list.' % transition_fct) 

    required_keys = ["label","init_file","init_content","target_file","target_content","momdip_key"]
    control_common.check_keys(required_keys,transitions_list,"Transitions list returned by the %s transition function" % transition_fct)
 
    # Console screen notification (we need to temporarily switch the standard outputs to show this message on screen and not in the log file)
    
    sys.stdout = original_stdout                                 
    print('%12s' % "[ DONE ]")
    sys.stdout = data_log    

    # =================================================================== #
    # =================================================================== #
    #                         DATA FILES CREATION                         #
    # =================================================================== #
    # =================================================================== #

    section_title = "2bis. Creating the data files"

    print("")
    print(''.center(len(section_title)+10, '*'))
    print(section_title.center(len(section_title)+10))
    print(''.center(len(section_title)+10, '*'))

    # Console screen notification (we need to temporarily switch the standard outputs to show this message on screen and not in the log file)

    sys.stdout = original_stdout                                 
    print ("{:<80}".format('\nCreating the data files ...'), end="")
    sys.stdout = data_log    

    # ========================================================= #
    # Creating the molecule directory                           #
    # ========================================================= #

    mol_dir = os.path.join(out_dir,source_name)

    os.makedirs(mol_dir,exist_ok=True)

    print ("{:<20} {:<100}".format('\nMolecule directory:',mol_dir))

    # ========================================================= #
    # Creating the profile subdirectory                         #
    # ========================================================= #

    pro_dir = os.path.join(mol_dir,profile)

    os.makedirs(pro_dir,exist_ok=True)

    print ("{:<20} {:<100}".format('\nProfile subdirectory:',pro_dir))

    # ========================================================= #
    # Creating the data subdirectory                            #
    # ========================================================= #

    data_dir = os.path.join(pro_dir,"data")

    if os.path.exists(data_dir) and not overwrite:
      raise control_common.ControlError ("ERROR: A data directory for the %s source file already exists in %s !" % (source_name, pro_dir))
    
    os.makedirs(data_dir,exist_ok=True)

    print ("{:<20} {:<100}".format('\nData directory:',data_dir))

    # ========================================================= #
    # Creating the system data files                            #
    # ========================================================= #

    # Energies

    energies_file = "energies"
    with open(os.path.join(data_dir,energies_file), "w") as f:
      for state in system['states_list']:
        print("{:1.10e}".format(state['energy']), file = f)
    print("    ├── The energies file ('%s') has been created into the directory" % energies_file)

    # Dipole moments matrices

    for momdip_key in system['momdip_mtx']:

      momdip_mtx_file = 'momdip_mtx_' + momdip_key
      np.savetxt(os.path.join(data_dir,momdip_mtx_file),system['momdip_mtx'][momdip_key],fmt='% 18.10e')
      print("    ├── The transition dipole moment matrix file corresponding to the '%s' key ('%s') has been created into the directory" % (momdip_key, momdip_mtx_file))

    # Conversion matrix from zero order states to eigenstates (if it has been provided in the modelling function, otherwise it will just be a unitary matrix)

    if 'conv_mtx' in system.keys():

      conv_mtx_file = 'conv_mtx'
      np.savetxt(os.path.join(data_dir,conv_mtx_file),system['conv_mtx'],fmt='% 18.10e')
      print("    ├── The conversion matrix file ('%s') has been created into the directory" % conv_mtx_file)

    else:

      conv_mtx_file = 'conv_mtx'
      conv_mtx = np.zeros((len(system['states_list']), len(system['states_list'])), dtype=float)
      np.fill_diagonal(conv_mtx, 1.0)
      np.savetxt(os.path.join(data_dir,conv_mtx_file),conv_mtx,fmt='% 18.10e')
      print("    ├── A matrix of ones ('%s') has been created into the directory, acting as a conversion matrix" % conv_mtx_file)

    # ========================================================= #
    # Creating the transition files                             #
    # ========================================================= #

    for transition in transitions_list:

      init_filename = transition["init_file"] + "1"
      if not os.path.exists(os.path.join(data_dir, init_filename)):
        with open(os.path.join(data_dir, init_filename), "w") as f:
          for line in transition["init_content"]:
            for val in line:
              print('( {0.real:.10e} , {0.imag:.10e} )'.format(val), end = " ", file = f)
            print('', file = f)
        print("    ├── The %s initial states file has been created into the directory" % init_filename)
      
      target_filename = transition["target_file"] + "1"
      if not os.path.exists(os.path.join(data_dir, target_filename)):
        with open(os.path.join(data_dir, target_filename), "w") as f:
          for line in transition["target_content"]:
            for val in line:
              print('( {0.real:.10e} , {0.imag:.10e} )'.format(val), end = " ", file = f)
            print('', file = f)
        print("    ├── The %s target states file has been created into the directory" % target_filename)
    
      if 'projector' in transition.keys():
        projector_filename = 'projector_' + transition["target_file"] + "1"
        if not os.path.exists(os.path.join(data_dir, projector_filename)):
          with open(os.path.join(data_dir, projector_filename), "w") as f:
            for line in transition["projector"]:
              for val in line:
                print('( {0.real:.10e} , {0.imag:.10e} )'.format(val), end = " ", file = f)
              print('', file = f)
          print("    ├── The projector for the %s target states file ('%s') has been created into the directory" % (target_filename,projector_filename))

    # ========================================================= #
    # Other files                                               #
    # ========================================================= #

    # Copying the source file into the data subdirectory
    
    shutil.copy(os.path.join(source_path,source_filename), data_dir)
    print("    └── The source file (%s) has been successfully copied into the directory." % source_filename)

    # Console screen notification (we need to temporarily switch the standard outputs to show this message on screen and not in the log file)
    
    sys.stdout = original_stdout                                 
    print('%12s' % "[ DONE ]")
    sys.stdout = data_log    

    # =================================================================== #
    # =================================================================== #
    #                             Job Scaling                             #
    # =================================================================== #
    # =================================================================== #

    section_title = "3. Job scaling"

    print("")
    print(''.center(len(section_title)+10, '*'))
    print(section_title.center(len(section_title)+10))
    print(''.center(len(section_title)+10, '*'))

    sys.stdout = original_stdout                                 
    print ("{:<80}".format('\nDetermining the job scale ...'), end="")
    sys.stdout = data_log  

    subsection_title = "A. Available job scales"

    print("")
    print("")
    print(subsection_title)
    print(''.center(len(subsection_title), '='))

    print("")
    print(''.center(136, '-'))
    print ("{:<15} {:<20} {:<20} {:<20} {:<20} {:<40}".format('Scale Limit','Label','Partition Name','Time','Memory (MB)','Delay Command'))
    print(''.center(136, '-'))
    for scale_limit, scale in job_scales.items():
      print ("{:<15} {:<20} {:<20} {:<20} {:<20} {:<40}".format(scale_limit, scale['label'], scale.get('partition_name', "not specified"), scale['time'], scale['memory'], scale.get('delay_command', "not specified")))
    print(''.center(136, '-'))

    subsection_title = "B. Calculation requirements"

    print("")
    print("")
    print(subsection_title)
    print(''.center(len(subsection_title), '='))

    # Use the number of states to determine the job scale

    scale_index = len(system['states_list'])

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
      raise control_common.ControlError ("ERROR: The number of states is too big for this cluster (%s). Please change cluster." % cluster_name)

    # Obtaining the information associated to our job scale

    job_partition = jobscale.get('partition_name')
    job_walltime = jobscale['time']
    job_memory = jobscale['memory']
    delay_command = jobscale.get("delay_command", '')

    print(''.center(50, '-'))
    print("{:<20} {:<30}".format("Cluster: ", cluster_name))
    print("{:<20} {:<30}".format("Job scale: ", jobscale["label"]))
    print("{:<20} {:<30}".format("Job scale limit: ", jobscale_limit))
    print(''.center(50, '-'))
    print("{:<20} {:<30}".format("Job partition: ", (job_partition or "not specified")))
    print("{:<20} {:<30}".format("Job walltime: ", job_walltime))
    print("{:<20} {:<30}".format("Memory (MB): ", job_memory))
    print("{:<20} {:<30}".format("Delay command: ", ("not specified" if delay_command == '' else delay_command)))
    print(''.center(50, '-'))

    # ========================================================= #
    # End of logging for the data files generation              #
    # ========================================================= #

    sys.stdout = original_stdout                                  # Reset the standard output to its original value
    print('%12s' % "[ DONE ]")
    data_log.close()
    if os.path.exists(os.path.join(data_dir,data_log_name)):
      os.remove(os.path.join(data_dir,data_log_name))
    shutil.move(os.path.join(out_dir,data_log_name), data_dir)    # Archive the log file in the data directory

    console_message = "End of the source file treatment"
    print("")
    print(''.center(len(console_message)+11, '*'))
    print(console_message.center(len(console_message)+10))
    print(''.center(len(console_message)+11, '*'))

  # ========================================================= #
  # Exception handling for the data files generation          #
  # ========================================================= #

  except control_common.ControlError as error:
    sys.stdout = original_stdout                        # Reset the standard output to its original value
    print(error)
    os.remove(os.path.join(out_dir,data_log_name))      # Remove the log file since there was a problem
    exit(-1)

  # =================================================================== #
  # =================================================================== #
  #           RENDERING THE TEMPLATES AND SUBMITTING THE JOBS           #
  # =================================================================== #
  # =================================================================== #

  problem_cf = [] # Empty list that will contain the names of the configuration files for which a problem has occurred (those configuration files will not be archived even if arch_cf was set)

  # For each transition-configuration combination, render the parameters file and run the corresponding job

  for transition in transitions_list:

    console_message = "Start procedure for the transition " + transition["label"]
    print("")
    print(''.center(len(console_message)+11, '*'))
    print(console_message.center(len(console_message)+10))
    print(''.center(len(console_message)+11, '*'))

    for config_filename in config_inp_list: 

      # For more information on try/except structures, see https://www.tutorialsteacher.com/python/exception-handling-in-python
      try:

        # Getting rid of the format extension to get the name of the configuration

        config_name = str(config_filename.split('.')[0])
        print("{:<80}".format("\nTreating '%s' transition with '%s' configuration ..." % (transition["label"], config_name)), end="")

        # Create an output log file for each transition - config combination

        log_name = transition["label"] + "_" + config_name + ".log"
        log = open(os.path.join(pro_dir,log_name), 'w', encoding='utf-8')

        # Redirect standard output to the log file (see https://stackabuse.com/writing-to-a-file-with-pythons-print-function/ for reference)

        sys.stdout = log

        # Define the name and path of the job directory for that specific transition-configuration combination

        job_dirname = transition["label"] + "_" + config_name
        job_dir = os.path.join(pro_dir,job_dirname)

        if os.path.exists(job_dir) and not overwrite:
          raise control_common.ControlError ("ERROR: A directory for the %s transition with the '%s' configuration already exists in %s !" % (transition["label"], config_name, pro_dir))

        # ========================================================= #
        # Rendering the templates                                   #
        # ========================================================= #

        section_title = "1. Rendering the templates"

        print("")
        print(''.center(len(section_title)+10, '*'))
        print(section_title.center(len(section_title)+10))
        print(''.center(len(section_title)+10, '*'))

        print ("{:<50} {:<100}".format('\nRendering function:',render_fct))

        # Load config file

        print ("{:<51}".format('\nLoading %s file ...' % config_filename), end="")
        with open(os.path.join(config_inp_path,config_filename), 'r', encoding='utf-8') as f_config:
          config = yaml.load(f_config, Loader=yaml.FullLoader)
        print("[ DONE ]")

        subsection_title = "A. Rendering function"

        print("")
        print("")
        print(subsection_title)
        print(''.center(len(subsection_title), '='))

        # Get the path to the jinja templates directory (a directory named "templates" in the same directory as this script)
        
        templates_dir = os.path.join(code_dir,"templates")

        # Build a dictionary that will contain all information related to the data directory

        data = {
          # Path of the generic data files
          "main_path" : data_dir,
          "energies_path" : os.path.join(data_dir,energies_file),
          "conv_mtx_path" : os.path.join(data_dir,conv_mtx_file),
          # Path of the data files specific to this transition          
          "init_path" : os.path.join(data_dir,transition['init_file']),
          "target_path" : os.path.join(data_dir,transition['target_file']),
          "momdip_mtx_path" : os.path.join(data_dir,'momdip_mtx_' + transition['momdip_key'])
        }

        if 'projector' in transition.keys():
          data.update({
            "projector_path" : os.path.join(data_dir,'projector_' + transition["target_file"])
          })

        # Build a dictionary that will contain all information related to the job

        job_specs = {
          "profile" : profile,
          "scale_index" : scale_index,
          "cluster_name" : cluster_name,
          "scale_label" : jobscale["label"],
          "scale_limit" : jobscale_limit,
          "partition" : job_partition,
          "walltime" : job_walltime,
          "memory" : job_memory
        }

        # Build a dictionary containing the additional variables that might be needed by the rendering function
        # If your rendering function needs anything else, you can add it in this dictionary

        misc = {  
            "code_dir" : code_dir,
            "templates_dir" : templates_dir,
            "source_name" : source_name,
            "source_content" : source_content,
            "mol_dir" : mol_dir,
            "pro_dir" : pro_dir,
            "config_name" : config_name,
            "job_dirname" : job_dirname,
            "transition" : transition,
            "transitions_list" : transitions_list
        }

        # Call the rendering function (defined in control_renderer.py, see the documentation for more information)

        try:
          rendered_content, rendered_script = eval("control_renderer." + render_fct)(clusters_cfg, config, system, data, job_specs, misc)
        except KeyError as error:
          raise control_common.ControlError ("ERROR: The '%s' rendering function tried to access an unknown key (%s). \nCheck your clusters configuration file ('clusters.yml') and the '%s' configuration file, as well as the spelling and definition of your variables in the rendering function." % (render_fct,error,config_filename))

        # Check the rendered_content dictionary

        if not isinstance(rendered_content, dict):
          raise control_common.ControlError ('ERROR: The "rendered_content" returned by the %s rendering function is not a dictionary.' % render_fct) 

        print("\nAll the templates have been succesfully rendered.")

        # ========================================================= #
        # Creating the job directory and its content                #
        # ========================================================= #

        subsection_title = "B. Creating the files"

        print("")
        print("")
        print(subsection_title)
        print(''.center(len(subsection_title), '='))

        print ("{:<20} {:<100}".format('\nJob directory:',job_dir))

        if os.path.exists(job_dir): # Overwrite was already checked previously, no need to check it again
          print("    /!\ Deleting the old %s directory ..." % job_dir, end="")
          shutil.rmtree(job_dir)
          print('%12s' % "[ DONE ]")

        os.makedirs(job_dir)

        # Write the content of each rendered file into its own file with the corresponding filename

        for filename, file_content in rendered_content.items():
          rendered_file_path = os.path.join(job_dir, filename)
          with open(rendered_file_path, "w", encoding='utf-8') as result_file:
            result_file.write(file_content)
          print("    ├── The %s file has been created into the directory" % filename)
        
        # Copying the config file into the job directory
        
        shutil.copy(os.path.join(config_inp_path,config_filename), job_dir)

        print("    └── The configuration file (%s) has been successfully copied into the directory." % config_filename)

        # ========================================================= #
        # Submitting the job                                        #
        # ========================================================= #

        section_title = "2. Submitting the job"

        print("")
        print("")
        print(''.center(len(section_title)+10, '*'))
        print(section_title.center(len(section_title)+10))
        print(''.center(len(section_title)+10, '*'))

        # Launch the job
        
        if not dry_run:

          print("{:<51}".format("\nLaunching the job ..."), end="")
          os.chdir(job_dir)

          # Define the launch command

          launch_command = submit_command + " " + delay_command + " " + rendered_script

          # Execute the command and get the command status

          retcode = os.system(launch_command)
          
          # If something went wrong when submitting the job, do not raise an exception and just quit the execution. It is likely a problem linked to the cluster.

          if retcode != 0 :
            sys.stdout = original_stdout                             # Reset the standard output to its original value
            print("Job submit encountered an issue")
            print("Aborting ...")
            exit(5)
        
          job_count += 1

          print("[ DONE ]")
        
        else:

          print("\nThe dry run option has been enabled, this job will not be submitted to the job scheduler.")
          sys.stdout = original_stdout                            # Reset the standard output to its original value
          print('%12s' % "[ DONE ]")

        # ================================================================== #
        # End of logging for that transition - config combination            #
        # ================================================================== #

        sys.stdout = original_stdout                            # Reset the standard output to its original value
        log.close()                                             # End of logging
        shutil.move(os.path.join(pro_dir,log_name), job_dir)    # Archive the log file in the job directory

      # ========================================================= #
      # Exception handling for the rendering and submitting steps #
      # ========================================================= #

      # In case of an error specific to the configuration file, skip it and do not archive the source file even if arch_src was set.

      except control_common.ControlError as error:
        sys.stdout = original_stdout                            # Reset the standard output to its original value
        print(error)
        print("Skipping configuration '%s'" % config_name)
        os.remove(os.path.join(pro_dir,log_name))               # Remove the log file since there was a problem
        problem_cf.append(config_filename)                      # Add the name of this configuration file to the list as to enforce not archiving it
        arch_src = False                                        # Flag to notify that a problem has occurred with this configuration file and to enforce not archiving the source file
        continue        

    console_message = "End of procedure for the transition " + transition["label"]
    print("")
    print(''.center(len(console_message)+10, '*'))
    print(console_message.center(len(console_message)+10))
    print(''.center(len(console_message)+10, '*'))

  # After all the transition-configuration combinations have been treated, archive the source file if arch_src has been set and there was no problem.

  if arch_src:
    launched_dir = os.path.join(source_path,"launched")      # Directory where the source file will be put after having been treated by this script, it will be created inside the directory where it was.
    os.makedirs(launched_dir, exist_ok=True)
    if os.path.exists(os.path.join(launched_dir,source_filename)):
      os.remove(os.path.join(launched_dir,source_filename))
    shutil.move(os.path.join(source_path,source_filename), launched_dir)
    print("\nSource file archived to %s" % launched_dir)

  # After all the transition-configuration combinations have been treated, archive the configuration files if arch_cf has been set and there was no problem.

  if arch_cf:
    for config_filename in config_inp_list:
      if config_filename not in problem_cf:
        launched_dir = os.path.join(config_inp_path,"launched") # Directory where the configuration files will be put after having been treated by this script, it will be created inside the directory where were all the configuration files.
        os.makedirs(launched_dir, exist_ok=True)
        if os.path.exists(os.path.join(launched_dir,config_filename)):
          os.remove(os.path.join(launched_dir,config_filename))
        shutil.move(os.path.join(config_inp_path,config_filename), launched_dir)
        print("\nThe '%s' config file has been archived to %s" % (config_filename,launched_dir))

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
