#!/usr/bin/env python3

################################################################################################################################################
##                                                 Ab Initio Input Builder & Job Launcher                                                     ##
##                                                                                                                                            ##
##                       For one or more geometry files, one or more configuration files and a given ab initio program,                       ##
##            this script prepares the input files needed for each calculation and launches the corresponding jobs on the cluster.            ##
##                                Extended documentation is available at https://chains-ulb.readthedocs.io/                                   ##
##                                                                                                                                            ##
##   /!\ In order to run, this script requires Python 3.5+ as well as YAML and Jinja2. Ask your cluster(s) administrator(s) if needed. /!\    ##
################################################################################################################################################

import argparse
import fnmatch
import os
import re
import shutil
import sys
from collections import OrderedDict
from inspect import getsourcefile

import jinja2  # Only needed in the renderer subscript, it is loaded here to check if your python installation does support jinja2
import yaml

# Subscripts (files that end with .py and must be placed in the same directory as this script)

import abin_errors
import geom_scan
import renderer
import scaling_fcts

# =================================================================== #
# =================================================================== #
#                       COMMAND LINE ARGUMENTS                        #
# =================================================================== #
# =================================================================== #

# Define the arguments needed for the script (here they are defined as named arguments rather than positional arguments, check https://stackoverflow.com/questions/24180527/argparse-required-arguments-listed-under-optional-arguments for more info).

parser = argparse.ArgumentParser(add_help=False, description="For one or more geometry files, one or more configuration files and a given ab initio program, this script prepares the input files needed for each calculation and launches the corresponding jobs on the cluster. Extended documentation is available at https://chains-ulb.readthedocs.io/. In order to run, this script requires Python 3.5+ as well as YAML and Jinja2.")

required = parser.add_argument_group('Required arguments')
required.add_argument("-m","--mol_inp", type=str, help="Path to either a geometry file or a directory containing multiple geometry files.", required=True)
required.add_argument('-cf', '--config', type=str, help="Path to either a YAML configuration file or a directory containing multiple YAML configuration files, extension must be .yml or .yaml.", required=True)
required.add_argument("-p","--program", type=str, help="Name of the program you wish to run jobs with, as defined in the YAML clusters configuration file.", required=True)
required.add_argument("-o","--out_dir", type=str, help="Path to the directory where you want to create the subdirectories for each job.", required=True)
required.add_argument('-cl', '--cluster_name', type=str, help="Name of the cluster where this script is running, as defined in the YAML clusters configuration file.", required=True)
#required.add_argument('-f', '--format', type=str, help="Format of the geometry files that need to be read.", required=True) #Uncomment this line if you add new scanning functions to geom_scan.py

optional = parser.add_argument_group('Optional arguments')
optional.add_argument('-h','--help',action='help',default=argparse.SUPPRESS,help='Show this help message and exit.')
optional.add_argument("-ow","--overwrite",action="store_true",help="If a job subdirectory for a geometry-configuration combination already exists, remove it before creating a new one.")
optional.add_argument("-d","--dry_run",action="store_true",help="Do not launch the jobs, just create the files and directories.")
optional.add_argument('--max_mol', type=int, help="Maximum number of geometry files that will be succesfully processed, limiting the number of jobs that will then be launched.")
optional.add_argument('--max_cf', type=int, help="Maximum number of configuration files that will be succesfully processed, limiting the number of jobs that will then be launched.")
optional.add_argument("-km","--keep_mol",action="store_true",help="Do not archive the geometry files after they have been processed and leave them where they are.")
optional.add_argument("-kc","--keep_cf",action="store_true",help="Do not archive the configuration files after they have been processed and leave them where they are.")

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

    # Save a reference to the original standard output as it will be modified later on (see https://stackabuse.com/writing-to-a-file-with-pythons-print-function/ for reference)

    original_stdout = sys.stdout

    # Get the size of the terminal in order to have a prettier console output, if you need something more robust, go check http://granitosaurus.rocks/getting-terminal-size.html

    columns, rows = shutil.get_terminal_size()

    # Output Header

    print("".center(columns,"*"))
    print("")
    print("EXECUTION OF THE AB INITIO INPUT BUILDER & JOB LAUNCHER BEGINS NOW".center(columns))
    print("")
    print("".center(columns,"*"))

    # ========================================================= #
    # Read command line arguments                               #
    # ========================================================= #

    args = parser.parse_args()

    # Required arguments

    mol_inp = args.mol_inp                   # Geometry file or directory containing the geometry files
    config_inp = args.config                 # YAML configuration file or directory containing the YAML configuration files

    prog = args.program                      # Name of the program for which files need to be created
    out_dir = args.out_dir                   # Directory where all jobs subdirectories will be created
    cluster_name = args.cluster_name         # Name of the cluster where this script is running, as defined in the clusters configuration YAML file

    # Optional arguments

    overwrite = args.overwrite               # Flag for removing job subdirectories before creating a new one, if they have the same name
    dry_run = args.dry_run                   # Flag to not launch the jobs and just create the files

    max_mol = args.max_mol                   # Maximum number of geometry files that will be treated
    max_cf = args.max_cf                     # Maximum number of configuration files that will be treated
    
    keep_mol = args.keep_mol                 # Flag for keeping the geometry files where they are
    keep_cf = args.keep_cf                   # Flag for keeping the configuration files where they are
    
    # Format of the molecule files

    mol_fmt = "xyz"                          # If you decide to add format as a command line argument, replace "xyz" by args.format
    mol_ext = "." + mol_fmt                  # Extension of the geometry files we're looking for

    # ========================================================= #
    # Define codes directory                                    #
    # ========================================================= #

    # Determined by getting the path to the directory of this script

    code_dir = os.path.dirname(os.path.realpath(os.path.abspath(getsourcefile(lambda:0))))

    print ("{:<40} {:<100}".format('\nCodes directory:',code_dir))

    # ========================================================= #
    # Load Mendeleev's periodic table                           #
    # ========================================================= #

    # Loading AlexGustafsson's Mendeleev Table (found at https://github.com/AlexGustafsson/molecular-data) which will be used notably by the scaling process.

    mendeleev_file = abin_errors.check_abspath(os.path.join(code_dir,"mendeleev.yml"),"Mendeleev periodic table YAML file","file")
    print ("{:<141}".format("\nLoading AlexGustafsson's Mendeleev Table ..."), end="")
    with open(mendeleev_file, 'r') as periodic_table:
      mendeleev = yaml.load(periodic_table, Loader=yaml.FullLoader)
    print('%12s' % "[ DONE ]")

    # ========================================================= #
    # Check and load the YAML clusters configuration file       #
    # ========================================================= #

    # Loading the clusters_file 

    clusters_file = abin_errors.check_abspath(os.path.join(code_dir,"clusters.yml"),"YAML clusters configuration file","file")
    print ("{:<141}".format('\nLoading the clusters configuration file "clusters.yml" ...'), end="")
    with open(clusters_file, 'r') as f_clusters:
      clusters_cfg = yaml.load(f_clusters, Loader=yaml.FullLoader)
    print('%12s' % "[ DONE ]")

    # Check the name of the cluster

    if cluster_name not in clusters_cfg:
      raise abin_errors.AbinError ("ERROR: There is no information about the %s cluster in the clusters configuration file. Please add relevant information or change the cluster before proceeding further." % cluster_name.upper())

    print("\nThis script is running on the %s cluster" % cluster_name.upper())

    # Check if the submit_command key has been defined

    submit_command = clusters_cfg[cluster_name].get("submit_command")

    if submit_command is None:
      raise abin_errors.AbinError ("ERROR: There is no defined submit_command for the %s cluster in the clusters configuration file." % cluster_name.upper()) 

    # Check if the program exists 

    if "progs" not in clusters_cfg[cluster_name]:
      raise abin_errors.AbinError ('ERROR: There is no "progs" key defined for the %s cluster in the clusters configuration file. Consult official documentation for details.' % cluster_name.upper())    

    if prog not in clusters_cfg[cluster_name]["progs"]:
      raise abin_errors.AbinError ("ERROR: The specified program (%s) is unknown on this cluster. Possible programs include: %s \nPlease use one of those, change cluster or add information for this program to the clusters configuration file." % (prog, ', '.join(program for program in clusters_cfg[cluster_name]["progs"].keys())))
    
    # Get the scaling function that will determine the scale_index of the molecule (necessary for determining the job scale) - defined in scaling_fcts.py and specified in the clusters configuration file

    scaling_fct = clusters_cfg[cluster_name]["progs"][prog].get("scaling_function")

    if scaling_fct is None:
      raise abin_errors.AbinError ("ERROR: There is no defined scaling function for the %s program in the %s cluster in the clusters configuration file." % (prog, cluster_name.upper()))
    if (scaling_fct) not in dir(scaling_fcts) or not callable(getattr(scaling_fcts, scaling_fct)):
      raise abin_errors.AbinError ("ERROR: There is no scaling function named %s defined in scaling_fcts.py." % scaling_fct)

    # ========================================================= #
    # Establishing the different job scales                     #
    # ========================================================= #

    # Gather all the different job scales from the clusters configuration file in a temporary dictionary

    job_scales_tmp = clusters_cfg[cluster_name]['progs'][prog].get('job_scales')

    if job_scales_tmp is None:
      raise abin_errors.AbinError ("ERROR: There is no defined job_scales for the %s program in the %s cluster in the clusters configuration file." % (prog, cluster_name.upper())) 

    # Defined the required keys in our job scales

    required_keys = frozenset({"label", "scale_limit", "time", "cores", "mem_per_cpu" })

    # Define an array for correct English spelling during printing

    special_numbers = {1:"st", 2:"nd", 3:"rd"}

    # Initialize the final dictionary where the job scales will be sorted by their upper limit

    job_scales = {}

    # Check the job scales

    for scale in job_scales_tmp:

      # Check if all the required keys are present

      for key in required_keys:
        if key not in scale:
          raise abin_errors.AbinError ('ERROR: There is no defined "%s" key for the %s%s job scale of the %s program in the %s cluster in the clusters configuration file.' % (key, job_scales_tmp.index(scale), ("th" if not job_scales_tmp.index(scale) in special_numbers else special_numbers[job_scales_tmp.index(scale)]), prog, cluster_name.upper()))           

      # Extract the scale upper limit from the job scales

      scale_limit = scale.pop('scale_limit')
      job_scales[scale_limit] = scale

    # Sort the different job scales by their upper limit and store them in the job_scales dictionary

    job_scales = OrderedDict(sorted(job_scales.items()))

    print("\nJob scales for %s on %s:" % (prog,cluster_name.upper()))
    print("")
    print(''.center(146, '-'))
    print ("{:<15} {:<20} {:<20} {:<20} {:<10} {:<20} {:<40}".format('Scale Limit','Label','Partition Name','Time','Cores','Mem per CPU (MB)','Delay Command'))
    print(''.center(146, '-'))
    for scale_limit, scale in job_scales.items():
      print ("{:<15} {:<20} {:<20} {:<20} {:<10} {:<20} {:<40}".format(scale_limit, scale['label'], scale.get('partition_name', "not specified"), scale['time'], scale['cores'], scale['mem_per_cpu'], scale.get('delay_command', "not specified")))
    print(''.center(146, '-'))

    # ========================================================= #
    # Check other important functions                           #
    # ========================================================= #

    # Define the scanning function that will extract information about the molecule from the molecule file (depends on the file format) - defined in geom_scan.py

    scan_fct = mol_fmt + "_scan"

    if (scan_fct) not in dir(geom_scan) or not callable(getattr(geom_scan, scan_fct)):
      raise abin_errors.AbinError ("ERROR: There is no function defined for the %s format in geom_scan.py." % mol_fmt)

    # Define the rendering function that will render the job script and the input file (depends on the program)  - defined in renderer.py

    render_fct = prog + "_render"

    if (render_fct) not in dir(renderer) or not callable(getattr(renderer, render_fct)):
      raise abin_errors.AbinError ("ERROR: There is no function defined for the %s program in renderer.py." % prog)

    # ========================================================= #
    # Check other arguments                                     #
    # ========================================================= #

    out_dir = abin_errors.check_abspath(out_dir,"Command line argument -o / --out_dir","directory")
    print ("{:<40} {:<100}".format('\nJobs main directory:',out_dir))

    if max_mol != None and max_mol <= 0:
      raise abin_errors.AbinError ("ERROR: The specified max_mol value (%s) must be a non-zero positive integer" % max_mol)

    if max_cf != None and max_cf <= 0:
      raise abin_errors.AbinError ("ERROR: The specified max_cf value (%s) must be a non-zero positive integer" % max_cf)

    # Check geometry file(s)
    # ======================

    mol_inp = abin_errors.check_abspath(mol_inp,"Command line argument -m / --mol_inp")

    # If the argument mol_inp is a directory, we need to look for every geometry file with the given format in that directory.

    if os.path.isdir(mol_inp):

      print("{:<40} {:<100}".format("\nLooking for %s geometry files in" % mol_ext, mol_inp + " ..."), end="")
      
      mol_inp_path = mol_inp

      # Define which type of file we are looking for (*.mol_ext) in a case-insensitive way (see https://gist.github.com/techtonik/5694830)

      rule = re.compile(fnmatch.translate("*" + mol_ext), re.IGNORECASE)

      # Find all matching files in mol_inp directory

      mol_inp_list = [mol_file for mol_file in os.listdir(mol_inp) if rule.match(mol_file)]

      if mol_inp_list == []:
        raise abin_errors.AbinError ("ERROR: Can't find any geometry of the %s format in %s" % (mol_ext,mol_inp_path))

      print('%12s' % "[ DONE ]")

    # If given a single geometry file as argument, check its extension.

    else:

      print ("{:<40} {:<100}".format('\nGeometry file:',mol_inp))

      if os.path.isfile(mol_inp) and os.path.splitext(mol_inp)[-1].lower() != (mol_ext.lower()):
        raise abin_errors.AbinError ("  ^ ERROR: This is not an %s file." % mol_fmt)

      mol_inp_path = os.path.dirname(mol_inp)
      mol_inp_file = os.path.basename(mol_inp)
      mol_inp_list = [mol_inp_file]

    # Check config file(s)
    # ====================

    config_inp = abin_errors.check_abspath(config_inp,"Command line argument -cf / --config")

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
        raise abin_errors.AbinError ("ERROR: Can't find any YAML config file with the .yml or .yaml extension in %s" % config_inp_path)

      print('%12s' % "[ DONE ]")

    # If given a single config file as argument, check its extension.

    else:

      print ("{:<40} {:<100}".format('\nConfiguration file:',config_inp))

      if os.path.isfile(config_inp) and os.path.splitext(config_inp)[-1].lower() != (".yml") and os.path.splitext(config_inp)[-1].lower() != (".yaml"):
        raise abin_errors.AbinError ("  ^ ERROR: This is not a YAML file (YAML file extension is either .yml or .yaml).")

      config_inp_path = os.path.dirname(config_inp)
      config_inp_file = os.path.basename(config_inp)
      config_inp_list = [config_inp_file]

  # ========================================================= #
  # Exception handling for the preparation step               #
  # ========================================================= #

  except abin_errors.AbinError as error:
    print("")
    print(error)
    exit(-1)

  # =================================================================== #
  # =================================================================== #
  #                   FILES MANIPULATION & GENERATION                   #
  # =================================================================== #
  # =================================================================== #

  problem_cf = [] # Empty list that will contain the names of the configuration files for which a problem has occurred (those configuration files will not be archived)

  if not dry_run:
    job_count = 0   # Launched jobs counter, this number will be showed on the console screen at the end of the execution

  # If a maximum number of geometry files has been given, initialize a geometry files counter

  if max_mol:
    mol_count = 0

  for mol_filename in mol_inp_list:

    # ========================================================= #
    # ========================================================= #
    #                Geometry file treatment                    #
    # ========================================================= #
    # ========================================================= #

    # For more information on try/except structures, see https://www.tutorialsteacher.com/python/exception-handling-in-python
    try:

      # Getting rid of the format extension to get the name of the molecule

      mol_name = str(mol_filename.split('.')[0])
      console_message = "Start procedure for the molecule " + mol_name
      print("")
      print(''.center(len(console_message)+11, '*'))
      print(console_message.center(len(console_message)+10))
      print(''.center(len(console_message)+11, '*'))

      # Create an output log file containing all the information about the geometry file treatment

      mol_log_name = mol_name + ".log"
      mol_log = open(os.path.join(out_dir,mol_log_name), 'w', encoding='utf-8')

      # Redirect standard output to the mol_log file (see https://stackabuse.com/writing-to-a-file-with-pythons-print-function/ for reference)

      sys.stdout = mol_log
      
      # ========================================================= #
      # Reading the content of the geometry file                  #
      # ========================================================= #
    
      print("{:<80}".format("\nScanning %s file ..." % mol_filename), end="")
      with open(os.path.join(mol_inp_path,mol_filename), 'r') as mol_file:
        mol_content = mol_file.read().splitlines()

      # Call the scanning function (defined in geom_scan.py, see the documentation for more information)

      file_data = eval("geom_scan." + scan_fct)(mol_content)

      print('%12s' % "[ DONE ]")

      # ========================================================= #
      # Check if all atom types do exist in Mendeleev's table     #
      # ========================================================= #

      for atom in file_data['chemical_formula'].keys():

        # Scan mendeleev looking for the atom symbol. If there is no match, return None thus raise an exception (see https://stackoverflow.com/questions/8653516/python-list-of-dictionaries-search for more details)

        if not next((element for element in mendeleev if element["symbol"] == atom), None):
          raise abin_errors.AbinError ("ERROR: Element %s is not defined in AlexGustafsson's Mendeleev Table YAML file (mendeleev.yml)" % atom)
      
      # ========================================================= #
      # Determining the scale index                               #
      # ========================================================= #
      
      section_title = "1. Scale index determination"

      print("")
      print("")
      print(''.center(len(section_title)+10, '*'))
      print(section_title.center(len(section_title)+10))
      print(''.center(len(section_title)+10, '*'))
      
      # Call the scaling function (defined in scaling_fcts.py, see the documentation for more information)
      
      scale_index = eval("scaling_fcts." + scaling_fct)(mendeleev, file_data)

      print("\nScale index: ", scale_index)
      
      # Pick the adequate job scale

      jobscale = None

      for scale_limit in job_scales:
        if scale_index > scale_limit:
          continue
        else:
          jobscale = job_scales[scale_limit]
          jobscale_limit = scale_limit
          break

      if not jobscale:
        raise abin_errors.AbinError("ERROR: The job scale of this molecule is too big for this cluster (%s). Please change cluster." % cluster_name.upper())
      
      # ========================================================= #
      # Determining the ressources needed for the job             #
      # ========================================================= #
      
      section_title = "2. Calculation requirements"

      print("")
      print("")
      print(''.center(len(section_title)+10, '*'))
      print(section_title.center(len(section_title)+10))
      print(''.center(len(section_title)+10, '*'))
      
      # Obtaining the information associated to our job scale
      
      job_partition = jobscale.get('partition_name')
      job_walltime = jobscale['time']
      job_cores = jobscale['cores']
      job_mem_per_cpu = jobscale['mem_per_cpu']
      delay_command = jobscale.get("delay_command", '')

      print("")
      print(''.center(50, '-'))
      print("{:<20} {:<30}".format("Scale index: ", scale_index))
      print(''.center(50, '-'))
      print("{:<20} {:<30}".format("Cluster: ", cluster_name))
      print("{:<20} {:<30}".format("Job scale: ", jobscale["label"]))
      print("{:<20} {:<30}".format("Job scale limit: ", jobscale_limit))
      print(''.center(50, '-'))
      print("{:<20} {:<30}".format("Job partition: ", (job_partition or "not specified")))
      print("{:<20} {:<30}".format("Job walltime: ", job_walltime))
      print("{:<20} {:<30}".format("Number of cores: ", job_cores))
      print("{:<20} {:<30}".format("Mem per CPU (MB): ", job_mem_per_cpu))
      print("{:<20} {:<30}".format("Delay command: ", ("not specified" if delay_command == '' else delay_command)))
      print(''.center(50, '-'))

      # ========================================================= #
      # End of logging for the geometry file                      #
      # ========================================================= #
  
      sys.stdout = original_stdout                       # Reset the standard output to its original value
      mol_log.close()

    # ========================================================= #
    # Exception handling for the geometry files                 #
    # ========================================================= #

    # In case of an error specific to ABIN LAUNCHER, skip the geometry file (and do not archive it)

    except abin_errors.AbinError as error:
      sys.stdout = original_stdout                       # Reset the standard output to its original value
      print(error)
      print("Skipping %s molecule" % mol_name)
      os.remove(os.path.join(out_dir,mol_log_name))      # Remove the log file since there was a problem
      keep_cf = True                                     # Flag to notify that a problem has occurred with this geometry and to not archive the configuration files
      continue

    # ========================================================= #
    # ========================================================= #
    #               Configuration file treatment                #
    # ========================================================= #
    # ========================================================= #

    keep_mol = args.keep_mol                 # Redefine keep_mol to its original value (this flag will be defined to True if there is a problem with one of the configuration files, as to not archive the geometry file)

    # If a maximum number of configuration files has been given, initialize a configuration files counter

    if max_cf:
      cf_count = 0

    # We are still inside the geometry files "for" loop and we are going to iterate over each configuration file with that geometry

    for config_filename in config_inp_list: 

      try:

        # Getting rid of the format extension to get the name of the configuration

        config_name = str(config_filename.split('.')[0])
        print("{:<80}".format("\nTreating '%s' geometry with '%s' configuration ..." % (mol_name, config_name)), end="")
        
        # Create an output log file for each geometry - config combination, using the mol_log file as a basis

        log_name = mol_name + "_" + config_name + ".log"
        shutil.copy(os.path.join(out_dir,mol_log_name),os.path.join(out_dir,log_name))
        log = open(os.path.join(out_dir,log_name), 'a', encoding='utf-8')

        # Redirect standard output to the log file (see https://stackabuse.com/writing-to-a-file-with-pythons-print-function/ for reference)

        sys.stdout = log
        
        # Check if a directory already exists for that geometry - config combination

        if os.path.exists(os.path.join(out_dir,mol_name + "_" + config_name)) and not overwrite:
          raise abin_errors.AbinError ("ERROR: A directory for the %s geometry with the '%s' configuration already exists in %s !" % (mol_name, config_name, out_dir))

        # ========================================================= #
        # Rendering the needed input files                          #
        # ========================================================= #

        section_title = "3. Generation of the job script and input files"

        print("")
        print("")
        print(''.center(len(section_title)+10, '*'))
        print(section_title.center(len(section_title)+10))
        print(''.center(len(section_title)+10, '*'))

        # Load config file

        print ("{:<80}".format('\nLoading the configuration file %s ...' % config_filename), end="")
        with open(os.path.join(config_inp_path,config_filename), 'r') as f_config:
          config = yaml.load(f_config, Loader=yaml.FullLoader)
        print('%12s' % "[ DONE ]")

        # Get the path to the jinja templates directory (a directory named "templates" in the same directory as this script)
        
        templates_dir = os.path.join(code_dir,"templates")

        # Build a dictionary that will contain all information related to the job

        job_specs = {
          "prog" : prog,
          "scaling_fct" : scaling_fct,
          "scale_index" : scale_index,
          "cluster_name" : cluster_name,
          "scale_label" : jobscale["label"],
          "scale_limit" : jobscale_limit,
          "partition" : job_partition,
          "walltime" : job_walltime,
          "cores" : job_cores,
          "mem_per_cpu" : job_mem_per_cpu
        }

        # Build a dictionary containing the additional variables that might be needed by the rendering function
        # If your rendering function needs anything else, you can add it in this dictionary

        misc = {  
            "code_dir" : code_dir,
            "templates_dir" : templates_dir,
            "mol_name" : mol_name,
            "config_name" : config_name
            }

        # Call the rendering function (defined in renderer.py, see the documentation for more information)

        rendered_content, rendered_script = eval("renderer." + render_fct)(mendeleev, clusters_cfg, config, file_data, job_specs, misc)
        
        # ========================================================= #
        # The end step                                              #
        # ========================================================= #

        section_title = "4. The end step"

        print("")
        print("")
        print(''.center(len(section_title)+10, '*'))
        print(section_title.center(len(section_title)+10))
        print(''.center(len(section_title)+10, '*'))

        # Creating the job subdirectory where the job will be launched and creating all the relevant files in it.
      
        job_dir = os.path.join(out_dir,mol_name + "_" + config_name)

        if os.path.exists(job_dir): # Overwrite was already checked previously, no need to check it again
          shutil.rmtree(job_dir)

        os.makedirs(job_dir)

        print("\nThe %s subdirectory has been created at %s" % (mol_name + "_" + config_name, out_dir))
        
        # Copying the config file and the geometry file into the job subdirectory
        
        shutil.copy(os.path.join(mol_inp_path,mol_filename), job_dir)
        shutil.copy(os.path.join(config_inp_path,config_filename), job_dir)

        print("\nThe files %s and %s have been successfully copied into the subdirectory." % (config_filename, mol_filename))
        print("")

        # Writing the content of each rendered file into its own file with the corresponding filename

        for filename, file_content in rendered_content.items():
          rendered_file_path = os.path.join(job_dir, filename)
          with open(rendered_file_path, "w", encoding='utf-8') as result_file:
            result_file.write(file_content)
          print("The %s file has been created into the subdirectory" % filename)

        # Launch the job
        
        if not dry_run:

          print("{:<80}".format("\nLaunching the job ..."), end="")
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

          print('%12s' % "[ DONE ]")

        # ================================================================== #
        # End of logging for that particular geometry - config combination   #
        # ================================================================== #

        sys.stdout = original_stdout                            # Reset the standard output to its original value
        log.close()                                             # End of logging for the config file
        shutil.move(os.path.join(out_dir,log_name), job_dir)    # Archive the log file in the job subdirectory

        # If a maximum number of configuration files has been given, increase the configuration counter and check if we have reached the max

        if max_cf:
          cf_count += 1
          if cf_count == max_cf:
            print("\nThe maximum number of configuration files has been reached, no more configuration files will be processed.")
            break

      # ========================================================= #
      # Exception handling for the configuration files            #
      # ========================================================= #

      # In case of an error specific to the configuration file, skip it (and do not archive it)

      except abin_errors.AbinError as error:
        sys.stdout = original_stdout                            # Reset the standard output to its original value
        print(error)
        print("Skipping configuration '%s'" % config_name)
        os.remove(os.path.join(out_dir,log_name))               # Remove the log file since there was a problem
        problem_cf.append(config_filename)                      # Add the name of this configuration file to the list as to not archive it
        keep_mol = True                                         # Flag to notify that a problem has occurred and to not archive the geometry file
        continue        

    # ========================================================= #
    # Archiving the geometry file                               #
    # ========================================================= #

    os.remove(os.path.join(out_dir,mol_log_name))               # Remove the geometry log file since we've finished treating this molecule

    # Archive the geometry file if keep_mol has not been set and there was no problem

    if not keep_mol:
      launched_dir = os.path.join(mol_inp_path,"launched")      # Directory where the geometry files will be put after having been treated by this script, it will be created inside the directory where were all the geometry files.
      os.makedirs(launched_dir, exist_ok=True)
      if os.path.exists(os.path.join(launched_dir,mol_filename)):
        os.remove(os.path.join(launched_dir,mol_filename))
      shutil.move(os.path.join(mol_inp_path,mol_filename), launched_dir)
      print("\nGeometry file archived to %s" % launched_dir)

    console_message = "End of procedure for the molecule " + mol_name
    print("")
    print(''.center(len(console_message)+10, '*'))
    print(console_message.center(len(console_message)+10))
    print(''.center(len(console_message)+10, '*'))

    # If a maximum number of geometry files has been given, increase the geometry counter and check if we have reached the max

    if max_mol:
      mol_count += 1
      if mol_count == max_mol:
        print("\nThe maximum number of geometry files has been reached, no more geometry files will be processed.")
        break

  # ========================================================= #
  # Archiving the configuration files                         #
  # ========================================================= #

  # After all geometries have been treated, archive the configuration files if keep_cf has not been set and there was no problem

  if not keep_cf:
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