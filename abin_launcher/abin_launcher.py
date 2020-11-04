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
required.add_argument("-m","--mol_inp", type=str, help="Relative or absolute path to either a geometry file or a directory containing multiple geometry files.", required=True)
required.add_argument('-cf', '--config', type=str, help="Relative or absolute path to either a YAML configuration file or a directory containing multiple YAML configuration files, extension must be .yml or .yaml.", required=True)
required.add_argument("-p","--program", type=str, help="Name of the program you wish to run jobs with, as defined in the YAML clusters configuration file.", required=True)
required.add_argument("-o","--out_dir", type=str, help="Relative or absolute path to the directory where you want to create the subdirectories for each job.", required=True)
required.add_argument('-cl', '--cluster_name', type=str, help="Name of the cluster where this script is running, as defined in the YAML clusters configuration file.", required=True)
#required.add_argument('-f', '--format', type=str, help="Format of the geometry files that need to be read.", required=True) #Uncomment this line if you add new scanning functions to geom_scan.py

optional = parser.add_argument_group('Optional arguments')
optional.add_argument('-h','--help',action='help',default=argparse.SUPPRESS,help='Show this help message and exit.')
optional.add_argument("-ow","--overwrite",action="store_true",help="If a job subdirectory for a geometry-configuration combination already exists, remove it before creating a new one.")
optional.add_argument('--max', type=int, help="Maximum number of geometry or configuration files that will be treated.")
optional.add_argument("-km","--keep_mol",action="store_true",help="Do not archive the geometry files after they have been processed and leave them where they are.")
optional.add_argument("-kc","--keep_cf",action="store_true",help="Do not archive the configuration files after they have been processed and leave them where they are.")
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

  # For more informations on try/except structures, see https://www.tutorialsteacher.com/python/exception-handling-in-python
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
    max = args.max                           # Maximum number of geometry or configuration files that will be treated
    keep_mol = args.keep_mol                 # Flag for keeping the geometry files where they are
    keep_cf = args.keep_cf                   # Flag for keeping the configuration files where they are
    dry_run = args.dry_run                   # Flag to not launch the jobs and just create the files

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
    # Check and load the YAML clusters configuration file       #
    # ========================================================= #

    # Loading the clusters_file for the information about the clusters

    clusters_file = abin_errors.check_abspath(os.path.join(code_dir,"clusters.yml"),"YAML clusters configuration file","file")
    print ("{:<40} {:<100}".format('\nLoading the clusters file',clusters_file + " ..."), end="")
    with open(clusters_file, 'r') as f_clusters:
      clusters_cfg = yaml.load(f_clusters, Loader=yaml.FullLoader)
    print('%12s' % "[ DONE ]")

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
    # Check the name of the cluster                             #
    # ========================================================= #

    if cluster_name not in clusters_cfg:
      raise abin_errors.AbinError ("ERROR: There is no information about the %s cluster in the %s file. Please add relevant information or change the cluster before proceeding further." % (cluster_name.upper(),clusters_file))

    print("\nThis script is running on the %s cluster" % cluster_name.upper())

    # ========================================================= #
    # Check program and subfunctions                            #
    # ========================================================= #

    # Check if the program exists in our clusters database. 

    if prog not in clusters_cfg[cluster_name]["progs"]:
      raise abin_errors.AbinError ("ERROR: The specified program (%s) is unknown on this cluster. Possible programs include: %s \nPlease use one of those, change cluster or add information for this program to the YAML cluster file." % (prog, ', '.join(program for program in clusters_cfg[cluster_name]["progs"].keys())))

    # Define the scanning function that will extract informations about the molecule from the molecule file (depends on the file format) - defined in geom_scan.py

    scan_fct = mol_fmt + "_scan"

    if (scan_fct) not in dir(geom_scan) or not callable(getattr(geom_scan, scan_fct)):
      raise abin_errors.AbinError ("ERROR: There is no function defined for the %s format in geom_scan.py." % mol_fmt)

    # Define the scaling function that will determine the scale_index of the molecule (necessary for determining the job scale) - defined in scaling_fcts.py

    scaling_fct = clusters_cfg[cluster_name]["progs"][prog]["scaling_function"]

    if (scaling_fct) not in dir(scaling_fcts) or not callable(getattr(scaling_fcts, scaling_fct)):
      raise abin_errors.AbinError ("ERROR: There is no scaling function named %s defined in scaling_fcts.py." % scaling_fct)

    # Define the rendering function that will render the job instructions file and the input file (depends on the program)  - defined in renderer.py

    render_fct = prog + "_render"

    if (render_fct) not in dir(renderer) or not callable(getattr(renderer, render_fct)):
      raise abin_errors.AbinError ("ERROR: There is no function defined for the %s program in renderer.py." % prog)

    # ========================================================= #
    # Establishing the different job scales                     #
    # ========================================================= #

    # Gather all the different job scales from the clusters configuration file in a temporary dictionary

    job_scales_tmp = clusters_cfg[cluster_name]['progs'][prog]['job_scales']

    # Initialize the final dictionary where the job scales will be sorted by their upper limit

    job_scales = {}

    # Sort the different job scales by their upper limit and store them in the job_scales dictionary

    for scale in job_scales_tmp:
      scale_limit = scale['scale_limit']
      del scale['scale_limit']
      job_scales[float(scale_limit)] = scale

    print("\nJob scales for %s on %s:" % (prog,cluster_name.upper()))
    job_scales = OrderedDict(sorted(job_scales.items()))

    print("")
    print(''.center(106, '-'))
    print ("{:<15} {:<20} {:<20} {:<20} {:<10} {:<20}".format('Scale Limit','Label','Partition Name','Time','Cores','Mem per CPU (MB)'))
    print(''.center(106, '-'))
    for scale_limit, scale in job_scales.items():
      print ("{:<15} {:<20} {:<20} {:<20} {:<10} {:<20}".format(scale_limit, scale['label'], scale['partition_name'], scale['time'], scale['cores'], scale['mem_per_cpu']))
    print(''.center(106, '-'))

    # ========================================================= #
    # Check other arguments                                     #
    # ========================================================= #

    out_dir = abin_errors.check_abspath(out_dir,"Command line argument -o / --out_dir","directory")
    print ("{:<40} {:<100}".format('\nJobs main directory:',out_dir))

    if max != None and max <= 0:
      raise abin_errors.AbinError ("ERROR: The specified max value (%s) must be a non-zero positive integer" % max)

    # Check geometry file(s)
    # ======================

    mol_inp = abin_errors.check_abspath(mol_inp,"Command line argument -m / --mol_inp")

    if os.path.isdir(mol_inp):

      # If the argument mol_inp is a directory, we need to look for every geometry file with the given format in that directory.

      print("{:<40} {:<100}".format("\nLooking for %s geometry files in" % mol_ext, mol_inp + " ..."), end="")
      mol_inp_path = mol_inp

      # Define which type of file we are looking for in a case-insensitive way (see https://gist.github.com/techtonik/5694830)

      rule = re.compile(fnmatch.translate("*." + mol_fmt), re.IGNORECASE)

      # Find all matching files in mol_inp directory

      mol_inp_list = [mol for mol in os.listdir(mol_inp) if rule.match(mol)]
      if mol_inp_list == []:
        raise abin_errors.AbinError ("ERROR: Can't find any geometry of the %s format in %s" % (mol_ext,mol_inp_path))
      if max:
        mol_inp_list = mol_inp_list[0:max]
      print('%12s' % "[ DONE ]")

    else:

      print ("{:<40} {:<100}".format('\nGeometry file:',mol_inp))

      # If given a single geometry file as argument, check its extension.

      if os.path.isfile(mol_inp) and os.path.splitext(mol_inp)[-1].lower() != (mol_ext.lower()):
        raise abin_errors.AbinError ("  ^ ERROR: This is not an %s file." % mol_fmt)
      mol_inp_path = os.path.dirname(mol_inp)
      mol_inp_file = os.path.basename(mol_inp)
      mol_inp_list = [mol_inp_file]

    # Check config file(s)
    # ====================

    config_inp = abin_errors.check_abspath(config_inp,"Command line argument -cf / --config")

    if os.path.isdir(config_inp):

      # If the argument config_inp is a directory, we need to look for every YAML configuration file in that directory.

      print("{:<40} {:<100}".format("\nLooking for .yml or .yaml files in", config_inp + " ..."), end="")
      config_inp_path = config_inp

      # Define which type of file we are looking for in a case-insensitive way (see https://gist.github.com/techtonik/5694830)

      rule = re.compile(fnmatch.translate("*.yml"), re.IGNORECASE)
      rule2 = re.compile(fnmatch.translate("*.yaml"), re.IGNORECASE)

      # Find all matching files in config_inp directory

      config_inp_list = [config for config in os.listdir(config_inp) if (rule.match(config) or rule2.match(config))]
      if config_inp_list == []:
        raise abin_errors.AbinError ("ERROR: Can't find any YAML config file with the .yml or .yaml extension in %s" % config_inp_path)
      if max:
        config_inp_list = config_inp_list[0:max]
      print('%12s' % "[ DONE ]")

    else:

      print ("{:<40} {:<100}".format('\nConfiguration file:',config_inp))

      # If given a single config file as argument, check its extension.

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

  for mol_filename in mol_inp_list:

    # ========================================================= #
    # ========================================================= #
    #                Geometry file treatment                    #
    # ========================================================= #
    # ========================================================= #

    # For more informations on try/except structures, see https://www.tutorialsteacher.com/python/exception-handling-in-python
    try:

      # Getting rid of the format extension to get the name of the molecule

      mol_name = str(mol_filename.split('.')[0])
      console_message = "Start procedure for the molecule " + mol_name
      print("")
      print(''.center(len(console_message)+11, '*'))
      print(console_message.center(len(console_message)+10))
      print(''.center(len(console_message)+11, '*'))

      # Create an output log file containing all the information about the molecule treatment

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

        # Scan mendeleev looking for the atom symbol. If there is no match, returns None thus raises an exception (see https://stackoverflow.com/questions/8653516/python-list-of-dictionaries-search for more details)

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
      
      # Obtaining the informations associated to our job scale
      
      job_partition = jobscale['partition_name']
      job_walltime = jobscale['time']
      job_cores = jobscale['cores']
      job_mem_per_cpu = jobscale['mem_per_cpu']

      print("")
      print(''.center(50, '-'))
      print("{:<20} {:<30}".format("Scale index: ", scale_index))
      print(''.center(50, '-'))
      print("{:<20} {:<30}".format("Cluster: ", cluster_name))
      print("{:<20} {:<30}".format("Job scale: ", jobscale["label"]))
      print("{:<20} {:<30}".format("Job scale limit: ", jobscale_limit))
      print(''.center(50, '-'))
      print("{:<20} {:<30}".format("Job partition: ", job_partition))
      print("{:<20} {:<30}".format("Job walltime: ", job_walltime))
      print("{:<20} {:<30}".format("Number of cores: ", job_cores))
      print("{:<20} {:<30}".format("Mem per CPU (MB): ", job_mem_per_cpu))
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

        section_title = "3. Generation of the job instructions file and input files"

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
        
        path_tpl_dir = os.path.join(code_dir,"templates")

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
            "path_tpl_dir" : path_tpl_dir,
            "mol_name" : mol_name,
            "config_name" : config_name
            }

        # Call the rendering function (defined in renderer.py, see the documentation for more information)

        rendered_content = eval("renderer." + render_fct)(mendeleev, clusters_cfg, config, file_data, job_specs, misc)
        
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

          print("\nLaunching the job ...", end="")
          os.chdir(job_dir)

          # Define the launch command

          subcommand = clusters_cfg[cluster_name]['subcommand']
          delay_command = str(jobscale.get("delay_command") or '')   # If delay_command is not defined in the job scales, an empty string is put in its place.
          job_inst = clusters_cfg[cluster_name]['progs'][prog]['job_instructions']
          launch_command = subcommand + " " + delay_command + " " + job_inst

          # Execute the command and get the command status

          retcode = os.system(launch_command)
          
          # If something went wrong when submitting the job, do not raise an exception and just quit the execution. It is likely a problem linked to the cluster.

          if retcode != 0 :
            sys.stdout = original_stdout                             # Reset the standard output to its original value
            print("Job submit encountered an issue")
            print("Aborting ...")
            exit(5)
        
          print('%12s' % "[ DONE ]")

        # ================================================================== #
        # End of logging for that particular geometry - config combination   #
        # ================================================================== #

        sys.stdout = original_stdout                            # Reset the standard output to its original value
        log.close()                                             # End of logging for the config file
        shutil.move(os.path.join(out_dir,log_name), job_dir)    # Archive the log file in the job subdirectory

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