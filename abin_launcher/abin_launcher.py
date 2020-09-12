#!/usr/bin/env python3

################################################################################################################################################
##                                                 Ab Initio Input Builder & Job Launcher                                                     ##
##                                                                                                                                            ##
##                             For one or more molecule files and a given ab initio program, this script prepares                             ##
##                       the input files needed for each calculation and launches the corresponding jobs on the cluster.                      ##
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

import mol_scan
import renderer
import scaling_fcts

# ===================================================================
# ===================================================================
# Function definitions
# ===================================================================
# ===================================================================

def check_abspath(path,type="either"):
    """Checks if a path towards a file or folder exists and makes sure it's absolute.

    Parameters
    ----------
    path : str
        The path towards the file or directory you want to test
    type : str (optional)
        The type of element for which you would like to test the path (file, folder or either)
        By default, checks if the path leads to either a file or a folder (type = either)
    
    Returns
    -------
    abspath : str
        Normalized absolutized version of the path
    """

    if type not in ["file","folder","either"]:
      # Not in try/except structure because the full error message will be need in this case
      raise ValueError ("The specified type for which the check_abspath function has been called is not one of 'file', 'folder or 'either'")

    # For more informations on try/except structures, see https://www.tutorialsteacher.com/python/exception-handling-in-python
    try:
      if not os.path.exists(path):
        raise IOError ("ERROR: The argument %s does not seem to exist." % path)
      elif type == "file":
        if not os.path.isfile(path):
          raise ValueError ("ERROR: The argument %s is not a file" % path)
      elif type == "folder":
        if not os.path.isdir(path):
          raise ValueError ("ERROR: The argument %s is not a directory" % path)
      elif type == "either":
        if not os.path.isdir(path) and not os.path.isfile(path):
          raise ValueError ("ERROR: The argument %s is neither a file nor a directory" % path)
    except Exception as error:
      print(error)
      exit(1)

    # If everything went well, get the normalized absolutized version of the path
    abspath = os.path.abspath(path)

    return abspath

# ===================================================================
# ===================================================================
# Command line arguments
# ===================================================================
# ===================================================================

# Define the arguments needed for the script (here they are defined as named arguments rather than positional arguments, check https://stackoverflow.com/questions/24180527/argparse-required-arguments-listed-under-optional-arguments for more info).

parser = argparse.ArgumentParser(add_help=False, description="For one or more molecule files and a given ab initio program, this script prepares the input files needed for each calculation and launches the corresponding job on the cluster.")

required = parser.add_argument_group('Required arguments')
required.add_argument("-p","--program", type=str, help="Program you wish to run jobs with, must be the same as the name mentioned in the clusters and configuration YAML files", required=True)
required.add_argument("-m","--mol_inp", type=str, help="Path to either a molecule file or a directory containing multiple molecule files", required=True)
required.add_argument('-cf', '--config', type=str, help="Path to either a YAML config file or a directory containing multiple YAML config files", required=True)
required.add_argument("-o","--out_dir", type=str, help="Path to the directory where you want to create the subdirectories for each job", required=True)

optional = parser.add_argument_group('Optional arguments')
optional.add_argument('-h','--help',action='help',default=argparse.SUPPRESS,help='Show this help message and exit')
optional.add_argument("-ow","--overwrite",action="store_true",help="Overwrite files if they already exists")
optional.add_argument("-k","--keep",action="store_true",help="Do not archive launched molecule files and leave them where they are")
optional.add_argument('-cl', '--clusters', type=str, help="Path to the YAML clusters file, default is this_script_directory/clusters.yml")

args = parser.parse_args()

# Define the variables corresponding to those arguments

prog = args.program                      # Name of the program for which files need to be created
mol_inp = args.mol_inp                   # Molecule file or folder containing the molecule files
config_inp = args.config                 # YAML configuration file or folder containing the YAML configuration files
out_dir = args.out_dir                   # Folder where all jobs subfolders will be created

overwrite = args.overwrite               # Flag for overwriting the files
keep = args.keep                         # Flag for keeping the molecule files where they are
clusters_file = args.clusters            # YAML file containing all informations about the clusters

# Other important variables that could become arguments if the need arises

mol_fmt = "xyz"                                 # Format of the molecule files we want to treat
mol_ext = "." + mol_fmt                         # Extension of the molecule files we're looking for

# ===================================================================
# ===================================================================
#                           PREPARATION STEP
# ===================================================================
# ===================================================================

# Save a reference to the original standard output as it will be modified later on (see https://stackabuse.com/writing-to-a-file-with-pythons-print-function/ for reference)

original_stdout = sys.stdout

# Get the size of the terminal in order to have a prettier output, if you need something more robust, go check http://granitosaurus.rocks/getting-terminal-size.html

columns, rows = shutil.get_terminal_size()

# Output Header

print("".center(columns,"*"))
print("")
print("EXECUTION OF THE AB INITIO INPUT BUILDER & JOB LAUNCHER BEGINS NOW".center(columns))
print("")
print("".center(columns,"*"))

# =========================================================
# Define the cluster
# =========================================================

cluster_name = os.environ['CLUSTER_NAME']
print("\nThis script is running on the %s cluster" % cluster_name.upper())

# =========================================================
# Define codes directory
# =========================================================

# Codes directory (determined by getting the path to the directory where this script is)
code_dir = os.path.dirname(os.path.realpath(os.path.abspath(getsourcefile(lambda:0))))
print ("{:<40} {:<100}".format('\nCodes directory:',code_dir))

# =========================================================
# Check molecule file(s)
# =========================================================

mol_inp = check_abspath(mol_inp)

if os.path.isdir(mol_inp):
  # If the argument mol_inp is a folder, we need to look for every molecule file with the given format in that folder.
  print("{:<40} {:<100}".format("\nLooking for %s molecule files in" % mol_ext, mol_inp + " ..."), end="")
  mol_inp_path = mol_inp
  # Define which type of file we are looking for in a case-insensitive way (see https://gist.github.com/techtonik/5694830)
  rule = re.compile(fnmatch.translate("*." + mol_fmt), re.IGNORECASE)
  # Find all matching files in mol_inp folder
  mol_inp_list = [mol for mol in os.listdir(mol_inp) if rule.match(mol)]
  if mol_inp_list == []:
    print("\nERROR: Can't find any molecule of the %s format in %s" % (mol_ext,mol_inp_path))
    exit(1)
  print('%12s' % "[ DONE ]")

else:
  print ("{:<40} {:<100}".format('\nMolecule file:',mol_inp))
  # If given a single molecule file as argument, check its extension.
  if os.path.isfile(mol_inp) and os.path.splitext(mol_inp)[-1].lower() != (mol_ext.lower()):
    print("  ^ ERROR: This is not an %s file." % mol_fmt)
    exit(1)
  mol_inp_path = os.path.dirname(mol_inp)
  mol_inp_file = os.path.basename(mol_inp)
  mol_inp_list = [mol_inp_file]
  
# =========================================================
# Check config file(s)
# =========================================================

config_inp = check_abspath(config_inp)

if os.path.isdir(config_inp):
  # If the argument config_inp is a folder, we need to look for every YAML configuration file in that folder.
  print("{:<40} {:<100}".format("\nLooking for .yml or .yaml files in", config_inp + " ..."), end="")
  config_inp_path = config_inp
  # Define which type of file we are looking for in a case-insensitive way (see https://gist.github.com/techtonik/5694830)
  rule = re.compile(fnmatch.translate("*.yml"), re.IGNORECASE)
  rule2 = re.compile(fnmatch.translate("*.yaml"), re.IGNORECASE)
  # Find all matching files in mol_inp folder
  config_inp_list = [config for config in os.listdir(config_inp) if (rule.match(config) or rule2.match(config))]
  if config_inp_list == []:
    print("\nERROR: Can't find any YAML config file with the .yml or .yaml extension in %s" % config_inp_path)
    exit(1)
  print('%12s' % "[ DONE ]")

else:
  print ("{:<40} {:<100}".format('\nConfiguration file:',config_inp))
  # If given a single config file as argument, check its extension.
  if os.path.isfile(config_inp) and os.path.splitext(config_inp)[-1].lower() != (".yml") and os.path.splitext(config_inp)[-1].lower() != (".yaml"):
    print("  ^ ERROR: This is not a YAML file (YAML file extension is either .yml or .yaml).")
    exit(1)
  config_inp_path = os.path.dirname(config_inp)
  config_inp_file = os.path.basename(config_inp)
  config_inp_list = [config_inp_file]

# =========================================================
# Check other arguments (program will be checked later)
# =========================================================

out_dir = check_abspath(out_dir,"folder")
print ("{:<40} {:<100}".format('\nJobs main directory:',out_dir))

if clusters_file: 
  clusters_file = check_abspath(clusters_file,"file")
else:
  # If no value has been provided through the command line, take the clusters.yml file in the same directory as this script 
  clusters_file = os.path.join(code_dir,"clusters.yml")

# =========================================================
# Load YAML files (except config)
# =========================================================

# Loading the clusters_file for the informations about the clusters

print ("{:<40} {:<100}".format('\nLoading the clusters file',clusters_file + " ..."), end="")
with open(clusters_file, 'r') as f_clusters:
  clusters_cfg = yaml.load(f_clusters, Loader=yaml.FullLoader)
print('%12s' % "[ DONE ]")

# Loading AlexGustafsson's Mendeleev Table (found at https://github.com/AlexGustafsson/molecular-data) which will be used for the scaling functions.

print ("{:<140}".format("\nLoading AlexGustafsson's Mendeleev Table ..."), end="")
with open(os.path.join(code_dir,'elements.yml'), 'r') as f_elements:
  elements = yaml.load(f_elements, Loader=yaml.FullLoader)
print('%12s' % "[ DONE ]")

# =========================================================
# Define the name of our subfunctions
# =========================================================

# Name of the scanning function that will extract informations about the molecule from the molecule file (depends on the file format) - defined in mol_scan.py
scan_fct = mol_fmt + "_scan"

# Name of the scaling function that will determine the scale_index of the molecule (necessary for determining the job scale) - defined in scaling_fcts.py
scaling_fct = clusters_cfg["progs"][prog]["scaling_function"]

# Name of the render function that will render the job manifest and the input file (depends on the program)  - defined in renderer.py
render_fct = prog + "_render"

# =========================================================
# Check program and subfunctions
# =========================================================

# Check if the program exists in our clusters database. 

if prog not in clusters_cfg[cluster_name]["progs"]:
  print("\nERROR: Program unknown on this cluster. Possible program(s) include" , ', '.join(program for program in clusters_cfg[cluster_name]["progs"].keys()))
  print("Please use one of those, change cluster or add informations for this program to the YAML cluster file.")
  exit(3) 

# Check if a rendering function has been defined for the program

if (render_fct) not in dir(renderer) or not callable(getattr(renderer, render_fct)):
  print("\nERROR: There is no function defined for the %s program in renderer.py." % prog)
  exit(3) 

# Check if a scanning function has been defined for the given format of the molecule file(s)

if (scan_fct) not in dir(mol_scan) or not callable(getattr(mol_scan, scan_fct)):
  print("\nERROR: There is no function defined for the %s format in mol_scan.py." % mol_fmt)
  exit(3) 

# Check if the chosen scaling function has been defined in scaling_fcts.py

if (scaling_fct) not in dir(scaling_fcts) or not callable(getattr(scaling_fcts, scaling_fct)):
  print("\nERROR: There is no scaling function named %s defined in scaling_fcts.py." % scaling_fct)
  exit(3) 

# =========================================================
# Establishing the different job scales
# =========================================================

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
print(''.center(105, '-'))
print ("{:<15} {:<20} {:<20} {:<20} {:<10} {:<20}".format('scale_limit','Label','Partition_name','Time','Cores','Mem per CPU (MB)'))
print(''.center(105, '-'))
for scale_limit, scale in job_scales.items():
  print ("{:<15} {:<20} {:<20} {:<20} {:<10} {:<20}".format(scale_limit, scale['label'], scale['partition_name'], scale['time'], scale['cores'], scale['mem_per_cpu']))
print(''.center(105, '-'))

# ===================================================================
# ===================================================================
#                   FILES MANIPULATION & GENERATION
# ===================================================================
# ===================================================================

for mol_filename in mol_inp_list:

  # Getting rid of the format extension to get the name of the molecule and the configuration
  mol_name = str(mol_filename.split('.')[0])
  print("Now treating %s molecule" % mol_name)

  # For more informations on try/except structures, see https://www.tutorialsteacher.com/python/exception-handling-in-python
  try:

    for config_filename in config_inp_list:

      # Getting rid of the format extension to get the name of the configuration
      config_name = str(config_inp_file.split('.')[0])
      print("{:<80}".format("    Treating %s molecule with %s configuration ..." % (mol_name, config_name)), end="")
      
      # Create an output log file for each molecule - config combination
      log_name = mol_name + "_" + config_name + ".log"
      log = open(os.path.join(out_dir,log_name), 'x', encoding='utf-8')

      # Redirect standard output to the log file (see https://stackabuse.com/writing-to-a-file-with-pythons-print-function/ for reference)
      sys.stdout = log

      # Check if a folder already exists for that molecule - config combination
      if os.path.exists(os.path.join(out_dir,mol_name + "_" + config_name)) and not overwrite:
        raise ValueError ("ERROR: A folder for the %s molecule with the %s configuration already exists in %s !" % (mol_name, config_name, out_dir))
      
      # =========================================================
      # Load config file
      # =========================================================

      print ("{:<40} {:<100}".format('\nLoading the configuration file',config_filename + " ..."), end="")
      #print("\nLoading the main configuration file  %s ..." % config_file, end="")
      with open(os.path.join(config_inp_path,config_filename), 'r') as f_config:
        config = yaml.load(f_config, Loader=yaml.FullLoader)
      print('%12s' % "[ DONE ]")

      # Check if the program exists in the config file. 

      if prog not in config:
        raise ValueError ("ERROR: No information provided for this program in the YAML config file")

      # =========================================================
      # Molecule file treatment
      # =========================================================
      
      # Reading the content of the molecule file
      
      print("{:<80}".format("\nScanning %s file ..." % mol_filename), end="")
      with open(os.path.join(mol_inp_path,mol_filename), 'r') as mol_file:
        mol_content = mol_file.read().splitlines()

      # Scanning the content of the file and extracting the relevant informations

      file_data = eval("mol_scan." + scan_fct)(mol_content)

      #! The file_data variable is a dictionary of the form {'chemical_formula':{}, 'atomic_coordinates':[]}
      #! The first key of file_data is a dictionary stating the chemical formula of the molecule in the form {'atom type 1':number of type 1 atoms, 'atom type 2':number of type 2 atoms, ...}, ex: {'Si':17, 'O':4, 'H':28}
      #! The second key is a list containing all atomic coordinates, as they will be used in the input file of the ab initio program
      #! If a problem arises when scanning the molecule file, a ValueError exception should be raised with a proper error message

      print('%12s' % "[ DONE ]")
      
      # =========================================================
      # Determining the scale_index
      # =========================================================
      
      section_title = "1. Scale index determination"

      print("")
      print("")
      print(''.center(len(section_title)+10, '*'))
      print(section_title.center(len(section_title)+10))
      print(''.center(len(section_title)+10, '*'))
      
      # scale_index determination
      
      scale_index = eval("scaling_fcts." + scaling_fct)(elements, file_data)

      print(("\nScale index: ", scale_index))
      
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
        raise ValueError("ERROR: This molecule job scale is too big for this cluster (%s). Please change cluster." % cluster_name.upper())
      
      # =========================================================
      # Determining the ressources needed for the job
      # =========================================================
      
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

      # =========================================================
      # Rendering the needed input files
      # =========================================================

      section_title = "3. Generation of the job manifest and input files"

      print("")
      print("")
      print(''.center(len(section_title)+10, '*'))
      print(section_title.center(len(section_title)+10))
      print(''.center(len(section_title)+10, '*'))
      
      # Get the path to jinja templates folder (a folder named "Templates" in the samed folder as this script)

      path_tpl_dir = os.path.join(code_dir,"Templates")

      # Dynamically call the inputs render function for the given program

      rendered_content = eval("renderer." + render_fct)(locals())  # Dictionary containing the text of all the rendered files in the form of <filename>: <rendered_content>
      
      # =========================================================
      # The end step
      # =========================================================

      section_title = "4. The end step"

      print("")
      print("")
      print(''.center(len(section_title)+10, '*'))
      print(section_title.center(len(section_title)+10))
      print(''.center(len(section_title)+10, '*'))

      # Creating the molecule subfolder where the job will be launched and creating all the relevant files in it.
    
      job_dir = os.path.join(out_dir,mol_name + "_" + config_name)

      if os.path.exists(job_dir):
        shutil.rmtree(job_dir)

      os.makedirs(job_dir)

      print("\nThe %s subfolder has been created at %s" % (mol_name + "_" + config_name, out_dir))
      
      # Copying the config file and the molecule file into the molecule subfolder
      
      shutil.copy(os.path.join(mol_inp_path,mol_filename), job_dir)
      shutil.copy(os.path.join(config_inp_path,config_filename), job_dir)

      print("\nThe files %s and %s have been successfully copied into the subfolder." % (config_filename, mol_filename))
      print("")

      # Writing the content of each rendered file into its own file with the corresponding filename

      for filename, file_content in rendered_content.items():
        rendered_file_path = os.path.join(job_dir, filename)
        with open(rendered_file_path, "w", encoding='utf-8') as result_file:
          result_file.write(file_content)
        print("The %s file has been created into the subfolder" % filename)

      # Launch the job
      
      print("\nLaunching the job ...")
      os.chdir(job_dir)
      subcommand = clusters_cfg[cluster_name]['subcommand']
      delay_command = str(jobscale["delay_command"] or '')
      manifest = config[prog]['rendered_files']['manifest']
      launch_command = subcommand + " " + delay_command + " " + manifest
      retcode = os.system(launch_command)
      if retcode != 0 :
        sys.stdout = original_stdout # Reset the standard output to its original value
        print("Job submit encountered an issue")
        print("Aborting ...")
        exit(5)
      
      # End of treatment for that particular molecule - config combination

      sys.stdout = original_stdout # Reset the standard output to its original value
      shutil.move(os.path.join(out_dir,log_name), job_dir) # Archive the log file in the job subfolder
      print('%12s' % "[ DONE ]")

  except ValueError as error:
    sys.stdout = original_stdout # Reset the standard output to its original value
    print(error)
    print("Skipping %s molecule" % mol_name)
    os.remove(os.path.join(out_dir,log_name))
    continue

  # Archive the molecule file in launched_dir if keep has not been set

  if not keep:
    launched_dir = os.path.join(mol_inp_path,"Launched")                             # Folder where the molecule files will be put after having been treated by this script, path is relative to the directory where are all the molecule files.
    os.makedirs(launched_dir, exist_ok=True)
    if os.path.exists(os.path.join(launched_dir,mol_filename)):
      os.remove(os.path.join(launched_dir,mol_filename))
    shutil.move(os.path.join(mol_inp_path,mol_filename), launched_dir)
    print("    Molecule original structure file archived to %s" % launched_dir)

  print("End of procedure for the molecule %s" % mol_name)
  print("")

print("")
print("".center(columns,"*"))
print("")
print("END OF EXECUTION".center(columns))
print("")
print("".center(columns,"*"))