################################################################################################################################################
##                                                                The Renderer                                                                ##
##                                                                                                                                            ##
##                                       This script contains the rendering functions for ABIN LAUNCHER,                                      ##
##                                consult the documentation at https://chains-ulb.readthedocs.io/ for details                                 ##
################################################################################################################################################

import os

import yaml
from jinja2 import Environment, FileSystemLoader

import abin_errors


def jinja_render(templates_dir:str, template_file:str, render_vars:dict):
    """Renders a file based on its Jinja template.

    Parameters
    ----------
    templates_dir : str
        The path towards the directory where the Jinja template is located.
    template_file : str
        The name of the Jinja template file.
    render_vars : dict
        Dictionary containing the definitions of all the variables present in the Jinja template.

    Returns
    -------
    output_text : str
        Content of the rendered file.
    """
   
    file_loader = FileSystemLoader(templates_dir)
    env = Environment(loader=file_loader)
    template = env.get_template(template_file)
    output_text = template.render(render_vars)
    
    return output_text

# =================================================================== #
# =================================================================== #
#                         Rendering functions                         #
# =================================================================== #
# =================================================================== #

def chains_orca_render(mendeleev:dict, clusters_cfg:dict, config:dict, file_data:dict, job_specs:dict, misc:dict):
    """Renders the job script and the input file associated with the ORCA program in CHAINS.

    Parameters
    ----------
    mendeleev : dict
        Content of AlexGustafsson's Mendeleev Table YAML file (found at https://github.com/AlexGustafsson/molecular-data).
        Unused in this function.
    clusters_cfg : dict
        Content of the YAML clusters configuration file.
    config : dict
        Content of the YAML configuration file.
    file_data : dict
        Information extracted by the scanning function from the geometry file.
    job_specs : dict
        Contains all information related to the job.
    misc : dict
        Contains all the additional variables that did not pertain to the other arguments.

    Returns
    -------
    rendered_content : dict
        Dictionary containing the text of all the rendered files in the form of <filename>: <rendered_content>.
    rendered_script : str
        Name of the rendered job script, necessary to launch the job.
    
    Notes
    -----
    Pay a particular attention to the render_vars dictionaries, they contain all the definitions of the variables appearing in your Jinja templates.
    """

    # ========================================================= #
    #                      Preparation step                     #
    # ========================================================= #

    # Check config file
    # =================

    # Check if a "general" block has been defined in the config file

    if not config.get('general'):
      raise abin_errors.AbinError ('ERROR: There is no "general" key defined in the "%s" configuration file.' % misc['config_name'])     

    # Check if an "orca" block has been defined in the config file

    if not config.get('orca'):
      raise abin_errors.AbinError ('ERROR: There is no "orca" key defined in the "%s" configuration file.' % misc['config_name'])      

    # Check the options defined in the config file

    pre_opt = config['orca'].get('pre_opt',False)

    if not isinstance(pre_opt, bool):
      raise abin_errors.AbinError ('ERROR: The "pre_opt" value given in the "orca" block of the "%s" configuration file is not a boolean (neither "True" nor "False").' % misc['config_name'])

    copy_files = config['orca'].get('copy_files',True)

    if not isinstance(copy_files, bool):
      raise abin_errors.AbinError ('ERROR: The "copy_files" value given in the "orca" block of the "%s" configuration file is not a boolean (neither "True" nor "False").' % misc['config_name'])

    benchmark = config['orca'].get('benchmark',False)

    if not isinstance(benchmark, bool):
      raise abin_errors.AbinError ('ERROR: The "benchmark" value given in the "orca" block of the "%s" configuration file is not a boolean (neither "True" nor "False").' % misc['config_name'])

    # Define the templates
    # ====================

    # Define the names of the templates.

    if pre_opt:
      template_input = "orca_preopt.inp.jinja"
    else:
      template_input = "orca.inp.jinja"

    template_script = "orca_job.sh.jinja"

    # Check if the specified templates exist in the "templates" directory of ABIN LAUNCHER.
    
    abin_errors.check_abspath(os.path.join(misc['templates_dir'],template_input),"Jinja template for the orca input file","file")
    abin_errors.check_abspath(os.path.join(misc['templates_dir'],template_script),"Jinja template for the orca job script","file")

    # Define rendered files
    # =====================

    # Define the names of the rendered files.

    rendered_input = misc['mol_name'] + ".inp"
    rendered_script = "orca_job.sh"

    # Initialize the dictionary that will be returned by the function

    rendered_content = {}

    # ========================================================= #
    #                  Rendering the input file                 #
    # ========================================================= #
  
    print("{:<80}".format("\nRendering the jinja template for the orca input file ...  "), end="")

    # Defining the Jinja variables
    # ============================

    # It is recommended to set the memory per CPU to 75% of the physical memory available (see https://sites.google.com/site/orcainputlibrary/orca-common-problems)

    orca_mem_per_cpu = int(0.75 * job_specs['mem_per_cpu']) # in MB

    # Variables not associated with the config file

    input_render_vars = {
      "orca_mem_per_cpu" : orca_mem_per_cpu,
      "job_cores" : job_specs['cores'],
      "coordinates" : file_data['atomic_coordinates']
    }

    # Variables associated with the "general" block of the config file

    try:
      input_render_vars.update({
        "charge" : config['general']['charge'],
        "multiplicity" : config['general']['multiplicity']
      })

    except KeyError as error:
      raise abin_errors.AbinError ('ERROR: The "%s" key is missing in the "general" block of the "%s" configuration file.' % (error,misc['config_name']))

    # Check if a "keywords" block has been defined in the "orca" block of the config file

    if not config['orca'].get('keywords'):
      raise abin_errors.AbinError ('ERROR: There is no "keywords" key in the "orca" block of the "%s" configuration file.' % misc['config_name'])    

    # Variables associated with the "keywords" block of the "orca" block in the config file

    try:
      input_render_vars.update({
        "method" : config['orca']['keywords']['method'],
        "basis_set" : config['orca']['keywords']['basis_set'],
        "aux_basis_set" : config['orca']['keywords']['aux_basis_set'],
        "other" : config['orca']['keywords']['other']
      })

    except KeyError as error:
      raise abin_errors.AbinError ('ERROR: The "%s" key is missing in the "keywords" block of the "orca" block in the "%s" configuration file.' % (error,misc['config_name']))

    # Variables specific to the pre_opt template

    if pre_opt:

      input_render_vars.update({
        "mol_name" : misc['mol_name']
      })

    # Rendering the file
    # ==================
     
    rendered_content[rendered_input] = jinja_render(misc['templates_dir'], template_input, input_render_vars)

    print('%12s' % "[ DONE ]")

    # ========================================================= #
    #                  Rendering the job script                 #
    # ========================================================= #

    # Get the path to the "check_scripts" directory because the job script needs to execute check_orca.py

    chains_path = os.path.dirname(misc['code_dir'])  
    check_script_path = os.path.join(chains_path,"check_scripts")

    # If we need to copy the output files to their respective results directory, load the CHAINS configuration file to get the necessary information

    if copy_files:

      chains_config_file = abin_errors.check_abspath(os.path.join(chains_path,"chains_config.yml"),"CHAINS configuration YAML file","file")
  
      print ("{:<80}".format("\nLoading CHAINS configuration YAML file ..."), end="")
      with open(chains_config_file, 'r') as chains:
        chains_config = yaml.load(chains, Loader=yaml.FullLoader)
      print('%12s' % "[ DONE ]")

    print("{:<80}".format("\nRendering the jinja template for the orca job script ..."), end="")

    # Defining the Jinja variables
    # ============================

    # Variables not associated with the config file

    script_render_vars = {  
      "mol_name" : misc['mol_name'],
      "job_walltime" : job_specs['walltime'],
      "job_cores" : job_specs['cores'],
      "job_mem_per_cpu" : job_specs['mem_per_cpu'], # in MB
      "cluster_name" : job_specs['cluster_name'],
      "partition" : job_specs['partition'],     
      "chains_dir" : chains_path,
      "check_dir" : check_script_path,
      "copy_files" : copy_files,
      "benchmark" : benchmark
    }

    # Variables associated with the "general" block of the config file

    try:
      script_render_vars.update({
        "user_email" : config['general']['user_email'],
        "mail_type" : config['general']['mail_type']
      })

    except KeyError as error:
      raise abin_errors.AbinError ('ERROR: The "%s" key is missing in the "general" block of the "%s" configuration file.' % (error,misc['config_name']))

    # Variables associated with the clusters configuration file

    try:
      script_render_vars.update({
        "set_env" : clusters_cfg[job_specs['cluster_name']]['profiles'][job_specs['profile']]['set_env'],       
        "command" : clusters_cfg[job_specs['cluster_name']]['profiles'][job_specs['profile']]['command']
      })

    except KeyError as error:
      raise abin_errors.AbinError ('ERROR: The "%s" key is missing in the "%s" profile of the clusters configuration file.' % (error,job_specs['profile']))

    # Variables specific to the copy_files portion of the template

    if copy_files:

      script_render_vars.update({
        "output_dir" : chains_config['output_dir']['orca'],
        "results_dir" : chains_config['results_dir'],
        "config_file" : misc['config_name'],
        "job_script" : rendered_script
      })

    # Variables specific to the benchmarking template
   
    if benchmark:

      script_render_vars.update({
        "benchmark_path" : "${CECIHOME}/BENCHMARK",
        "prefix": job_specs['profile'] + "_" + job_specs['cluster_name'],
        "profile" : job_specs['profile'],
        "cluster_name" : job_specs['cluster_name'],
        "jobscale_label" : job_specs['scale_label'],
        "job_walltime" : job_specs['walltime'],
        "job_mem_per_cpu" : job_specs['mem_per_cpu'], # in MB
        "scaling_function" : job_specs['scaling_fct'],
        "scale_index" : job_specs['scale_index']
      })
    
    # Rendering the file
    # ==================

    rendered_content[rendered_script] = jinja_render(misc['templates_dir'], template_script, script_render_vars)

    print('%12s' % "[ DONE ]")

    return rendered_content, rendered_script


######################################################################################################################################


def chains_qchem_render(mendeleev:dict, clusters_cfg:dict, config:dict, file_data:dict, job_specs:dict, misc:dict):
    """Renders the job script and the input file associated with the Q-CHEM program in CHAINS.

    Parameters
    ----------
    mendeleev : dict
        Content of AlexGustafsson's Mendeleev Table YAML file (found at https://github.com/AlexGustafsson/molecular-data).
        Unused in this function.
    clusters_cfg : dict
        Content of the YAML clusters configuration file.
    config : dict
        Content of the YAML configuration file.
    file_data : dict
        Information extracted by the scanning function from the geometry file.
    job_specs : dict
        Contains all information related to the job.
    misc : dict
        Contains all the additional variables that did not pertain to the other arguments.

    Returns
    -------
    rendered_content : dict
        Dictionary containing the text of all the rendered files in the form of <filename>: <rendered_content>.
    rendered_script : str
        Name of the rendered job script, necessary to launch the job.
    
    Notes
    -----
    Pay a particular attention to the render_vars dictionaries, they contain all the definitions of the variables appearing in your Jinja templates.
    """

    # ========================================================= #
    #                      Preparation step                     #
    # ========================================================= #

    # Check config file
    # =================

    # Check if a "general" block has been defined in the config file

    if not config.get('general'):
      raise abin_errors.AbinError ('ERROR: There is no "general" key defined in the "%s" configuration file.' % misc['config_name'])     

    # Check if a "qchem" block has been defined in the config file

    if not config.get('qchem'):
      raise abin_errors.AbinError ('ERROR: There is no "qchem" key defined in the "%s" configuration file.' % misc['config_name'])      

    # Check the options defined in the config file

    copy_files = config['qchem'].get('copy_files',True)

    if not isinstance(copy_files, bool):
      raise abin_errors.AbinError ('ERROR: The "copy_files" value given in the "qchem" block of the "%s" configuration file is not a boolean (neither "True" nor "False").' % misc['config_name'])

    benchmark = config['qchem'].get('benchmark',False)

    if not isinstance(benchmark, bool):
      raise abin_errors.AbinError ('ERROR: The "benchmark" value given in the "qchem" block of the "%s" configuration file is not a boolean (neither "True" nor "False").' % misc['config_name'])

    # Define the templates
    # ====================

    # Define the names of the templates.

    template_input = "qchem.in.jinja"
    template_script = "qchem_job.sh.jinja"

    # Check if the specified templates exist in the "templates" directory of ABIN LAUNCHER.
    
    abin_errors.check_abspath(os.path.join(misc['templates_dir'],template_input),"Jinja template for the qchem input file","file")
    abin_errors.check_abspath(os.path.join(misc['templates_dir'],template_script),"Jinja template for the qchem job script","file")

    # Define rendered files
    # =====================

    # Define the names of the rendered files.

    rendered_input = misc['mol_name'] + ".in"
    rendered_script = "qchem_job.sh"

    # Initialize the dictionary that will be returned by the function

    rendered_content = {}

    # ========================================================= #
    #                  Rendering the input file                 #
    # ========================================================= #

    print("{:<80}".format("\nRendering the jinja template for the qchem input file ...  "), end="")

    # Defining the Jinja variables
    # ============================

    # Variables not associated with the config file

    input_render_vars = {
      "mem_total" : job_specs['cores'] * job_specs['mem_per_cpu'],
      "coordinates" : file_data['atomic_coordinates']
    }

    # Variables associated with the "general" block of the config file

    try:
      input_render_vars.update({
        "charge" : config['general']['charge'],
        "multiplicity" : config['general']['multiplicity']
      })

    except KeyError as error:
      raise abin_errors.AbinError ('ERROR: The "%s" key is missing in the "general" block of the "%s" configuration file.' % (error,misc['config_name']))

    # Check if a "keywords" block has been defined in the "qchem" block of the config file

    if not config['qchem'].get('keywords'):
      raise abin_errors.AbinError ('ERROR: There is no "keywords" key in the "qchem" block of the "%s" configuration file.' % misc['config_name'])    

    # Variables associated with the "keywords" block of the "qchem" block in the config file

    try:
      input_render_vars.update({
        "job_type" : config['qchem']['keywords']['job_type'],
        "exchange" : config['qchem']['keywords']['exchange'],
        "basis_set" : config['qchem']['keywords']['basis_set'],
        "cis_n_roots" : config['qchem']['keywords']['cis_n_roots']
      })

    except KeyError as error:
      raise abin_errors.AbinError ('ERROR: The "%s" key is missing in the "keywords" block of the "qchem" block in the "%s" configuration file.' % (error,misc['config_name']))

    # Rendering the file
    # ==================

    rendered_content[rendered_input] = jinja_render(misc['templates_dir'], template_input, input_render_vars)

    print('%12s' % "[ DONE ]")

    # ========================================================= #
    #                  Rendering the job script                 #
    # ========================================================= #

    # Get the path to the "check_scripts" directory because the job script needs to execute check_qchem.py

    chains_path = os.path.dirname(misc['code_dir'])  
    check_script_path = os.path.join(chains_path,"check_scripts")

    # If we need to copy the output files to their respective results directory, load the CHAINS configuration file to get the necessary information

    if copy_files:

      chains_config_file = abin_errors.check_abspath(os.path.join(chains_path,"chains_config.yml"),"CHAINS configuration YAML file","file")
  
      print ("{:<80}".format("\nLoading CHAINS configuration YAML file ..."), end="")
      with open(chains_config_file, 'r') as chains:
        chains_config = yaml.load(chains, Loader=yaml.FullLoader)
      print('%12s' % "[ DONE ]")

    print("{:<80}".format("\nRendering the jinja template for the qchem job script ..."), end="")

    # Defining the Jinja variables
    # ============================

    # Variables not associated with the config file

    script_render_vars = {  
      "mol_name" : misc['mol_name'],
      "job_walltime" : job_specs['walltime'],
      "job_cores" : job_specs['cores'],
      "job_mem_per_cpu" : job_specs['mem_per_cpu'], # in MB
      "partition" : job_specs['partition'],     
      "chains_dir" : chains_path,
      "check_dir" : check_script_path,
      "copy_files" : copy_files,
      "benchmark" : benchmark
    }

    # Variables associated with the "general" block of the config file

    try:
      script_render_vars.update({
        "user_email" : config['general']['user_email'],
        "mail_type" : config['general']['mail_type']
      })

    except KeyError as error:
      raise abin_errors.AbinError ('ERROR: The "%s" key is missing in the "general" block of the "%s" configuration file.' % (error,misc['config_name']))

    # Variables associated with the clusters configuration file

    try:
      script_render_vars.update({
        "set_env" : clusters_cfg[job_specs['cluster_name']]['profiles'][job_specs['profile']]['set_env'],       
        "command" : clusters_cfg[job_specs['cluster_name']]['profiles'][job_specs['profile']]['command']
      })

    except KeyError as error:
      raise abin_errors.AbinError ('ERROR: The "%s" key is missing in the "%s" profile of the clusters configuration file.' % (error,job_specs['profile']))

    # Variables specific to the copy_files portion of the template

    if copy_files:

      script_render_vars.update({
        "output_dir" : chains_config['output_dir']['qchem'],
        "results_dir" : chains_config['results_dir'],
        "config_file" : misc['config_name'],
        "job_script" : rendered_script
      })

    # Variables specific to the benchmarking template
   
    if benchmark:

      script_render_vars.update({
        "benchmark_path" : "${CECIHOME}/BENCHMARK",
        "prefix": job_specs['profile'] + "_" + job_specs['cluster_name'],
        "profile" : job_specs['profile'],
        "cluster_name" : job_specs['cluster_name'],
        "jobscale_label" : job_specs['scale_label'],
        "job_walltime" : job_specs['walltime'],
        "job_mem_per_cpu" : job_specs['mem_per_cpu'], # in MB
        "scaling_function" : job_specs['scaling_fct'],
        "scale_index" : job_specs['scale_index']
      })
    
    # Rendering the file
    # ==================

    rendered_content[rendered_script] = jinja_render(misc['templates_dir'], template_script, script_render_vars)

    print('%12s' % "[ DONE ]")
   
    return rendered_content, rendered_script
