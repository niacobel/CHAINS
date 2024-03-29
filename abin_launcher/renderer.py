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

def basic_rendering(mendeleev:dict, clusters_cfg:dict, config:dict, file_data:dict, job_specs:dict, misc:dict):
    """Renders the job script and the input file associated with simple jobs. The templates filenames are specified in the YAML configuration file and the jinja variables are taken as is from it (outside of the variables related to the clusters YAML configuration file).

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

    # Check if a "templates" block has been defined in the config file

    if not config.get('templates'):
      raise abin_errors.AbinError ('ERROR: There is no "templates" key defined in the "%s" configuration file.' % misc['config_name'])     

    # Check if a "jinja_variables" block has been defined in the config file

    if not config.get('jinja_variables'):
      raise abin_errors.AbinError ('ERROR: There is no "jinja_variables" key defined in the "%s" configuration file.' % misc['config_name'])      

    # Define the templates
    # ====================

    # Define the names of the templates.

    try:
      template_input = config['templates']['input']
      template_script = config['templates']['job_script']

    except KeyError as error:
      raise abin_errors.AbinError ('ERROR: The "%s" key is missing in the "templates" block of the "%s" configuration file.' % (error,misc['config_name']))

    # Check if the specified templates exist in the "templates" directory of ABIN LAUNCHER.
    
    abin_errors.check_abspath(os.path.join(misc['templates_dir'],template_input),"Jinja template for the input file","file")
    abin_errors.check_abspath(os.path.join(misc['templates_dir'],template_script),"Jinja template for the job script","file")

    # Define rendered files
    # =====================

    # Define the names of the rendered files.

    rendered_input = misc['mol_name'] + os.path.splitext(template_input.rstrip(".jinja"))[1]
    rendered_script = template_script.rstrip(".jinja")

    # Initialize the dictionary that will be returned by the function

    rendered_content = {}

    # ========================================================= #
    #                  Rendering the input file                 #
    # ========================================================= #
  
    print("{:<80}".format("\nRendering the jinja template for the input file ...  "), end="")

    # Preparing the Jinja variables
    # =============================

    # Variables given by the config file

    input_render_vars = config['jinja_variables']

    # Other useful variables

    input_render_vars.update({
      "mol_name" : misc['mol_name'],
      "config_file" : misc['config_name'],
      "coordinates" : file_data['atomic_coordinates'],
      "job_cores" : job_specs['cores'],
      "job_mem_per_cpu" : job_specs['mem_per_cpu'] # in MB
    })

    # Rendering the file
    # ==================
     
    rendered_content[rendered_input] = jinja_render(misc['templates_dir'], template_input, input_render_vars)

    print('%12s' % "[ DONE ]")

    # ========================================================= #
    #                  Rendering the job script                 #
    # ========================================================= #

    # Preparing the Jinja variables
    # =============================

    # Variables given by the config file

    script_render_vars = config['jinja_variables']

    # Variables not associated with the config file

    script_render_vars.update({
      "mol_name" : misc['mol_name'],
      "config_file" : misc['config_name'],
      "job_script" : rendered_script,
      "job_walltime" : job_specs['walltime'],
      "job_cores" : job_specs['cores'],
      "job_mem_per_cpu" : job_specs['mem_per_cpu'], # in MB
      "cluster_name" : job_specs['cluster_name'],
      "partition" : job_specs['partition']
    })

    # Variables associated with the clusters configuration file

    try:
      script_render_vars.update({
        "set_env" : clusters_cfg[job_specs['cluster_name']]['profiles'][job_specs['profile']]['set_env'],       
        "command" : clusters_cfg[job_specs['cluster_name']]['profiles'][job_specs['profile']]['command']
      })

    except KeyError as error:
      raise abin_errors.AbinError ('ERROR: The "%s" key is missing in the "%s" profile of the clusters configuration file.' % (error,job_specs['profile']))
    
    # Rendering the file
    # ==================

    rendered_content[rendered_script] = jinja_render(misc['templates_dir'], template_script, script_render_vars)

    print('%12s' % "[ DONE ]")

    return rendered_content, rendered_script


######################################################################################################################################


def chains_gaussian_render(mendeleev:dict, clusters_cfg:dict, config:dict, file_data:dict, job_specs:dict, misc:dict):
    """Renders the job script and the input file associated with the GAUSSIAN program in CHAINS.

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

    # Check if a "gaussian" block has been defined in the config file

    if not config.get('gaussian'):
      raise abin_errors.AbinError ('ERROR: There is no "gaussian" key defined in the "%s" configuration file.' % misc['config_name'])      

    # Check the options defined in the config file

    auto_restart = config['gaussian'].get('auto_restart',False)

    if not isinstance(auto_restart, bool):
      raise abin_errors.AbinError ('ERROR: The "auto_restart" value given in the "gaussian" block of the "%s" configuration file is not a boolean (neither "True" nor "False").' % misc['config_name'])

    benchmark = config['gaussian'].get('benchmark',False)

    if not isinstance(benchmark, bool):
      raise abin_errors.AbinError ('ERROR: The "benchmark" value given in the "gaussian" block of the "%s" configuration file is not a boolean (neither "True" nor "False").' % misc['config_name'])

    copy_files = config['gaussian'].get('copy_files',False)

    if not isinstance(copy_files, bool):
      raise abin_errors.AbinError ('ERROR: The "copy_files" value given in the "gaussian" block of the "%s" configuration file is not a boolean (neither "True" nor "False").' % misc['config_name'])

    if copy_files:
      ip_calc = str(config['gaussian'].get('ip_calc',"no_key")).lower()
      if ip_calc == "no_key":
        # If there was no "ip_calc" key in the config file, use the scale index to define how the ionization potential will be calculated
        if job_specs['scale_index'] > 650:
          ip_calc = "vertical"
        else:
          ip_calc = "adiabatic"
      elif ip_calc not in ["none","vertical","adiabatic"]:
        raise abin_errors.AbinError ('ERROR: The "ip_calc" value given in the "gaussian" block of the "%s" configuration file is neither "None", "Vertical" nor "Adiabatic" (This is not case sensitive).' % misc['config_name'])
    else:
      ip_calc = None

    # Define the templates
    # ====================

    # Define the names of the templates.

    template_input = "gaussian.com.jinja"
    template_script = "gaussian_job.sh.jinja"

    # Check if the specified templates exist in the "templates" directory of ABIN LAUNCHER.
    
    abin_errors.check_abspath(os.path.join(misc['templates_dir'],template_input),"Jinja template for the gaussian input file","file")
    abin_errors.check_abspath(os.path.join(misc['templates_dir'],template_script),"Jinja template for the gaussian job script","file")

    # Define rendered files
    # =====================

    # Define the names of the rendered files.

    rendered_input = misc['mol_name'] + ".com"
    rendered_script = "gaussian_job.sh"

    # Initialize the dictionary that will be returned by the function

    rendered_content = {}

    # ========================================================= #
    #                  Rendering the input file                 #
    # ========================================================= #
  
    print("{:<80}".format("\nRendering the jinja template for the gaussian input file ...  "), end="")

    # Defining the mandatory Jinja variables
    # ======================================

    # Variables not associated with the config file

    input_render_vars = {
      "mol_name" : misc['mol_name'],
      "mem_total" : job_specs['cores'] * job_specs['mem_per_cpu'],
      "job_cores" : job_specs['cores'],
      "coordinates" : file_data['atomic_coordinates'],
      "ip_calc" : ip_calc # Associated with the config file, but it has already been verified
    }

    # Variables associated with the "general" block of the config file

    try:
      input_render_vars.update({
        "charge" : config['general']['charge'],
        "multiplicity" : config['general']['multiplicity']
      })

    except KeyError as error:
      raise abin_errors.AbinError ('ERROR: The "%s" key is missing in the "general" block of the "%s" configuration file.' % (error,misc['config_name']))

    # Check if a "keywords" block has been defined in the "gaussian" block of the config file

    if not config['gaussian'].get('keywords'):
      raise abin_errors.AbinError ('ERROR: There is no "keywords" key in the "gaussian" block of the "%s" configuration file.' % misc['config_name'])    

    # Variables associated with the "keywords" block of the "gaussian" block in the config file

    try:
      input_render_vars.update({
        "method" : config['gaussian']['keywords']['method'],
        "basis_set" : config['gaussian']['keywords']['basis_set'],
        "other" : config['gaussian']['keywords']['other']
      })

    except KeyError as error:
      raise abin_errors.AbinError ('ERROR: The "%s" key is missing in the "keywords" block of the "gaussian" block in the "%s" configuration file.' % (error,misc['config_name']))

    # Defining the specific Jinja variables
    # =====================================

    # Variables specific to the ip_calc == vertical portion of the template

    if ip_calc == "vertical":

      # Determine the new charge and multiplicity of the cation

      charge_cation = int(config['general']['charge']) + 1

      if int(config['general']['multiplicity']) == 1: # In the case of a singlet ground state, we break a pair of electrons and thus "create" a new unpaired electron
        multiplicity_cation = 2
      else:
        multiplicity_cation = int(config['general']['multiplicity']) - 1 # We removed one unpaired electron

      # Variables associated with the config file but they have already been verified

      input_render_vars.update({
        "charge_cation" : charge_cation,
        "multiplicity_cation" : multiplicity_cation
      })

    # Rendering the file
    # ==================
     
    rendered_content[rendered_input] = jinja_render(misc['templates_dir'], template_input, input_render_vars)

    print('%12s' % "[ DONE ]")

    # ========================================================= #
    #                  Rendering the job script                 #
    # ========================================================= #

    # Get the path to the "check_scripts" directory because the job script needs to execute gaussian_check.py

    chains_path = os.path.dirname(misc['code_dir'])  
    check_script_path = os.path.join(chains_path,"check_scripts")

    # If we need to copy the output files to their respective results directory (or compute the ionization potentials through the cation), load the CHAINS configuration file to get the necessary information

    if copy_files:

      chains_config_file = abin_errors.check_abspath(os.path.join(chains_path,"configs","chains_config.yml"),"CHAINS configuration YAML file","file")
  
      print ("{:<80}".format("\nLoading CHAINS configuration YAML file ..."), end="")
      with open(chains_config_file, 'r') as chains:
        chains_config = yaml.load(chains, Loader=yaml.FullLoader)
      print('%12s' % "[ DONE ]")

    print("{:<80}".format("\nRendering the jinja template for the gaussian job script ..."), end="")

    # Defining the mandatory Jinja variables
    # ======================================

    # Variables not associated with the config file

    script_render_vars = {  
      "mol_name" : misc['mol_name'],
      "config_file" : misc['config_name'],
      "job_walltime" : job_specs['walltime'],
      "job_cores" : job_specs['cores'],
      "job_mem_per_cpu" : job_specs['mem_per_cpu'], # in MB
      "cluster_name" : job_specs['cluster_name'],
      "partition" : job_specs['partition'],     
      "chains_dir" : chains_path,
      "check_dir" : check_script_path,
      "auto_restart" : auto_restart, # Associated with the config file, but it has already been verified
      "copy_files" : copy_files,     # Associated with the config file, but it has already been verified
      "benchmark" : benchmark        # Associated with the config file, but it has already been verified
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

    # Defining the specific Jinja variables
    # =====================================

    # Variables specific to the copy_files portion of the template

    if copy_files:

      # Variables not associated with the config file

      script_render_vars.update({
        "ip_calc" : ip_calc, # Associated with the config file, but it has already been verified
        "job_script" : rendered_script
      })

      # Variables associated with the CHAINS configuration file

      try:
        script_render_vars.update({
          "output_dir" : chains_config['output_gaussian'],
          "results_dir" : chains_config['results_dir']
        })

        if ip_calc == 'vertical' or ip_calc == 'adiabatic':
          script_render_vars.update({
            "ip_file" : chains_config['ip_file']
          })

      except KeyError as error:
        raise abin_errors.AbinError ('ERROR: The "%s" key is missing in the CHAINS configuration file (chains_config.yml).' % error)

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


def chains_gaussian_cation_render(mendeleev:dict, clusters_cfg:dict, config:dict, file_data:dict, job_specs:dict, misc:dict):
    """Renders the job script and the input file associated with cation calculations using the GAUSSIAN program in CHAINS.

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

    # Check if a "gaussian" block has been defined in the config file

    if not config.get('gaussian'):
      raise abin_errors.AbinError ('ERROR: There is no "gaussian" key defined in the "%s" configuration file.' % misc['config_name'])      

    # Check the options defined in the config file

    benchmark = config['gaussian'].get('benchmark',False)

    if not isinstance(benchmark, bool):
      raise abin_errors.AbinError ('ERROR: The "benchmark" value given in the "gaussian" block of the "%s" configuration file is not a boolean (neither "True" nor "False").' % misc['config_name'])

    # Define the templates
    # ====================

    # Define the names of the templates.

    template_input = "gaussian_cation.com.jinja"
    template_script = "gaussian_cation_job.sh.jinja"

    # Check if the specified templates exist in the "templates" directory of ABIN LAUNCHER.
    
    abin_errors.check_abspath(os.path.join(misc['templates_dir'],template_input),"Jinja template for the gaussian cation input file","file")
    abin_errors.check_abspath(os.path.join(misc['templates_dir'],template_script),"Jinja template for the gaussian cation job script","file")

    # Define rendered files
    # =====================

    # Define the names of the rendered files.

    rendered_input = misc['mol_name'] + ".com"
    rendered_script = "gaussian_cation_job.sh"

    # Initialize the dictionary that will be returned by the function

    rendered_content = {}

    # ========================================================= #
    #                  Rendering the input file                 #
    # ========================================================= #
  
    print("{:<80}".format("\nRendering the jinja template for the gaussian cation input file ...  "), end="")

    # Defining the mandatory Jinja variables
    # ======================================

    # Variables not associated with the config file

    input_render_vars = {
      "mol_name" : misc['mol_name'],
      "mem_total" : job_specs['cores'] * job_specs['mem_per_cpu'],
      "job_cores" : job_specs['cores'],
      "coordinates" : file_data['atomic_coordinates']
    }

    # Variables associated with the "general" block of the config file

    try:
      charge = int(config['general']['charge'])
      multiplicity = int(config['general']['multiplicity'])

    except KeyError as error:
      raise abin_errors.AbinError ('ERROR: The "%s" key is missing in the "general" block of the "%s" configuration file.' % (error,misc['config_name']))

    # Determine the new charge and multiplicity of the cation

    charge_cation = charge + 1

    if multiplicity == 1: # In the case of a singlet ground state, we break a pair of electrons and thus "create" a new unpaired electron
      multiplicity_cation = 2
    else:
      multiplicity_cation = multiplicity - 1 # We removed one unpaired electron

    # Add the new variables

    input_render_vars.update({
      "charge" : charge,
      "multiplicity" : multiplicity,
      "charge_cation" : charge_cation,
      "multiplicity_cation" : multiplicity_cation
    })

    # Check if a "keywords" block has been defined in the "gaussian" block of the config file

    if not config['gaussian'].get('keywords'):
      raise abin_errors.AbinError ('ERROR: There is no "keywords" key in the "gaussian" block of the "%s" configuration file.' % misc['config_name'])    

    # Variables associated with the "keywords" block of the "gaussian" block in the config file

    try:
      input_render_vars.update({
        "method" : config['gaussian']['keywords']['method'],
        "basis_set" : config['gaussian']['keywords']['basis_set'],
        "other" : config['gaussian']['keywords']['other']
      })

    except KeyError as error:
      raise abin_errors.AbinError ('ERROR: The "%s" key is missing in the "keywords" block of the "gaussian" block in the "%s" configuration file.' % (error,misc['config_name']))

    # Rendering the file
    # ==================
     
    rendered_content[rendered_input] = jinja_render(misc['templates_dir'], template_input, input_render_vars)

    print('%12s' % "[ DONE ]")

    # ========================================================= #
    #                  Rendering the job script                 #
    # ========================================================= #

    # Get the path to the "check_scripts" directory because the job script needs to execute gaussian_check.py

    chains_path = os.path.dirname(misc['code_dir'])  
    check_script_path = os.path.join(chains_path,"check_scripts")

    # Load the CHAINS configuration file to get the path towards the IP file

    chains_config_file = abin_errors.check_abspath(os.path.join(chains_path,"configs","chains_config.yml"),"CHAINS configuration YAML file","file")

    print ("{:<80}".format("\nLoading CHAINS configuration YAML file ..."), end="")
    with open(chains_config_file, 'r') as chains:
      chains_config = yaml.load(chains, Loader=yaml.FullLoader)
    print('%12s' % "[ DONE ]")

    print("{:<80}".format("\nRendering the jinja template for the gaussian job script ..."), end="")

    # Defining the mandatory Jinja variables
    # ======================================

    # Variables not associated with the config file

    script_render_vars = {  
      "mol_name" : misc['mol_name'],
      "config_file" : misc['config_name'],
      "job_walltime" : job_specs['walltime'],
      "job_cores" : job_specs['cores'],
      "job_mem_per_cpu" : job_specs['mem_per_cpu'], # in MB
      "cluster_name" : job_specs['cluster_name'],
      "partition" : job_specs['partition'],     
      "chains_dir" : chains_path,
      "check_dir" : check_script_path,
      "benchmark" : benchmark    # Associated with the config file, but it has already been verified
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

    # Variables associated with the CHAINS configuration file

    try:
      script_render_vars.update({
        "ip_file" : chains_config['ip_file']
      })

    except KeyError as error:
      raise abin_errors.AbinError ('ERROR: The "%s" key is missing in the CHAINS configuration file (chains_config.yml).' % error)

    # Defining the specific Jinja variables
    # =====================================

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

    copy_files = config['qchem'].get('copy_files',False)

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

    # Define the memory usage

    mem_total = job_specs['cores'] * job_specs['mem_per_cpu']
    mem_static = int(0.02*mem_total) if int(0.02*mem_total) > 200 else 200

    # Variables not associated with the config file

    input_render_vars = {
      "mem_total" : mem_total,
      "mem_static" : mem_static,
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
        "scf_algorithm" : config['qchem']['keywords']['scf_algorithm'],
        "max_scf_cycles" : config['qchem']['keywords']['max_scf_cycles'],
        "cis_n_roots" : config['qchem']['keywords']['cis_n_roots'],
        "iqmol_fchk":  config['qchem']['keywords']['iqmol_fchk'],
        "cis_ampl_anal":  config['qchem']['keywords']['cis_ampl_anal']
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

    # Get the path to the "check_scripts" directory because the job script needs to execute qchem_check.py

    chains_path = os.path.dirname(misc['code_dir'])  
    check_script_path = os.path.join(chains_path,"check_scripts")

    # If we need to copy the output files to their respective results directory, load the CHAINS configuration file to get the necessary information

    if copy_files:

      chains_config_file = abin_errors.check_abspath(os.path.join(chains_path,"configs","chains_config.yml"),"CHAINS configuration YAML file","file")
  
      print ("{:<80}".format("\nLoading CHAINS configuration YAML file ..."), end="")
      with open(chains_config_file, 'r') as chains:
        chains_config = yaml.load(chains, Loader=yaml.FullLoader)
      print('%12s' % "[ DONE ]")

    print("{:<80}".format("\nRendering the jinja template for the qchem job script ..."), end="")

    # Defining the mandatory Jinja variables
    # ======================================

    # Variables not associated with the config file

    script_render_vars = {  
      "mol_name" : misc['mol_name'],
      "config_file" : misc['config_name'],
      "job_walltime" : job_specs['walltime'],
      "job_cores" : job_specs['cores'],
      "job_mem_per_cpu" : job_specs['mem_per_cpu'], # in MB
      "partition" : job_specs['partition'],     
      "chains_dir" : chains_path,
      "check_dir" : check_script_path,
      "copy_files" : copy_files, # Associated with the config file, but it has already been verified
      "benchmark" : benchmark    # Associated with the config file, but it has already been verified
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

    # Defining the specific Jinja variables
    # =====================================

    # Variables specific to the copy_files portion of the template

    if copy_files:

      # Variables not associated with the config file

      script_render_vars.update({
        "job_script" : rendered_script
      })

      # Variables associated with the CHAINS configuration file

      try:

        if config['qchem']['keywords']['basis_set'].lower() == "def2-tzvp": # The output directories are different for the TZVP basis set

          script_render_vars.update({
            "output_dir" : chains_config['output_qchem_tzvp'],
            "results_dir" : chains_config['results_dir_tzvp']
          })

        else:

          script_render_vars.update({
            "output_dir" : chains_config['output_qchem'],
            "results_dir" : chains_config['results_dir']
          })

      except KeyError as error:
        raise abin_errors.AbinError ('ERROR: The "%s" key is missing in the CHAINS configuration file (chains_config.yml).' % error)

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

    copy_files = config['orca'].get('copy_files',True)

    if not isinstance(copy_files, bool):
      raise abin_errors.AbinError ('ERROR: The "copy_files" value given in the "orca" block of the "%s" configuration file is not a boolean (neither "True" nor "False").' % misc['config_name'])

    benchmark = config['orca'].get('benchmark',False)

    if not isinstance(benchmark, bool):
      raise abin_errors.AbinError ('ERROR: The "benchmark" value given in the "orca" block of the "%s" configuration file is not a boolean (neither "True" nor "False").' % misc['config_name'])

    # Define the templates
    # ====================

    # Define the names of the templates.

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
      "mol_name" : misc['mol_name'],
      "orca_mem_per_cpu" : orca_mem_per_cpu,
      "job_cores" : job_specs['cores']
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
        "other" : config['orca']['keywords']['other'],
        "nroots" : config['orca']['keywords']['nroots'],
        "printlevel_tdm" : config['orca']['keywords']['printlevel_tdm'],
        "printlevel_soc" : config['orca']['keywords']['printlevel_soc']
      })

    except KeyError as error:
      raise abin_errors.AbinError ('ERROR: The "%s" key is missing in the "keywords" block of the "orca" block in the "%s" configuration file.' % (error,misc['config_name']))

    # Rendering the file
    # ==================
     
    rendered_content[rendered_input] = jinja_render(misc['templates_dir'], template_input, input_render_vars)

    print('%12s' % "[ DONE ]")

    # ========================================================= #
    #                  Rendering the job script                 #
    # ========================================================= #

    # Get the path to the "check_scripts" directory because the job script needs to execute orca_check.py

    chains_path = os.path.dirname(misc['code_dir'])  
    check_script_path = os.path.join(chains_path,"check_scripts")

    # If we need to copy the output files to their respective results directory, load the CHAINS configuration file to get the necessary information

    if copy_files:

      chains_config_file = abin_errors.check_abspath(os.path.join(chains_path,"configs","chains_config.yml"),"CHAINS configuration YAML file","file")
  
      print ("{:<80}".format("\nLoading CHAINS configuration YAML file ..."), end="")
      with open(chains_config_file, 'r') as chains:
        chains_config = yaml.load(chains, Loader=yaml.FullLoader)
      print('%12s' % "[ DONE ]")

    print("{:<80}".format("\nRendering the jinja template for the orca job script ..."), end="")

    # Defining the mandatory Jinja variables
    # ======================================

    # Variables not associated with the config file

    script_render_vars = {  
      "mol_name" : misc['mol_name'],
      "config_file" : misc['config_name'],
      "job_walltime" : job_specs['walltime'],
      "job_cores" : job_specs['cores'],
      "job_mem_per_cpu" : job_specs['mem_per_cpu'], # in MB
      "cluster_name" : job_specs['cluster_name'],
      "partition" : job_specs['partition'],     
      "chains_dir" : chains_path,
      "check_dir" : check_script_path,
      "copy_files" : copy_files, # Associated with the config file, but it has already been verified
      "benchmark" : benchmark    # Associated with the config file, but it has already been verified
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

    # Defining the specific Jinja variables
    # =====================================

    # Variables specific to the copy_files portion of the template

    if copy_files:

      # Variables not associated with the config file

      script_render_vars.update({
        "config_file" : misc['config_name'],
        "job_script" : rendered_script
      })

      # Variables associated with the CHAINS configuration file

      try:
        script_render_vars.update({
          "output_dir" : chains_config['output_orca'],
          "results_dir" : chains_config['results_dir']
        })

      except KeyError as error:
        raise abin_errors.AbinError ('ERROR: The "%s" key is missing in the CHAINS configuration file (chains_config.yml).' % error)

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