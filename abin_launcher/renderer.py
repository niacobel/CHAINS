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

    # Load the CHAINS configuration file to get the additional information

    chains_path = os.path.dirname(misc['code_dir'])
    chains_config_file = abin_errors.check_abspath(os.path.join(chains_path,"chains_config.yml"),"CHAINS configuration YAML file","file")

    print ("{:<80}".format("\nLoading CHAINS configuration YAML file ..."), end="")
    with open(chains_config_file, 'r') as chains:
      chains_config = yaml.load(chains, Loader=yaml.FullLoader)
    print('%12s' % "[ DONE ]")

    # Check if we need to also perform a pre-optimization with ORCA (recommended for large molecules)

    pre_opt = config['orca'].get('pre_opt',False)

    # Define the names of the templates.

    if pre_opt:
        template_input = "orca_preopt.inp.jinja"
    else:
        template_input = "orca.inp.jinja"

    template_script = "orca_job.sh.jinja"

    # Check if the specified templates exist in the "templates" directory of ABIN LAUNCHER.
    
    abin_errors.check_abspath(os.path.join(misc['templates_dir'],template_input),"Jinja template for the orca input file","file")
    abin_errors.check_abspath(os.path.join(misc['templates_dir'],template_script),"Jinja template for the orca job script","file")

    # Define the names of the rendered files.

    rendered_input = misc['mol_name'] + ".inp"
    rendered_script = "orca_job.sh"

    # Initialize the dictionary that will be returned by the function

    rendered_content = {}

    # ========================================================= #
    #                  Rendering the input file                 #
    # ========================================================= #
  
    print("{:<80}".format("\nRendering the jinja template for the orca input file ...  "), end="")

    # It is recommended to set the memory per CPU to 75% of the physical memory available (see https://sites.google.com/site/orcainputlibrary/orca-common-problems)

    orca_mem_per_cpu = int(0.75 * job_specs['mem_per_cpu']) # in MB
    
    # Defining the Jinja variables

    input_render_vars = {
        "orca_mem_per_cpu" : orca_mem_per_cpu,
        "job_cores" : job_specs['cores'],
        "charge" : config['general']['charge'],
        "multiplicity" : config['general']['multiplicity'],
        "coordinates" : file_data['atomic_coordinates'],
        "method" : config['orca']['method'],
        "basis_set" : config['orca']['basis_set'],
        "aux_basis_set" : config['orca']['aux_basis_set'],
        "other" : config['orca']['other']
    }

    if pre_opt:

        input_render_vars.update({
            "mol_name" : misc['mol_name']
        })

    # Rendering the file
     
    rendered_content[rendered_input] = jinja_render(misc['templates_dir'], template_input, input_render_vars)

    print('%12s' % "[ DONE ]")

    # ========================================================= #
    #                  Rendering the job script                 #
    # ========================================================= #

    print("{:<80}".format("\nRendering the jinja template for the orca job script ..."), end="")

    # Get the path to the "check_scripts" directory because the job script needs to execute check_orca.py
    
    check_script_path = os.path.join(chains_path,"check_scripts")

    # Defining the Jinja variables
  
    script_render_vars = {  
        "mol_name" : misc['mol_name'],
        "user_email" : config['general']['user_email'],
        "mail_type" : config['general']['mail_type'],
        "job_walltime" : job_specs['walltime'],
        "job_cores" : job_specs['cores'],
        "job_mem_per_cpu" : job_specs['mem_per_cpu'], # in MB
        "partition" : job_specs['partition'],     
        "set_env" : clusters_cfg[job_specs['cluster_name']]['profiles'][job_specs['profile']]['set_env'],       
        "command" : clusters_cfg[job_specs['cluster_name']]['profiles'][job_specs['profile']]['command'],
        "output_dir" : chains_config['output_dir'][job_specs['profile']],
        "results_dir" : chains_config['results_dir'],
        "chains_dir" : chains_path,
        "check_dir" : check_script_path,
        "results_subdir" : job_specs['profile'].upper(),
        "job_script" : rendered_script,
        "config_file" : misc['config_name']
    }

    # Add variables specific to the benchmarking template
   
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

    # Load the CHAINS configuration file to get the additional information

    chains_path = os.path.dirname(misc['code_dir'])
    chains_config_file = abin_errors.check_abspath(os.path.join(chains_path,"chains_config.yml"),"CHAINS configuration YAML file","file")

    print ("{:<80}".format("\nLoading CHAINS configuration YAML file ..."), end="")
    with open(chains_config_file, 'r') as chains:
      chains_config = yaml.load(chains, Loader=yaml.FullLoader)
    print('%12s' % "[ DONE ]")

    # Define the names of the templates.

    template_input = "qchem.in.jinja"
    template_script = "qchem_job.sh.jinja"

    # Check if the specified templates exist in the "templates" directory of ABIN LAUNCHER.
    
    abin_errors.check_abspath(os.path.join(misc['templates_dir'],template_input),"Jinja template for the qchem input file","file")
    abin_errors.check_abspath(os.path.join(misc['templates_dir'],template_script),"Jinja template for the qchem job script","file")

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

    input_render_vars = {
        "mem_total" : job_specs['cores'] * job_specs['mem_per_cpu'],
        "job_type" : config['qchem']['job_type'],
        "exchange" : config['qchem']['exchange'],
        "basis_set" : config['qchem']['basis_set'],
        "cis_n_roots" : config['qchem']['cis_n_roots'],
        "charge" : config['general']['charge'],
        "multiplicity" : config['general']['multiplicity'],
        "coordinates" : file_data['atomic_coordinates']
    }

    # Rendering the file

    rendered_content[rendered_input] = jinja_render(misc['templates_dir'], template_input, input_render_vars)

    print('%12s' % "[ DONE ]")

    # ========================================================= #
    #                  Rendering the job script                 #
    # ========================================================= #

    print("{:<80}".format("\nRendering the jinja template for the qchem job script ..."), end="")

    # Get the path to the "check_scripts" directory because the job script needs to execute check_qchem.py
    
    check_script_path = os.path.join(chains_path,"check_scripts")

    # Defining the Jinja variables
    
    script_render_vars = {
        "mol_name" : misc['mol_name'],
        "user_email" : config['general']['user_email'],
        "mail_type" : config['general']['mail_type'],
        "job_walltime" : job_specs['walltime'],
        "job_cores" : job_specs['cores'],
        "job_mem_per_cpu" : job_specs['mem_per_cpu'], # in MB
        "partition" : job_specs['partition'],     
        "set_env" : clusters_cfg[job_specs['cluster_name']]['profiles'][job_specs['profile']]['set_env'],       
        "command" : clusters_cfg[job_specs['cluster_name']]['profiles'][job_specs['profile']]['command'],
        "output_dir" : chains_config['output_dir'][job_specs['profile']],
        "results_dir" : chains_config['results_dir'],
        "chains_dir" : chains_path,
        "check_dir" : check_script_path,
        "results_subdir" : job_specs['profile'].upper(),
        "job_script" : rendered_script,
        "config_file" : misc['config_name']
    }

    # Add variables specific to the benchmarking template
   
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

    rendered_content[rendered_script] = jinja_render(misc['templates_dir'], template_script, script_render_vars)

    print('%12s' % "[ DONE ]")
   
    return rendered_content, rendered_script
