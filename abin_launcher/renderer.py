################################################################################################################################################
##                                                                The Renderer                                                                ##
##                                                                                                                                            ##
##                                       This script contains the rendering functions for ABIN LAUNCHER,                                      ##
##                                consult the documentation at https://chains-ulb.readthedocs.io/ for details                                 ##
################################################################################################################################################

import os

from jinja2 import Environment, FileSystemLoader

import abin_errors


def jinja_render(templates_dir:str, template_file:str, render_vars:dict):
    """Renders a file based on its jinja template.

    Parameters
    ----------
    templates_dir : str
        The path towards the directory where the jinja template is located.
    template_file : str
        The name of the jinja template file.
    render_vars : dict
        Dictionary containing the definitions of all the variables present in the jinja template.

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

def orca_render(mendeleev:dict, clusters_cfg:dict, config:dict, file_data:dict, job_specs:dict, misc:dict):
    """Renders the job instructions file and the input file associated with the ORCA program.

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
    rendered_instructions : str
        Name of the rendered job instructions file, necessary to launch the job.
    
    Notes
    -----
    Pay a particular attention to the render_vars dictionary, it contains all the definitions of the variables appearing in your jinja template.
    """

    # Check if all the files specified in the clusters YAML file exists in the templates directory of abin_launcher.
    
    for filename in clusters_cfg[job_specs['cluster_name']]['progs'][job_specs['prog']]['jinja_templates'].values():    
        abin_errors.check_abspath(os.path.join(misc['templates_dir'],filename),"Jinja template","file")

    # Define the names of all the template and rendered files, given in the YAML clusters configuration file.

    template_input = clusters_cfg[job_specs['cluster_name']]['progs'][job_specs['prog']]['jinja_templates']['input']                     # Jinja template file for the orca input
    template_instructions = clusters_cfg[job_specs['cluster_name']]['progs'][job_specs['prog']]['jinja_templates']['job_instructions']         # Jinja template file for the orca job instructions file

    rendered_input = misc['mol_name'] + ".inp"                                                                                         # Name of the rendered input file (automatically named after the molecule and not defined in the clusters file)
    rendered_instructions = "orca_job.sh"                            # Name of the rendered job instructions file

    # Initialize our dictionary that will contain all the text of the rendered files

    rendered_content = {}

    # Get the path to the chains directory and the check_scripts directory because the job instructions file needs to execute check_orca.py and source load_modules.sh

    chains_path = os.path.dirname(misc['code_dir'])                      # Get the parent directory from the codes_dir (which is the abin_launcher directory)                         
    check_script_path = os.path.join(chains_path,"check_scripts")

    # Rendering the jinja template for the orca job instructions file
  
    print("{:<80}".format("\nRendering the jinja template for the orca job instructions file ..."), end="")

    render_vars = {  
        "mol_name" : misc['mol_name'],
        "user_email" : config['general']['user-email'],
        "mail_type" : config['general']['mail-type'],
        "job_walltime" : job_specs['walltime'],
        "job_cores" : job_specs['cores'],
        "job_mem_per_cpu" : job_specs['mem_per_cpu'], # in MB
        "partition" : job_specs['partition'],     
        "set_env" : clusters_cfg[job_specs['cluster_name']]['progs'][job_specs['prog']]['set_env'],       
        "command" : clusters_cfg[job_specs['cluster_name']]['progs'][job_specs['prog']]['command'],
        "output_folder" : config[job_specs['prog']]['output-folder'],
        "results_folder" : config['results']['main_folder'],
        "chains_folder" : chains_path,
        "check_folder" : check_script_path,
        "results_subfolder" : config['results'][job_specs['prog']]['folder_name'],
        "job_manifest" : rendered_instructions,
        "config_file" : misc['config_name'],
        "benchmark" : config['general']['benchmark'],
        "benchmark_folder" : config['general']['benchmark-folder'],
        "prog" : job_specs['prog'],
        "jobscale_label" : job_specs['scale_label'],
        "scaling_function" : job_specs['scaling_fct'],
        "scale_index" : job_specs['scale_index']
        }
    
    rendered_content[rendered_instructions] = jinja_render(misc['templates_dir'], template_instructions, render_vars)

    print('%12s' % "[ DONE ]")
   
    # Rendering the jinja template for the orca input file
  
    print("{:<80}".format("\nRendering the jinja template for the orca input file ...  "), end="")

    orca_mem_per_cpu = int(0.75 * job_specs['mem_per_cpu']) # in MB
    
    render_vars = {
        "method" : config[job_specs['prog']]['method'],
        "basis_set" : config[job_specs['prog']]['basis-set'],
        "aux_basis_set" : config[job_specs['prog']]['aux-basis-set'],
        "job_type" : config[job_specs['prog']]['job-type'],
        "other" : config[job_specs['prog']]['other'],
        "job_cores" : job_specs['cores'],
        "orca_mem_per_cpu" : orca_mem_per_cpu,
        "charge" : config['general']['charge'],
        "multiplicity" : config['general']['multiplicity'],
        "coordinates" : file_data['atomic_coordinates']
        }
      
    rendered_content[rendered_input] = jinja_render(misc['templates_dir'], template_input, render_vars)

    print('%12s' % "[ DONE ]")

    return rendered_content, rendered_instructions

def qchem_render(mendeleev:dict, clusters_cfg:dict, config:dict, file_data:dict, job_specs:dict, misc:dict):
    """Renders the job instructions file and the input file associated with the Q-CHEM program.

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
    rendered_instructions : str
        Name of the rendered job instructions file, necessary to launch the job.
    
    Notes
    -----
    Pay a particular attention to the render_vars dictionary, it contains all the definitions of the variables appearing in your jinja template.
    """

    # Check if all the files specified in the clusters YAML file exists in the templates directory of abin_launcher.
    
    for filename in clusters_cfg[job_specs['cluster_name']]['progs'][job_specs['prog']]['jinja_templates'].values():    
        abin_errors.check_abspath(os.path.join(misc['templates_dir'],filename),"Jinja template","file")

    # Define the names of all the template and rendered files, given in the YAML clusters configuration file.

    template_input = clusters_cfg[job_specs['cluster_name']]['progs'][job_specs['prog']]['jinja_templates']['input']                     # Jinja template file for the qchem input
    template_instructions = clusters_cfg[job_specs['cluster_name']]['progs'][job_specs['prog']]['jinja_templates']['job_instructions']         # Jinja template file for the qchem job instructions file

    rendered_input = misc['mol_name'] + ".in"                                                                                          # Name of the rendered input file (automatically named after the molecule and not defined in the clusters file)
    rendered_instructions = "qchem_job.sh"                            # Name of the rendered job instructions file

    # Initialize our dictionary that will contain all the text of the rendered files

    rendered_content = {}

    # Get the path to the chains directory and the check_scripts directory because the job instructions file needs to execute check_qchem.py and source load_modules.sh

    chains_path = os.path.dirname(misc['code_dir'])                      # Get the parent directory from the codes_dir (which is the abin_launcher directory)                         
    check_script_path = os.path.join(chains_path,"check_scripts")

    # Rendering the jinja template for the qchem job instructions file
  
    print("{:<80}".format("\nRendering the jinja template for the qchem job instructions file ..."), end="")
  
    render_vars = {
        "mol_name" : misc['mol_name'],
        "user_email" : config['general']['user-email'],
        "mail_type" : config['general']['mail-type'],
        "job_walltime" : job_specs['walltime'],
        "job_cores" : job_specs['cores'],
        "job_mem_per_cpu" : job_specs['mem_per_cpu'], # in MB
        "partition" : job_specs['partition'],     
        "set_env" : clusters_cfg[job_specs['cluster_name']]['progs'][job_specs['prog']]['set_env'],       
        "command" : clusters_cfg[job_specs['cluster_name']]['progs'][job_specs['prog']]['command'],
        "output_folder" : config[job_specs['prog']]['output-folder'],
        "results_folder" : config['results']['main_folder'],
        "chains_folder" : chains_path,
        "check_folder" : check_script_path,
        "results_subfolder" : config['results'][job_specs['prog']]['folder_name'],
        "job_manifest" : rendered_instructions,
        "config_file" : misc['config_name'],
        "benchmark" : config['general']['benchmark'],
        "benchmark_folder" : config['general']['benchmark-folder'],
        "prog" : job_specs['prog'],
        "jobscale_label" : job_specs['scale_label'],
        "scaling_function" : job_specs['scaling_fct'],
        "scale_index" : job_specs['scale_index']
        }
    
    rendered_content[rendered_instructions] = jinja_render(misc['templates_dir'], template_instructions, render_vars)

    print('%12s' % "[ DONE ]")
   
    # Rendering the jinja template for the qchem input file
  
    print("{:<80}".format("\nRendering the jinja template for the qchem input file ...  "), end="")
    
    render_vars = {
        "job_type" : config[job_specs['prog']]['job-type'],
        "exchange" : config[job_specs['prog']]['exchange'],
        "basis_set" : config[job_specs['prog']]['basis-set'],
        "cis_n_roots" : config[job_specs['prog']]['cis-n-roots'],
        "charge" : config['general']['charge'],
        "multiplicity" : config['general']['multiplicity'],
        "coordinates" : file_data['atomic_coordinates']
        }
      
    rendered_content[rendered_input] = jinja_render(misc['templates_dir'], template_input, render_vars)

    print('%12s' % "[ DONE ]")

    return rendered_content, rendered_instructions
