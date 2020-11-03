################################################################################################################################################
##                                                                The Renderer                                                                ##
##                                                                                                                                            ##
##                                       This script contains the rendering functions for ABIN LAUNCHER,                                      ##
##                                consult the documentation at https://chains-ulb.readthedocs.io/ for details                                 ##
################################################################################################################################################

import os

from jinja2 import Environment, FileSystemLoader

import abin_errors


def jinja_render(path_tpl_dir:str, tpl:str, render_vars:dict):
    """Renders a file based on its jinja template.

    Parameters
    ----------
    path_tpl_dir : str
        The path towards the directory where the jinja template is located.
    tpl : str
        The name of the jinja template file.
    render_vars : dict
        Dictionary containing the definitions of all the variables present in the jinja template.

    Returns
    -------
    output_text : str
        Content of the rendered file.
    """
   
    file_loader = FileSystemLoader(path_tpl_dir)
    env = Environment(loader=file_loader)
    template = env.get_template(tpl)
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
    
    Notes
    -----
    Pay a particular attention to the render_vars dictionary, it contains all the definitions of the variables appearing in your jinja template.
    """
    
    # Define the names of the templates
    
    tpl_inp = "sample_orca.inp.jinja"
    tpl_inst = "sample_orca_job.sh.jinja"
    
    # Define the names of the rendered files
    
    rnd_input = misc['mol_name'] + ".inp"
    rnd_inst = clusters_cfg[job_specs['cluster_name']]['progs'][job_specs['prog']]['job_instructions']
    
    # Initialize the dictionary that will be returned by the function
    
    rendered_content = {}
    
    # Render the template for the input file

    print("{:<80}".format("\nRendering the jinja template for the orca input file ..."), end="")
    
    render_vars = {
      "method" : config['method'],
      "basis_set" : config['basis-set'],
      "job_type" : config['job-type'],
      "charge" : config['charge'],
      "multiplicity" : config['multiplicity'],
      "coordinates" : file_data['atomic_coordinates']
    }
    
    rendered_content[rnd_input] = jinja_render(misc['path_tpl_dir'], tpl_inp, render_vars)
    
    print('%12s' % "[ DONE ]")

    # Render the template for the job instructions file

    print("{:<80}".format("\nRendering the jinja template for the orca job instructions file ..."), end="")
    
    render_vars = {
      "mol_name" : misc['mol_name'],
      "config_name" : misc['config_name'],
      "user_email" : config['user-email'],
      "mail_type" : config['mail-type'],
      "job_walltime" : job_specs['walltime'],
      "job_cores" : job_specs['cores'],
      "job_mem_per_cpu" : job_specs['mem_per_cpu'],
      "partition" : job_specs['partition'],
      "set_env" : clusters_cfg[job_specs['cluster_name']]['progs'][job_specs['prog']]['set_env'],
      "command" : clusters_cfg[job_specs['cluster_name']]['progs'][job_specs['prog']]['command'],
      "prog" : job_specs['prog']
    }
    
    rendered_content[rnd_inst] = jinja_render(misc['path_tpl_dir'], tpl_inst, render_vars)
    
    print('%12s' % "[ DONE ]")

    # Return the content of the rendered files
    
    return rendered_content
