################################################################################################################################################
##                                                                The Renderer                                                                ##
##                                                                                                                                            ##
##                                     This script contains the rendering function for CONTROL LAUNCHER,                                      ##
##                                consult the documentation at https://chains-ulb.readthedocs.io/ for details                                 ##
################################################################################################################################################

import os

from jinja2 import Environment, FileSystemLoader

import control_errors


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

def sample_qoctra_render(clusters_cfg:dict, config:dict, system:dict, data:dict, job_specs:dict, misc:dict):
    """Renders the job script and the parameters file associated with the QOCT-GRAD program.

    Parameters
    ----------
    clusters_cfg : dict
        Content of the YAML clusters configuration file.
    config : dict
        Content of the YAML configuration file.
    system : dict
        Information extracted by the parsing function and derived from it.
    data : dict
        Data directory path and the name of its files.
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

    # Define the names of the templates

    template_param = "sample_param.nml.jinja"
    template_script = "sample_qoctra_job.sh.jinja"

    # Define the names of the rendered files

    rendered_param = "param.nml"
    rendered_script = "qoctra_job.sh"

    # Initialize the dictionary that will be returned by the function

    rendered_content = {}

    # Render the template for the parameters file

    print("{:<80}".format("\nRendering the jinja template for the parameters file ...  "), end="")

    param_render_vars = {
      "source_name" : misc['source_name'],
      "energies_file_path" : os.path.join(data['path'],data['eigenvalues_file']),
      "momdip_e_path" : os.path.join(data['path'],data['momdip_es_mtx_file']),
      "init_file_path" : os.path.join(data['path'],data['transition']['init_file']),
      "final_file_path" : os.path.join(data['path'],"final"),
      "proj_file_path" : os.path.join(data['path'],data['transition']['target_file']),
      "transition" : data['transition']['label'],
      "nstep" : config[job_specs['profile']]['control']['nstep'],
      "dt" : config[job_specs['profile']]['control']['dt'],
      "processus" : "OPM",
      "source" : " ",
      "niter" : config[job_specs['profile']]['control']['niter'],
      "threshold" : config[job_specs['profile']]['control']['threshold'],
      "alpha0" : config[job_specs['profile']]['control']['alpha0'],
      "ndump" : config[job_specs['profile']]['control']['ndump'],
      "numericincrements" : config[job_specs['profile']]['guess_pulse']['numericincrements'],
      "numberofpixels" : config[job_specs['profile']]['guess_pulse']['numberofpixels'],
      "inputenergy" : config[job_specs['profile']]['guess_pulse']['inputenergy'],
      "widthhalfmax" : config[job_specs['profile']]['guess_pulse']['widthhalfmax'],
      "omegazero" : sum(system['eigenvalues']) / (len(system['eigenvalues']) - 1)
    }
     
    rendered_content[rendered_param] = jinja_render(misc['templates_dir'], template_param, param_render_vars)

    print('%12s' % "[ DONE ]")

    # Render the template for the job script

    print("{:<80}".format("\nRendering the jinja template for the qoctra job script ..."), end="")

    script_render_vars = {
      "source_name" : misc['source_name'],
      "transition" : data['transition']['label'],
      "config_name" : misc['config_name'],
      "user_email" : config['user_email'],
      "mail_type" : config['mail_type'],
      "job_walltime" : job_specs['walltime'],
      "job_memory" : job_specs['memory'], # in MB
      "partition" : job_specs['partition'],
      "set_env" : clusters_cfg[job_specs['cluster_name']]['profiles'][job_specs['profile']]['set_env'],
      "rendered_param" : rendered_param
    }
    
    rendered_content[rendered_script] = jinja_render(misc['templates_dir'], template_script, script_render_vars)

    print('%12s' % "[ DONE ]")

    # Return the content of the rendered files and the name of the rendered job script

    return rendered_content, rendered_script