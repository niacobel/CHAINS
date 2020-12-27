################################################################################################################################################
##                                                                The Renderer                                                                ##
##                                                                                                                                            ##
##                                     This script contains the rendering function for CONTROL LAUNCHER,                                      ##
##                                consult the documentation at https://chains-ulb.readthedocs.io/ for details                                 ##
################################################################################################################################################

import os

import yaml
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

def chains_qoctra_render(clusters_cfg:dict, config:dict, system:dict, data:dict, job_specs:dict, misc:dict):
    """Renders the job script and the two parameters files (OPM & PCP) associated with the QOCT-RA program in CHAINS.

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
    # ========================================================= #
    #                      Preparation step                     #
    # ========================================================= #

    # Check config file
    # =================

    # Check if a "general" block has been defined in the config file

    if not config.get('general'):
      raise control_errors.ControlError ('ERROR: There is no "general" key defined in the "%s" configuration file.' % misc['config_name'])     

    # Check if a "qoctra" block has been defined in the config file

    if not config.get('qoctra'):
      raise control_errors.ControlError ('ERROR: There is no "qoctra" key defined in the "%s" configuration file.' % misc['config_name'])      

    # Check the options defined in the config file

    copy_files = config['qoctra'].get('copy_files',True)

    if not isinstance(copy_files, bool):
      raise control_errors.ControlError ('ERROR: The "copy_files" value given in the "qoctra" block of the "%s" configuration file is not a boolean (neither "True" nor "False").' % misc['config_name'])

    # Define the templates
    # ====================

    # Define the names of the templates.

    template_param = "param.nml.jinja"
    template_script = "qoctra_job.sh.jinja"

    # Check if the specified templates exist in the "templates" directory of CONTROL LAUNCHER.
    
    control_errors.check_abspath(os.path.join(misc['templates_dir'],template_param),"Jinja template for the qoctra parameters files","file")
    control_errors.check_abspath(os.path.join(misc['templates_dir'],template_script),"Jinja template for the qoctra job script","file")

    # Define rendered files
    # =====================

    # Define the names of the rendered files.

    rendered_param_opm = "param_" + misc['transition_label'] + "_OPM.nml"
    rendered_param_pcp = "param_" + misc['transition_label'] + "_PCP.nml"
    rendered_script = "qoctra_job.sh"

    # Initialize the dictionary that will be returned by the function

    rendered_content = {}

    # Load the CHAINS configuration file to get the additional information

    chains_path = os.path.dirname(misc['code_dir'])
    chains_config_file = control_errors.check_abspath(os.path.join(chains_path,"chains_config.yml"),"CHAINS configuration YAML file","file")

    print ("{:<80}".format("\nLoading CHAINS configuration YAML file ..."), end="")
    with open(chains_config_file, 'r') as chains:
      chains_config = yaml.load(chains, Loader=yaml.FullLoader)
    print('%12s' % "[ DONE ]")

    # ========================================================= #
    #             Rendering the OPM parameters file             #
    # ========================================================= #

    print("{:<80}".format("\nRendering the jinja template for the OPM parameters file ...  "), end="")

    # Determine the central frequency of the guess pulse in cm-1 (here defined as the average of the eigenvalues)

    central_frequency = sum(system['eigenvalues']) / (len(system['eigenvalues']) - 1) # -1 because the ground state doesn't count
    
    # Define here the number of iterations for QOCT-RA, as it will be used multiple times later on

    niter = config['qoctra']['control']['niter']

    # Defining the Jinja variables

    param_render_vars = {
      "source_name" : misc['source_name'],
      "energies_file_path" : data['eigenvalues_path'],
      "momdip_e_path" : data['momdip_es_mtx_path'],
      "init_file_path" : data['init_path'],
      "final_file_path" : os.path.join(data['main_path'],"final_"),
      "proj_file_path" : data['target_path'],
      "transition" : misc['transition_label'],
      "nstep" : config['qoctra']['control']['nstep'],
      "dt" : config['qoctra']['control']['dt'],
      "processus" : "OPM",
      "source" : " ",
      "niter" : niter,
      "threshold" : config['qoctra']['control']['threshold'],
      "alpha0" : config['qoctra']['control']['alpha0'],
      "ndump" : config['qoctra']['control']['ndump'],
      "ndump2" : config['qoctra']['post_control']['ndump2'],
      "mat_et0_path" : data['eigenvectors_path'],
      "numericincrements" : config['qoctra']['guess_pulse']['numericincrements'],
      "numberofpixels" : config['qoctra']['guess_pulse']['numberofpixels'],
      "inputenergy" : config['qoctra']['guess_pulse']['inputenergy'],
      "widthhalfmax" : config['qoctra']['guess_pulse']['widthhalfmax'],
      "omegazero" : central_frequency
    }

    # Rendering the file
     
    rendered_content[rendered_param_opm] = jinja_render(misc['templates_dir'], template_param, param_render_vars)

    print('%12s' % "[ DONE ]")

    # ========================================================= #
    #             Rendering the PCP parameters file             #
    # ========================================================= #

    print("{:<80}".format("\nRendering the jinja template for the PCP parameters file ...  "), end="")

    # Defining the Jinja variables by updating the dictionary defined for the OPM parameters file

    param_render_vars.update({
      "processus" : "PCP",
      "source" : "../Pulse/Pulse_iter" + str(niter)
    })

    # Rendering the file
     
    rendered_content[rendered_param_pcp] = jinja_render(misc['templates_dir'], template_param, param_render_vars)

    print('%12s' % "[ DONE ]")

    # ========================================================= #
    #                  Rendering the job script                 #
    # ========================================================= #

    print("{:<80}".format("\nRendering the jinja template for the qoctra job script ..."), end="")

    # Defining the Jinja variables

    script_render_vars = {
      "source_name" : misc['source_name'],
      "transition" : misc['transition_label'],
      "user_email" : config['general']['user_email'],
      "mail_type" : config['general']['mail_type'],
      "job_walltime" : job_specs['walltime'],
      "job_memory" : job_specs['memory'], # in MB
      "partition" : job_specs['partition'],
      "rendered_param" : rendered_param_opm,
      "rendered_param_PCP" : rendered_param_pcp,
      "set_env" : clusters_cfg[job_specs['cluster_name']]['profiles'][job_specs['profile']]['set_env'],
      "mol_dir" : misc['mol_dir'],
      "nb_transitions" : len(misc['transitions_list']),
      "output_dir" : chains_config['output_dir']['qoctra'],
      "results_dir" : config['results']['main_dir'],
      "results_subdir" : config['results']['qoctra']['dir_name'],
      "data_dir" : data['main_path'],
      "job_dirname" : misc['transition_label'] + "_" + misc['config_name'],
      "job_script" : rendered_script,
      "niter" : niter
    }
    
    # Rendering the file

    rendered_content[rendered_script] = jinja_render(misc['templates_dir'], template_script, script_render_vars)

    print('%12s' % "[ DONE ]")

    return rendered_content, rendered_script