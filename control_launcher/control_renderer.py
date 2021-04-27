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
import re

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

    # Check the process keyword

    process = config['qoctra'].get('process')

    if not process:
      raise control_errors.ControlError ('ERROR: There is no "process" key in the "qoctra" block of the "%s" configuration file.' % misc['config_name'])

    # Define the templates
    # ====================

    # Define the names of the default templates.

    template_param = "param.nml.jinja"
    template_script = "qoctra_job.sh.jinja"

    # Check if the specified templates exist in the "templates" directory of CONTROL LAUNCHER.
    
    control_errors.check_abspath(os.path.join(misc['templates_dir'],template_param),"Jinja template for the qoctra parameters files","file")
    control_errors.check_abspath(os.path.join(misc['templates_dir'],template_script),"Jinja template for the qoctra job script","file")

    # Other specific templates

    if process == 'OPC':
      template_pulse = "guess_pulse_OPC.jinja"
      control_errors.check_abspath(os.path.join(misc['templates_dir'],template_pulse),"Jinja template for the OPC guess pulse","file")
    if process == 'OPM':
      template_pulse = "guess_pulse_OPM.jinja"
      control_errors.check_abspath(os.path.join(misc['templates_dir'],template_pulse),"Jinja template for the OPM guess pulse","file")

    # Define rendered files
    # =====================

    # Define the names of the rendered files.

    rendered_script = "qoctra_job.sh"
    rendered_pulse = "guess_pulse"
    rendered_param = "param_" + misc['transition']['label'] + ".nml"
    rendered_param_pcp = "param_" + misc['transition']['label'] + "_PCP.nml"

    # Initialize the dictionary that will be returned by the function

    rendered_content = {}

    # ========================================================= #
    #           Rendering the control parameters file           #
    # ========================================================= #

    print("{:<80}".format("\nRendering the jinja template for the control parameters file ...  "), end="")

    # Defining the mandatory Jinja variables
    # ======================================

    # Variables not associated with the config file

    param_render_vars = {
      # GENERAL
      "source_name" : misc['source_name'],
      "process" : process, # Associated with the config file, but it has already been verified
      # DATA FILES
      "energies_file_path" : data['eigenvalues_path'],
      "momdip_e_path" : data['momdip_es_mtx_path'],
      "init_file_path" : data['init_path'],
      "target_file_path" : data['target_path']
    }

    # Check if a "control" block has been defined in the "qoctra" block of the config file

    if not config['qoctra'].get('control'):
      raise control_errors.ControlError ('ERROR: There is no "control" key in the "qoctra" block of the "%s" configuration file.' % misc['config_name'])    

    # Variables associated with the "control" block of the "qoctra" block in the config file

    try:
      
      time_step = config['qoctra']['control']['time_step'] # Stored into a variable since it will be reused
      time_step = re.compile(r'(\d*\.\d*)[dD]([-+]?\d+)').sub(r'\1E\2', time_step) # Replace the possible d/D from Fortran double precision float format with an "E", understandable by Python

      param_render_vars.update({
        # CONTROL
        "max_iter" : config['qoctra']['control']['max_iter'],
        "threshold" : config['qoctra']['control']['threshold'],
        "time_step" : time_step,
        "start_pulse" : " ",
        "guess_pulse" : rendered_pulse
      })

    except KeyError as error:
      raise control_errors.ControlError ('ERROR: The "%s" key is missing in the "control" block of the "qoctra" block in the "%s" configuration file.' % (error,misc['config_name']))

    # Defining the specific Jinja variables
    # =====================================

    if process == 'OPC':

      # Check if a "opc" block has been defined in the "qoctra" block of the config file

      if not config['qoctra'].get('opc'):
        raise control_errors.ControlError ('ERROR: There is no "opc" key in the "qoctra" block of the "%s" configuration file.' % misc['config_name'])    

      # Variables associated with the "opc" block of the "qoctra" block in the config file

      try:

        nb_steps = config['qoctra']['opc']['nb_steps'] # Stored into a variable since it will be reused

        param_render_vars.update({
          # OPC
          "nb_steps" : nb_steps,
          "alpha" : config['qoctra']['opc']['alpha'],
          "write_freq" : config['qoctra']['opc']['write_freq']
        })

      except KeyError as error:
        raise control_errors.ControlError ('ERROR: The "%s" key is missing in the "opc" block of the "qoctra" block in the "%s" configuration file.' % (error,misc['config_name']))

    #! Temporary

    param_render_vars.update({
      # POST CONTROL
      "mat_et0_path" : data['eigenvectors_path']
    })

    # Rendering the file
    # ==================
     
    rendered_content[rendered_param] = jinja_render(misc['templates_dir'], template_param, param_render_vars)

    print('%12s' % "[ DONE ]")

    # ========================================================= #
    #             Rendering the PCP parameters file             #
    # ========================================================= #

    print("{:<80}".format("\nRendering the jinja template for the PCP parameters file ...  "), end="")

    # Defining the Jinja variables (via updating the dictionary from the control parameters file)
    # ============================

    # Variables not associated with the config file

    param_render_vars.update({
      # GENERAL
      "process" : "PCP",
      # CONTROL
      "start_pulse" : "../Pulse/Pulse",
      # POST CONTROL
      "mat_et0_path" : data['eigenvectors_path']
    })

    # Check if a "post_control" block has been defined in the "qoctra" block of the config file

    if not config['qoctra'].get('post_control'):
      raise control_errors.ControlError ('ERROR: There is no "post_control" key in the "qoctra" block of the "%s" configuration file.' % misc['config_name'])    

    # Variables associated with the "post_control" block of the "qoctra" block in the config file

    try:
      param_render_vars.update({
        # POST CONTROL
        "analy_freq" : config['qoctra']['post_control']['analy_freq']
      })

    except KeyError as error:
      raise control_errors.ControlError ('ERROR: The "%s" key is missing in the "post_control" block of the "qoctra" block in the "%s" configuration file.' % (error,misc['config_name']))

    # Rendering the file
    # ==================
     
    rendered_content[rendered_param_pcp] = jinja_render(misc['templates_dir'], template_param, param_render_vars)

    print('%12s' % "[ DONE ]")

    # ========================================================= #
    #              Rendering the guess pulse file               #
    # ========================================================= #

    print("{:<80}".format("\nRendering the jinja template for the guess pulse file ...  "), end="")

    # Defining the Jinja variables for OPC
    # ====================================

    if process == 'OPC':

      # Set the variables associated with the config file
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      # Check if a "guess_pulse" block has been defined in the "qoctra" block of the config file

      if not config['qoctra'].get('guess_pulse'):
        raise control_errors.ControlError ('ERROR: There is no "guess_pulse" key in the "qoctra" block of the "%s" configuration file.' % misc['config_name'])    

      # Variables associated with the "guess_pulse" block of the "qoctra" block in the config file

      try:

        tbp = config['qoctra']['guess_pulse']['tbp']

        pulse_render_vars = {
          "pulse_type" : config['qoctra']['guess_pulse']['pulse_type'],
          "subpulse_type" : config['qoctra']['guess_pulse']['subpulse_type'],
          "amplitude" : config['qoctra']['guess_pulse']['amplitude'],
          "phase_change" : config['qoctra']['guess_pulse']['phase_change'],
          "sign" : config['qoctra']['guess_pulse']['sign']       
        }

      except KeyError as error:
        raise control_errors.ControlError ('ERROR: The "%s" key is missing in the "guess_pulse" block of the "qoctra" block in the "%s" configuration file.' % (error,misc['config_name']))

      # Determine the subpulses constituting the guess pulse
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      # Determine the spectral bandwidth of the pulse

      duration = int(nb_steps) * float(time_step) 
      duration_fwhm = duration / 6 * 2.35 # Using 2.35 as the standard deviation for a Gaussian
      bandwidth_fwhm = float(tbp) / duration_fwhm
      bandwidth = bandwidth_fwhm * 6 / 2.35

      # Determine the frequency range of subpulses (delimited by the spectral bandwidth, centered around the transition energy)

      upper_limit = misc['transition']['energy'] + bandwidth/2
      lower_limit = misc['transition']['energy'] - bandwidth/2

      # Initialize the list of subpulses

      subpulses = []

      # Consider each pair of excited states

      for state_1 in range(1, len(system['eigenstates_list'])): # Starting at 1 to exclude the ground state
        for state_2 in range(state_1 + 1, len(system['eigenstates_list'])): # Starting at "state_1 + 1" to exclude the cases where both states are the same

          # If the transition energy of this pair is comprised in our frequency range, add it to the subpulses list
         
          energy_diff = abs(system['eigenstates_list'][state_1]['energy'] - system['eigenstates_list'][state_2]['energy'])

          if lower_limit <= energy_diff <= upper_limit:
            subpulses.append(str(state_1 + 1) + " \t " + str(state_2 + 1)) # +1 because Fortran starts numbering at 1 while Python starts at 0.

      # Set the variables not associated with the config file
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      pulse_render_vars.update({
        "basis" : data['eigenvalues_path'],
        "nb_subpulses" : len(subpulses),
        "subpulses" : subpulses
      })

    # Rendering the file
    # ==================
     
    rendered_content[rendered_pulse] = jinja_render(misc['templates_dir'], template_pulse, pulse_render_vars)

    print('%12s' % "[ DONE ]")

    # ========================================================= #
    #                  Rendering the job script                 #
    # ========================================================= #

    # If we need to copy the output files to their respective results directory, load the CHAINS configuration file to get the necessary information

    if copy_files:

      chains_path = os.path.dirname(misc['code_dir']) 
      chains_config_file = control_errors.check_abspath(os.path.join(chains_path,"configs","chains_config.yml"),"CHAINS configuration YAML file","file")
  
      print ("{:<80}".format("\nLoading CHAINS configuration YAML file ..."), end="")
      with open(chains_config_file, 'r') as chains:
        chains_config = yaml.load(chains, Loader=yaml.FullLoader)
      print('%12s' % "[ DONE ]")

    print("{:<80}".format("\nRendering the jinja template for the qoctra job script ..."), end="")

    # Defining the mandatory Jinja variables
    # ======================================

    # Variables not associated with the config file

    script_render_vars = {
      "source_name" : misc['source_name'],
      "transition" : misc['transition']['label'],
      "config_name" : misc['config_name'],
      "job_walltime" : job_specs['walltime'],
      "job_memory" : job_specs['memory'], # in MB
      "partition" : job_specs['partition'],
      "rendered_param" : rendered_param,
      "rendered_param_PCP" : rendered_param_pcp,
      "guess_pulse" : rendered_pulse,
      "copy_files" : copy_files # Associated with the config file, but it has already been verified
    }

    # Variables associated with the "general" block of the config file

    try:
      script_render_vars.update({
        "user_email" : config['general']['user_email'],
        "mail_type" : config['general']['mail_type']
      })

    except KeyError as error:
      raise control_errors.ControlError ('ERROR: The "%s" key is missing in the "general" block of the "%s" configuration file.' % (error,misc['config_name']))

    # Variables associated with the clusters configuration file

    try:
      script_render_vars.update({
        "set_env" : clusters_cfg[job_specs['cluster_name']]['profiles'][job_specs['profile']]['set_env']     
      })

    except KeyError as error:
      raise control_errors.ControlError ('ERROR: The "%s" key is missing in the "%s" profile of the clusters configuration file.' % (error,job_specs['profile']))

    # Defining the specific Jinja variables
    # =====================================

    # Variables specific to the copy_files portion of the template

    if copy_files:

      # Variables not associated with the config file

      script_render_vars.update({
        "mol_dir" : misc['mol_dir'],
        "nb_transitions" : len(misc['transitions_list']),
        "data_dir" : data['main_path'],
        "job_script" : rendered_script
      })

      # Variables associated with the CHAINS configuration file

      try:
        script_render_vars.update({
          "output_dir" : chains_config['output_qoctra'],
          "results_dir" : chains_config['results_dir']
        })

      except KeyError as error:
        raise control_errors.ControlError ('ERROR: The "%s" key is missing in the CHAINS configuration file (chains_config.yml).' % error)

    # Rendering the file
    # ==================

    rendered_content[rendered_script] = jinja_render(misc['templates_dir'], template_script, script_render_vars)

    print('%12s' % "[ DONE ]")

    return rendered_content, rendered_script