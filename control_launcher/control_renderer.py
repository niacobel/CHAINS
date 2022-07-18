################################################################################################################################################
##                                                                The Renderer                                                                ##
##                                                                                                                                            ##
##                                     This script contains the rendering function for CONTROL LAUNCHER,                                      ##
##                                consult the documentation at https://chains-ulb.readthedocs.io/ for details                                 ##
################################################################################################################################################

import csv
import math
import os
import re

import numpy as np
import yaml
from jinja2 import Environment, FileSystemLoader
from scipy.spatial import ConvexHull, distance
from scipy import constants

import control_common


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

def alpha_duration_param_search(clusters_cfg:dict, config:dict, system:dict, data:dict, job_specs:dict, misc:dict):
    """Renders the job script and the parameters files for each combination of two parameters values: the penalty factor (alpha) and the duration of the pulse (given by the number of steps for a specified time step). The job script will launch a job array where each job is a unique combination of an alpha value and a duration value. A PCP parameter file and a guess pulse file are also rendered but are independant of the chosen alpha and duration values.

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
      raise control_common.ControlError ('ERROR: There is no "general" key defined in the "%s" configuration file.' % misc['config_name'])     

    # Check if a "qoctra" block has been defined in the config file

    if not config.get('qoctra'):
      raise control_common.ControlError ('ERROR: There is no "qoctra" key defined in the "%s" configuration file.' % misc['config_name'])      

    # Check if a "parameters" block has been defined in the config file

    if not config.get('parameters'):
      raise control_common.ControlError ('ERROR: There is no "parameters" key defined in the "%s" configuration file.' % misc['config_name'])  

    # Check the options defined in the config file

    copy_files = config['qoctra'].get('copy_files',True)

    if not isinstance(copy_files, bool):
      raise control_common.ControlError ('ERROR: The "copy_files" value given in the "qoctra" block of the "%s" configuration file is not a boolean (neither "True" nor "False").' % misc['config_name'])

    # Define the templates
    # ====================

    # Define the names of the default templates.

    template_param = "param.nml.jinja"
    template_script = "aldu_param_search_job.sh.jinja"
    template_pulse = "guess_pulse_OPC.jinja"
    template_treatment = "aldu_param_treatment_job.sh.jinja"

    # Check if the specified templates exist in the "templates" directory of CONTROL LAUNCHER.
    
    control_common.check_abspath(os.path.join(misc['templates_dir'],template_param),"Jinja template for the qoctra parameters files","file")
    control_common.check_abspath(os.path.join(misc['templates_dir'],template_script),"Jinja template for the qoctra job script","file")
    control_common.check_abspath(os.path.join(misc['templates_dir'],template_pulse),"Jinja template for the OPC guess pulse","file")
    control_common.check_abspath(os.path.join(misc['templates_dir'],template_treatment),"Jinja template for the treatment job script","file")

    # Define rendered files
    # =====================

    # Define the names of the rendered files.

    rendered_script = "aldu_param_search_job.sh"
    rendered_treatment = "aldu_param_treatment_job.sh"
    rendered_pulse = "guess_pulse"
    rendered_param_pcp = "PCP_param.nml"

    # Define the prefix for the parameters files.

    prefix_param = "param_"

    # Define the name of the text file containing all the parameters filenames

    input_names = "input_filenames.txt"

    # Initialize the dictionary that will be returned by the function

    rendered_content = {}

    # Load CHAINS configuration file if needed
    # ========================================

    chains_path = os.path.dirname(misc['code_dir']) 

    if copy_files:
      
      chains_config_file = control_common.check_abspath(os.path.join(chains_path,"configs","chains_config.yml"),"CHAINS configuration YAML file","file")

      print ("{:<80}".format("\nLoading CHAINS configuration YAML file ..."), end="")
      with open(chains_config_file, 'r') as chains:
        chains_config = yaml.load(chains, Loader=yaml.FullLoader)
      print('%12s' % "[ DONE ]")

    # ========================================================= #
    #           Rendering the control parameters files          #
    # ========================================================= #

    print("{:<80}".format("\nRendering the jinja template for the control parameters files ...  "), end="")

    # Defining the main Jinja variables
    # =================================

    # Variables not associated with the config file

    param_render_vars = {
      # GENERAL
      "source_name" : misc['source_name'],
      "process" : "OPC",
      "energies_file_path" : data['energies_path'],
      "momdip_path" : data['momdip_mtx_path'],
      "init_file_path" : data['init_path'],
      "target_file_path" : data['target_path'],
      "projector_path" : data.get('projector_path',"no")
    }

    # Check if a "control" block has been defined in the "qoctra" block of the config file

    if not config['qoctra'].get('control'):
      raise control_common.ControlError ('ERROR: There is no "control" key in the "qoctra" block of the "%s" configuration file.' % misc['config_name'])    

    # Variables associated with the "control" block of the "qoctra" block in the config file

    try:
      
      time_step = float(config['qoctra']['control']['time_step'])
      time_step_for = "{:.5e}".format(time_step).replace('e','d') # Replace the 'e' from Python with the 'd' from Fortran double precision

      threshold = float(config['qoctra']['control']['threshold'])
      threshold_for = "{:.5e}".format(threshold).replace('e','d')

      conv_thresh = float(config['qoctra']['control']['conv_thresh'])
      conv_thresh_for = "{:.5e}".format(conv_thresh).replace('e','d')

      param_render_vars.update({
        # CONTROL
        "max_iter" : config['qoctra']['control']['max_iter'],
        "threshold" : threshold_for,
        "conv_thresh" : conv_thresh_for,
        "time_step" : time_step_for,
        "start_pulse" : " ",
        "guess_pulse" : rendered_pulse,
        "write_freq" : config['qoctra']['control']['write_freq']
      })

    except KeyError as error:
      raise control_common.ControlError ('ERROR: The "%s" key is missing in the "control" block of the "qoctra" block in the "%s" configuration file.' % (error,misc['config_name']))

    #! Temporary

    param_render_vars.update({
      # POST CONTROL
      "conv_mtx_path" : data['conv_mtx_path']
    })

    # Handle the parameters
    # =====================

    # Check the parameters

    param_keys = ["alpha","duration"]

    for key in param_keys:
      if not config['parameters'].get(key):
        raise control_common.ControlError ('ERROR: There is no "%s" key in the "parameters" block of the "%s" configuration file.' % (key,misc['config_name']))
      if not isinstance(config['parameters'][key], list):
        raise control_common.ControlError ('ERROR: The value of the "%s" key in the "parameters" block of the "%s" configuration file is not a list.' % (key,misc['config_name']))
      if not all(isinstance(value,(int,float)) for value in config['parameters'][key]):
        raise control_common.ControlError ('ERROR: Some or all of the values in the list of the "%s" key in the "parameters" block of the "%s" configuration file are neither integer nor floats.' % (key,misc['config_name']))

    # Iterate over the parameters value unique combinations

    alpha_list = config['parameters']['alpha']
    duration_list = config['parameters']['duration'] # in ps
    filenames = []

    for alpha in alpha_list:
      for duration in duration_list:

        #! Check if the total duration of the pulse is not too long compared to the lifetime of our most populated target state

        alpha_for = "{:.5e}".format(alpha).replace('e','d')
        nb_steps = round((duration * 1e-12/constants.value('atomic unit of time'))/time_step)

        param_render_vars.update({
          # OPC
          "alpha" : alpha_for,
          "nb_steps" : nb_steps
        })

        # Rendering the file for those parameters
        # =======================================

        rendered_param = prefix_param + "alpha" + str(alpha) + "_dur" + str(duration) + ".nml"
        rendered_content[rendered_param] = jinja_render(misc['templates_dir'], template_param, param_render_vars)

        # Add the name of this file to the list

        filenames.append(rendered_param)

    # Add to the dictionary the content of the text file containing all the parameters filenames

    rendered_content[input_names] = "\n".join(filenames)

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
      "conv_mtx_path" : data['conv_mtx_path']
    })

    # Check if a "post_control" block has been defined in the "qoctra" block of the config file

    if not config['qoctra'].get('post_control'):
      raise control_common.ControlError ('ERROR: There is no "post_control" key in the "qoctra" block of the "%s" configuration file.' % misc['config_name'])    

    # Variables associated with the "post_control" block of the "qoctra" block in the config file

    try:
      param_render_vars.update({
        # POST CONTROL
        "analy_freq" : config['qoctra']['post_control']['analy_freq']
      })

    except KeyError as error:
      raise control_common.ControlError ('ERROR: The "%s" key is missing in the "post_control" block of the "qoctra" block in the "%s" configuration file.' % (error,misc['config_name']))

    # Rendering the files (one for each orientation)
    # ==============================================

    for momdip_key in system['momdip_mtx']:

      param_render_vars.update({
        # GENERAL
        "momdip_path" : os.path.join(data['main_path'],'momdip_mtx_' + momdip_key)
      })

      rendered_content[momdip_key + "_" + rendered_param_pcp] = jinja_render(misc['templates_dir'], template_param, param_render_vars)

    print('%12s' % "[ DONE ]")

    # ========================================================= #
    #              Rendering the guess pulse file               #
    # ========================================================= #

    print("{:<80}".format("\nRendering the jinja template for the guess pulse file ...  "), end="")

    # Defining the Jinja variables
    # ============================

    # Check if a "guess_pulse" block has been defined in the "qoctra" block of the config file

    if not config['qoctra'].get('guess_pulse'):
      raise control_common.ControlError ('ERROR: There is no "guess_pulse" key in the "qoctra" block of the "%s" configuration file.' % misc['config_name'])    

    # Variables associated with the "guess_pulse" block of the "qoctra" block in the config file

    try:

      amplitude = float(config['qoctra']['guess_pulse']['amplitude']) / constants.value('atomic unit of electric field')
      amplitude_for = "{:.5e}".format(amplitude).replace('e','d') # Replace the 'e' from Python with the 'd' from Fortran double precision

      phase_change = float(config['qoctra']['guess_pulse']['phase_change'])
      phase_change_for = "{:.5e}".format(phase_change).replace('e','d')

      pulse_render_vars = {
        "amplitude" : amplitude_for,
        "pulse_type" : config['qoctra']['guess_pulse']['pulse_type'],
        "subpulse_type" : config['qoctra']['guess_pulse']['subpulse_type'],
        "phase_change" : phase_change_for,
        "sign" : config['qoctra']['guess_pulse']['sign']
      }

    except KeyError as error:
      raise control_common.ControlError ('ERROR: The "%s" key is missing in the "guess_pulse" block of the "qoctra" block in the "%s" configuration file.' % (error,misc['config_name']))

    # Determine the subpulses constituting the guess pulse
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subpulses = []
    max_energy_diff = 0

    # Consider each pair of excited states

    for state_1 in range(len(system['states_list'])):
      if state_1 == 0:
        continue # Skip the ground state
      for state_2 in range(state_1 + 1, len(system['states_list'])): # Starting at "state_1 + 1" to exclude the cases where both states are the same

        # Add it to the subpulses list

        subpulses.append(str(state_1 + 1) + " \t " + str(state_2 + 1)) # +1 because Fortran starts numbering at 1 while Python starts at 0.

        # Compute the energy difference

        energy_diff = abs(system['states_list'][state_1]['energy'] - system['states_list'][state_2]['energy'])

        # Store the maximum transition energy to apply Nyquist–Shannon sampling theorem

        if energy_diff > max_energy_diff:
          max_energy_diff = energy_diff

    # Check the time step
    # ~~~~~~~~~~~~~~~~~~~

    # Nyquist–Shannon sampling theorem: check if the sampling rate is bigger than the double of the highest frequency

    sampling_freq = 1 / (time_step * constants.value('atomic unit of time'))
    max_freq = control_common.energy_unit_conversion(max_energy_diff,'ha','hz')
    
    if not sampling_freq > (2 * max_freq):
      raise control_common.ControlError ('ERROR: The time step value (%s a.u.) is too big for this transition. The sampling rate (%e Hz) is not bigger than the double of the highest frequency of the signal (%e Hz).' % (time_step,sampling_freq, max_freq))

    # Set the variables not associated with the config file
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    pulse_render_vars.update({
      "basis" : data['energies_path'],
      "nb_subpulses" : len(subpulses),
      "subpulses" : subpulses
    })

    # Rendering the file
    # ==================
     
    rendered_content[rendered_pulse] = jinja_render(misc['templates_dir'], template_pulse, pulse_render_vars)

    print('%12s' % "[ DONE ]")

    # ========================================================= #
    #               Rendering the main job script               #
    # ========================================================= #

    print("{:<80}".format("\nRendering the jinja template for the qoctra job script ..."), end="")

    # Defining the mandatory Jinja variables
    # ======================================

    # Variables not associated with the config file

    array_size = (len(alpha_list) * len(duration_list)) - 1 # array begins at 0 in the template

    script_render_vars = {
      "source_name" : misc['source_name'],
      "transition" : misc['transition']['label'],
      "config_name" : misc['config_name'],
      "job_walltime" : job_specs['walltime'],
      "job_memory" : job_specs['memory'], # in MB
      "partition" : job_specs['partition'],
      "array_size" : array_size,
      "cluster_name": job_specs['cluster_name'],
      "treatment_script" : rendered_treatment,
      "input_names": input_names,
      "prefix_param" : prefix_param,
      "momdip_keys" : list(system['momdip_mtx'].keys()),
      "rendered_param_PCP" : rendered_param_pcp,
      "guess_pulse" : rendered_pulse,
      "profile" : job_specs['profile']
    }

    # Variables associated with the "general" block of the config file

    try:
      script_render_vars.update({
        "user_email" : config['general']['user_email'],
        "mail_type" : config['general']['mail_type']
      })

    except KeyError as error:
      raise control_common.ControlError ('ERROR: The "%s" key is missing in the "general" block of the "%s" configuration file.' % (error,misc['config_name']))

    # Variables associated with the clusters configuration file

    try:
      script_render_vars.update({
        "set_env" : clusters_cfg[job_specs['cluster_name']]['profiles'][job_specs['profile']]['set_env']     
      })

    except KeyError as error:
      raise control_common.ControlError ('ERROR: The "%s" key is missing in the "%s" profile of the clusters configuration file.' % (error,job_specs['profile']))

    # Rendering the file
    # ==================

    rendered_content[rendered_script] = jinja_render(misc['templates_dir'], template_script, script_render_vars)

    print('%12s' % "[ DONE ]")

    # ========================================================= #
    #             Rendering the treatment job script            #
    # ========================================================= #

    print("{:<80}".format("\nRendering the jinja template for the treatment job script ..."), end="")

    # Defining the mandatory Jinja variables
    # ======================================

    # Variables not associated with the config file

    treatment_render_vars = {
      "source_name" : misc['source_name'],
      "transition" : misc['transition']['label'],
      "config_name" : misc['config_name'],
      "cluster_name": job_specs['cluster_name'],
      "profile" : job_specs['profile'],
      "chains_dir" : chains_path,
      "copy_files" : copy_files # Associated with the config file, but it has already been verified
    }

    # Variables associated with the "general" block of the config file

    try:
      treatment_render_vars.update({
        "user_email" : config['general']['user_email'],
        "mail_type" : config['general']['mail_type']
      })

    except KeyError as error:
      raise control_common.ControlError ('ERROR: The "%s" key is missing in the "general" block of the "%s" configuration file.' % (error,misc['config_name']))

    # Defining the specific Jinja variables
    # =====================================

    # Variables specific to the copy_files portion of the template

    if copy_files:

      # Variables not associated with the config file

      treatment_render_vars.update({
        "data_dir" : data['main_path'],
        "momdip_keys" : list(system['momdip_mtx'].keys()),
        "rendered_param_PCP" : rendered_param_pcp,
        "guess_pulse" : rendered_pulse,
        "job_script" : rendered_script,
        "treatment_script" : rendered_treatment,
        "momdip_key" : misc['transition']['momdip_key'],
        "pro_dir" : misc['pro_dir']
      })

      # Variables associated with the CHAINS configuration file

      try:
        treatment_render_vars.update({
          "output_dir" : chains_config['output_aldu'],
          "results_dir" : chains_config['results_dir']
        })

      except KeyError as error:
        raise control_common.ControlError ('ERROR: The "%s" key is missing in the CHAINS configuration file (chains_config.yml).' % error)

    # Rendering the file
    # ==================

    rendered_content[rendered_treatment] = jinja_render(misc['templates_dir'], template_treatment, treatment_render_vars)

    print('%12s' % "[ DONE ]")

    return rendered_content, rendered_script


######################################################################################################################################


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
      raise control_common.ControlError ('ERROR: There is no "general" key defined in the "%s" configuration file.' % misc['config_name'])     

    # Check if a "qoctra" block has been defined in the config file

    if not config.get('qoctra'):
      raise control_common.ControlError ('ERROR: There is no "qoctra" key defined in the "%s" configuration file.' % misc['config_name'])      

    # Check the options defined in the config file

    copy_files = config['qoctra'].get('copy_files',True)

    if not isinstance(copy_files, bool):
      raise control_common.ControlError ('ERROR: The "copy_files" value given in the "qoctra" block of the "%s" configuration file is not a boolean (neither "True" nor "False").' % misc['config_name'])

    # Check the process keyword

    process = config['qoctra'].get('process')

    if not process:
      raise control_common.ControlError ('ERROR: There is no "process" key in the "qoctra" block of the "%s" configuration file.' % misc['config_name'])

    # Define the templates
    # ====================

    # Define the names of the default templates.

    template_param = "param.nml.jinja"
    template_script = "qoctra_job.sh.jinja"

    # Check if the specified templates exist in the "templates" directory of CONTROL LAUNCHER.
    
    control_common.check_abspath(os.path.join(misc['templates_dir'],template_param),"Jinja template for the qoctra parameters files","file")
    control_common.check_abspath(os.path.join(misc['templates_dir'],template_script),"Jinja template for the qoctra job script","file")

    # Other specific templates

    if process == 'OPC':
      template_pulse = "guess_pulse_OPC.jinja"
      control_common.check_abspath(os.path.join(misc['templates_dir'],template_pulse),"Jinja template for the OPC guess pulse","file")
    if process == 'OPM':
      template_pulse = "guess_pulse_OPM.jinja"
      control_common.check_abspath(os.path.join(misc['templates_dir'],template_pulse),"Jinja template for the OPM guess pulse","file")

    # Define rendered files
    # =====================

    # Define the names of the rendered files.

    rendered_script = "qoctra_job.sh"
    rendered_pulse = "guess_pulse"
    rendered_param = "param_" + misc['transition']['label'] + ".nml"
    rendered_param_pcp = "param_" + misc['transition']['label'] + "_PCP.nml"

    # Initialize the dictionary that will be returned by the function

    rendered_content = {}

    # Load CHAINS configuration file
    # ==============================

    chains_path = os.path.dirname(misc['code_dir']) 
    chains_config_file = control_common.check_abspath(os.path.join(chains_path,"configs","chains_config.yml"),"CHAINS configuration YAML file","file")

    print ("{:<80}".format("\nLoading CHAINS configuration YAML file ..."), end="")
    with open(chains_config_file, 'r') as chains:
      chains_config = yaml.load(chains, Loader=yaml.FullLoader)
    print('%12s' % "[ DONE ]")

    # Load ionization potentials CSV file
    # ===================================

    ip_file = control_common.check_abspath(chains_config['ip_file'],"Ionization potentials CSV file","file")

    print ("{:<80}".format("\nLoading ionization potentials CSV file ..."), end="")
    with open(ip_file, 'r', newline='') as csv_file:
      ip_content = csv.DictReader(csv_file, delimiter=';')
      ip_list = list(ip_content)
    print('%12s' % "[ DONE ]")

    # Load the pulse shapers file
    # ===========================

    shapers_file = control_common.check_abspath(os.path.join(chains_path,"shapers.yml"),"Pulse shapers YAML file","file")

    print ("{:<80}".format("\nLoading the pulse shapers YAML file ..."), end="")
    with open(shapers_file, 'r', encoding='utf-8') as shape:
      shapers = yaml.load(shape, Loader=yaml.FullLoader)
    print('%12s' % "[ DONE ]")

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
      "energies_file_path" : data['energies_path'],
      "momdip_path" : data['momdip_mtx_path'],
      "init_file_path" : data['init_path'],
      "target_file_path" : data['target_path']
    }

    # Check if a "control" block has been defined in the "qoctra" block of the config file

    if not config['qoctra'].get('control'):
      raise control_common.ControlError ('ERROR: There is no "control" key in the "qoctra" block of the "%s" configuration file.' % misc['config_name'])    

    # Variables associated with the "control" block of the "qoctra" block in the config file

    try:
      
      time_step_raw = config['qoctra']['control']['time_step'] # Stored into a variable since it will be reused
      time_step = float(re.compile(r'(\d*\.\d*)[dD]([-+]?\d+)').sub(r'\1E\2', time_step_raw)) # Replace the possible d/D from Fortran double precision float format with an "E", understandable by Python

      param_render_vars.update({
        # CONTROL
        "max_iter" : config['qoctra']['control']['max_iter'],
        "threshold" : config['qoctra']['control']['threshold'],
        "time_step" : time_step_raw,
        "start_pulse" : " ",
        "guess_pulse" : rendered_pulse
      })

    except KeyError as error:
      raise control_common.ControlError ('ERROR: The "%s" key is missing in the "control" block of the "qoctra" block in the "%s" configuration file.' % (error,misc['config_name']))

    # Central frequency
    # ~~~~~~~~~~~~~~~~~

    # Compute the transition energy that will act as the central frequency

    init_number = np.argmax(np.diagonal(misc['transition']['init_content']))
    target_number = np.argmax(np.diagonal(misc['transition']['target_content']))
    omegazero = abs(system['states_list'][init_number]['energy'] - system['states_list'][target_number]['energy'])

    # Convert the central frequency to nanometers and to wavenumbers

    omegazero_nm = control_common.energy_unit_conversion(omegazero,'ha','nm')
    omegazero_cm = control_common.energy_unit_conversion(omegazero,'ha','cm-1')

    # Choosing the pulse shaper
    # ~~~~~~~~~~~~~~~~~~~~~~~~~

    min_diff = float('inf')

    # Find the shaper with the closest frequency to our central frequency

    for shaper in shapers:
      for profile in shaper['delay']:
        diff = abs(profile['frequency'] - omegazero_nm)
        if diff < min_diff:
          min_diff = diff
          ref_shaper = shaper
          ref_frequency = profile['frequency']
          ref_duration = profile['duration']
    
    # Pulse duration
    # ~~~~~~~~~~~~~~
    
    # Fix the duration of the pulse using the reference shaper

    duration_s = ref_duration * 1e-12
    duration = duration_s / constants.value('atomic unit of time')

    # Defining the specific Jinja variables
    # =====================================

    if process == 'OPC':

      # Compute the number of steps based on the total duration of the pulse and the time step

      nb_steps = round(duration/time_step)

      # Check if a "opc" block has been defined in the "qoctra" block of the config file

      if not config['qoctra'].get('opc'):
        raise control_common.ControlError ('ERROR: There is no "opc" key in the "qoctra" block of the "%s" configuration file.' % misc['config_name'])    

      # Variables associated with the "opc" block of the "qoctra" block in the config file

      try:

        bandwidth = config['qoctra']['opc']['bandwidth']

        spectral_filter = config['qoctra']['opc'].get('spectral_filter', None)
        max_fluence = config['qoctra']['opc'].get('max_fluence', False)

        param_render_vars.update({

          'spectral_filter' : spectral_filter.upper() if spectral_filter else "None",
          'max_fluence' : max_fluence,

          # OPC
          "nb_steps" : nb_steps,
          "write_freq" : config['qoctra']['opc']['write_freq']

        })

        if spectral_filter:

          # Check the form of the filter function and act accordingly

          if spectral_filter.lower() == 'sgw':
            param_render_vars.update({
              # OPC
              "spectral_filter_center" : "{:.5e}".format(omegazero_cm).replace('e','d'),
              "spectral_filter_fwhm" : "{:.5e}".format(bandwidth).replace('e','d')
            })

          elif spectral_filter.lower() == 'spgw':

            full_bandwidth = ( bandwidth * 6 ) / 2.35 # Conversion from FWHM to full width for a gaussian

            param_render_vars.update({
              # OPC
              "spectral_filter_center" : "{:.5e}".format(omegazero_cm).replace('e','d'),
              "spectral_filter_fwhm" : "{:.5e}".format(full_bandwidth).replace('e','d')
            })

          elif spectral_filter.lower() == 'mgw':
            param_render_vars.update({
              # OPC
              "spectral_filter_fwhm" : "30.d0"
            })

          else:
            raise control_common.ControlError ('ERROR: The given value for the "spectral_filter" key (%s) of the "opc" block in the "qoctra" block from the "%s" configuration file is not supported. Supported values include: SGW, SPGW, MGW and None (This is not case sensitive).' % (spectral_filter,misc['config_name']))

        if not max_fluence:
          param_render_vars.update({
            # OPC
            "alpha" : config['qoctra']['opc']['alpha']
          })
        
      except KeyError as error:
        raise control_common.ControlError ('ERROR: The "%s" key is missing in the "opc" block of the "qoctra" block in the "%s" configuration file.' % (error,misc['config_name']))

      # Maximum fluence (in J/m²)
      # ~~~~~~~~~~~~~~~

      if max_fluence:

        # Ionization potential
        # --------------------

        ip = -1

        for line in ip_list:
          if line['Molecule'] == misc['source_name']:
            if line['IP (adiabatic)'] != "N/A":
              ip = float(line['IP (adiabatic)'])
            elif line['IP (vertical)'] != "N/A":
              ip = float(line['IP (vertical)'])
            else:
              ip = float(line['IP (Koopmans)'])
            break

        if ip == -1:
          raise control_common.ControlError ("ERROR: Unable to find the ionization potential of this molecule in the %s file." % ip_file)

        # Convert the IP from Ha to Joules and divide by 1000 to get the maximum pulse energy (must remain a perturbation to the electrons in the system)

        energy = control_common.energy_unit_conversion(ip,'ha','j') / 1000

        # Affected area of the molecule
        # -----------------------------

        # Initialize some variables

        coord_list = []
        section_found = False

        # Define the expression patterns for the lines containing the atomic coordinates of the molecule in the QCHEM output file
      
        coord_rx = {

          # Pattern for finding the "Standard Nuclear Orientation (Angstroms)" line (which marks the start of the section)
          'start': re.compile(r'^\s*Standard Nuclear Orientation \(Angstroms\)\s*$'),

          # Pattern for finding lines looking like '   16      Si     -2.7647071137    -0.0043786180     2.7647071137'
          'coordinates': re.compile(r'^\s*\d+\s+[a-zA-Z]{1,3}\s+(?P<coord_x>-?\d+\.\d+)\s+(?P<coord_y>-?\d+\.\d+)\s+(?P<coord_z>-?\d+\.\d+)\s*$'),

          # Pattern for finding the "Total QAlloc Memory Limit" line (which marks the end of the section, no need to get further)
          'end': re.compile(r'^\s*Total QAlloc Memory Limit')

        }

        # Parse the qchem output file to get the information

        for line in misc['source_content']:

          # Define when the section begins and ends (ignore beta orbitals)

          if not section_found:
            if coord_rx['start'].match(line):
              section_found = True
        
          elif section_found and coord_rx['end'].match(line):
            break

          # If the line matches our coordinates pattern, extract the values and store them (after converting them from Angstroms to meters)

          elif coord_rx['coordinates'].match(line):

            x1 = float(coord_rx['coordinates'].match(line).group('coord_x'))*1e-10
            y1 = float(coord_rx['coordinates'].match(line).group('coord_y'))*1e-10
            z1 = float(coord_rx['coordinates'].match(line).group('coord_z'))*1e-10

            coord_list.append([x1,y1,z1])

        # Raise an exception if the section has not been found

        if not section_found:
          raise control_common.ControlError ("ERROR: Unable to find the 'Standard Nuclear Orientation (Angstroms)' section in the QCHEM output file")

        # Raise an exception if the atomic coordinates have not been found

        if coord_list == []:
          raise control_common.ControlError ("ERROR: Unable to find the atomic coordinates of the molecule in the QCHEM output file")

        # Compute the convex hull of the molecule and get the list of points consituting it

        points = np.array(coord_list)
        hull = ConvexHull(points)
        pts_hull = points[hull.vertices,:]

        # Compute the average distance between the points of the hull and their centroid ("radius")

        centroid = np.array(np.mean(pts_hull,axis=0),ndmin=2)
        dist = distance.cdist(pts_hull,centroid)
        radius = np.mean(dist)

        # Compute the area of the molecule affected by the laser (consider the shape of a sphere, and compute the surface of the hemisphere)

        area = 2 * np.pi *(radius**2)

        # Maximum fluence of the molecule
        # -------------------------------

        fluence = energy / area

        # Maximum fluence of the shaper
        # -----------------------------

        # Convert and compute the values using the shaper parameters to compute the shaper max fluence

        shaper_energy = float(ref_shaper['input_beam']['energy']) * 1e-6
        shaper_area = np.pi * ( (float(ref_shaper['input_beam']['diameter']) * 1e-3) ** 2 )

        shaper_fluence = shaper_energy / shaper_area

        # Store the lowest fluence
        # ------------------------

        if shaper_fluence < fluence:
          fluence = shaper_fluence

        param_render_vars.update({
            # OPC
            "fluence" : "{:.5e}".format(fluence).replace('e','d')
          })

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
      "start_pulse" : "../Pulse/Pulse_best",
      # POST CONTROL
      "mat_et0_path" : data['eigenvectors_path']
    })

    # Check if a "post_control" block has been defined in the "qoctra" block of the config file

    if not config['qoctra'].get('post_control'):
      raise control_common.ControlError ('ERROR: There is no "post_control" key in the "qoctra" block of the "%s" configuration file.' % misc['config_name'])    

    # Variables associated with the "post_control" block of the "qoctra" block in the config file

    try:
      param_render_vars.update({
        # POST CONTROL
        "analy_freq" : config['qoctra']['post_control']['analy_freq']
      })

    except KeyError as error:
      raise control_common.ControlError ('ERROR: The "%s" key is missing in the "post_control" block of the "qoctra" block in the "%s" configuration file.' % (error,misc['config_name']))

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
        raise control_common.ControlError ('ERROR: There is no "guess_pulse" key in the "qoctra" block of the "%s" configuration file.' % misc['config_name'])    

      # Variables associated with the "guess_pulse" block of the "qoctra" block in the config file

      try:

        amplitude = config['qoctra']['guess_pulse'].get('amplitude')

        pulse_render_vars = {
          "pulse_type" : config['qoctra']['guess_pulse']['pulse_type'],
          "subpulse_type" : config['qoctra']['guess_pulse']['subpulse_type'],
          "phase_change" : config['qoctra']['guess_pulse']['phase_change'],
          "sign" : config['qoctra']['guess_pulse']['sign']       
        }

      except KeyError as error:
        raise control_common.ControlError ('ERROR: The "%s" key is missing in the "guess_pulse" block of the "qoctra" block in the "%s" configuration file.' % (error,misc['config_name']))

      # Determine the subpulses constituting the guess pulse
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      # Convert the bandwidth from cm-1 to atomic units

      bandwidth_au = control_common.energy_unit_conversion(bandwidth,'cm-1','ha')
      full_bandwidth = ( bandwidth_au * 6 ) / 2.35 # Conversion from FWHM to full width for a gaussian

      # Determine the frequency range of subpulses (delimited by the spectral bandwidth, centered around the transition energy)

      upper_limit = omegazero + full_bandwidth/2
      lower_limit = omegazero - full_bandwidth/2

      # Initialize the variables 

      subpulses = []
      max_energy_diff = 0

      # Consider each pair of excited states

      for state_1 in range(len(system['states_list'])):
        for state_2 in range(state_1 + 1, len(system['states_list'])): # Starting at "state_1 + 1" to exclude the cases where both states are the same

          # Compute the energy difference and compare it to our frequency range

          energy_diff = abs(system['states_list'][state_1]['energy'] - system['states_list'][state_2]['energy'])

          if lower_limit <= energy_diff <= upper_limit:

            # Add it to the subpulses list

            subpulses.append(str(state_1 + 1) + " \t " + str(state_2 + 1)) # +1 because Fortran starts numbering at 1 while Python starts at 0.

            # Store the maximum transition energy to apply Nyquist–Shannon sampling theorem

            if energy_diff > max_energy_diff:
              max_energy_diff = energy_diff

      # Check the duration and time step
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      # Nyquist–Shannon sampling theorem: check if the sampling rate is bigger than the double of the highest frequency

      sampling_freq = 1 / (time_step * constants.value('atomic unit of time'))
      max_freq = control_common.energy_unit_conversion(max_energy_diff,'ha','hz')
      
      if not sampling_freq > 2*max_freq:
        raise control_common.ControlError ('ERROR: The time step value (%s a.u.) is too big for this transition. The sampling rate (%e Hz) is not bigger than the double of the highest frequency of the signal (%e Hz).' % (time_step,sampling_freq, max_freq))

      # Check if the total duration of the pulse is not too long compared to the lifetime of our most populated target state

      time_limit = 0.01 * ( system['states_list'][target_number]['lifetime'] / constants.value('atomic unit of time') )

      if duration > time_limit:
        raise control_common.ControlError ('ERROR: The duration of the pulse (%s a.u.) is too long compared to the radiative lifetime of this excited state (%s a.u.).' % (duration,system['states_list'][target_number]['lifetime']))

      #! Correct the values rather than raise an exception? Or calculate the values based on nstep (user-supplied) and time_limit (duration = time_limit and time_step = time_limit / nstep, then check against NYQ-SHA and decrease time_step and duration if necessary.)

      # Define the amplitude of the subpulses
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      if amplitude:

        amplitude_for = amplitude
        amplitude_au = float(re.compile(r'(\d*\.\d*)[dD]([-+]?\d+)').sub(r'\1E\2', amplitude_for)) # Replace the possible d/D from Fortran double precision float format with an "E", understandable by Python
        amplitude = amplitude_au * constants.value('atomic unit of electric field')

      elif not amplitude and max_fluence:
        
        # Compute the amplitude using the maximum fluence as reference

        intensity = fluence / duration_s
        amplitude = math.sqrt( intensity / (0.5 * constants.value('speed of light in vacuum') * constants.value('vacuum electric permittivity')) ) / 100
        amplitude_au = amplitude / constants.value('atomic unit of electric field')
        amplitude_for = "{:.5e}".format(amplitude_au).replace('e','d') # Replace the 'e' from Python with the 'd' from Fortran double precision      

      elif not amplitude and not max_fluence:
        raise control_common.ControlError ('ERROR: The "amplitude" key is missing in the "guess_pulse" block of the "qoctra" block in the "%s" configuration file while the maximum fluence has not been set to True.' % misc['config_name'])

      # Set the variables not associated with the config file
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      pulse_render_vars.update({
        "basis" : data['energies_path'],
        "nb_subpulses" : len(subpulses),
        "subpulses" : subpulses,
        "amplitude" : amplitude_for
      })

    # Rendering the file
    # ==================
     
    rendered_content[rendered_pulse] = jinja_render(misc['templates_dir'], template_pulse, pulse_render_vars)

    print('%12s' % "[ DONE ]")

    # ========================================================= #
    #                  Rendering the job script                 #
    # ========================================================= #

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
      raise control_common.ControlError ('ERROR: The "%s" key is missing in the "general" block of the "%s" configuration file.' % (error,misc['config_name']))

    # Variables associated with the clusters configuration file

    try:
      script_render_vars.update({
        "set_env" : clusters_cfg[job_specs['cluster_name']]['profiles'][job_specs['profile']]['set_env']     
      })

    except KeyError as error:
      raise control_common.ControlError ('ERROR: The "%s" key is missing in the "%s" profile of the clusters configuration file.' % (error,job_specs['profile']))

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
        raise control_common.ControlError ('ERROR: The "%s" key is missing in the CHAINS configuration file (chains_config.yml).' % error)

    # Rendering the file
    # ==================

    rendered_content[rendered_script] = jinja_render(misc['templates_dir'], template_script, script_render_vars)

    print('%12s' % "[ DONE ]")

    # ========================================================= #
    #                  Show the computed values                 #
    # ========================================================= #

    print("")
    print(''.center(60, '-'))
    print('Shaper values'.center(60, ' '))
    print(''.center(60, '-'))
    print("{:<30} {:<30}".format("Label: ", ref_shaper['label']))
    print("{:<30} {:<30}".format("Input beam energy (µJ): ", ref_shaper['input_beam']['energy']))
    print("{:<30} {:<30}".format("Input beam diameter (mm): ", ref_shaper['input_beam']['diameter']))
    if max_fluence:
      print("{:<30} {:<30}".format("Shaper max fluence (J/m^2): ", "{:.4e}".format(shaper_fluence)))
    print("{:<30} {:<30}".format("Reference frequency (nm): ", ref_frequency))
    print("{:<30} {:<30}".format("Reference duration (ps): ", ref_duration))
    print(''.center(60, '-'))
    print("")
    print(''.center(60, '-'))
    print('Pulse characteristics'.center(60, ' '))
    print(''.center(60, '-'))
    print("{:<30} {:<30}".format("Duration (a.u.): ", "{:.4e}".format(duration)))
    print("{:<30} {:<30}".format("Duration (s): ", "{:.4e}".format(duration_s)))
    print("{:<30} {:<30}".format("Bandwidth FWHM (a.u.): ", "{:.4e}".format(bandwidth_au)))
    print("{:<30} {:<30}".format("Bandwidth FWHM (cm^-1): ", "{:.4e}".format(bandwidth)))    
    print("{:<30} {:<30}".format("Central frequency (a.u.): ", "{:.4e}".format(omegazero))) 
    print("{:<30} {:<30}".format("Central frequency (cm^-1): ", "{:.4e}".format(omegazero_cm)))
    if max_fluence:
      print("{:<30} {:<30}".format("Maximum fluence (J/m^2): ", "{:.4e}".format(fluence)))
    print(''.center(60, '-'))
    print("")
    print(''.center(60, '-'))
    print('Guess pulse characteristics'.center(60, ' '))
    print(''.center(60, '-'))
    print("{:<30} {:<30}".format("Nb subpulses: ", len(subpulses)))
    print("{:<30} {:<30}".format("Amplitude (V/m): ", "{:.4e}".format(amplitude)))
    print("{:<30} {:<30}".format("Amplitude (a.u.): ", "{:.4e}".format(amplitude_au)))
    print(''.center(60, '-'))

    return rendered_content, rendered_script


######################################################################################################################################


def basic_opc_pcp_render(clusters_cfg:dict, config:dict, system:dict, data:dict, job_specs:dict, misc:dict):
    """Renders the job script and the two parameters files (OPC & PCP) associated with the QOCT-RA program in CHAINS.

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

    # Define the names of the templates

    template_param = "param.nml.jinja"
    template_script = "basic_opc_pcp_job.sh.jinja"
    template_pulse = "guess_pulse_OPC.jinja"

    # Check if the specified templates exist in the "templates" directory of CONTROL LAUNCHER.
    
    control_common.check_abspath(os.path.join(misc['templates_dir'],template_param),"Jinja template for the qoctra parameters files","file")
    control_common.check_abspath(os.path.join(misc['templates_dir'],template_script),"Jinja template for the qoctra job script","file")
    control_common.check_abspath(os.path.join(misc['templates_dir'],template_pulse),"Jinja template for the OPC guess pulse","file")

    # Define the names of the rendered files.

    rendered_script = "basic_opc_pcp_job.sh"
    rendered_pulse = "guess_pulse"
    rendered_param = "param_" + misc['transition']['label'] + ".nml"
    rendered_param_pcp = "param_" + misc['transition']['label'] + "_PCP.nml"

    # Initialize the dictionary that will be returned by the function

    rendered_content = {}

    # ===================================================== #
    #           Rendering the OPC parameters file           #
    # ===================================================== #

    print("{:<80}".format("\nRendering the jinja template for the OPC parameters file ...  "), end="")

    # Defining the Jinja variables
    # ============================

    param_render_vars = {
      # GENERAL
      "source_name" : misc['source_name'],
      "process" : "OPC",
      # DATA FILES
      "energies_file_path" : data['energies_path'],
      "momdip_path" : data['momdip_mtx_path'],
      "init_file_path" : data['init_path'],
      "target_file_path" : data['target_path'],
      # CONTROL
      "max_iter" : config['control']['max_iter'],
      "threshold" : config['control']['threshold'],
      "time_step" : config['control']['time_step'],
      "start_pulse" : " ",
      "guess_pulse" : rendered_pulse,
      # OPC
      "nb_steps" : config['opc']['nb_steps'],
      "alpha" : config['opc']['alpha'],
      "write_freq" : config['opc']['write_freq'],
      #! Temporary
      # POST CONTROL
      "mat_et0_path" : data['eigenvectors_path']
    }

    # Rendering the file
    # ==================
     
    rendered_content[rendered_param] = jinja_render(misc['templates_dir'], template_param, param_render_vars)

    print('%12s' % "[ DONE ]")

    # ========================================================= #
    #             Rendering the PCP parameters file             #
    # ========================================================= #

    print("{:<80}".format("\nRendering the jinja template for the PCP parameters file ...  "), end="")

    # Defining the Jinja variables (via updating the dictionary from the OPC parameters file)
    # ============================

    param_render_vars.update({
      # GENERAL
      "process" : "PCP",
      # CONTROL
      "start_pulse" : "../Pulse/Pulse",
      # POST CONTROL
      "analy_freq" : config['post_control']['analy_freq'],
      "mat_et0_path" : data['eigenvectors_path']
    })

    # Rendering the file
    # ==================
     
    rendered_content[rendered_param_pcp] = jinja_render(misc['templates_dir'], template_param, param_render_vars)

    print('%12s' % "[ DONE ]")

    # ========================================================= #
    #              Rendering the guess pulse file               #
    # ========================================================= #

    print("{:<80}".format("\nRendering the jinja template for the guess pulse file ...  "), end="")

    # Defining the generic Jinja variables
    # ====================================

    pulse_render_vars = {
      "basis" : data['energies_path'],
      "pulse_type" : config['guess_pulse']['pulse_type'],
      "subpulse_type" : config['guess_pulse']['subpulse_type'],
      "init_strength" : config['guess_pulse']['amplitude'],
      "phase_change" : config['guess_pulse']['phase_change'],
      "sign" : config['guess_pulse']['sign']       
    }

    # Determine the subpulses constituting the guess pulse
    # ====================================================

    # Initialize the variable

    subpulses = []

    # Add all transitions from the ground state to each of the excited states

    for state in range(1, len(system['states_list'])): # Starting at 1 to exclude the ground state 

      subpulses.append("1 \t " + str(state + 1)) # +1 because Fortran starts numbering at 1 while Python starts at 0.

    # Set the other variables associated with the subpulses
    # =====================================================

    pulse_render_vars.update({
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

    print("{:<80}".format("\nRendering the jinja template for the qoctra job script ..."), end="")

    # Defining the Jinja variables
    # ============================

    script_render_vars = {
      "source_name" : misc['source_name'],
      "transition" : misc['transition']['label'],
      "config_name" : misc['config_name'],
      "user_email" : config['user_email'],
      "mail_type" : config['mail_type'],
      "job_walltime" : job_specs['walltime'],
      "job_memory" : job_specs['memory'], # in MB
      "partition" : job_specs['partition'],
      "rendered_param" : rendered_param,
      "rendered_param_PCP" : rendered_param_pcp,
      "guess_pulse" : rendered_pulse,
      "set_env" : clusters_cfg[job_specs['cluster_name']]['profiles'][job_specs['profile']]['set_env']
    }

    # Rendering the file
    # ==================

    rendered_content[rendered_script] = jinja_render(misc['templates_dir'], template_script, script_render_vars)

    print('%12s' % "[ DONE ]")

    return rendered_content, rendered_script