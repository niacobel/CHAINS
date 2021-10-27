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

    # Check if a "ref_duration" block has been defined in the config file

    if not config.get('ref_duration'):
      raise control_common.ControlError ('ERROR: There is no "ref_duration" key defined in the "%s" configuration file.' % misc['config_name'])

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

    # Pulse duration
    # ~~~~~~~~~~~~~~
    
    # Compute the transition energy that will act as the central frequency

    init_number = misc['transition']['init_states'].index(max(misc['transition']['init_states']))
    target_number = misc['transition']['target_states'].index(max(misc['transition']['target_states']))
    omegazero = abs(system['eigenstates_list'][init_number]['energy'] - system['eigenstates_list'][target_number]['energy'])

    # Convert the central frequency to nanometers

    omegazero_nm = control_common.energy_unit_conversion(omegazero,'ha','nm')

    # Find the closest reference to our central frequency (see https://stackoverflow.com/questions/12141150/from-list-of-integers-get-number-closest-to-a-given-value)

    closest_ref = min(config['ref_duration'].keys(), key=lambda x : abs(x - omegazero_nm))

    # Fix the duration of the pulse using that reference

    duration_s = config['ref_duration'][closest_ref] * 1e-12
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

        bandwidth = control_common.energy_unit_conversion(config['qoctra']['opc']['bandwidth'],'cm-1','ha')

        amplitude_filter = config['qoctra']['opc'].get('amplitude_filter', False)
        bandwidth_filter = config['qoctra']['opc'].get('bandwidth_filter', False)

        param_render_vars.update({
          # OPC
          "nb_steps" : nb_steps,
          "alpha" : config['qoctra']['opc']['alpha'],
          "write_freq" : config['qoctra']['opc']['write_freq']
        })

        if bandwidth_filter:
          param_render_vars.update({
            # OPC
            "omegazero" : omegazero,
            "bandwidth" : bandwidth
          })

      except KeyError as error:
        raise control_common.ControlError ('ERROR: The "%s" key is missing in the "opc" block of the "qoctra" block in the "%s" configuration file.' % (error,misc['config_name']))

      # Amplitude filtering
      # ~~~~~~~~~~~~~~~~~~~

      if amplitude_filter:

        # Affected area of the molecule (necessary to compute the maximum field strength)
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

        # Determine the maximum field strength
        # ------------------------------------

        # Get the ionization potential

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

        # Convert the IP from Ha to Joules and divide by 1000 to get the maximul pulse energy (must remain a perturbation to the electrons in the system)

        energy = control_common.energy_unit_conversion(ip,'ha','j') / 1000

        # Compute the maximum field strength then convert it to atomic units

        max_strength = math.sqrt( (2 * energy) / (constants.value('speed of light in vacuum') * constants.value('vacuum electric permittivity') * area * duration_s) )
        max_strength_au = max_strength / constants.value('atomic unit of electric field')

        # Store the value

        param_render_vars.update({
            # OPC
            "max_strength" : "{:.5e}".format(max_strength_au).replace('e','d')
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
      "start_pulse" : "../Pulse/Pulse",
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

        if amplitude_filter:
          init_strength = max_strength_au / 100
          init_strength_for = "{:.5e}".format(init_strength).replace('e','d') # Replace the 'e' from Python with the 'd' from Fortran double precision
        else:
          init_strength_for = config['qoctra']['guess_pulse']['amplitude']
          init_strength = float(re.compile(r'(\d*\.\d*)[dD]([-+]?\d+)').sub(r'\1E\2', init_strength_for)) # Replace the possible d/D from Fortran double precision float format with an "E", understandable by Python

        pulse_render_vars = {
          "pulse_type" : config['qoctra']['guess_pulse']['pulse_type'],
          "subpulse_type" : config['qoctra']['guess_pulse']['subpulse_type'],
          "init_strength" : init_strength_for,
          "phase_change" : config['qoctra']['guess_pulse']['phase_change'],
          "sign" : config['qoctra']['guess_pulse']['sign']       
        }

      except KeyError as error:
        raise control_common.ControlError ('ERROR: The "%s" key is missing in the "guess_pulse" block of the "qoctra" block in the "%s" configuration file.' % (error,misc['config_name']))

      # Determine the subpulses constituting the guess pulse
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      # Determine the frequency range of subpulses (delimited by the spectral bandwidth, centered around the transition energy)

      upper_limit = omegazero + bandwidth/2
      lower_limit = omegazero - bandwidth/2

      # Initialize the variables 

      subpulses = []
      max_energy_diff = 0

      # Consider each pair of excited states

      for state_1 in range(1, len(system['eigenstates_list'])): # Starting at 1 to exclude the ground state
        for state_2 in range(state_1 + 1, len(system['eigenstates_list'])): # Starting at "state_1 + 1" to exclude the cases where both states are the same

          # Compute the energy difference and compare it to our frequency range

          energy_diff = abs(system['eigenstates_list'][state_1]['energy'] - system['eigenstates_list'][state_2]['energy'])

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

      time_limit = 0.01 * system['eigenstates_list'][target_number]['lifetime']

      if duration > time_limit:
        raise control_common.ControlError ('ERROR: The duration of the pulse (%s a.u.) is too long compared to the radiative lifetime of this excited state (%s a.u.).' % (duration,system['eigenstates_list'][target_number]['lifetime']))

      #! Correct the values rather than raise an exception? Or calculate the values based on nstep (user-supplied) and time_limit (duration = time_limit and time_step = time_limit / nstep, then check against NYQ-SHA and decrease time_step and duration if necessary.)

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
    print('Pulse values'.center(60, ' '))
    print(''.center(60, '-'))

    print("{:<30} {:<30}".format("Duration (a.u.): ", "{:.4e}".format(duration)))
    print("{:<30} {:<30}".format("Duration (s): ", "{:.4e}".format(duration_s)))
    print("{:<30} {:<30}".format("Initial strength (V/m): ", "{:.4e}".format(init_strength * constants.value('atomic unit of electric field'))))
    print("{:<30} {:<30}".format("Initial strength (a.u.): ", "{:.4e}".format(init_strength)))
    print("{:<30} {:<30}".format("Bandwidth (a.u.): ", "{:.4e}".format(bandwidth)))
    print("{:<30} {:<30}".format("Bandwidth (cm^-1): ", "{:.4e}".format(control_common.energy_unit_conversion(bandwidth,'ha','cm-1'))))    
    print("{:<30} {:<30}".format("Central frequency (a.u.): ", "{:.4e}".format(omegazero))) 
    print("{:<30} {:<30}".format("Central frequency (cm^-1): ", "{:.4e}".format(control_common.energy_unit_conversion(omegazero,'ha','cm-1')))) 
    print("{:<30} {:<30}".format("Nb subpulses: ", len(subpulses)))

    if amplitude_filter:
      print("{:<30} {:<30}".format("Max energy (J): ", "{:.4e}".format(energy)))
      print("{:<30} {:<30}".format("Molecular surface (m^2): ", "{:.4e}".format(area)))
      print("{:<30} {:<30}".format("Maximum strentgh (V/m): ", "{:.4e}".format(max_strength)))
      print("{:<30} {:<30}".format("Maximum strentgh (a.u.): ", "{:.4e}".format(max_strength_au)))

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
      "energies_file_path" : data['eigenvalues_path'],
      "momdip_e_path" : data['momdip_es_mtx_path'],
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
      "basis" : data['eigenvalues_path'],
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

    # Add all transitions from the ground state to each of the excited eigenstates

    for state in range(1, len(system['eigenstates_list'])): # Starting at 1 to exclude the ground state 

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