################################################################################################################################################
##                                                             Modelling functions                                                            ##
##                                                                                                                                            ##
##                                      This script contains the modelling functions for CONTROL LAUNCHER,                                    ##
##                                consult the documentation at https://chains-ulb.readthedocs.io/ for details                                 ##
################################################################################################################################################

import math
import re

import numpy as np
from scipy import constants, linalg

import control_common

# =================================================================== #
# =================================================================== #
#                    ORCA TD-DFT Modelling Function                   #
# =================================================================== #
# =================================================================== #

def orca_tddft(source_content:list):

    # Initialize the system dictionary that will be returned by the function

    system = {}

    # Define an array for correct English spelling during printing

    special_numbers = {1:"st", 2:"nd", 3:"rd"}

    # =================================================================== #
    # =================================================================== #
    #                  Non-Relativistic States Treatment                  #
    # =================================================================== #
    # =================================================================== #

    section_title = "Non-Relativistic States Treatment"

    print("")
    print(''.center(len(section_title)+30, '='))
    print(section_title.center(len(section_title)+30))
    print(''.center(len(section_title)+30, '='))

    # ========================================================= #
    #                      Number of roots                      #
    # ========================================================= #

    nb_roots = False
    section_found = False

    # Define the expression patterns for the lines containing information about the number of roots, that number will be used to check if all the needed values have been collected.

    nb_roots_rx = {

      # Pattern for finding the "INPUT FILE" line (which marks the start of the section)
      'start': re.compile(r'^\s*INPUT\s+FILE\s*$'),
      
      # Pattern for finding lines looking like "|  8> nroots 5"
      'value': re.compile(r'^\|\s*\d+>\s+nroots\s+(?P<n_roots>\d+)$'),

      # Pattern for finding the "****END OF INPUT****" line (which marks the end of the section)
      'end': re.compile(r'^\|\s*\d+>\s+\*+END OF INPUT\*+\s*$')

    }
  
    # Parse the source file to get the information

    for line in source_content:

      # Define when the section begins and ends

      if not section_found:
        if nb_roots_rx['start'].match(line):
          section_found = True
    
      elif section_found and nb_roots_rx['end'].match(line):
        break
  
      # Get the normal number of roots

      elif nb_roots_rx['value'].match(line):
        if nb_roots:
          new_nb_roots = int(nb_roots_rx['value'].match(line).group('n_roots'))
          if new_nb_roots != nb_roots:
            raise control_common.ControlError ("ERROR: There is a discrepancy for the number of roots in the source file, the last value found (%s) does not match the previous value (%s)" % (new_nb_roots, nb_roots))
        else:
          nb_roots = int(nb_roots_rx['value'].match(line).group('n_roots'))

    # Raise an exception if the number of roots has not been found

    if not nb_roots:
      raise control_common.ControlError ("ERROR: Unable to find the number of roots in the source file")
    if not section_found:
      raise control_common.ControlError ("ERROR: Unable to find the 'INPUT FILE' section in the source file")

    print("{:<50} {:<10}".format("\nNumber of roots: ",nb_roots))

    # ========================================================= #
    #                Non-relativistic states List               #
    # ========================================================= #

    print("{:<50}".format("\nParsing the excited states ...  "), end="")

    # Initialization of the variables

    section_found = False
    cnt_state = 0
    cnt_triplet = 0
    cnt_singlet = 0
    state_number = 0

    # Define the ground state of our molecule (which is the first state and has a zero energy)
    
    system['zero_states_list'] = [{'number': 0, 'orca_number': 0, 'label': 'S0', 'energy' : 0.0}]

    # Define the expression patterns for the lines containing information about the states

    states_rx = {

      # Pattern for finding the "$$$$$$$$$$$$$$$$  JOB NUMBER  1 $$$$$$$$$$$$$$" line (which marks the start of the section)
      'start': re.compile(r'^\s*\$+\s*JOB\s+NUMBER\s+2\s+\$*\s*$'),
      
      # Pattern for finding lines looking like 'STATE  2:  E=   0.160763 au      4.375 eV    35283.3 cm**-1 <S**2> =   2.000000'
      'state_line': re.compile(r'^\s*STATE\s+(?P<number>\d+):\s+E=\s+(?P<energy>\d+\.\d+)\s+au\s+\d+\.\d+\s+eV\s+\d+\.\d+\s+cm\*\*-1\s+<S\*\*2>\s+=\s+(?P<spin>\d+\.\d+)\s*$'),

      # Pattern for finding the "$$$$$$$$$$$$$$$$  JOB NUMBER  2 $$$$$$$$$$$$$$" line (which marks the end of the section)
      'end': re.compile(r'^\s*TD-DFT\/TDA\s+SPIN-ORBIT\s+COUPLING\s*$')

    }

    # Parse the source file to get the information and build the states list

    for line in source_content:

      # Define when the section begins and ends

      if not section_found:
        if states_rx['start'].match(line):
          section_found = True
    
      elif section_found and states_rx['end'].match(line):
        break
        
      # Get the state number and multiplicity and its associated excitation energy

      elif states_rx['state_line'].match(line):

        cnt_state += 1
        init_number = int(states_rx['state_line'].match(line).group("number"))
        exc_energy = float(states_rx['state_line'].match(line).group("energy"))
        spin = round(float(states_rx['state_line'].match(line).group("spin")))

        if spin == 2:
          multiplicity = "Triplet"
        elif spin == 0:
          multiplicity = "Singlet"

        # Check the multiplicity and increase the counter (cnt) for that multiplicity (needed to check if we've found everything)

        if multiplicity == "Triplet":
          first_letter = "T"
          cnt_triplet += 1

        elif multiplicity == "Singlet":
          first_letter = "S"
          cnt_singlet += 1

        else:
          raise control_common.ControlError ("ERROR: Multiplicity of the %s%s state is of unknown value (%s)" % (init_number, ("th" if not init_number in special_numbers else special_numbers[init_number]),multiplicity))

        # Append information about the current state to the zero_states_list key

        if multiplicity == "Singlet":
          state_number += 1
          system['zero_states_list'].append({'number': state_number, 'orca_number': cnt_state, 'label': (first_letter + str(init_number)), 'energy': exc_energy})

        elif multiplicity == "Triplet":

          # Add the substate ms=0
          state_number += 1
          system['zero_states_list'].append({'number': state_number, 'orca_number': cnt_state, 'label': (first_letter + str(init_number) + "(ms=0)"), 'energy': exc_energy})
          # Add the substate ms=1
          state_number += 1
          system['zero_states_list'].append({'number': state_number, 'orca_number': cnt_state, 'label': (first_letter + str(init_number) + "(ms=1)"), 'energy': exc_energy})
          # Add the substate ms=-1
          state_number += 1
          system['zero_states_list'].append({'number': state_number, 'orca_number': cnt_state, 'label': (first_letter + str(init_number) + "(ms=-1)"), 'energy': exc_energy})
          
    # Raise an exception if the section has not been found

    if not section_found:
      raise control_common.ControlError ("ERROR: Unable to find the 'JOB NUMBER  1' section in the source file")

    # Raise an exception if not all the values have been found

    if cnt_state != 2*nb_roots:
      raise control_common.ControlError ("ERROR: The modelling function could not find the right number of excited states in the source file (%s of the %s expected states have been found)" % (cnt_state,2*nb_roots))
    if cnt_triplet != nb_roots:
      raise control_common.ControlError ("ERROR: The modelling function could not find the right number of excited triplet states in the source file (%s of the %s expected triplet states have been found)" % (cnt_triplet,nb_roots))
    if cnt_singlet != nb_roots:
      raise control_common.ControlError ("ERROR: The modelling function could not find the right number of excited singlet states in the source file (%s of the %s expected singlet states have been found)" % (cnt_singlet,nb_roots))

    # Raise an exception if the state numbers are not consecutive and starting at 0

    control_common.is_consecutive(list(dict.fromkeys([state['orca_number'] for state in system['zero_states_list']])),"Excited state numbers from the source file")

    print("[ DONE ]")
  
    # ========================================================= #
    #                    Dipole Moments List                    #
    # ========================================================= #

    print("{:<50}".format("\nParsing the transition dipole moments ...  "), end="")

    # Initialization of the variables

    momdip_section_found = False
    abs_section_found = False
    trans_section_found = False
    sub_trans_section_found = False
    momdip_list = []

    # Define the start and end expression patterns for the lines containing information about the dipole moments

    moment_rx = {

      # Pattern for finding the "$$$$$$$$$$$$$$$$  JOB NUMBER  1 $$$$$$$$$$$$$$" line (which marks the start of the section)
      'start': re.compile(r'^\s*\$+\s*JOB\s+NUMBER\s+1\s+\$*\s*$'),

      # Pattern for finding the "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS" line (which marks the start of the absorption spectrum section)
      'abs_start': re.compile(r'^\s*ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS\s*$'),

      # Pattern for finding lines looking like '   1   31006.3    322.5   0.022098171   0.23463   0.00009  -0.00010   0.48439'
      'moment': re.compile(r'^\s*(?P<state_2>\d+)\s+\d+\.\d+\s+\d+\.\d+\s+\d+\.\d+\s+\d+\.\d+\s+(?P<mom_x>-?\d+\.\d+)\s+(?P<mom_y>-?\d+\.\d+)\s+(?P<mom_z>-?\d+\.\d+)\s*$'),

      # Pattern for finding lines looking like '   6   29671.3    337.0   spin forbidden '
      'forb_moment': re.compile(r'^\s*(?P<state_2>\d+)\s+\d+\.\d+\s+\d+\.\d+\s+spin forbidden\s+'),
      
      # Pattern for finding the "ABSORPTION SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS" line (which marks the end of the absorption spectrum section)
      'abs_end': re.compile(r'^\s*ABSORPTION SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS\s*$'),
      
      # Pattern for finding the "TRANSIENT TD-DFT/TDA-EXCITATION SPECTRA" line (which marks the start of the transient spectra section)
      'trans_start': re.compile(r'^^\s*TRANSIENT TD-DFT/TDA-EXCITATION SPECTRA\s*$'),

      # Pattern for finding lines looking like 'Transitions starting from IROOT 1:'
      'trans_state_1': re.compile(r'^\s*Transitions starting from IROOT (?P<state_1>\d+):\s*$'),

      # Pattern for finding the "TRANSIENT ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS" line (which marks the start of the transient absorption spectrum section for a specific IROOT)
      'sub_trans_start': re.compile(r'^\s*TRANSIENT ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS\s*$'),

      # Pattern for finding the "TRANSIENT ABSORPTION SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS" line (which marks the end of the athe transient absorption spectrum section for a specific IROOT)
      'sub_trans_end': re.compile(r'^\s*TRANSIENT ABSORPTION SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS\s*$'),

      # Pattern for finding the "*** ORCA-CIS/TD-DFT FINISHED WITHOUT ERROR ***" line (which marks the end of the transient spectra section)
      'trans_end': re.compile(r'^\s*\**\s*ORCA-CIS/TD-DFT FINISHED WITHOUT ERROR\s*\**\s*$'),
      
      # Pattern for finding the "$$$$$$$$$$$$$$$$  JOB NUMBER  2 $$$$$$$$$$$$$$" line (which marks the end of the section)
      'end': re.compile(r'^\s*\$+\s*JOB\s+NUMBER\s+2\s+\$*\s*$')

    }
  
    # Parse the source file to get the information and build the dipole moments list

    for line in source_content:

      # Define when the global section begins and ends

      if not momdip_section_found:
        if moment_rx['start'].match(line):
          momdip_section_found = True
    
      elif momdip_section_found and moment_rx['end'].match(line):
        break

      else: 

        # Extract the relevant information from the absorption spectrum and add it to the momdip_list
        
        if not abs_section_found:
          if moment_rx['abs_start'].match(line):
            abs_section_found = True
  
        elif abs_section_found and moment_rx['abs_end'].match(line):
          abs_section_found = False
  
        elif abs_section_found:
          
          if moment_rx['moment'].match(line):
          
            matching_line = moment_rx['moment'].match(line)
    
            state_1 = 0
            state_2 = int(matching_line.group('state_2'))
            value_x = float(matching_line.group('mom_x'))
            value_y = float(matching_line.group('mom_y'))
            value_z = float(matching_line.group('mom_z'))
            momdip = (state_1, state_2, value_x, value_y, value_z)
          
            # Add the new line to the momdip_list
            momdip_list.append(momdip)
  
          elif moment_rx['forb_moment'].match(line):
          
            matching_line = moment_rx['forb_moment'].match(line)
    
            state_1 = 0
            state_2 = int(matching_line.group('state_2'))
            value_x = float(0)
            value_y = float(0)
            value_z = float(0)
            momdip = (state_1, state_2, value_x, value_y, value_z)
          
            # Add the new line to the momdip_list
            momdip_list.append(momdip)
 
        # Extract the relevant information from the transient spectra and add it to the momdip_list
  
        if not trans_section_found:
          if moment_rx['trans_start'].match(line):
            trans_section_found = True
  
        elif trans_section_found and moment_rx['trans_end'].match(line):
          trans_section_found = False
  
        elif trans_section_found:

          if moment_rx['trans_state_1'].match(line):
            state_1 = int(moment_rx['trans_state_1'].match(line).group('state_1'))

          elif not sub_trans_section_found:
            if moment_rx['sub_trans_start'].match(line):
              sub_trans_section_found = True
    
          elif sub_trans_section_found and moment_rx['sub_trans_end'].match(line):
            sub_trans_section_found = False

          elif sub_trans_section_found and moment_rx['moment'].match(line):
            
            matching_line = moment_rx['moment'].match(line)

            state_2 = int(matching_line.group('state_2'))
            value_x = float(matching_line.group('mom_x'))
            value_y = float(matching_line.group('mom_y'))
            value_z = float(matching_line.group('mom_z'))
            momdip = (state_1, state_2, value_x, value_y, value_z)
          
            # Add the new line to the momdip_list
            momdip_list.append(momdip)
        
    # Raise an exception if the section has not been found

    if not section_found:
      raise control_common.ControlError ("ERROR: Unable to find the 'STATE-TO-STATE TRANSITION MOMENTS' section in the source file")

    # Raise an exception if not all the values have been found

    nb_momdip = nb_roots + nb_roots + (nb_roots*(nb_roots-1)/2) # Ground to Singlet + Ground to Triplet + Singlet to Singlet
    if len(momdip_list) != nb_momdip:
      raise control_common.ControlError ("ERROR: The parsing function could not find the right number of state-to-state transition moments in the source file (%s of the %s expected values have been found)" % (len(momdip_list),nb_momdip))

    print("[ DONE ]")

    # ========================================================= #
    #                   Dipole Moments Matrices                 #
    # ========================================================= #

    print("{:<50}".format("\nBuilding transition dipole moments matrices ... "), end="")

    # Intialize the momdip_o_mtx dictionary

    system['momdip_o_mtx'] = {}

    # Initialize the matrices as zero-filled matrices

    system['momdip_o_mtx']['X'] = np.zeros((len(system['zero_states_list']), len(system['zero_states_list'])), dtype=float)
    system['momdip_o_mtx']['Y'] = np.zeros((len(system['zero_states_list']), len(system['zero_states_list'])), dtype=float)
    system['momdip_o_mtx']['Z'] = np.zeros((len(system['zero_states_list']), len(system['zero_states_list'])), dtype=float)

    # Filling the matrices

    for momdip in momdip_list:

      states_1 = [state['number'] for state in system['zero_states_list'] if state['orca_number'] == momdip[0]] # For triplets, there are three states for each 'qchem_number'
      states_2 = [state['number'] for state in system['zero_states_list'] if state['orca_number'] == momdip[1]]

      for k1 in states_1:
        for k2 in states_2:
          
          system['momdip_o_mtx']['X'][k1][k2] = momdip[2]
          system['momdip_o_mtx']['X'][k2][k1] = momdip[2]    # For symetry purposes
    
          system['momdip_o_mtx']['Y'][k1][k2] = momdip[3]
          system['momdip_o_mtx']['Y'][k2][k1] = momdip[3]    # For symetry purposes
    
          system['momdip_o_mtx']['Z'][k1][k2] = momdip[4]
          system['momdip_o_mtx']['Z'][k2][k1] = momdip[4]    # For symetry purposes

    print("[ DONE ]")

    # ========================================================= #
    #                            MIME                           #
    # ========================================================= #

    print("{:<50}".format("\nParsing the full SOC matrix ...  "), end="")

    # Initialization of the variables

    section_found = False
    real_section_found = False
    im_section_found = False
    nb_values = 0
    mime_real = np.zeros((len(system['zero_states_list']), len(system['zero_states_list'])), dtype=float)
    mime_imag = np.zeros((len(system['zero_states_list']), len(system['zero_states_list'])), dtype=float)
    
    # Define the expression patterns for the lines containing information about the SOC
    
    matrix_rx = {

      # Pattern for finding the "The full SOC matrix:" line (which marks the start of the section)
      'start': re.compile(r'^\s*The\s+full\s+SOC\s+matrix:\s*$'),

      # Pattern for finding the 'Real part:' line
      'real_start': re.compile(r'^\s*Real\s+part:\s*$'),

      # Pattern for finding the 'Image part:' line
      'im_start': re.compile(r'^\s*Image\s+part:\s*$'),

      # Pattern for finding the '... done' line
      'im_end': re.compile(r'^\s*\.\.\.\s*done\s*$'),   

      # Pattern for finding lines looking like '0          1          2          3          4          5'
      'state_2_line': re.compile(r'^\s*(?:\d+\s+)*\d+\s*$'),

      # Pattern for finding lines looking like '3    2.946908e-05  4.460731e-05  -1.435899e-04  2.718156e-06  -4.520891e-06  5.879732e-05'
      'matrix_line': re.compile(r'^\s*\d+\s+(?:-?\d+\.\d+e[\+-]?\d+\s*)+$'),
      
      # Pattern for finding the "Diagonalizing the SOC matrix:   ... done" line (which marks the end of the section)
      'end': re.compile(r'^\s*Diagonalizing\s+the\s+SOC\s+matrix:\s+\.\.\.\s+done\s*$')

    }

    # Parse the source file to get the information and build the SOC list

    for line in source_content:

      # Define when the section begins and ends

      if not section_found:
        if matrix_rx['start'].match(line):
          section_found = True
    
      elif section_found and matrix_rx['end'].match(line):
        break

      else:
                
        # Fetch the real part of the SOC matrix
        
        if not real_section_found:
          if matrix_rx['real_start'].match(line):
            real_section_found = True
  
        elif real_section_found and matrix_rx['im_start'].match(line):
          real_section_found = False
  
        elif real_section_found:
          
          if matrix_rx['state_2_line'].match(line):
            states_2 = [int(nb) for nb in line.split(" ") if nb != '']

          elif matrix_rx['matrix_line'].match(line):
            line_list = [value for value in line.split(" ") if value != '']
            state_1 = int(line_list.pop(0))
            for idx, value in enumerate(line_list):
              state_2 = states_2[idx]
              mime_real[state_1][state_2] = value
              nb_values += 1
            
        # Fetch the imaginary part of the SOC matrix
  
        if not im_section_found:
          if matrix_rx['im_start'].match(line):
            im_section_found = True
  
        elif im_section_found and matrix_rx['im_end'].match(line):
          im_section_found = False
  
        elif im_section_found:

          if matrix_rx['state_2_line'].match(line):
            states_2 = [int(nb) for nb in line.split(" ") if nb != '']

          elif matrix_rx['matrix_line'].match(line):
            line_list = [value for value in line.split(" ") if value != '']
            state_1 = int(line_list.pop(0))
            for idx, value in enumerate(line_list):
              state_2 = states_2[idx]
              mime_imag[state_1][state_2] = value
              nb_values += 1
        
    # Raise an exception if the section has not been found

    if not section_found:
      raise control_common.ControlError ("ERROR: Unable to find the full SOC matrix in the source file")

    # Raise an exception if not all the values have been found

    nb_soc = ((nb_roots * 4) + 1) ** 2
    if nb_values != 2*nb_soc:
      raise control_common.ControlError ("ERROR: The modelling function could not find the right number of spin-orbit couplings in the source file (%s of the %s expected values have been found)" % (nb_values,nb_soc))

    # Form the total MIME by merging the two arrays

    system['mime'] = np.zeros((len(system['zero_states_list']), len(system['zero_states_list'])), dtype=complex)
    system['mime'].real = mime_real
    system['mime'].imag = mime_imag
  
    print("[ DONE ]")

    # ========================================================= #
    #                          SOC List                         #
    # ========================================================= #

    print("{:<50}".format("\nEstablishing the spin-orbit couplings list ..."), end="")

    # Initialization of the variables

    soc_list = []
    
    # Iterate over the MIME to get the information and build the SOC list

    it = np.nditer(system['mime'], flags=['multi_index'])
  
    for soc in it:

      state_1 = it.multi_index[0]
      label_1 = [state['label'] for state in system['zero_states_list'] if state_1 == state['number']][0]
      state_2 = it.multi_index[1]
      label_2 = [state['label'] for state in system['zero_states_list'] if state_2 == state['number']][0]
      value = complex(soc)

      # Skip diagonal values and singlet-singlet values
      if state_1 == state_2 or (label_1.startswith('S') and label_2.startswith('S')):
        continue

      # Add the information to the soc_list
      soc_line = (state_1, label_1, state_2, label_2, value)
      eq_soc_line = (state_2, label_2, state_1, label_1, value.conjugate())
      if eq_soc_line not in soc_list:
        soc_list.append(soc_line)
  
    print("[ DONE ]")
  
    """
    # ========================================================= #
    #           Radiative lifetime of excited states            #
    # ========================================================= #

    print("{:<50}".format("\nComputing radiative lifetimes ... "), end="")

    # This calculation is based on the A_mn Einstein Coefficients and their link with the transition dipole moment
    # See https://aapt.scitation.org/doi/pdf/10.1119/1.12937 for reference
    # Note that this calculation is performed using atomic units, which means the Planck constant equals 2*pi and the vacuum permittivity equals 1/(4*pi)

    # Constants

    light_speed_au = constants.value('speed of light in vacuum') / constants.value('atomic unit of velocity')

    # Calculate the radiative lifetime of each excited state

    for state in system['zero_states_list']:

      sum_einstein_coeffs = 0

      for other_state in system['zero_states_list']:

        if other_state['energy'] < state['energy']:

          # Compute and convert the energy difference

          energy_diff = state['energy'] - other_state['energy']

          # Compute the square of the transition dipole moment

          square_dipole = 0
          for momdip_key in system['momdip_o_mtx']:
            square_dipole += system['momdip_o_mtx'][momdip_key][state['number']][other_state['number']] ** 2

          # Calculate the A Einstein Coefficient   
         
          einstein_coeff = (4/3) * square_dipole * (energy_diff**3) / (light_speed_au**3)
          sum_einstein_coeffs += einstein_coeff

      if sum_einstein_coeffs == 0:
        state['lifetime'] = float('inf')
      else:
        state['lifetime'] = 1 / sum_einstein_coeffs
        state['lifetime'] = state['lifetime'] * constants.value('atomic unit of time')

    print("[ DONE ]")
    """
  
    # ========================================================= #
    #              Printing values in the log file              #
    # ========================================================= #

    # Print the non-relativistic states list
    # ======================================

    table_width = 73
    print("")
    print(''.center(table_width, '-'))
    print('Non-relativistic states list'.center(table_width, ' '))
    print(''.center(table_width, '-'))
    print("{:<10} {:<15} {:<15} {:<15} {:<15}".format('Number','Label','Energy (cm-1)','Energy (Ha)','Energy (nm)'))
    print(''.center(table_width, '-'))
    for state in system['zero_states_list']:
      print("{:<10} {:<15} {:<15.2f} {:<15.5f} {:<15.2f}".format(state['number'],state['label'],control_common.energy_unit_conversion(state['energy'],"ha","cm-1"),state['energy'],control_common.energy_unit_conversion(state['energy'],"ha","nm")))
    print(''.center(table_width, '-'))

    # Print the SOC list
    # ==================

    table_width = 73
    print("")
    print(''.center(table_width, '-'))
    print('Spin-orbit couplings'.center(table_width, ' '))
    print(''.center(table_width, '-'))
    print("{:<15} {:<15} {:<20} {:<20}".format('State 1','State 2','Real value (Ha)','Imag value (Ha)'))
    print(''.center(table_width, '-'))
    for soc in soc_list:
      column_1 = soc[1]
      column_2 = soc[3]
      print("{:<15} {:<15} {:<20.5g} {:<20.5g}".format(column_1,column_2,soc[4].real,soc[4].imag))
    print(''.center(table_width, '-'))    

    # Print the transition dipole moments list
    # ========================================

    table_width = 73
    print("")
    print(''.center(table_width, '-'))
    print('Transition dipole moments'.center(table_width, ' '))
    print(''.center(table_width, '-'))
    print("{:<10} {:<10} {:<12} {:<12} {:<12} {:<12}".format('State 1','State 2','X (a.u.)','Y (a.u.)','Z (a.u.)','Tot (a.u.)'))
    print(''.center(table_width, '-'))
    for momdip in momdip_list:
      column_1 = [state['label'] for state in system['zero_states_list'] if state['orca_number'] == momdip[0]][0].partition("(")[0]
      column_2 = [state['label'] for state in system['zero_states_list'] if state['orca_number'] == momdip[1]][0].partition("(")[0]
      column_6 = math.sqrt(momdip[2]**2 + momdip[3]**2 + momdip[4]**2)
      print("{:<10} {:<10} {:<12.4g} {:<12.4g} {:<12.4g} {:<12.4g}".format(column_1,column_2,momdip[2],momdip[3],momdip[4],column_6))
    print(''.center(table_width, '-'))

    # =================================================================== #
    # =================================================================== #
    #                    Relativistic States Treatment                    #
    # =================================================================== #
    # =================================================================== #

    section_title = "Relativistic States Treatment"

    print("\n\n")
    print(''.center(len(section_title)+30, '='))
    print(section_title.center(len(section_title)+30))
    print(''.center(len(section_title)+30, '='))

    # ========================================================= #
    #                    MIME diagonalization                   #
    # ========================================================= #

    print("{:<50}".format("\nDiagonalizing the MIME ..."), end="")

    # Diagonalization
    # ===============

    # Use SciPy to diagonalize the matrix (see https://personal.math.ubc.ca/~pwalls/math-python/linear-algebra/eigenvalues-eigenvectors/ for reference)
    # Documentation page for the function used here : https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.eigh.html  
     
    eigenvalues, system['eigenvectors'] = linalg.eigh(system['mime'])
    eigenvalues = eigenvalues.tolist()

    # Build the relativistic states list
    # ==================================

    # Initialize the variables

    system['states_list'] = []

    # Build the list

    for eigenvalue in eigenvalues:

      eigenstate = {
        'energy' : eigenvalue,
        'number' : eigenvalues.index(eigenvalue),
        'label' : "E" + str(eigenvalues.index(eigenvalue))
      }

      system['states_list'].append(eigenstate)

    # Transpose the eigenvectors list
    # ===============================

    # Using SciPy to invert the eigenvectors matrix
    # Note that the inverse of an orthonormal matrix is equal to its transpose, so each line of this matrix corresponds to an eigenvector. (see https://math.stackexchange.com/questions/156735/in-which-cases-is-the-inverse-of-a-matrix-equal-to-its-transpose)

    system['eigenvectors_inv'] = linalg.inv(system['eigenvectors'])

    print("[ DONE ]")

    # Evaluate the diagonalization
    # ============================

    # Using NumPy to convert the MIME from the zero order basis set to the eigenstates basis set through a matrix product (see https://numpy.org/doc/stable/reference/generated/numpy.matmul.html#numpy.matmul for reference)
    mime_diag = np.matmul(np.matmul(system['eigenvectors_inv'],system['mime']),system['eigenvectors'])

    # Get the absolute value of the matrix
    abs_mime_diag = np.abs(mime_diag)

    # Average the diagonal elements
    diag_mean = np.mean(np.trace(abs_mime_diag))

    # Average the non-diagonal elements (sum everything minus the trace and divide by n*(n-1), where n is the number of states)
    nondiag_mean = (np.sum(abs_mime_diag) - np.trace(abs_mime_diag)) / (len(system['eigenvectors']) * (len(system['eigenvectors']) - 1))

    # Evaluate the ratio of the diagonalization
    ratio = nondiag_mean / diag_mean
        
    print ("{:<50} {:<.2e}".format('\nDiagonalization ratio (non-diag/diag): ',ratio))

    # ========================================================= #
    #          Relativistic transition dipole moments           #
    # ========================================================= #

    print("{:<50}".format("\nComputing new transition dipole moments ..."), end="")

    # Convert the matrices
    # ====================    

    # Intialize the momdip_mtx dictionary

    system['momdip_mtx'] = {}

    # Convert each matrix from the non-relativistic basis set to the relativistic basis set through a matrix product (see https://numpy.org/doc/stable/reference/generated/numpy.matmul.html#numpy.matmul for reference)

    for momdip_key in system["momdip_o_mtx"]:
      system['momdip_mtx'][momdip_key] = np.absolute(np.matmul(np.matmul(system['eigenvectors_inv'],system['momdip_o_mtx'][momdip_key]),system['eigenvectors']))

    # Build the new list
    # ==================

    # Initialization of the variables

    rel_momdip_list = []
    
    # Iterate over the MIME to get the information and build the SOC list

    it = np.nditer(system['momdip_mtx']['X'], flags=['multi_index'])
  
    for momdip in it:

      state_1 = it.multi_index[0]
      state_2 = it.multi_index[1]
      value_x = momdip
      value_y = system['momdip_mtx']['Y'][state_1][state_2]
      value_z = system['momdip_mtx']['Z'][state_1][state_2]

      # Add the information to the soc_list
      momdip_line = (state_1, state_2, value_x, value_y, value_z)
      eq_line = False
      for line in rel_momdip_list:
        if line[0] == state_2 and line[1] == state_1:
          eq_line = True
          break
      if not eq_line:
        rel_momdip_list.append(momdip_line)

    print("[ DONE ]")

    # ========================================================= #
    #                     Mixing percentages                    #
    # ========================================================= #

    print("{:<50}".format("\nComputing mixing percentages ..."), end="")

    for state in system['states_list']:
      state['sing_percent'] = 0.0
      state['trip_percent'] = 0.0
      weights_list = [val**2 for val in np.absolute(system['eigenvectors'][state['number']])]
      for idx, weight in enumerate(weights_list):
        if system['zero_states_list'][idx]['label'].startswith('S'):
          state['sing_percent'] += weight
        elif system['zero_states_list'][idx]['label'].startswith('T'):
          state['trip_percent'] += weight        

    print("[ DONE ]")
    
    """
    
    # Radiative lifetime of excited states
    # ====================================

    #! Temporary: add the degeneracy key to the list of states

    for state in system['states_list']:
      state['degeneracy'] = 1

    # This calculation is based on the A_mn Einstein Coefficients and their link with the transition dipole moment
    # See https://aapt.scitation.org/doi/pdf/10.1119/1.12937 for reference
    # Note that this calculation is performed using atomic units, which means the Planck constant equals 2*pi and the vacuum permittivity equals 1/(4*pi)

    # Constants

    light_speed_au = constants.value('speed of light in vacuum') / constants.value('atomic unit of velocity')

    # Iterate over each excited state

    for state_m in system['states_list']:

      sum_einstein_coeffs = 0

      # Get the index of the state (to locate it in the matrices)

      m_index = system['states_list'].index(state_m)

      # Iterate over each state with an energy lower than the current one

      for state_n in [state_n for state_n in system['states_list'] if state_n['energy'] < state_m['energy']]:

        # Get the index of the state (to locate it in the matrices)

        n_index = system['states_list'].index(state_n)

        # Compute the energy difference

        energy_diff = state_m['energy'] - state_n['energy']

        # Compute the square of the transition dipole moment

        square_dipole = 0
        
        for momdip_key in system['momdip_mtx']:
          square_dipole += system['momdip_mtx'][momdip_key][m_index][n_index] ** 2

        # Calculate the A Einstein Coefficient          

        einstein_coeff = (state_n['degeneracy']/state_m['degeneracy']) * (4/3) * square_dipole * (energy_diff**3) / (light_speed_au**3)
        sum_einstein_coeffs += einstein_coeff

      # Compute the radiative lifetime

      if sum_einstein_coeffs == 0:
        state_m['lifetime'] = float('inf')
      else:
        state_m['lifetime'] = 1 / sum_einstein_coeffs
        state_m['lifetime'] = state_m['lifetime'] * constants.value('atomic unit of time')
    """

    # ========================================================= #
    #              Printing values in the log file              #
    # ========================================================= #

    # Print the relativistic states list
    # ==================================

    table_width = 73
    print("")
    print(''.center(table_width, '-'))
    print('Relativistic states list'.center(table_width, ' '))
    print(''.center(table_width, '-'))
    print("{:<9} {:<9} {:<15} {:<15} {:<10} {:<10}".format('Number','Label','Energy (Ha)','Energy (cm-1)','% Singlet','% Triplet'))
    print(''.center(table_width, '-'))
    for state in system['states_list']:
      print("{:<9} {:<9} {:<15.5e} {:<15.4f} {:<10.2f} {:<10.2f}".format(state['number'],state['label'],state['energy'],control_common.energy_unit_conversion(state['energy'],"ha","cm-1"),state['sing_percent']*100,state['trip_percent']*100))
    print(''.center(table_width, '-'))

    # Print the relativistic transition dipole moments list
    # =====================================================

    table_width = 73
    print("")
    print(''.center(table_width, '-'))
    print('Relativistic transition dipole moments'.center(table_width, ' '))
    print(''.center(table_width, '-'))
    print("{:<10} {:<10} {:<12} {:<12} {:<12} {:<12}".format('State 1','State 2','X (a.u.)','Y (a.u.)','Z (a.u.)','Tot (a.u.)'))
    print(''.center(table_width, '-'))
    for momdip in rel_momdip_list:
      column_6 = math.sqrt(momdip[2]**2 + momdip[3]**2 + momdip[4]**2)
      print("{:<10} {:<10} {:<12.4g} {:<12.4g} {:<12.4g} {:<12.4g}".format(momdip[0],momdip[1],momdip[2],momdip[3],momdip[4],column_6))
    print(''.center(table_width, '-'))

    # ========================================================= #
    #                    End of the function                    #
    # ========================================================= #
  
    #! Temporary: get the absolute values of the eigenvectors

    system['eigenvectors'] = np.absolute(system['eigenvectors'])
    system['eigenvectors_inv'] = np.absolute(system['eigenvectors_inv'])

    return system

# =================================================================== #
# =================================================================== #
#                   Q-Chem TD-DFT Modelling Function                  #
# =================================================================== #
# =================================================================== #

def qchem_tddft(source_content:list):
    """Parses the content of a Q-CHEM SOC TD-DFT calculation output file, looking to build the non-relavistic MIME and the non-relavistic transition dipole moments matrices of the molecule. It then diagonalizes the MIME to build the eigenstates basis set (relativistic states) and convert the dipole moments matrices into this new basis set.

    Parameters
    ----------
    source_content : list
        Content of the Q-CHEM output file. Each element of the list is a line of the file.
    
    Returns
    -------
    system : dict
        The extracted information of the source file. It contains the two mandatory keys and their associated values: ``states_list`` and ``momdip_mtx`` where
        
        - ``states_list`` is a list of dictionaries containing three keys each: ``energy``, ``number`` and ``label`` where

          - ``energy`` is the excitation energy of the relativistic state, in Ha.
          - ``number`` is the number of the state in an increasing order of energy, starting at 0 (which is the ground state).
          - ``label`` is the label of the state, in the form of "E" + number of that state.

        - ``momdip_mtx`` is a dictionary containing three keys and their associated values:
        
          - ``X`` is a NumPy array representing the relativistic transition dipole moments matrix along the X axis (in atomic units).
          - ``Y`` is a NumPy array representing the relativistic transition dipole moments matrix along the Y axis (in atomic units).
          - ``Z`` is a NumPy array representing the relativistic transition dipole moments matrix along the Z axis (in atomic units).

        It also contains two additional keys and their associated values that will be used by the rendering function: ``eigenvectors``, ``eigenvectors_inv`` where

        - ``eigenvectors`` is a NumPy array, containing the eigenvectors of the eigenstates in its columns.
        - ``eigenvectors_inv`` is the inverse of the eigenvectors matrix, containing the eigenvectors of the eigenstates in its lines.

        It also contains three additional keys and their associated values that will be used by another script: ``zero_states_list``, ``mime`` and ``momdip_o_mtx`` where
        
        - ``zero_states_list`` is a list of dictionaries containing four keys each: ``number``, ``type``, ``label`` and ``energy`` where

          - ``number`` is the number of the non-relativistic state, starting at 0 (which is the ground state).
          - ``type`` reflects the selection rule for transitions going from the ground state to this state, i.e. either "Bright" or "Dark".
          - ``energy`` is the excitation energy of the state, in Ha.
          - ``label`` is the label of the state, in the form of the first letter of its multiplicity + number of that state of this multiplicity (e.g. T1 for the first triplet, S3 for the third singlet, ...).

        - ``mime`` is a NumPy array, representing the Matrix Image of the MoleculE which acts as an effective Hamiltonian. It contains the excitation energies on the diagonal elements, and the spin-orbit couplings on the non-diagonal elements (in Hartree).

        - ``momdip_o_mtx`` is a dictionary containing three keys and their associated values:
        
          - ``X`` is a NumPy array representing the transition dipole moments matrix along the X axis (in atomic units).
          - ``Y`` is a NumPy array representing the transition dipole moments matrix along the Y axis (in atomic units).
          - ``Z`` is a NumPy array representing the transition dipole moments matrix along the Z axis (in atomic units).

    Raises
    ------
    ControlError
        If some of the needed values are missing or unknown.

    """

    # Initialize the system dictionary that will be returned by the function

    system = {}

    # Define an array for correct English spelling during printing

    special_numbers = {1:"st", 2:"nd", 3:"rd"}

    # =================================================================== #
    # =================================================================== #
    #                  Non-Relativistic States Treatment                  #
    # =================================================================== #
    # =================================================================== #

    section_title = "Non-Relativistic States Treatment"

    print("")
    print(''.center(len(section_title)+30, '='))
    print(section_title.center(len(section_title)+30))
    print(''.center(len(section_title)+30, '='))

    # ========================================================= #
    #                      Number of roots                      #
    # ========================================================= #

    nb_roots = False
    old_nb_roots = False

    # Define the expression patterns for the lines containing information about the number of roots, that number will be used to check if all the needed values have been collected.

    nb_roots_rx = {

      # Pattern for finding lines looking like "CIS_N_ROOTS          4"
      'normal': re.compile(r'^\s*CIS_N_ROOTS\s*(?P<n_roots>\d+)\s*$'),

      # Pattern for finding lines looking like "NRoots was altered as:  4 --> 12"
      'altered': re.compile(r'^\s*NRoots was altered as:\s*\d+\s*-->\s*(?P<new_n_roots>\d+)\s*$'),

      # Pattern for finding the "TDDFT/TDA Excitation Energies" line (which marks the end of the section)
      'end': re.compile(r'TDDFT/TDA Excitation Energies')

    }
  
    # Parse the source file to get the information

    for line in source_content:

      # Define when the section ends

      if nb_roots_rx['end'].match(line): # The line matches the 'end' pattern
        break  

      # Get the normal number of roots

      elif nb_roots_rx['normal'].match(line):
        nb_roots = int(nb_roots_rx['normal'].match(line).group('n_roots'))

      # Get the altered number of roots

      elif nb_roots_rx['altered'].match(line):
        old_nb_roots = nb_roots
        nb_roots = int(nb_roots_rx['altered'].match(line).group('new_n_roots'))

    # Raise an exception if the number of roots has not been found

    if not nb_roots:
      raise control_common.ControlError ("ERROR: Unable to find the number of roots in the source file")

    print("{:<50} {:<10}".format("\nNumber of roots: ",nb_roots))
    if old_nb_roots:
      print("{:<50} {:<10}".format("\nInitial number of roots: ",old_nb_roots))

    # ========================================================= #
    #                        Zero states List                   #
    # ========================================================= #

    print("{:<50}".format("\nParsing the excited states ...  "), end="")

    # Initialization of the variables

    section_found = False
    cnt_state = 0
    cnt_triplet = 0
    cnt_singlet = 0
    state_number = 0

    # Define the ground state of our molecule (which is the first state and has a zero energy)
    
    system['zero_states_list'] = [{'number': 0, 'qchem_number': 0, 'label': 'S0', 'energy' : 0.0}]

    # Define the expression patterns for the lines containing information about the states

    states_rx = {

      # Pattern for finding the "TDDFT/TDA Excitation Energies" line (which marks the start of the section)
      'start': re.compile(r'TDDFT/TDA Excitation Energies'),
      
      # Pattern for finding lines looking like 'Excited state   1: excitation energy (eV) =    4.6445'
      'state_energy': re.compile(r'^\s*Excited state\s+(?P<state>\d+): excitation energy \(eV\) =\s+(?P<energy>[+]?\d*\.\d+|\d+)$'),

      # Pattern for finding lines looking like '    Multiplicity: Triplet'
      'state_mp': re.compile(r'^\s*Multiplicity: (?P<mp>\w+)$'),

      # Pattern for finding the "SETman timing summary" line (which marks the end of the section)
      'end': re.compile(r'SETman timing summary')

    }

    # Parse the source file to get the information and build the states list

    for line in source_content:

      # Define when the section begins and ends

      if not section_found:
        if states_rx['start'].match(line):
          section_found = True
    
      elif section_found and states_rx['end'].match(line):
        break
        
      # Get the state number and its associated excitation energy

      elif states_rx['state_energy'].match(line):

        exc_state = int(states_rx['state_energy'].match(line).group("state"))
        exc_energy = float(states_rx['state_energy'].match(line).group("energy"))
        cnt_state += 1
            
      # Get the corresponding state multiplicity and add the data to the zero_states_list

      elif states_rx['state_mp'].match(line):

        multiplicity = states_rx['state_mp'].match(line).group("mp")

        # Check the multiplicity and increase the counter (cnt) for that multiplicity (needed for the state label)

        cnt = -1

        if multiplicity == "Triplet":
          first_letter = "T"
          cnt_triplet += 1
          cnt = cnt_triplet

        elif multiplicity == "Singlet":
          first_letter = "S"
          cnt_singlet += 1
          cnt = cnt_singlet

        else:
          raise control_common.ControlError ("ERROR: Multiplicity of the %s%s state is of unknown value (%s)" % (exc_state, ("th" if not exc_state in special_numbers else special_numbers[exc_state]),multiplicity))

        # Append information about the current state to the zero_states_list key

        if multiplicity == "Singlet":
          state_number += 1
          system['zero_states_list'].append({'number': state_number, 'qchem_number': exc_state, 'label': (first_letter + str(cnt)), 'energy': control_common.energy_unit_conversion(exc_energy,"ev","cm-1")})

        elif multiplicity == "Triplet":

          # Add the substate ms=0
          state_number += 1
          system['zero_states_list'].append({'number': state_number, 'qchem_number': exc_state, 'label': (first_letter + str(cnt) + "(ms=0)"), 'energy': control_common.energy_unit_conversion(exc_energy,"ev","cm-1")})
          # Add the substate ms=1
          state_number += 1
          system['zero_states_list'].append({'number': state_number, 'qchem_number': exc_state, 'label': (first_letter + str(cnt) + "(ms=1)"), 'energy': control_common.energy_unit_conversion(exc_energy,"ev","cm-1")})
          # Add the substate ms=-1
          state_number += 1
          system['zero_states_list'].append({'number': state_number, 'qchem_number': exc_state, 'label': (first_letter + str(cnt) + "(ms=-1)"), 'energy': control_common.energy_unit_conversion(exc_energy,"ev","cm-1")})
          
    # Raise an exception if the section has not been found

    if not section_found:
      raise control_common.ControlError ("ERROR: Unable to find the 'TDDFT/TDA Excitation Energies' section in the source file")

    # Raise an exception if not all the values have been found

    if cnt_state != 2*nb_roots:
      raise control_common.ControlError ("ERROR: The modelling function could not find the right number of excited states in the source file (%s of the %s expected states have been found)" % (cnt_state,2*nb_roots))
    if cnt_triplet != nb_roots:
      raise control_common.ControlError ("ERROR: The modelling function could not find the right number of excited triplet states in the source file (%s of the %s expected triplet states have been found)" % (cnt_triplet,nb_roots))
    if cnt_singlet != nb_roots:
      raise control_common.ControlError ("ERROR: The modelling function could not find the right number of excited singlet states in the source file (%s of the %s expected singlet states have been found)" % (cnt_singlet,nb_roots))

    # Raise an exception if the state numbers are not consecutive and starting at 0

    control_common.is_consecutive(list(dict.fromkeys([state['qchem_number'] for state in system['zero_states_list']])),"Excited state numbers from the source file")

    print("[ DONE ]")
  
    # ========================================================= #
    #                          SOC List                         #
    # ========================================================= #

    print("{:<50}".format("\nParsing the spin-orbit couplings ...  "), end="")

    # Initialization of the variables

    section_found = False
    soc_list = []
    
    # Define the expression patterns for the lines containing information about the SOC
    
    soc_rx = {

      # Pattern for finding the "SPIN-ORBIT COUPLING JOB BEGINS HERE" line (which marks the start of the section)
      'start': re.compile(r'^\*+SPIN-ORBIT COUPLING JOB BEGINS HERE\*+$'),

      # Pattern for finding lines looking like 'SOC between the singlet ground state and excited triplet states (ms=0):'
      'ground_to_triplets': re.compile(r'^\s*SOC between the singlet ground state and excited triplet states \(ms=-?\d\):$'),

      # Pattern for finding lines looking like 'SOC between the S9 state and excited triplet states (ms=1):'
      'singlet_to_triplets': re.compile(r'^\s*SOC between the (?P<state_1>[A-Z]\d+) state and excited triplet states \(ms=-?\d\):$'),

      # Pattern for finding lines looking like 'SOC between the T7 (ms=0) state and excited triplet states (ms=1):'
      'triplet_to_triplets': re.compile(r'^\s*SOC between the (?P<state_1>[A-Z]\d+) \(ms=(?P<ms>-?\d)\) state and excited triplet states \(ms=-?\d\):$'),      

      # Pattern for finding lines looking like 'T5(ms=-1)     0.000000  + (0.002867i)    cm-1'
      'soc_value': re.compile(r'^\s*(?P<state_2>[A-Z]\d+)\(ms=(?P<ms>-?\d)\)\s+(?P<real_value>-?\(?-?\d+\.?\d*\)?)\s+(?P<im_value>[\+-]\s+\(?-?\d+\.?\d*i\)?)\s+cm-1$'),

      # Pattern for finding the "SOC CODE ENDS HERE" line (which marks the end of the section)
      'end': re.compile(r'^\*+SOC CODE ENDS HERE\*+$')

    }

    # Parse the source file to get the information and build the SOC list

    for line in source_content:

      # Define when the section begins and ends

      if not section_found:
        if soc_rx['start'].match(line):
          section_found = True
    
      elif section_found and soc_rx['end'].match(line):
        break

      # Get the number and label of the first state

      elif soc_rx['ground_to_triplets'].match(line):
        state_1 = 0
        label_1 = "S0"

      elif soc_rx['singlet_to_triplets'].match(line):

        label_1 = soc_rx['singlet_to_triplets'].match(line).group('state_1')

      elif soc_rx['triplet_to_triplets'].match(line):

        label_1 = soc_rx['triplet_to_triplets'].match(line).group('state_1') + "(ms=%s)" % soc_rx['triplet_to_triplets'].match(line).group('ms')
        
      # Get the number and label of the second state and the corresponding SOC value, before adding the data to the soc_list

      elif soc_rx['soc_value'].match(line):

        label_2 = soc_rx['soc_value'].match(line).group('state_2') + "(ms=%s)" % soc_rx['soc_value'].match(line).group('ms')

        # Convert the labels to the state numbers

        state_1 = -1
        for state in system['zero_states_list']:
          if label_1 == state['label']:
            state_1 = state['number']
            break
        if state_1 == -1:
          raise control_common.ControlError ("ERROR: Unknown excited state (%s) has been catched during the SOC parsing." % label_1)

        state_2 = -1
        for state in system['zero_states_list']:
          if label_2 == state['label']:
            state_2 = state['number']
            break
        if state_2 == -1:
          raise control_common.ControlError ("ERROR: Unknown excited state (%s) has been catched during the SOC parsing." % label_2)
          
        # Get the complex value of the SOC
        
        real_raw = re.sub(r'[\+\s]', '', soc_rx['soc_value'].match(line).group('real_value'))
        
        if real_raw.startswith('-(') and real_raw.endswith(')'):
          real_raw = (real_raw.replace('-(','')).replace(')','')
          real_value = -float(real_raw)
        else:
          real_raw = (real_raw.replace('(','')).replace(')','')
          real_value = float(real_raw)
        
        im_raw = re.sub(r'[i\+\s]', '', soc_rx['soc_value'].match(line).group('im_value'))
        
        if im_raw.startswith('-(') and im_raw.endswith(')'):
          im_raw = (im_raw.replace('-(','')).replace(')','')
          im_value = -float(im_raw)
        else:
          im_raw = (im_raw.replace('(','')).replace(')','')
          im_value = float(im_raw)
  
        value = complex(real_value,im_value)
        
        # Add the information to the soc_list

        soc_line = (state_1, label_1, state_2, label_2, value)
        soc_list.append(soc_line)

    # Raise an exception if the section has not been found

    if not section_found:
      raise control_common.ControlError ("ERROR: Unable to find the 'SPIN-ORBIT COUPLING' section in the source file")

    # Raise an exception if not all the values have been found

    nb_soc = 3*nb_roots*(nb_roots+1) + ((6*nb_roots*(nb_roots-1))/2) # Singlet to Triplet + Triplet to Triplet
    if len(soc_list) != nb_soc:
      raise control_common.ControlError ("ERROR: The modelling function could not find the right number of spin-orbit couplings in the source file (%s of the %s expected values have been found)" % (len(soc_list),nb_soc))

    print("[ DONE ]")

    # ========================================================= #
    #                       MIME Creation                       #
    # ========================================================= #

    print("{:<50}".format("\nBuilding the MIME ... "), end="")   

    # Initialize the MIME as a zero-filled matrix

    system['mime'] = np.zeros((len(system['zero_states_list']), len(system['zero_states_list'])), dtype=complex)

    # Creation of the MIME - Non-diagonal values (SOC)

    for soc in soc_list:
      k1 = soc[0]
      k2 = soc[2]
      val = soc[4]
      system['mime'][k1][k2] = val
      system['mime'][k2][k1] = val.conjugate()

    # Creation of the MIME - Diagonal values (Excitation energies)

    for state in system['zero_states_list']:
      system['mime'][state['number']][state['number']] = complex(state['energy'])

    # Convert the MIME to Hartree units

    system['mime'] = system['mime'] * control_common.energy_unit_conversion(1,"cm-1","ha")

    print("[ DONE ]")

    # ========================================================= #
    #                    Dipole Moments List                    #
    # ========================================================= #

    print("{:<50}".format("\nParsing the transition dipole moments ...  "), end="")

    # Initialization of the variables

    section_found = False
    momdip_list = []

    # Define the expression patterns for the lines containing information about the dipole moments

    moment_rx = {

      # Pattern for finding the "STATE-TO-STATE TRANSITION MOMENTS" line (which marks the start of the section)
      'start': re.compile(r'^STATE-TO-STATE TRANSITION MOMENTS$'),

      # Pattern for finding lines looking like '    1    2   0.001414  -0.001456   0.004860   1.240659E-10'
      'moment': re.compile(r'^\s*(?P<state_1>\d+)\s+(?P<state_2>\d+)\s+(?P<mom_x>-?\d+\.\d+)\s+(?P<mom_y>-?\d+\.\d+)\s+(?P<mom_z>-?\d+\.\d+)\s+(?P<strength>\d|\d\.\d+|\d\.\d+E[-+]\d+)$'),

      # Pattern for finding the "END OF TRANSITION MOMENT CALCULATION" line (which marks the end of the section)
      'end': re.compile(r'^END OF TRANSITION MOMENT CALCULATION$')

    }

    # Parse the source file to get the information and build the dipole moments list

    for line in source_content:

      # Define when the section begins and ends

      if not section_found:
        if moment_rx['start'].match(line):
          section_found = True
    
      elif section_found and moment_rx['end'].match(line):
        break

      # Extract the relevant information and add it to the momdip_list
  
      elif moment_rx['moment'].match(line):
  
        matching_line = moment_rx['moment'].match(line)

        state_1 = int(matching_line.group('state_1'))
        state_2 = int(matching_line.group('state_2'))
        value_x = float(matching_line.group('mom_x'))
        value_y = float(matching_line.group('mom_y'))
        value_z = float(matching_line.group('mom_z'))
        momdip = (state_1, state_2, value_x, value_y, value_z)
      
        # Add the new line to the momdip_list
        momdip_list.append(momdip)

    # Raise an exception if the section has not been found

    if not section_found:
      raise control_common.ControlError ("ERROR: Unable to find the 'STATE-TO-STATE TRANSITION MOMENTS' section in the source file")

    # Raise an exception if not all the values have been found

    nb_momdip = nb_roots * (nb_roots+1) # Ground to Singlet (nb_roots) + Ground to Triplet (nb_roots) + Singlet to Singlet (nb_roots*(nb_roots-1)/2) + Triplet to Triplet (nb_roots*(nb_roots-1)/2)
    if len(momdip_list) != nb_momdip:
      raise control_common.ControlError ("ERROR: The modelling function could not find the right number of state-to-state transition moments in the source file (%s of the %s expected values have been found)" % (len(momdip_list),nb_momdip))

    print("[ DONE ]")

    # ========================================================= #
    #                   Dipole Moments Matrices                 #
    # ========================================================= #

    print("{:<50}".format("\nBuilding transition dipole moments matrices ... "), end="")

    # Intialize the momdip_o_mtx dictionary

    system['momdip_o_mtx'] = {}

    # Initialize the matrices as zero-filled matrices

    system['momdip_o_mtx']['X'] = np.zeros((len(system['zero_states_list']), len(system['zero_states_list'])), dtype=float)
    system['momdip_o_mtx']['Y'] = np.zeros((len(system['zero_states_list']), len(system['zero_states_list'])), dtype=float)
    system['momdip_o_mtx']['Z'] = np.zeros((len(system['zero_states_list']), len(system['zero_states_list'])), dtype=float)

    # Filling the matrices
  
    for momdip in momdip_list:

      states_1 = [state['number'] for state in system['zero_states_list'] if state['qchem_number'] == momdip[0]] # For triplets, there are three states for each 'qchem_number'
      states_2 = [state['number'] for state in system['zero_states_list'] if state['qchem_number'] == momdip[1]]

      for k1 in states_1:
        for k2 in states_2:
          
          system['momdip_o_mtx']['X'][k1][k2] = momdip[2]
          system['momdip_o_mtx']['X'][k2][k1] = momdip[2]    # For symetry purposes
    
          system['momdip_o_mtx']['Y'][k1][k2] = momdip[3]
          system['momdip_o_mtx']['Y'][k2][k1] = momdip[3]    # For symetry purposes
    
          system['momdip_o_mtx']['Z'][k1][k2] = momdip[4]
          system['momdip_o_mtx']['Z'][k2][k1] = momdip[4]    # For symetry purposes

    print("[ DONE ]")

    # ========================================================= #
    #           Reducing number of states (if needed)           #
    # ========================================================= #

    if old_nb_roots:

      # Define the states that need to be kept
      
      states_to_keep = []
      qchem_states_to_keep = []
      for state in system['zero_states_list']:
        label = state['label'].partition("(")[0]
        label_number = int(re.sub(r'[a-zA-Z]','',label))
        if label_number <= old_nb_roots:
          states_to_keep.append(state['number'])
          qchem_states_to_keep.append(state['qchem_number'])

      # Remove the rest from the MIME and the transition dipole moment matrices

      system['mime'] = system['mime'][np.ix_(states_to_keep,states_to_keep)]
      for momdip_key in system['momdip_o_mtx']:
        system['momdip_o_mtx'][momdip_key] = system['momdip_o_mtx'][momdip_key][np.ix_(states_to_keep,states_to_keep)]

      # Remove the rest from the states list
        
      system['zero_states_list'] = [state for state in system['zero_states_list'] if state['number'] in states_to_keep]
      for state in system['zero_states_list']:
        state['number'] = system['zero_states_list'].index(state)

      # Remove the rest from the SOC list

      soc_list = [soc for soc in soc_list if soc[0] in states_to_keep and soc[2] in states_to_keep]

      # Remove the rest from the transition dipole moments list
    
      momdip_list = [momdip for momdip in momdip_list if momdip[0] in qchem_states_to_keep and momdip[1] in qchem_states_to_keep]

    """        
    # ========================================================= #
    #           Radiative lifetime of excited states            #
    # ========================================================= #

    print("{:<50}".format("\nComputing radiative lifetimes ... "), end="")

    # This calculation is based on the A_mn Einstein Coefficients and their link with the transition dipole moment
    # See https://aapt.scitation.org/doi/pdf/10.1119/1.12937 for reference
    # Note that this calculation is performed using atomic units, which means the Planck constant equals 2*pi and the vacuum permittivity equals 1/(4*pi)

    # Constants

    light_speed_au = constants.value('speed of light in vacuum') / constants.value('atomic unit of velocity')

    # Calculate the radiative lifetime of each excited state

    for state in system['zero_states_list']:

      sum_einstein_coeffs = 0

      for other_state in system['zero_states_list']:

        if other_state['energy'] < state['energy']:

          # Compute and convert the energy difference

          energy_diff = control_common.energy_unit_conversion(state['energy'],"cm-1","ha") - control_common.energy_unit_conversion(other_state['energy'],"cm-1","ha")

          # Compute the square of the transition dipole moment

          square_dipole = 0
          for momdip_key in system['momdip_o_mtx']:
            square_dipole += system['momdip_o_mtx'][momdip_key][state['number']][other_state['number']] ** 2

          # Calculate the A Einstein Coefficient   
         
          einstein_coeff = (4/3) * square_dipole * (energy_diff**3) / (light_speed_au**3)
          sum_einstein_coeffs += einstein_coeff

      if sum_einstein_coeffs == 0:
        state['lifetime'] = float('inf')
      else:
        state['lifetime'] = 1 / sum_einstein_coeffs
        state['lifetime'] = state['lifetime'] * constants.value('atomic unit of time')

    print("[ DONE ]")
    """

    # ========================================================= #
    #              Printing values in the log file              #
    # ========================================================= #

    # Print the non-relativistic states list
    # ======================================

    table_width = 73
    print("")
    print(''.center(table_width, '-'))
    print('Non-relativistic states list'.center(table_width, ' '))
    print(''.center(table_width, '-'))
    print("{:<10} {:<15} {:<15} {:<15} {:<15}".format('Number','Label','Energy (cm-1)','Energy (Ha)','Energy (nm)'))
    print(''.center(table_width, '-'))
    for state in system['zero_states_list']:
      print("{:<10} {:<15} {:<15.2f} {:<15.5f} {:<15.2f}".format(state['number'],state['label'],state['energy'],control_common.energy_unit_conversion(state['energy'],"cm-1","ha"),control_common.energy_unit_conversion(state['energy'],"cm-1","nm")))
    print(''.center(table_width, '-'))

    # Print the SOC list
    # ==================

    table_width = 73
    print("")
    print(''.center(table_width, '-'))
    print('Spin-orbit couplings'.center(table_width, ' '))
    print(''.center(table_width, '-'))
    print("{:<15} {:<15} {:<20} {:<20}".format('State 1','State 2','Real value (Ha)','Imag value (Ha)'))
    print(''.center(table_width, '-'))
    for soc in soc_list:
      column_1 = soc[1]
      column_2 = soc[3]
      column_3 = control_common.energy_unit_conversion(soc[4].real,"cm-1","ha")
      column_4 = control_common.energy_unit_conversion(soc[4].imag,"cm-1","ha")
      print("{:<15} {:<15} {:<20.5g} {:<20.5g}".format(column_1,column_2,column_3,column_4))
    print(''.center(table_width, '-'))  

    # Print the transition dipole moments list
    # ========================================

    table_width = 73
    print("")
    print(''.center(table_width, '-'))
    print('Transition dipole moments'.center(table_width, ' '))
    print(''.center(table_width, '-'))
    print("{:<10} {:<10} {:<12} {:<12} {:<12} {:<12}".format('State 1','State 2','X (a.u.)','Y (a.u.)','Z (a.u.)','Tot (a.u.)'))
    print(''.center(table_width, '-'))
    for momdip in momdip_list:
      column_1 = [state['label'] for state in system['zero_states_list'] if state['qchem_number'] == momdip[0]][0].partition("(")[0]
      column_2 = [state['label'] for state in system['zero_states_list'] if state['qchem_number'] == momdip[1]][0].partition("(")[0]
      column_6 = math.sqrt(momdip[2]**2 + momdip[3]**2 + momdip[4]**2)
      print("{:<10} {:<10} {:<12.4g} {:<12.4g} {:<12.4g} {:<12.4g}".format(column_1,column_2,momdip[2],momdip[3],momdip[4],column_6))
    print(''.center(table_width, '-'))

    # =================================================================== #
    # =================================================================== #
    #                    Relativistic States Treatment                    #
    # =================================================================== #
    # =================================================================== #

    section_title = "Relativistic States Treatment"

    print("\n\n")
    print(''.center(len(section_title)+30, '='))
    print(section_title.center(len(section_title)+30))
    print(''.center(len(section_title)+30, '='))

    # ========================================================= #
    #                    MIME diagonalization                   #
    # ========================================================= #

    print("{:<50}".format("\nDiagonalizing the MIME ..."), end="")

    # Diagonalization
    # ===============

    # Use SciPy to diagonalize the matrix (see https://personal.math.ubc.ca/~pwalls/math-python/linear-algebra/eigenvalues-eigenvectors/ for reference)
    # Documentation page for the function used here : https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.eigh.html  
     
    eigenvalues, system['eigenvectors'] = linalg.eigh(system['mime'])
    eigenvalues = eigenvalues.tolist()

    # Build the eigenstates list
    # ==========================

    # Initialize the variables

    system['states_list'] = []

    # Build the list

    for eigenvalue in eigenvalues:

      eigenstate = {
        'energy' : eigenvalue,
        'number' : eigenvalues.index(eigenvalue),
        'label' : "E" + str(eigenvalues.index(eigenvalue))
      }

      system['states_list'].append(eigenstate)

    # Transpose the eigenvectors list
    # ===============================

    # Using SciPy to invert the eigenvectors matrix
    # Note that the inverse of an orthonormal matrix is equal to its transpose, so each line of this matrix corresponds to an eigenvector. (see https://math.stackexchange.com/questions/156735/in-which-cases-is-the-inverse-of-a-matrix-equal-to-its-transpose)

    system['eigenvectors_inv'] = linalg.inv(system['eigenvectors'])

    print("[ DONE ]")

    # Evaluate the diagonalization
    # ============================

    # Using NumPy to convert the MIME from the zero order basis set to the eigenstates basis set through a matrix product (see https://numpy.org/doc/stable/reference/generated/numpy.matmul.html#numpy.matmul for reference)
    mime_diag = np.matmul(np.matmul(system['eigenvectors_inv'],system['mime']),system['eigenvectors'])

    # Get the absolute value of the matrix
    abs_mime_diag = np.abs(mime_diag)

    # Average the diagonal elements
    diag_mean = np.mean(np.trace(abs_mime_diag))

    # Average the non-diagonal elements (sum everything minus the trace and divide by n*(n-1), where n is the number of states)
    nondiag_mean = (np.sum(abs_mime_diag) - np.trace(abs_mime_diag)) / (len(system['eigenvectors']) * (len(system['eigenvectors']) - 1))

    # Evaluate the ratio of the diagonalization
    ratio = nondiag_mean / diag_mean
        
    print ("{:<50} {:<.2e}".format('\nDiagonalization ratio (non-diag/diag): ',ratio))

    # ========================================================= #
    #          Relativistic transition dipole moments           #
    # ========================================================= #

    print("{:<50}".format("\nComputing new transition dipole moments ..."), end="")

    # Convert the matrices
    # ==================== 

    # Intialize the momdip_mtx dictionary

    system['momdip_mtx'] = {}

    # Convert each matrix from the zero order basis set to the eigenstates basis set through a matrix product (see https://numpy.org/doc/stable/reference/generated/numpy.matmul.html#numpy.matmul for reference)

    for momdip_key in system["momdip_o_mtx"]:
      system['momdip_mtx'][momdip_key] = np.absolute(np.matmul(np.matmul(system['eigenvectors_inv'],system['momdip_o_mtx'][momdip_key]),system['eigenvectors']))

    """
    # ========================================================= #
    # Handling states degeneracies                              #
    # ========================================================= #

    deg_threshold = 1e-5
    print ("{:<50} {:<.2e}".format('\nDegeneracy threshold: ',deg_threshold))

    # Add the degeneracy key to the list of states

    for state in system['states_list']:
      state['degeneracy'] = 1

    # Initialize the list of degeneracies (each item of this list will be a group of degenerated states)

    deg_list = []

    # Look for every group of degenerated states

    for state in system['states_list']:

      for other_state in [other_state for other_state in system['states_list'] if other_state['number'] < state['number']]:

        if abs(state['energy'] - other_state['energy']) < deg_threshold:

          # Define if and how to add this pair of degenerated states to the list of degeneracies

          added = False

          for group in deg_list:

            if other_state['number'] in group and state['number'] not in group:
              group.append(state['number'])
              added = True

            elif state['number'] in group and other_state['number'] not in group:
              group.append(other_state['number'])
              added = True
              
            elif other_state['number'] in group and state['number'] in group:
              added = True

          if not added:
              deg_list.append([state['number'],other_state['number']])

    # Update the transition dipole moments matrices
    # =============================================

    # Build a version of the list of degeneracies that uses the indices of the states rather than their number (to locate them in the matrices)

    deg_list_ind = []

    for group in deg_list:

      group_ind = []

      for number in group:

        # Fetch the index corresponding to a particular state through its number and add it to the group_ind
        index = system['states_list'].index(next(state for state in system['states_list'] if state['number'] == number))
        group_ind.append(index)
      
      deg_list_ind.append(group_ind)

    # Iterate over each matrix separately

    for momdip_key in system["momdip_mtx"]:

      # Initialize the new matrix that will replace the old one

      new_mtx = np.copy(system["momdip_mtx"][momdip_key])

      # Intialize the list that will contain the indices of the lines that need to be removed

      to_remove = []

      # Iterate over each group separately

      for group_ind in deg_list_ind:

        # Determine where the new line will be placed

        spot = min(group_ind)

        # Determine the new line by taking the maximum (disregarding the sign) of each dipole moments from each line of the group
        # See https://stackoverflow.com/questions/51209928/get-maximum-of-absolute-along-axis for details

        max_idx = np.argmax(np.absolute(new_mtx[group_ind]), axis=0)
        values = new_mtx[group_ind][tuple([max_idx,np.arange(new_mtx[group_ind].shape[1])])]

        #! Other possibility: determine the new line by summing the dipole moments from each line of the group using the reduce method from NumPy (see https://numpy.org/doc/stable/reference/generated/numpy.ufunc.reduce.html)
        #! values = np.sum.reduce(new_mtx[group_ind])

        # Insert the new line and prepare to remove the old ones (do not remove them immediately to not mess with the other groups)

        for index in group_ind:
          if index == spot:
            new_mtx[spot] = values
            new_mtx[:,spot] = values
          else:
            to_remove.append(index)

      # Remove the now useless lines and columns

      new_mtx = np.delete(np.delete(new_mtx,to_remove,0), to_remove, 1)

      # Replace the old matrix with the new one

      system["momdip_mtx"][momdip_key] = new_mtx

      #! Temporary (print it as a table rather than a matrix)

      print("\nDipole moments matrix with the '%s' key (atomic units)" % momdip_key)
      print('')
      for row in system['momdip_mtx'][momdip_key]:
        for val in row:
          print(np.format_float_scientific(val,precision=3,unique=False,pad_left=2), end = " ")
        print('')

    # Update the states list
    # ===========================

    for group in deg_list:

      # Initialize the new state resulting from the combination

      new_state = {}

      # Determine the new values for the new degenerated state

      new_state['number'] = min(group)
      new_state['label'] = "E" + "-".join(map(str,sorted(group))) # e.g. for a group including the states number 2, 3 and 4, its label will be E2-3-4
      new_state['energy'] = np.mean([state['energy'] for state in system['states_list'] if state['number'] in group])
      new_state['degeneracy'] = len(group)

      # Get the index of the state that has the same number as the new one

      index = system['states_list'].index(next(state for state in system['states_list'] if state['number'] == new_state['number']))

      # Remove the degenerated states from the states_list (by keeping only the states not included in the group)

      system['states_list'] = [state for state in system['states_list'] if state['number'] not in group]

      # Add the new state to the list at the position occupied by the state that had the same number

      system['states_list'].insert(index,new_state)

    # Once all the list has been updated, correct the state numbers

    energies = sorted([state['energy'] for state in system['states_list']])

    for state in system['states_list']:
      state['number'] = energies.index(state['energy'])
    """

    # Build the new list
    # ==================

    # Initialization of the variables

    rel_momdip_list = []
    
    # Iterate over the MIME to get the information and build the SOC list

    it = np.nditer(system['momdip_mtx']['X'], flags=['multi_index'])
  
    for momdip in it:

      state_1 = it.multi_index[0]
      state_2 = it.multi_index[1]
      value_x = momdip
      value_y = system['momdip_mtx']['Y'][state_1][state_2]
      value_z = system['momdip_mtx']['Z'][state_1][state_2]

      # Add the information to the soc_list
      momdip_line = (state_1, state_2, value_x, value_y, value_z)
      eq_line = False
      for line in rel_momdip_list:
        if line[0] == state_2 and line[1] == state_1:
          eq_line = True
          break
      if not eq_line:
        rel_momdip_list.append(momdip_line)

    print("[ DONE ]")

    # ========================================================= #
    #                     Mixing percentages                    #
    # ========================================================= #

    print("{:<50}".format("\nComputing mixing percentages ..."), end="")

    for state in system['states_list']:
      state['sing_percent'] = 0.0
      state['trip_percent'] = 0.0
      weights_list = [val**2 for val in np.absolute(system['eigenvectors'][state['number']])]
      for idx, weight in enumerate(weights_list):
        if system['zero_states_list'][idx]['label'].startswith('S'):
          state['sing_percent'] += weight
        elif system['zero_states_list'][idx]['label'].startswith('T'):
          state['trip_percent'] += weight        

    print("[ DONE ]")

    """
    # Radiative lifetime of excited states
    # ====================================

    #! Temporary: add the degeneracy key to the list of states

    for state in system['states_list']:
      state['degeneracy'] = 1

    # This calculation is based on the A_mn Einstein Coefficients and their link with the transition dipole moment
    # See https://aapt.scitation.org/doi/pdf/10.1119/1.12937 for reference
    # Note that this calculation is performed using atomic units, which means the Planck constant equals 2*pi and the vacuum permittivity equals 1/(4*pi)

    # Constants

    light_speed_au = constants.value('speed of light in vacuum') / constants.value('atomic unit of velocity')

    # Iterate over each excited state

    for state_m in system['states_list']:

      sum_einstein_coeffs = 0

      # Get the index of the state (to locate it in the matrices)

      m_index = system['states_list'].index(state_m)

      # Iterate over each state with an energy lower than the current one

      for state_n in [state_n for state_n in system['states_list'] if state_n['energy'] < state_m['energy']]:

        # Get the index of the state (to locate it in the matrices)

        n_index = system['states_list'].index(state_n)

        # Compute the energy difference

        energy_diff = state_m['energy'] - state_n['energy']

        # Compute the square of the transition dipole moment

        square_dipole = 0
        
        for momdip_key in system['momdip_mtx']:
          square_dipole += system['momdip_mtx'][momdip_key][m_index][n_index] ** 2

        # Calculate the A Einstein Coefficient          

        einstein_coeff = (state_n['degeneracy']/state_m['degeneracy']) * (4/3) * square_dipole * (energy_diff**3) / (light_speed_au**3)
        sum_einstein_coeffs += einstein_coeff

      # Compute the radiative lifetime

      if sum_einstein_coeffs == 0:
        state_m['lifetime'] = float('inf')
      else:
        state_m['lifetime'] = 1 / sum_einstein_coeffs
        state_m['lifetime'] = state_m['lifetime'] * constants.value('atomic unit of time')
    """

    # ========================================================= #
    #              Printing values in the log file              #
    # ========================================================= #

    # Print the relativistic states list
    # ==================================

    table_width = 73
    print("")
    print(''.center(table_width, '-'))
    print('Relativistic states list'.center(table_width, ' '))
    print(''.center(table_width, '-'))
    print("{:<9} {:<9} {:<15} {:<15} {:<10} {:<10}".format('Number','Label','Energy (Ha)','Energy (cm-1)','% Singlet','% Triplet'))
    print(''.center(table_width, '-'))
    for state in system['states_list']:
      print("{:<9} {:<9} {:<15.5e} {:<15.4f} {:<10.2f} {:<10.2f}".format(state['number'],state['label'],state['energy'],control_common.energy_unit_conversion(state['energy'],"ha","cm-1"),state['sing_percent']*100,state['trip_percent']*100))
    print(''.center(table_width, '-'))

    # Print the relativistic transition dipole moments list
    # =====================================================

    table_width = 73
    print("")
    print(''.center(table_width, '-'))
    print('Relativistic transition dipole moments'.center(table_width, ' '))
    print(''.center(table_width, '-'))
    print("{:<10} {:<10} {:<12} {:<12} {:<12} {:<12}".format('State 1','State 2','X (a.u.)','Y (a.u.)','Z (a.u.)','Tot (a.u.)'))
    print(''.center(table_width, '-'))
    for momdip in rel_momdip_list:
      column_6 = math.sqrt(momdip[2]**2 + momdip[3]**2 + momdip[4]**2)
      print("{:<10} {:<10} {:<12.4g} {:<12.4g} {:<12.4g} {:<12.4g}".format(momdip[0],momdip[1],momdip[2],momdip[3],momdip[4],column_6))
    print(''.center(table_width, '-'))

    # ========================================================= #
    #                    End of the function                    #
    # ========================================================= #

    # Converting the states energy from cm-1 to Ha

    for state in system['zero_states_list']:
      state['energy'] = control_common.energy_unit_conversion(state['energy'],"cm-1","ha")

    #! Temporary: get the absolute values of the eigenvectors

    system['eigenvectors'] = np.absolute(system['eigenvectors'])
    system['eigenvectors_inv'] = np.absolute(system['eigenvectors_inv'])

    return system

# =================================================================== #
# =================================================================== #
#                    Custom File Modelling Function                   #
# =================================================================== #
# =================================================================== #

def custom_file(source_content:list):
    """Parses the content of a custom source file, looking to build the MIME and the transition dipole moments matrices of the molecule. Consult the official documentation for more information on the format of the custom file.

    Parameters
    ----------
    source_content : list
        Content of the custom file. Each element of the list is a line of the file.
    
    Returns
    -------
    system : dict
        The extracted information of the source file. It contains three keys and their associated values: ``states_list``, ``mime`` and ``momdip_mtx`` where
        
        - ``states_list`` is a list of dictionaries containing four keys each: ``number``, ``type``, ``label`` and ``energy`` where

          - ``number`` is the number of the state, starting at 0 (which is the ground state).
          - ``type`` reflects the selection rule for transitions going from the ground state to this state, i.e. either "Bright" or "Dark".
          - ``energy`` is the excitation energy of the state, in Ha.
          - ``label`` is the label of the state, as specified in the custom file.

        - ``mime`` is a NumPy array, representing the Matrix Image of the MoleculE which acts as an effective Hamiltonian. It contains the state energies on the diagonal elements, and the state couplings on the non-diagonal elements (in Hartree).

        - ``momdip_mtx`` is a dictionary containing the transition dipole moments matrix (in atomic units).

    Raises
    ------
    ControlError
        If some of the needed values are missing or unknown.

    """

    # Initialize the system dictionary that will be returned by the function

    system = {}

    # ========================================================= #
    #                        States List                        #
    # ========================================================= #

    print("{:<50}".format("\nParsing the states list ...  "), end="")

    # Initialization of the variables

    section_found = False
    system['states_list'] = []

    # Define the expression patterns for the lines containing information about the states

    states_rx = {

      # Pattern for finding the "  1. States list" line (which marks the start of the section)
      'start': re.compile(r'^\s*1\.\s+States\s+list\s*$'),

      # Pattern for finding lines looking like 'Number      Type        Label       Energy (cm-1)'
      'unit_line': re.compile(r'^\s*Number\s+Type\s+Label\s+Energy\s+\((?P<unit>[a-zA-Z_0-9\-]+)\)\s*$'),

      # Pattern for finding lines looking like '1           Bright      EA          4375.800918'
      'state_line': re.compile(r'^\s*(?P<number>\d+)\s+(?P<type>[a-zA-Z_0-9\-]+)\s+(?P<label>[a-zA-Z_0-9\-]+)\s+(?P<energy>\d+|\d+\.\d+|\d+\.\d+E[-+]\d+)\s*$'),

      # Pattern for finding the "  2. Couplings" line (which marks the end of the section)
      'end': re.compile(r'^\s*2\.\s+Couplings\s*$')

    }

    # Parse the source file to get the information and build the states list

    for line in source_content:

      # Define when the section begins and ends

      if not section_found:
        if states_rx['start'].match(line):
          section_found = True
    
      elif section_found and states_rx['end'].match(line):
        break
        
      # Get the energy unit

      elif states_rx['unit_line'].match(line):
        energy_unit = states_rx['unit_line'].match(line).group("unit")
            
      # Extract the values of each state line

      elif states_rx['state_line'].match(line):

        state_number = int(states_rx['state_line'].match(line).group("number"))
        state_type = states_rx['state_line'].match(line).group("type")
        state_label = states_rx['state_line'].match(line).group("label")
        state_energy = float(states_rx['state_line'].match(line).group("energy"))

        # Convert the energy value to Ha

        state_energy = control_common.energy_unit_conversion(state_energy,energy_unit,'Ha')

        # Append information about the current state to the states_list variable

        system['states_list'].append({'number': state_number, 'type': state_type, 'label': state_label, 'energy': state_energy})

    # Raise an exception if the section has not been found

    if not section_found:
      raise control_common.ControlError ("ERROR: Unable to find the 'States list' section in the source file")

    # Raise an exception if the state numbers are not consecutive and starting at 0

    control_common.is_consecutive([state['number'] for state in system['states_list']],"Excited state numbers from the source file")

    print("[ DONE ]")

    # ========================================================= #
    #                       Couplings List                      #
    # ========================================================= #

    print("{:<50}".format("\nParsing the coupling values ...  "), end="")

    # Initialization of the variables

    section_found = False
    coupling_list = []
    
    # Define the expression patterns for the lines containing information about the couplings
    
    coupling_rx = {

      # Pattern for finding the "  2. Couplings" line (which marks the start of the section)
      'start': re.compile(r'^\s*2\.\s+Couplings\s*$'),

      # Pattern for finding lines looking like 'State 1     State 2     Value (cm-1)'
      'unit_line': re.compile(r'^\s*State\s+1\s+State\s+2\s+Value\s+\((?P<unit>[a-zA-Z_0-9\-]+)\)\s*$'),

      # Pattern for finding lines looking like '2           3          -4.669548'
      'coupling_line': re.compile(r'^\s*(?P<state_1>\d+)\s+(?P<state_2>\d+)\s+(?P<coupling>[-]?\d+|[-]?\d+\.\d+|[-]?\d+\.\d+E[-+]\d+)\s*$'),

      # Pattern for finding the "  3. Transition dipole moments" line (which marks the end of the section)
      'end': re.compile(r'^\s*3\.\s+Transition\s+dipole\s+moments\s*$')

    }

    # Parse the source file to get the information and build the SOC list

    for line in source_content:

      # Define when the section begins and ends

      if not section_found:
        if coupling_rx['start'].match(line):
          section_found = True
    
      elif section_found and coupling_rx['end'].match(line):
        break

      # Get the energy unit

      elif coupling_rx['unit_line'].match(line):
        coupling_unit = coupling_rx['unit_line'].match(line).group("unit")
            
      # Extract the values of each coupling line

      elif coupling_rx['coupling_line'].match(line):

        state_1 = int(coupling_rx['coupling_line'].match(line).group("state_1"))
        state_2 = int(coupling_rx['coupling_line'].match(line).group("state_2"))
        coupling_energy = float(coupling_rx['coupling_line'].match(line).group("coupling"))

        # Raise an exception if one of the state numbers has not been registered previously

        if state_1 not in [state['number'] for state in system['states_list']]:
          raise control_common.ControlError ("ERROR: One of the 'State 1' numbers in the 'Couplings' section of the source file does not appear in the 'States list' section")

        if state_2 not in [state['number'] for state in system['states_list']]:
          raise control_common.ControlError ("ERROR: One of the 'State 2' numbers in the 'Couplings' section of the source file does not appear in the 'States list' section")

        # Convert the coupling value to Ha

        coupling_energy = control_common.energy_unit_conversion(coupling_energy,coupling_unit,'Ha')

        # Append information about the current coupling to the coupling_list variable

        coupling_list.append((state_1,state_2,coupling_energy))

    # Raise an exception if the section has not been found

    if not section_found:
      raise control_common.ControlError ("ERROR: Unable to find the 'Couplings' section in the source file")

    print("[ DONE ]")

    # ========================================================= #
    #                       MIME Creation                       #
    # ========================================================= #

    print("{:<50}".format("\nBuilding the MIME ... "), end="")   

    # Initialize the MIME as a zero-filled matrix

    system['mime'] = np.zeros((len(system['states_list']), len(system['states_list'])))

    # Creation of the MIME - Non-diagonal values (couplings)

    for coupling in coupling_list:
      k1 = coupling[0]
      k2 = coupling[1]
      val = coupling[2]
      system['mime'][k1][k2] = val
      system['mime'][k2][k1] = val    # For symetry purposes

    # Creation of the MIME - Diagonal values (state energies)

    for state in system['states_list']:
      system['mime'][state['number']][state['number']] = state['energy']

    print("[ DONE ]")

    # ========================================================= #
    #                    Dipole Moments List                    #
    # ========================================================= #

    print("{:<50}".format("\nParsing the transition dipole moments ...  "), end="")

    # Initialization of the variables

    section_found = False
    momdip_list = []

    # Define the expression patterns for the lines containing information about the dipole moments

    moment_rx = {

      # Pattern for finding the "  3. Transition dipole moments" line (which marks the start of the section)
      'start': re.compile(r'^\s*3\.\s+Transition\s+dipole\s+moments\s*$'),

      # Pattern for finding lines looking like 'State 1     State 2     X value (D)   Y value (D)   Z value (D)'
      'unit_line': re.compile(r'^\s*State\s+1\s+State\s+2\s+X\s+value\s+\((?P<unit>[a-zA-Z_0-9\-]+)\)\s+Y\s+value\s+Z\s+value\s*$'),

      # Pattern for finding lines looking like '2           3          -4.669548'
      'momdip_line': re.compile(r'^\s*(?P<state_1>\d+)\s+(?P<state_2>\d+)\s+(?P<mom_x>[-]?\d+|[-]?\d+\.\d+|[-]?\d+\.\d+E[-+]\d+)\s+(?P<mom_y>[-]?\d+|[-]?\d+\.\d+|[-]?\d+\.\d+E[-+]\d+)\s+(?P<mom_z>[-]?\d+|[-]?\d+\.\d+|[-]?\d+\.\d+E[-+]\d+)\s*$')

    }

    # Parse the source file to get the information and build the dipole moments list

    for line in source_content:

      # Define when the section begins

      if not section_found:
        if moment_rx['start'].match(line):
          section_found = True

      # Get the unit

      elif moment_rx['unit_line'].match(line):

        moment_unit = moment_rx['unit_line'].match(line).group("unit").lower()

        # Check the unit

        supported_units = ['d','au']

        if moment_unit not in supported_units:
          raise control_common.ControlError ("ERROR: The unit of the transition dipole moment (%s) is currently not supported. Supported values include: %s" % (moment_unit,supported_units))

      # Extract the values of each moment line
  
      elif moment_rx['momdip_line'].match(line):
  
        state_1 = int(moment_rx['momdip_line'].match(line).group("state_1"))
        state_2 = int(moment_rx['momdip_line'].match(line).group("state_2"))
        momdip_x = float(moment_rx['momdip_line'].match(line).group("mom_x"))
        momdip_y = float(moment_rx['momdip_line'].match(line).group("mom_y"))
        momdip_z = float(moment_rx['momdip_line'].match(line).group("mom_z"))

        # Raise an exception if one of the state numbers has not been registered previously

        if state_1 not in [state['number'] for state in system['states_list']]:
          raise control_common.ControlError ("ERROR: One of the 'State 1' numbers in the 'Transition dipole moments' section of the source file does not appear in the 'States list' section")

        if state_2 not in [state['number'] for state in system['states_list']]:
          raise control_common.ControlError ("ERROR: One of the 'State 2' numbers in the 'Transition dipole moments' section of the source file does not appear in the 'States list' section")

        # If necessary, convert the moment values to atomic units (conversion factor from https://link.springer.com/content/pdf/bbm%3A978-3-319-89972-5%2F1.pdf)

        if moment_unit == 'd':
          conv_factor = 0.3934303
          momdip_x = momdip_x * conv_factor
          momdip_y = momdip_y * conv_factor
          momdip_z = momdip_z * conv_factor

        # Append information about the current moment to the momdip_list variable

        momdip_list.append((state_1,state_2,momdip_x,momdip_y,momdip_z))

    # Raise an exception if the section has not been found

    if not section_found:
      raise control_common.ControlError ("ERROR: Unable to find the 'Transition dipole moments' section in the source file")

    print("[ DONE ]")

    # ========================================================= #
    #                   Dipole Moments Matrices                 #
    # ========================================================= #

    print("{:<50}".format("\nBuilding transition dipole moments matrices ... "), end="")

    # Intialize the momdip_mtx dictionary

    system['momdip_mtx'] = {}

    # Initialize the matrices as zero-filled matrices

    system['momdip_mtx']['X'] = np.zeros((len(system['states_list']), len(system['states_list'])), dtype=float)
    system['momdip_mtx']['Y'] = np.zeros((len(system['states_list']), len(system['states_list'])), dtype=float)
    system['momdip_mtx']['Z'] = np.zeros((len(system['states_list']), len(system['states_list'])), dtype=float)

    # Filling the matrices

    for momdip in momdip_list:

      k1 = momdip[0]
      k2 = momdip[1]

      system['momdip_mtx']['X'][k1][k2] = momdip[2]
      system['momdip_mtx']['X'][k2][k1] = momdip[2]    # For symetry purposes

      system['momdip_mtx']['Y'][k1][k2] = momdip[3]
      system['momdip_mtx']['Y'][k2][k1] = momdip[3]    # For symetry purposes

      system['momdip_mtx']['Z'][k1][k2] = momdip[4]
      system['momdip_mtx']['Z'][k2][k1] = momdip[4]    # For symetry purposes

    print("[ DONE ]")

    # ========================================================= #
    #           Radiative lifetime of excited states            #
    # ========================================================= #

    print("{:<50}".format("\nComputing radiative lifetimes ... "), end="")

    # This calculation is based on the A_mn Einstein Coefficients and their link with the transition dipole moment
    # See https://aapt.scitation.org/doi/pdf/10.1119/1.12937 for reference
    # Note that this calculation is performed using atomic units, which means the Planck constant equals 2*pi and the vacuum permittivity equals 1/(4*pi)

    # Constants

    light_speed_au = constants.value('speed of light in vacuum') / constants.value('atomic unit of velocity')

    # Calculate the radiative lifetime of each excited state

    for state in system['states_list']:

      sum_einstein_coeffs = 0

      for other_state in system['states_list']:

        if other_state['energy'] < state['energy']:

          # Compute and convert the energy difference

          energy_diff = state['energy']- other_state['energy']

          # Compute the square of the transition dipole moment

          square_dipole = 0
          for momdip_key in system['momdip_mtx']:
            square_dipole += system['momdip_mtx'][momdip_key][state['number']][other_state['number']] ** 2

          # Calculate the A Einstein Coefficient   
         
          einstein_coeff = (4/3) * square_dipole * (energy_diff**3) / (light_speed_au**3)
          sum_einstein_coeffs += einstein_coeff

      if sum_einstein_coeffs == 0:
        state['lifetime'] = float('inf')
      else:
        state['lifetime'] = 1 / sum_einstein_coeffs

    print("[ DONE ]")

    # ========================================================= #
    #              Printing values in the log file              #
    # ========================================================= #

    # Printing the states list
    # ========================

    print("")
    print(''.center(70, '-'))
    print('States list'.center(70, ' '))
    print(''.center(70, '-'))
    print("{:<10} {:<15} {:<10} {:<15} {:<15}".format('Number','Type','Label','Energy (Ha)','Lifetime (a.u.)'))
    print(''.center(70, '-'))
    for state in system['states_list']:
      print("{:<10} {:<15} {:<10} {:<15.5f} {:<15.5e}".format(state['number'],state['type'],state['label'],state['energy'],state['lifetime']))
    print(''.center(70, '-'))

    # Printing the couplings list
    # ===========================

    print("")
    print(''.center(40, '-'))
    print('Couplings'.center(40, ' '))
    print(''.center(40, '-'))
    print("{:<10} {:<10} {:<20}".format('State 1','State 2','Value (Ha)'))
    print(''.center(40, '-'))
    for coupling in coupling_list:
      print("{:<10} {:<10} {:<.6g}".format(coupling[0],coupling[1],coupling[2]))
    print(''.center(40, '-'))    

    # Printing the momdip list
    # ========================

    print("")
    print(''.center(85, '-'))
    print('Transition dipole moments'.center(85, ' '))
    print(''.center(85, '-'))
    print("{:<10} {:<10} {:<20} {:<20} {:<20}".format('State 1','State 2','X value (a.u.)','Y value (a.u.)','Z value (a.u.)'))
    print(''.center(85, '-'))
    for momdip in momdip_list:
      print("{:<10} {:<10} {:<20} {:<20} {:<20}".format(momdip[0],momdip[1],momdip[2],momdip[3],momdip[4]))
    print(''.center(85, '-'))    

    # ========================================================= #
    #                    End of the function                    #
    # ========================================================= #

    # Remove dipole moment matrices filled with zeroes as they are useless

    system['momdip_mtx'] = { momdip_key : matrix for momdip_key, matrix in system['momdip_mtx'].items() if not np.all((matrix == 0)) }

    return system
