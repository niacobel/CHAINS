################################################################################################################################################
##                                                               Source File Parser                                                           ##
##                                                                                                                                            ##
##                                       This script contains the parsing functions for CONTROL LAUNCHER,                                     ##
##                                consult the documentation at https://chains-ulb.readthedocs.io/ for details                                 ##
################################################################################################################################################

import math
import re

import numpy

import control_errors


def energy_unit_conversion(value:float,init:str,target:str) -> float:
    """|  Converts an energy value from an initial unit to a target unit by using atomic units of energy (Hartree) as an intermediary.
    |  Currently supported units: Hartree, cm\ :sup:`-1`\ , eV, nm, Hz and Joules

    Parameters
    ----------
    value : float
        The energy value we need to convert.
    init : str
        The unit of the value we need to convert.
    target : str
        The unit we must convert the value to.
    
    Returns
    -------
    conv_value : float
        The converted energy value.
    """

    # Define the dictionary of conversion factors, from atomic units (Hartree) to any unit you want. - Taken from the NIST website (https://physics.nist.gov/)

    conv_factors = {
      # 1 Hartree equals:
      "Ha" : 1,
      "cm-1" : 2.1947463136320e+05,
      "eV" : 27.211386245988,
      "nm" : 2.1947463136320e+05 / 1e+07,
      "Hz" : 6.579683920502e+15,
      "J" : 4.3597447222071e-18
      }

    # Put everything in lower cases, to make it case insensitive

    init_low = init.lower()
    target_low = target.lower()
    conv_factors_low = dict((key.lower(), value) for key, value in conv_factors.items())

    # Check if the desired units are supported

    if init_low not in conv_factors_low.keys():
      raise control_errors.ControlError ("ERROR: The unit of the value you want to convert (%s) is currently not supported. Supported values include: %s" % (init, ', '.join(unit for unit in conv_factors.keys())))
    elif target_low not in conv_factors_low.keys():
      raise control_errors.ControlError ("ERROR: The unit you want to convert the value to (%s) is currently not supported. Supported values include: %s" % (target, ', '.join(unit for unit in conv_factors.keys())))
    
    # Convert the value

    conv_value = (value / conv_factors_low[init_low]) * conv_factors_low[target_low]

    return conv_value

# =================================================================== #
# =================================================================== #
#                    Q-Chem TD-DFT Parsing Function                   #
# =================================================================== #
# =================================================================== #

def qchem_tddft(source_content:list):
    """Parses the content of a Q-CHEM TD-DFT calculation output file, looking to build the MIME and the transition dipole moments matrix of the molecule. The ground state of the calculation must be a singlet, otherwise this function will need to be slightly modified.

    Parameters
    ----------
    source_content : list
        Content of the Q-CHEM output file. Each element of the list is a line of the file.
    
    Returns
    -------
    system : dict
        The extracted information of the source file. It contains three keys and their associated values: ``states_list``, ``mime`` and ``momdip_mtx`` where
        
        - ``states_list`` is a list of dictionaries containing four keys each: ``number``, ``multiplicity``, ``energy`` and ``label`` where

          - ``number`` is the number of the state, starting at 0 (which is the ground state).
          - ``multiplicity`` is the multiplicity of the state (Singlet, Triplet, ...).
          - ``energy`` is the excitation energy of the state, in cm\ :sup:`-1`\ .
          - ``label`` is the label of the state, in the form of the first letter of multiplicity + number of that state of this multiplicity (e.g. T1 for the first triplet, S3 for the third singlet, ...).

        - ``mime`` is a NumPy array, representing the Matrix Image of the MoleculE which acts as an effective Hamiltonian. It contains the excitation energies on the diagonal elements, and the spin-orbit couplings on the non-diagonal elements (in Hartree).
        - ``momdip_mtx`` is a NumPy array representing the transition dipole moments (in atomic units).

    Raises
    ------
    ControlError
        If some of the needed values are missing or unknown.

    """

    # Initialize the system dictionary that will be returned by the function

    system = {}

    # Define an array for correct English spelling during printing

    special_numbers = {1:"st", 2:"nd", 3:"rd"}

    # ========================================================= #
    #                      Number of roots                      #
    # ========================================================= #

    nb_roots = False

    # Define the expression patterns for the lines containing information about the number of roots, that number will be used to check if all the needed values have been collected.

    nb_roots_rx = {

      # Look for the "CIS_N_ROOTS          4" type of line
      'normal': re.compile(
          r'^\s*CIS_N_ROOTS\s*(?P<n_roots>\d+)\s*$'),

      # Look for the "NRoots was altered as:  4 --> 12" type of line
      'altered': re.compile(
          r'^\s*NRoots was altered as:\s*\d+\s*-->\s*(?P<new_n_roots>\d+)\s*$'),

      # No need to go further than the "TDDFT/TDA Excitation Energies" line
        'end': re.compile(
          r'TDDFT/TDA Excitation Energies')
    }

    # Parse the source file to get the information

    for line in source_content:
        
        # Compare each line to all the nb_roots_rx expressions

        for key, rx in nb_roots_rx.items():

          matching_line = rx.match(line)

          # If the line matches the rx pattern

          if matching_line is not None:

            # Get the info corresponding to the pattern

            if key == 'normal':
              nb_roots = int(matching_line.group('n_roots'))

            elif key == 'altered':
              nb_roots = int(matching_line.group('new_n_roots'))
              break          

            elif key == 'end':
              break         

    # Raise an exception if the number of roots has not been found

    if not nb_roots:
      raise control_errors.ControlError ("ERROR: Unable to find the number of roots in the source file")

    # ========================================================= #
    #                        States List                        #
    # ========================================================= #

    print("{:<50}".format("\nParsing the excited states ...  "), end="")

    # Initialization of the variables

    section_found = False
    search_energy = True

    cnt_triplet = 0
    cnt_singlet = 0

    # Define the ground state of our molecule (which is the first state and has a null energy)
    
    system['states_list'] = [{'number': 0, 'multiplicity': 'Singlet', 'energy': 0.0, 'label': 'S0'}]

    # Define the START and END expression patterns of the "TDDFT/TDA Excitation Energies" section of the output file

    section_rx = {
      'section_start': re.compile(
          r'TDDFT/TDA Excitation Energies'),
      'section_end': re.compile(
          r'SETman timing summary')
    }

    # Define the expression patterns for the lines containing information about the states

    states_rx = {

      # Pattern for finding lines looking like 'Excited state   1: excitation energy (eV) =    4.6445'
      'state_energy': re.compile(
          r'^\s*Excited state\s+(?P<state>\d+): excitation energy \(eV\) =\s+(?P<energy>[+]?\d*\.\d+|\d+)$'),

      # Pattern for finding lines looking like '    Multiplicity: Triplet'
      'state_mp': re.compile(
          r'^\s*Multiplicity: (?P<mplicity>\w+)$')
    }

    # Parse the source file to get the information and build the states list

    for line in source_content:

      # Define when the section begins and ends

      if not section_found:
        if section_rx['section_start'].match(line): # The line matches the section_start pattern
          section_found = True
          continue
    
      if section_found and section_rx['section_end'].match(line):
        break
        
      # Process the section lines

      else:

        # First, look for the excited energy of the state

        if search_energy:

          matching_line = states_rx['state_energy'].match(line)

          # If the line matches our pattern, get the state number and its associad excitation energy

          if matching_line is not None:

            exc_state = int(matching_line.group("state"))
            exc_energy = float(matching_line.group("energy"))
            search_energy = False

            continue
            
        # Second, look for the corresponding state multiplicity

        else:

          matching_line = states_rx['state_mp'].match(line)

          # If the line matches our pattern, get the mutiplicity and increase the counter (cnt) for that multiplicity (needed for the state label)

          if matching_line is not None:

            multiplicity = matching_line.group("mplicity")

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
              raise control_errors.ControlError ("ERROR: Multiplicity of the %s%s state is of unknown value (%s)" % (exc_state, ("th" if not exc_state in special_numbers else special_numbers[exc_state]),multiplicity))

            search_energy = True # Resume searching for the energy of the next state

            # Append information about the current state to the states_list variable

            system['states_list'].append({'number': exc_state, 'multiplicity': multiplicity, 'energy': energy_unit_conversion(exc_energy,"ev","cm-1"), 'label': (first_letter + str(cnt))})

            continue

    # Raise an exception if the section has not been found

    if not section_found:
      raise control_errors.ControlError ("ERROR: Unable to find the 'TDDFT/TDA Excitation Energies' section in the source file")

    # Raise an exception if not all the values have been found

    nb_states = (2 * nb_roots) + 1  # Ground state (1) + Triplet states (nb_roots) + Singlet states (nb_roots)
    if len(system['states_list']) != nb_states:
      raise control_errors.ControlError ("ERROR: The parsing function could not find the right number of excited states in the source file (%s of the %s expected states have been found)" % (len(system['states_list']),nb_states))

    print("[ DONE ]")

    # Print the states list

    print("")
    print(''.center(50, '-'))
    print('States List'.center(50, ' '))
    print(''.center(50, '-'))
    print("{:<10} {:<15} {:<15} {:<10}".format('Number','Multiplicity','Energy (cm-1)','Label'))
    print(''.center(50, '-'))
    for state in system['states_list']:
      print("{:<10} {:<15} {:<15.3f} {:<10}".format(state['number'],state['multiplicity'],state['energy'],state['label']))
    print(''.center(50, '-'))

    # ========================================================= #
    #                          SOC List                         #
    # ========================================================= #

    print("{:<50}".format("\nParsing the spin-orbit couplings ...  "), end="")

    # Initialization of the variables

    section_found = False
    soc_list = []
    
    # Define the START and END expression patterns of the "SPIN-ORBIT COUPLING" section of the output file

    section_rx = {
      'section_start': re.compile(
        r'^\*+SPIN-ORBIT COUPLING JOB BEGINS HERE\*+$'),
      'section_end': re.compile(
        r'^\*+SOC CODE ENDS HERE\*+$')
    }

    # Define the expression patterns for the lines containing information about the SOC
    
    soc_rx = {

      # Pattern for finding lines looking like 'Total SOC between the singlet ground state and excited triplet states:'
      'ground_to_triplets': re.compile(
        r'^\s*Total SOC between the singlet ground state and excited triplet states:$'),

      # Pattern for finding lines looking like 'Total SOC between the S1 state and excited triplet states:'
      'between_excited_states': re.compile(
        r'^\s*Total SOC between the (?P<s_key>[A-Z]\d+) state and excited triplet states:$'),

      # Pattern for finding lines looking like 'T2      76.018382    cm-1'
      'soc_value': re.compile(
        r'^\s*(?P<soc_key>[A-Z]\d+)\s+(?P<soc_value>\d+\.?\d*)\s+cm-1$')
    }

    # Parse the source file to get the information and build the SOC list

    for line in source_content:

      # Define when the section begins and ends

      if not section_found:
        if section_rx['section_start'].match(line): # The line matches the section_start pattern
          section_found = True
          continue
    
      if section_found and section_rx['section_end'].match(line):
        break

      # Process the section lines

      else:

        # Compare each line to all the soc_rx expressions
        
        for key, rx in soc_rx.items():

          matching_line = rx.match(line)

          # If the line matches the rx pattern

          if matching_line is not None:

            # Get the info corresponding to the pattern

            if key == 'ground_to_triplets':
              state_1 = 0

            elif key == 'between_excited_states':
              state_1 = matching_line.group('s_key')

              # state_1 is the state label, we need to convert it to the state number
              match = False
              for state in system['states_list']:
                if state_1 == state['label']:
                  state_1 = state['number']
                  match = True
                  break
              if not match:
                raise control_errors.ControlError ("ERROR: Unknown excited state (%s) has been catched during the SOC parsing." % state_1)

            elif key == 'soc_value':

              state_2 = matching_line.group('soc_key')

              # state_2 is the state label, we need to convert it to the state number
              match = False
              for state in system['states_list']:
                if state_2 == state['label']:
                  state_2 = state['number']
                  match = True
                  break
              if not match:
                raise control_errors.ControlError ("ERROR: Unknown excited state (%s) has been catched during the SOC parsing." % state_2)

              value = float(matching_line.group('soc_value'))

              soc_line = (state_1, state_2, value)
        
              # Add the new line to the soc_list
              soc_list.append(soc_line)

    # Raise an exception if the section has not been found

    if not section_found:
      raise control_errors.ControlError ("ERROR: Unable to find the 'SPIN-ORBIT COUPLING' section in the source file")

    # Raise an exception if not all the values have been found

    nb_soc = nb_roots/2 + (3/2 * (nb_roots**2)) # Ground to Triplet (nb_roots) + Triplet to Triplet (nb_roots*(nb_roots-1)/2) + Singlet to Triplet (nb_roots**2)
    if len(soc_list) != nb_soc:
      raise control_errors.ControlError ("ERROR: The parsing function could not find the right number of spin-orbit couplings in the source file (%s of the %s expected values have been found)" % (len(soc_list),nb_soc))

    print("[ DONE ]")

    # ========================================================= #
    #                       MIME Creation                       #
    # ========================================================= #

    print("{:<50}".format("\nBuilding the MIME ... "), end="")   

    # Initialize the MIME as a zero-filled matrix

    system['mime'] = numpy.zeros((len(system['states_list']), len(system['states_list'])))

    # Creation of the MIME - Non-diagonal values (SOC)

    for soc in soc_list:
      k1 = soc[0]
      k2 = soc[1]
      val = soc[2]
      system['mime'][k1][k2] = val
      system['mime'][k2][k1] = val    # For symetry purposes

    # Creation of the MIME - Diagonal values (Excitation energies)

    for state in system['states_list']:
      system['mime'][state['number']][state['number']] = state['energy']

    print("[ DONE ]")

    print("\nMIME (cm-1)")
    print('')
    for row in system['mime']:
      for val in row:
        print(numpy.format_float_positional(val,precision=3,unique=False,pad_left=8,trim="0",pad_right=3), end = " ")
      print('')

    # Convert the MIME to Hartree units

    system['mime'] = system['mime'] * energy_unit_conversion(1,"cm-1","ha")

    print("\nMIME (Ha)")
    print('')
    for row in system['mime']:
      for val in row:
        print(numpy.format_float_scientific(val,precision=5,unique=False,pad_left=2), end = " ")
      print('')    

    # ========================================================= #
    #                    Dipole Moments List                    #
    # ========================================================= #

    print("{:<50}".format("\nParsing the transition dipole moments ...  "), end="")

    # Initialization of the variables

    section_found = False
    momdip_list = []

    # Define the START and END expression patterns of the "STATE-TO-STATE TRANSITION MOMENTS" section of the output file

    section_rx = {
        'section_start': re.compile(
            r'^STATE-TO-STATE TRANSITION MOMENTS$'),
        'section_end': re.compile(
            r'^END OF TRANSITION MOMENT CALCULATION$')
    }

    # Define the expression patterns for the lines containing information about the dipole moments

    moment_rx = {
        # Pattern for finding lines looking like '    1    2   0.001414  -0.001456   0.004860   1.240659E-10'
        'moment': re.compile(
            r'^\s*(?P<mom_key1>\d+)\s+(?P<mom_key2>\d+)\s+(?P<mom_x>-?\d+\.\d+)\s+(?P<mom_y>-?\d+\.\d+)\s+(?P<mom_z>-?\d+\.\d+)\s+(?P<strength>\d|\d\.\d+|\d\.\d+E[-+]\d+)$')
    }

    # Parse the source file to get the information and build the dipole moments list

    for line in source_content:

      # Define when the section begins and ends

      if not section_found:
        if section_rx['section_start'].match(line): # The line matches the section_start pattern
          section_found = True
          continue
    
      if section_found and section_rx['section_end'].match(line):
        break

      # Process the section lines
  
      else:
  
        matching_line = moment_rx['moment'].match(line)

        # If the line matches our pattern, get the two states number and their associated strength
  
        if matching_line is not None:

          state_1 = matching_line.group('mom_key1')
          state_2 = matching_line.group('mom_key2')
          value = math.sqrt((float(matching_line.group('mom_x'))**2)+(float(matching_line.group('mom_y'))**2)+(float(matching_line.group('mom_z'))**2)) # Calculate the module of the transition dipole moment vector
          mom_line = (state_1, state_2, value)
        
          # Add the new line to the momdip_list
          momdip_list.append(mom_line)

    # Raise an exception if the section has not been found

    if not section_found:
      raise control_errors.ControlError ("ERROR: Unable to find the 'STATE-TO-STATE TRANSITION MOMENTS' section in the source file")

    # Raise an exception if not all the values have been found

    nb_momdip = nb_roots * (nb_roots+1) # Ground to Singlet (nb_roots) + Ground to Triplet (nb_roots) + Singlet to Singlet (nb_roots*(nb_roots-1)/2) + Triplet to Triplet (nb_roots*(nb_roots-1)/2)
    if len(momdip_list) != nb_momdip:
      raise control_errors.ControlError ("ERROR: The parsing function could not find the right number of state-to-state transition moments in the source file (%s of the %s expected values have been found)" % (len(momdip_list),nb_momdip))

    print("[ DONE ]")

    # ========================================================= #
    #                   Dipole Moments Matrix                   #
    # ========================================================= #

    # Initialize the matrix as a zero-filled matrix

    system['momdip_mtx'] = numpy.zeros((len(system['states_list']), len(system['states_list'])), dtype=float)

    # Filling the matrix

    for momdip in momdip_list:
      k1 = int(momdip[0])
      k2 = int(momdip[1])
      val = float(momdip[2])
      system['momdip_mtx'][k1][k2] = val
      system['momdip_mtx'][k2][k1] = val    # For symetry purposes

    print("\nDipole moments matrix (atomic units)")
    print('')
    for row in system['momdip_mtx']:
      for val in row:
        print(numpy.format_float_positional(val,precision=6,unique=False,pad_left=2,trim="0",pad_right=6), end = " ")
      print('')

    # End of the function, return the needed variables
    
    return system
