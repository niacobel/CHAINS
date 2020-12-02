################################################################################################################################################
##                                                               Source File Parser                                                           ##
##                                                                                                                                            ##
##                                       This script contains the parsing functions for CONTROL LAUNCHER,                                     ##
##                                consult the documentation at https://chains-ulb.readthedocs.io/ for details                                 ##
################################################################################################################################################

import re
import control_errors

def ev_to_cm1(ev: float) -> float:
    """Converts a value from eV to cm-1.

    Parameters
    ----------
    ev : float
        The value in eV we need to convert
    
    Returns
    -------
    cm : float
        The eV value converted to cm-1
    """

    cm = ev * 8065.6

    return cm

# =================================================================== #
# =================================================================== #
#                    Q-Chem TD-DFT Parsing Function                   #
# =================================================================== #
# =================================================================== #

def qchem_tddft_parser(file_content:list):
    """Parses the content of a Q-Chem TD-DFT calculation output file, looking to extract the list of the excited electronic states, the spin-orbit couplings and the transition dipole moments between the electronic states of the molecule. The ground state of the calculation must be a singlet, otherwise this function will need to be slightly modified.

    Parameters
    ----------
    file_content : list
        Content of the qchem output file
    
    Returns
    -------
    states_list : list
        A list of tuples of the form [[0, Multiplicity, Energy, Label], [1, Multiplicity, Energy, Label], [2, Multiplicity, Energy, Label], ...]
        The first element of each tuple is the state number, starting at 0
        Multiplicity is the multiplicity of the state (ex : singlet, triplet)
        Energy is the energy of the state, in cm-1
        Label is the label of the state, in the form of first letter of multiplicity + number of that state of this multiplicity (ex : T1 for the first triplet, S3 for the third singlet)
    coupling_list : list
        List of tuples of the form [[State0, State1, SOC_0-1], [State0, State2, SOC_0-2], [State1, State2, SOC_1-2], ...]
        The first two elements of each tuple are the number of the two states and the third one is the value of the spin-orbit coupling between them (in cm-1)
    momdip_list : list
        List of tuples of the form [[State0, State1, MomDip0-1], [State0, State2, MomDip0-2], [State1, State2, MomDip1-2], ...]
        The first two elements of each tuple are the number of the two states and the third one is the value of the transition dipole moment associated with the transition between them (in atomic units)

    Raises
    ------
    ControlError
        If one of the informations that need to be extracted presents some kind of problem.

    Notes
    -----
    This function makes heavy use of regexes (Regular Expressions). If you're unfamiliar with them, please consult the "What You Need To Know" section of the documentation for information.
    """

    # Define an array for correct English spelling during printing

    special_numbers = {1:"st", 2:"nd", 3:"rd"}

    # ========================================================= #
    #                        States List                        #
    # ========================================================= #

    # Initialization of the variables

    section_found = False
    search_energy = True

    cpt_triplet = 0
    cpt_singlet = 0

    # Define the first state of our molecule - the ground state (wich is the first state and has a null energy)
    
    states_list = [(0, 'Singlet', 0.0, 'S0')]

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
          r'^Excited state\s+(?P<state>\d+): excitation energy \(eV\) =\s+(?P<energy>[-+]?\d*\.\d+|\d+)$'),

      # Pattern for finding lines looking like '    Multiplicity: Triplet'
      'state_mp': re.compile(
          r'^\s*Multiplicity: (?P<mplicity>\w+)$')
    }

    # Parse the source file to get the information and build the states list

    for line in file_content:

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
            if exc_state <= 0:
              raise control_errors.ControlError ("ERROR: The number of one of the found excited states (%s) is not a non zero positive integer." % exc_state)

            exc_energy = float(matching_line.group("energy"))
            if exc_energy < 0:
              raise control_errors.ControlError ("ERROR: The excitation energy of the %s%s excited state is negative (%s)" % (exc_state, ("th" if not exc_state in special_numbers else special_numbers[exc_state]),exc_energy))

            search_energy = False

            continue
            
        # Second, look for the corresponding state multiplicity

        else:

          matching_line = states_rx['state_mp'].match(line)

          # If the line matches our pattern, get the mutiplicity and increase the counter for that multiplicity (needed for the state label)

          if matching_line is not None:

            multiplicity = matching_line.group("mplicity")

            cpt = -1

            if multiplicity == "Triplet":
              first_letter = "T"
              cpt_triplet += 1
              cpt = cpt_triplet

            elif multiplicity == "Singlet":
              first_letter = "S"
              cpt_singlet += 1
              cpt = cpt_singlet

            else:
              raise control_errors.ControlError ("ERROR: Multiplicity of the %s%s state is of unknown value (%s)" % (exc_state, ("th" if not exc_state in special_numbers else special_numbers[exc_state]),multiplicity))

            search_energy = True # Resume searching for the energy of the next state

            # Append information about the current state to the states_list variable with the format: (state_number, state_multiplicity, energy value (cm-1), label)
            states_list.append((exc_state, multiplicity, ev_to_cm1(exc_energy), (first_letter + str(cpt))))

            continue

    # ========================================================= #
    #                          SOC List                         #
    # ========================================================= #

    # Initialization of the variables

    section_found = False
    coupling_list = []
    tpl = None   # Tuple that will look like "(state_1, state_2, value)"
    
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
        r'^Total SOC between the singlet ground state and excited triplet states:$'),

      # Pattern for finding lines looking like 'Total SOC between the S1 state and excited triplet states:'
      'between_excited_states': re.compile(
        r'^Total SOC between the (?P<s_key>[A-Z]\d) state and excited triplet states:$'),

      # Pattern for finding lines looking like 'T2      76.018382    cm-1'
      'soc_value': re.compile(
        r'^(?P<soc_key>T\d)\s+(?P<soc_value>\d+\.?\d*)\s+cm-1$')
    }

    # Parse the source file to get the information and build the SOC list

    for line in file_content:

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
              state_1 = 'S0' # The ground state label is S0 by convention

            elif key == 'between_excited_states':
              state_1 = matching_line.group('s_key')
              if state_1 not in [x[3] for x in states_list]:
                raise control_errors.ControlError ("ERROR: Unknown excited state (%s) has been catched during parsing." % state_1)

            elif key == 'soc_value':

              state_2 = matching_line.group('soc_key')
              if state_2 not in [x[3] for x in states_list]:
                raise control_errors.ControlError ("ERROR: Unknown excited state (%s) has been catched during parsing." % state_2)

              value = float(matching_line.group('soc_value'))
              if value < 0:
                raise control_errors.ControlError ("ERROR: The SOC value between the %s and %s states is negative (%s)" % (state_1,state_2,value))

              tpl = (state_1, state_2, value)
       
        # When all the relevant information has been found, append a new line to the coupling_list

        if tpl is not None:
          coupling_list.append(tpl)
          tpl = None
          continue

    # Rewrite coupling_list by replacing all state labels by their state number (needs the states_list that was determined in the previous section)
  
    for index, soc in enumerate(coupling_list):
      coupling_list[index] = (
        # Filter the states_list, looking for the state corresponding to the label (see https://stackoverflow.com/questions/3013449/list-comprehension-vs-lambda-filter for reference)
        int([ x for x in states_list if x[3] == soc[0] ][0][0]),
        int([ x for x in states_list if x[3] == soc[1] ][0][0]),
        soc[2]
      )

    # ========================================================= #
    #                    Dipole Moments List                    #
    # ========================================================= #

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
            r'^(?P<mom_key1>\d+)\s+(?P<mom_key2>\d+)(\s+-?\d\.\d+){3}\s+(?P<strength>(\d|\d\.\d+|\d\.\d+E[-+]\d{2}))$')
    }

    # Parse the source file to get the information and build the dipole moments list

    for line in file_content:

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
          value = float(matching_line.group('strength'))
          line = (state_1, state_2, value)
        
          # Add the new line to the momdip_list
          momdip_list.append(line)

    
    return states_list,coupling_list,momdip_list