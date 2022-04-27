################################################################################################################################################
##                                                             Modelling functions                                                            ##
##                                                                                                                                            ##
##                                      This script contains the modelling functions for CONTROL LAUNCHER,                                    ##
##                                consult the documentation at https://chains-ulb.readthedocs.io/ for details                                 ##
################################################################################################################################################

import math
import re

import numpy as np
import scipy
from scipy import constants

import control_common

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

    # ========================================================= #
    #                      Number of roots                      #
    # ========================================================= #

    nb_roots = False

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
        nb_roots = int(nb_roots_rx['altered'].match(line).group('new_n_roots'))

    # Raise an exception if the number of roots has not been found

    if not nb_roots:
      raise control_common.ControlError ("ERROR: Unable to find the number of roots in the source file")

    print("{:<50} {:<10}".format("\nNumber of roots: ",nb_roots))

    # ========================================================= #
    #                        Zero states List                        #
    # ========================================================= #

    print("{:<50}".format("\nParsing the excited zero order states ...  "), end="")

    # Initialization of the variables

    section_found = False
    cnt_state = 0
    cnt_triplet = 0
    cnt_singlet = 0

    # Define the ground state of our molecule (which is the first state and has a zero energy)
    
    system['zero_states_list'] = [{'number': 0, 'type': '-', 'label': 'S0', 'energy' : 0.0}]

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
          state_type = "Dark"
          first_letter = "T"
          cnt_triplet += 1
          cnt = cnt_triplet

        elif multiplicity == "Singlet":
          state_type = "Bright"
          first_letter = "S"
          cnt_singlet += 1
          cnt = cnt_singlet

        else:
          raise control_common.ControlError ("ERROR: Multiplicity of the %s%s state is of unknown value (%s)" % (exc_state, ("th" if not exc_state in special_numbers else special_numbers[exc_state]),multiplicity))

        # Append information about the current state to the zero_states_list key

        system['zero_states_list'].append({'number': exc_state, 'type': state_type, 'label': (first_letter + str(cnt)), 'energy': control_common.energy_unit_conversion(exc_energy,"ev","cm-1")})

    # Raise an exception if the section has not been found

    if not section_found:
      raise control_common.ControlError ("ERROR: Unable to find the 'TDDFT/TDA Excitation Energies' section in the source file")

    # Raise an exception if not all the values have been found

    if cnt_state != 2*nb_roots:
      raise control_common.ControlError ("ERROR: The parsing function could not find the right number of excited states in the source file (%s of the %s expected states have been found)" % (cnt_state,2*nb_roots))
    if cnt_triplet != nb_roots:
      raise control_common.ControlError ("ERROR: The parsing function could not find the right number of excited triplet states in the source file (%s of the %s expected triplet states have been found)" % (cnt_triplet,nb_roots))
    if cnt_singlet != nb_roots:
      raise control_common.ControlError ("ERROR: The parsing function could not find the right number of excited singlet states in the source file (%s of the %s expected singlet states have been found)" % (cnt_singlet,nb_roots))

    # Raise an exception if the state numbers are not consecutive and starting at 0

    control_common.is_consecutive([state['number'] for state in system['zero_states_list']],"Excited state numbers from the QCHEM output file")

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

      # Pattern for finding lines looking like 'Total SOC between the singlet ground state and excited triplet states:'
      'ground_to_triplets': re.compile(r'^\s*Total SOC between the singlet ground state and excited triplet states:$'),

      # Pattern for finding lines looking like 'Total SOC between the S1 state and excited triplet states:'
      'between_excited_states': re.compile(r'^\s*Total SOC between the (?P<state_1>[A-Z]\d+) state and excited triplet states:$'),

      # Pattern for finding lines looking like 'T2      76.018382    cm-1'
      'soc_value': re.compile(r'^\s*(?P<state_2>[A-Z]\d+)\s+(?P<soc_value>\d+\.?\d*)\s+cm-1$'),

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

      elif soc_rx['between_excited_states'].match(line):

        label_1 = soc_rx['between_excited_states'].match(line).group('state_1')

        # Convert the label to the state number

        state_1 = False
        for state in system['zero_states_list']:
          if label_1 == state['label']:
            state_1 = state['number']
            break
        if not state_1:
          raise control_common.ControlError ("ERROR: Unknown excited state (%s) has been catched during the SOC parsing." % label_1)

      # Get the number and label of the second state and the corresponding SOC value, before adding the data to the soc_list

      elif soc_rx['soc_value'].match(line):

        label_2 = soc_rx['soc_value'].match(line).group('state_2')
        value = float(soc_rx['soc_value'].match(line).group('soc_value'))

        # Convert the label to the state number

        state_2 = False
        for state in system['zero_states_list']:
          if label_2 == state['label']:
            state_2 = state['number']
            break
        if not state_2:
          raise control_common.ControlError ("ERROR: Unknown excited state (%s) has been catched during the SOC parsing." % label_2)

        # Add the information to the soc_list

        soc_line = (state_1, label_1, state_2, label_2, value)
        soc_list.append(soc_line)

    # Raise an exception if the section has not been found

    if not section_found:
      raise control_common.ControlError ("ERROR: Unable to find the 'SPIN-ORBIT COUPLING' section in the source file")

    # Raise an exception if not all the values have been found

    nb_soc = nb_roots/2 + (3/2 * (nb_roots**2)) # Ground to Triplet (nb_roots) + Triplet to Triplet (nb_roots*(nb_roots-1)/2) + Singlet to Triplet (nb_roots**2)
    if len(soc_list) != nb_soc:
      raise control_common.ControlError ("ERROR: The parsing function could not find the right number of spin-orbit couplings in the source file (%s of the %s expected values have been found)" % (len(soc_list),nb_soc))

    print("[ DONE ]")

    # ========================================================= #
    #                       MIME Creation                       #
    # ========================================================= #

    print("{:<50}".format("\nBuilding the MIME ... "), end="")   

    # Initialize the MIME as a zero-filled matrix

    system['mime'] = np.zeros((len(system['zero_states_list']), len(system['zero_states_list'])))

    # Creation of the MIME - Non-diagonal values (SOC)

    for soc in soc_list:
      k1 = soc[0]
      k2 = soc[2]
      val = soc[4]
      system['mime'][k1][k2] = val
      system['mime'][k2][k1] = val    # For symetry purposes

    # Creation of the MIME - Diagonal values (Excitation energies)

    for state in system['zero_states_list']:
      system['mime'][state['number']][state['number']] = state['energy']

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

      k1 = int(momdip[0])
      k2 = int(momdip[1])

      system['momdip_o_mtx']['X'][k1][k2] = momdip[2]
      system['momdip_o_mtx']['X'][k2][k1] = momdip[2]    # For symetry purposes

      system['momdip_o_mtx']['Y'][k1][k2] = momdip[3]
      system['momdip_o_mtx']['Y'][k2][k1] = momdip[3]    # For symetry purposes

      system['momdip_o_mtx']['Z'][k1][k2] = momdip[4]
      system['momdip_o_mtx']['Z'][k2][k1] = momdip[4]    # For symetry purposes

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

    # ========================================================= #
    #              Printing values in the log file              #
    # ========================================================= #

    # Printing the states list
    # ========================

    print("")
    print(''.center(85, '-'))
    print('States List'.center(85, ' '))
    print(''.center(85, '-'))
    print("{:<10} {:<15} {:<10} {:<15} {:<15} {:<15}".format('Number','Type','Label','Energy (cm-1)','Energy (Ha)','Lifetime (s)'))
    print(''.center(85, '-'))
    for state in system['zero_states_list']:
      print("{:<10} {:<15} {:<10} {:<15.2f} {:<15.5f} {:<15.5e}".format(state['number'],state['type'],state['label'],state['energy'],control_common.energy_unit_conversion(state['energy'],"cm-1","ha"),state['lifetime']))
    print(''.center(85, '-'))

    # Printing the SOC list
    # =====================

    print("")
    print(''.center(40, '-'))
    print('Spin-orbit couplings'.center(40, ' '))
    print(''.center(40, '-'))
    print("{:<10} {:<10} {:<20}".format('State 1','State 2','Value (cm-1)'))
    print(''.center(40, '-'))
    for soc in soc_list:
      column_1 = soc[1] + " (" + str(soc[0]) + ")"
      column_2 = soc[3] + " (" + str(soc[2]) + ")"
      print("{:<10} {:<10} {:<.6g}".format(column_1,column_2,soc[4]))
    print(''.center(40, '-'))    

    # Printing the momdip list
    # ========================

    print("")
    print(''.center(105, '-'))
    print('Transition dipole moments'.center(105, ' '))
    print(''.center(105, '-'))
    print("{:<10} {:<10} {:<20} {:<20} {:<20} {:<20}".format('State 1','State 2','X value (a.u.)','Y value (a.u.)','Z value (a.u.)','Tot value (a.u.)'))
    print(''.center(105, '-'))
    for momdip in momdip_list:
      column_1 = system['zero_states_list'][momdip[0]]["label"] + " (" + str(momdip[0]) + ")"
      column_2 = system['zero_states_list'][momdip[1]]["label"] + " (" + str(momdip[1]) + ")"
      column_6 = math.sqrt(momdip[2]**2 + momdip[3]**2 + momdip[4]**2)
      print("{:<10} {:<10} {:<20} {:<20} {:<20} {:<.6g}".format(column_1,column_2,momdip[2],momdip[3],momdip[4],column_6))
    print(''.center(105, '-'))    

    # ========================================================= #
    # MIME diagonalization                                      #
    # ========================================================= #

    # Converting the states energy from cm-1 to Ha

    for state in system['zero_states_list']:
      state['energy'] = control_common.energy_unit_conversion(state['energy'],"cm-1","ha")

    # Diagonalization
    # ===============

    # Use SciPy to diagonalize the matrix (see https://personal.math.ubc.ca/~pwalls/math-python/linear-algebra/eigenvalues-eigenvectors/ for reference)
    # Documentation page for the function used here : https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.eigh.html  
     
    print("{:<50}".format("\nDiagonalizing the MIME ..."), end="")
    eigenvalues, system['eigenvectors'] = scipy.linalg.eigh(system['mime'])
    eigenvalues = eigenvalues.tolist()
    print("[ DONE ]")

    # Build the eigenstates list
    # ==========================

    print("{:<50}".format("\nBuilding the eigenstates list ..."), end="")

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

    print("[ DONE ]")

    # Transpose the eigenvectors list
    # ===============================

    # Using SciPy to invert the eigenvectors matrix
    # Note that the inverse of an orthonormal matrix is equal to its transpose, so each line of this matrix corresponds to an eigenvector. (see https://math.stackexchange.com/questions/156735/in-which-cases-is-the-inverse-of-a-matrix-equal-to-its-transpose)

    print("{:<50}".format("\nInverting eigenvectors matrix ..."), end="")
    system['eigenvectors_inv'] = scipy.linalg.inv(system['eigenvectors'])
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
    # Dipole moment matrices in the eigenstates basis set       #
    # ========================================================= #

    # Intialize the momdip_mtx dictionary

    system['momdip_mtx'] = {}

    # Convert each matrix from the zero order basis set to the eigenstates basis set through a matrix product (see https://numpy.org/doc/stable/reference/generated/numpy.matmul.html#numpy.matmul for reference)

    for momdip_key in system["momdip_o_mtx"]:
      system['momdip_mtx'][momdip_key] = np.matmul(np.matmul(system['eigenvectors_inv'],system['momdip_o_mtx'][momdip_key]),system['eigenvectors'])

    # ========================================================= #
    #                    End of the function                    #
    # ========================================================= #

    return system

# =================================================================== #
# =================================================================== #
#                     Custom File Parsing Function                    #
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
    print('States List'.center(70, ' '))
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
