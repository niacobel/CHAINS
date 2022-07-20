################################################################################################################################################
##                                                            Transition functions                                                            ##
##                                                                                                                                            ##
##                                     This script contains the transition functions for CONTROL LAUNCHER,                                    ##
##                                consult the documentation at https://chains-ulb.readthedocs.io/ for details                                 ##
################################################################################################################################################

import numpy as np
import math

import control_common

# =================================================================== #
# =================================================================== #
#                          Support functions                          #
# =================================================================== #
# =================================================================== #

def is_indistinguishable(energy1, energy2):
    """Check if two states can be considered indistinguishable by comparing their energy differences to the spectral resolution of the pulse shapers associated with their energy domain. The spectral resolutions are fetched from:
    - the specifications of the DAZZLER ultrafast pulse shapers provided by the FASTLITE company (https://fastlite.com/produits/dazzler-ultrafast-pulse-shaper/)
    - the Quickshaper IR provided by the PhaseTech company (http://phasetechspectroscopy.com/)
    - the compact optical parametric amplifier platform from Jakob et al. (https://doi.org/10.1364/oe.27.026979)

    Parameters
    ----------
    energy1 : float
        Energy of the first state in cm-1.
    energy2 : float
        Energy of the second state in cm-1.

    Returns
    -------
    bool
        "True" if the two states are indistinguishable, "False" otherwise.
    """  

    # Define the different energy domains and their associated spectral resolution

    resolutions = {
      40000: 15.99, # Dazzler Qz-250-400
      25000: 12.49, # Dazzler Qz-250-400
      20000:	8.00, # Dazzler WR-460-740
      14286:	6.12, # Dazzler WR-460-740
      11111:	6.17, # Dazzler WR-510-950
      9709:  	4.71, # Dazzler UHR-900-1700
      6452:	  4.16, # Dazzler UHR-900-1700
      4762:	  4.53, # Dazzler WB45-1450-3000
      3846:	  4.43, # Dazzler WB45-1450-3000
      2857:	  5.70, # Dazzler WB45-2000-3700
      1818:	  5.00, # Phasetech Quickshape IR
      926:    1.97  # Jakob et al.
    }

    # Define the domain of our values (see https://www.geeksforgeeks.org/python-find-closest-number-to-k-in-given-list/)

    domains = np.asarray(list(resolutions.keys()))

    domain1 = domains[(np.abs(domains - energy1)).argmin()]
    domain2 = domains[(np.abs(domains - energy2)).argmin()]

    # If the domains of the values differ, take the average resolution

    if domain1 != domain2:
      res = (resolutions[domain1] + resolutions[domain2]) / 2
    else:
      res = resolutions[domain1]

    # Check the energy difference of the states and compare it to the spectral resolution

    if math.isclose(energy1, energy2, abs_tol=res):
      return True
    else:
      return False

# =================================================================== #
# =================================================================== #
#                        Transition functions                         #
# =================================================================== #
# =================================================================== #

def brightest_to_darkest(system:dict):
    """Determines the content of the files needed by QOCT-RA for the transition going from the brightest state to the darkest state using the transition dipole moments matrix presenting the highest values on average. The brightest and darkest states are identified by their transition dipole moment with the ground state. This function ensures a stimulated emmission controle scheme so the chosen dark state will have a lower energy compared to the chosen bright state (other cases will be skipped).

    Parameters
    ----------
    system : dict
        Information extracted by the parsing function and derived from it.

    Returns
    -------
    transitions_list : list
        List of a single dictionary containing six keys: 

          - ``label`` is the label of the transition, which will be used for the name of the job directories.
          - ``init_file`` is the name of the initial state file, minus the number at the end.
          - ``init_content`` is the content of the initial state file.
          - ``target_file`` is the name of the target state file, minus the number at the end.
          - ``target_content`` is the content of the target state file.
          - ``momdip_key`` is the key of the transition dipole moments matrix used for this transition.
    """  

    # Initialize the list of dictionaries that will be returned by the function

    transitions_list = []

    """
    # ========================================================= #
    #              Identifying degenerated states               #
    # ========================================================= #

    # Initialize some variables

    deg_threshold = control_common.energy_unit_conversion(5,'cm-1','ha') # Degeneracy threshold (in Ha)
    deg_groups = []

    # Iterate over the states and identify groups of degenerated states

    for state in system['states_list']:

      degenerated = False

      for other_state in [other_state for other_state in system['states_list'] if other_state['number'] < state['number']]:

        if math.isclose(state['energy'],other_state['energy'],abs_tol=deg_threshold):

          degenerated = True

          # Define if and how to add this pair of degenerated states to the list of degeneracies

          added = False

          for group in deg_groups:

            if other_state['number'] in group and state['number'] not in group:
              group.append(state['number'])
              added = True

            elif state['number'] in group and other_state['number'] not in group:
              group.append(other_state['number'])
              added = True
              
            elif other_state['number'] in group and state['number'] in group:
              added = True

          if not added:
            deg_groups.append([state['number'],other_state['number']])
      
      if not degenerated:
        deg_groups.append([state['number']])

    # Pretty recap for the log file

    print("")
    print(''.center(70, '-'))
    print("{:<70}".format("Degenerated groups"))
    print(''.center(70, '-'))
    for idx, group in enumerate(deg_groups):
      print("{:<30} {:<40}".format("Group %s: " % str(idx), " - ".join(map(str,group)) if len(group) > 1 else str(group[0])))
    print(''.center(70, '-'))
    """

    # ========================================================= #
    #                   Defining orientation                    #
    # ========================================================= #

    # Determining the best transition dipole moments matrix

    highest_mean = -float('inf')
    for momdip_key in system['momdip_mtx']:
      mtx = np.absolute(system['momdip_mtx'][momdip_key])
      if mtx.mean() > highest_mean:
        highest_mean = mtx.mean()
        best_mtx = momdip_key

    # ========================================================= #
    #                Defining the pair of states                #
    # ========================================================= #

    # Sorting the bright states
    # =========================

    # Consider the dipole moments for transitions involving the ground state (and make sure it's a list and not a NumPy array)
    gs_line = system['momdip_mtx'][best_mtx][0]
    if isinstance(gs_line, np.ndarray):
      gs_line = gs_line.tolist()

    # Sort the indices of the states by decreasing order of the absolute value of their transition dipole moments
    abs_gs_line = [abs(mom) for mom in gs_line]
    bright_max_indices = [abs_gs_line.index(mom) for mom in sorted(abs_gs_line, reverse=True) if abs_gs_line.index(mom) != 0]

    # Sorting the dark states
    # =======================

    # Sort the indices of the states by decreasing order of their radiative lifetime (not taking into account the lowest energy level which has an infinite lifetime)
    # dark_max_indices = [system['states_list'].index(state) for state in sorted(system['states_list'], key = lambda i: i['lifetime'], reverse=True) if state['lifetime'] != float('inf')]

    # Sort the indices of the states by increasing order of the absolute value of their transition dipole moments 
    dark_max_indices = [abs_gs_line.index(mom) for mom in sorted(abs_gs_line) if abs_gs_line.index(mom) != 0]

    # Iterate over the states
    # =======================

    pair_found = False

    # Start iterating over the brightest states
    for b_index in bright_max_indices:

      if pair_found:
        break

      # Skip if the bright state is in fact not a bright state
      if system['states_list'][b_index]["sing_percent"] < 0.5:
        continue

      # Start iterating over the darkest states
      for d_index in dark_max_indices:
  
        # Skip if both states are indistinguishable
        b_energy = control_common.energy_unit_conversion(system['states_list'][b_index]["energy"], "ha", "cm-1")
        d_energy = control_common.energy_unit_conversion(system['states_list'][d_index]["energy"], "ha", "cm-1")
        if is_indistinguishable(b_energy, d_energy):
          continue

        # Skip if the dark state energy is higher than the one of the brightest state (stimulated emission scheme)
        if system['states_list'][b_index]["energy"] < system['states_list'][d_index]["energy"]:
          continue

        # Skip if the dark state is in fact not a dark state
        if system['states_list'][d_index]["trip_percent"] < 0.5:
          continue

        # Identify the pair

        bright_index = b_index
        bright_label = system['states_list'][bright_index]["label"]

        dark_index = d_index
        dark_label = system['states_list'][dark_index]["label"]

        pair_found = True

        break

    # ========================================================= #
    #               Defining the density matrices               #
    # ========================================================= #

    # Initial states
    # ~~~~~~~~~~~~~~

    # Initialize the density matrix

    init_content = np.zeros((len(system['states_list']), len(system['states_list'])),dtype=complex)

    # Define the initial density matrix file name (will be modified if there are multiple states)

    init_file = bright_label

    # Check if there are multiple states in the initial state

    bright_group = []
    bright_energy = control_common.energy_unit_conversion(system['states_list'][bright_index]["energy"], "ha", "cm-1")

    for state in system['states_list']:
      state_energy = control_common.energy_unit_conversion(state["energy"], "ha", "cm-1")
      if is_indistinguishable(state_energy, bright_energy):
        bright_group.append(state['number'])

    # The starting populations are proportional to the transition dipole moments between each state of the "bright group" and the ground state
    # If there is only one state, its starting population will be 1

    total_momdip = sum([system['momdip_mtx'][best_mtx][0][state_number] for state_number in bright_group])
    for state_number in bright_group:
      if state_number != bright_index:
        init_file += "-" + str(state_number) # Add the number of this state to the init filename
      init_content[state_number][state_number] = complex(system['momdip_mtx'][best_mtx][0][state_number] / total_momdip)

    # Target states
    # ~~~~~~~~~~~~~

    target_content = np.zeros((len(system['states_list']), len(system['states_list'])),dtype=complex)
    target_content[dark_index][dark_index] = complex(1)

    # ========================================================= #
    #                  Defining the projector                   #
    # ========================================================= #

    # Initialize the projector matrix

    projector_content = np.zeros((len(system['states_list']), len(system['states_list'])),dtype=complex)

    # Define the target density matrix file name (will be modified if there are multiple states)

    target_file = dark_label

    # Check if there are multiple states in the target state

    dark_group = []
    dark_energy = control_common.energy_unit_conversion(system['states_list'][dark_index]["energy"], "ha", "cm-1")

    for state in system['states_list']:
      state_energy = control_common.energy_unit_conversion(state["energy"], "ha", "cm-1")
      if state_energy >= bright_energy:
        continue
      if is_indistinguishable(bright_energy - dark_energy, bright_energy - state_energy):
        dark_group.append(state['number'])

    # We need to 'open the channel' (put a 1) for each state belonging to the "dark group"

    for state_number in dark_group:
      if state_number != dark_index:
          target_file += "-" + str(state_number) # Add the number of this state to the target filename
      projector_content[state_number][state_number] = complex(1)

    # ========================================================= #
    #           Building the transition dictionary              #
    # ========================================================= #

    # Building the transition dictionary

    transition_label = best_mtx + "_" + init_file + "_" + target_file

    transition = {
      "label" : transition_label,
      "init_file" : init_file + "_",
      "init_content" : init_content,
      "target_file" : target_file + "_",
      "target_content" : target_content,
      "momdip_key" : best_mtx,
      "projector" : projector_content
      }

    # Pretty recap for the log file

    print("")
    print(''.center(70, '-'))
    print("{:<30} {:<40}".format("Label: ", transition["label"]))
    print(''.center(70, '-'))
    print(''.center(70, '-'))
    print("{:<30} {:<40}".format("Transition dipole matrix: ", best_mtx))
    print(''.center(70, '-'))
    print("{:<30} {:<40}".format("Initial state: ", "%s" % init_file))
    print(''.center(70, '-'))
    print("{:<30} {:<40}".format("Target state: ", "%s" % target_file))
    print(''.center(70, '-'))

    # Add the transition to the transitions list

    transitions_list.append(transition)

    return transitions_list

#######################################################################

def gs_or_brightests_to_darkests_and_reverse(system:dict):
    """Determines the content of the files needed by QOCT-RA for the transitions going from either the ground state or one of the first B brightest states to one of the first D darkest states, and vice versa. The brightest states are identified by their transition dipole moment with the ground state while the darkest states are identified by their radiative lifetime. B and D are fixed parameters that can easily be changed at the beginning of the function definition.

    Parameters
    ----------
    system : dict
        Information extracted by the parsing function and derived from it.

    Returns
    -------
    transitions_list : list
        List of dictionaries containing six keys each: 

          - ``label`` is the label of the transition, which will be used for the name of the job directories.
          - ``init_file`` is the name of the initial state file, minus the number at the end.
          - ``init_content`` is the content of the initial state file.
          - ``target_file`` is the name of the target state file, minus the number at the end.
          - ``target_content`` is the content of the target state file.
          - ``momdip_key`` is the key of the transition dipole moments matrix used for this transition.
    """  

    # Initialize the list of dictionaries that will be returned by the function

    transitions_list = []

    # Set the number of bright and dark states to consider
    #! Keep in mind that you will end up with a total number of transitions corresponding to bmax+1 * dmax * the number of transition dipole moments matrices!

    bmax = 2 # Maximum number of bright states to consider (e.g. if bmax = 2, then only the two brightest states will be considered)
    dmax = 2 # Maximum number of dark states to consider (e.g. if dmax = 2, then only the two darkest states will be considered)

    # ========================================================= #
    #                   Defining transitions                    #
    # ========================================================= #

    # We need to consider each transition dipole moments matrix separately
    for momdip_key in system['momdip_mtx']:

      # ========================================================= #
      #                Defining the pair of states                #
      # ========================================================= #

      # Sorting the bright states
      # =========================

      # Consider the dipole moments for transitions involving the ground state (and make sure it's a list and not a NumPy array)
      gs_line = system['momdip_mtx'][momdip_key][0]
      if isinstance(gs_line, np.ndarray):
        gs_line = gs_line.tolist()

      # Sort the indices of the states by decreasing order of the absolute value of their transition dipole moments
      # e.g. if gs_line = [0.0, 0.0, 0.0, 0.0, 0.0, -1.2, 0.8, 1.0, -0.4] then bright_max_indices = [5, 7, 6, 8]
      abs_gs_line = [abs(mom) for mom in gs_line]
      bright_max_indices = [abs_gs_line.index(mom) for mom in sorted(abs_gs_line, reverse=True)]

      # Add the ground state
      # ====================

      # Add the ground state (index 0) to the beginning of the list of bright states
      bright_max_indices.insert(0,0)

      # Sorting the dark states
      # =======================

      # Sort the indices of the states by decreasing order of their radiative lifetime (not taking into account the lowest energy level which has an infinite lifetime)
      dark_max_indices = [system['states_list'].index(state) for state in sorted(system['states_list'], key = lambda i: i['lifetime'], reverse=True) if state['lifetime'] != float('inf')]

      # Iterate over the dark states
      # ============================

      # Start iterating over the first Nth darkest states, where N is the dmax argument (see https://stackoverflow.com/questions/36106712/how-can-i-limit-iterations-of-a-loop-in-python for details)
      for iter_dark, dark_index in zip(range(dmax), dark_max_indices):

      # iter_dark is the number of the current iteration of this loop, e.g if dmax = 2, then iter_dark will be 0, then 1. (useful to label the transition)
      # dark_index is the number of the state currently considered, e.g. if dark_max_indices = [4, 3, 2, 1] and iter_bright = 0, then dark_index = 4.

        # Define the target density matrix file name and content

        dark_label = system['states_list'][dark_index]["label"]

        target_file = dark_label + "_"

        target_content = np.zeros((len(system['states_list']), len(system['states_list'])),dtype=complex)  # Quick init of a zero-filled matrix
        target_content[dark_index][dark_index] = complex(1)

        # Iterate over the bright states
        # ==============================

        # Start iterating over the first Nth brightest states, where N is the bmax argument* (see https://stackoverflow.com/questions/36106712/how-can-i-limit-iterations-of-a-loop-in-python for details)
        # *Since the first "bright" state is the ground state, we need to increase bmax by 1.
        for iter_bright, bright_index in zip(range(bmax+1), bright_max_indices):

        # iter_bright is the number of the current iteration of this loop, e.g if bmax = 3, then iter_bright will be 0, then 1, then 2. (useful to label the transition)
        # bright_index is the number of the state currently considered, e.g. if bright_max_indices = [0, 5, 7, 6, 8] and iter_bright = 1, then bright_index = 5.

          # Exit if both states are the same (that transition will be skipped)
          if bright_index == dark_index:
            state_label = system['states_list'][bright_index]["label"]
            transition_label = momdip_key + "_" + str(iter_bright) + "B" + str(iter_dark+1) + "D"
            print("\nNOTICE: The transition %s will be skipped since both states are the same (%s)." % (transition_label,state_label))
            continue

          # Define the initial density matrix file name and content

          bright_label = system['states_list'][bright_index]["label"]

          init_file = bright_label + "_"

          init_content = np.zeros((len(system['states_list']), len(system['states_list'])),dtype=complex)  # Quick init of a zero-filled matrix
          init_content[bright_index][bright_index] = complex(1)

          # ========================================================= #
          #           Building the transition dictionary              #
          # ========================================================= #

          # Building the transition dictionary

          if bright_index == 0:
            transition_label = momdip_key + "_GS" + str(iter_dark+1) + "D_" + bright_label + "-" + dark_label
          else:
            transition_label = momdip_key + "_" + str(iter_bright) + "B" + str(iter_dark+1) + "D_" + bright_label + "-" + dark_label

          transition = {
            "label" : transition_label,
            "init_file" : init_file,
            "init_content" : init_content,
            "target_file" : target_file,
            "target_content" : target_content,
            "momdip_key" : momdip_key
            }

          # Pretty recap for the log file

          print("")
          print(''.center(70, '-'))
          print("{:<30} {:<40}".format("Label: ", transition["label"]))
          print(''.center(70, '-'))
          print(''.center(70, '-'))
          print("{:<30} {:<40}".format("Transition dipole matrix: ", momdip_key))
          print(''.center(70, '-'))
          print("{:<30} {:<40}".format("Initial state: ", "%s (%s)" % (bright_label,bright_index)))
          print(''.center(70, '-'))
          print("{:<30} {:<40}".format("Target state: ", "%s (%s)" % (dark_label,dark_index)))
          print(''.center(70, '-'))

          # Add the transition to the transitions list

          transitions_list.append(transition)

          # ========================================================= #
          #       Building the reverse transition dictionary          #
          # ========================================================= #

          # Building the reverse transition dictionary

          if bright_index == 0:
            transition_label = "R_" + momdip_key + "_" + str(iter_dark+1) + "D" + "GS_" + dark_label + "-" +  bright_label
          else:
            transition_label = "R_" + momdip_key + "_" + str(iter_dark+1) + "D" + str(iter_bright) + "B_" + dark_label + "-" +  bright_label

          transition = {
            "label" : transition_label,
            "init_file" : target_file,
            "init_content" : target_content,
            "target_file" : init_file,
            "target_content" : init_content,
            "momdip_key" : momdip_key
            }

          # Pretty recap for the log file

          print("")
          print(''.center(70, '-'))
          print("{:<30} {:<40}".format("Label: ", transition["label"]))
          print(''.center(70, '-'))
          print(''.center(70, '-'))
          print("{:<30} {:<40}".format("Transition dipole matrix: ", momdip_key))
          print(''.center(70, '-'))
          print("{:<30} {:<40}".format("Initial state: ", "%s (%s)" % (dark_label,dark_index)))
          print(''.center(70, '-'))
          print("{:<30} {:<40}".format("Target state: ", "%s (%s)" % (bright_label,bright_index)))
          print(''.center(70, '-'))

          # Add the reverse transition to the transitions list, before proceeding with the next darkest state (or the next brightest state, if this was the last one)

          transitions_list.append(transition)

    return transitions_list

#######################################################################

def dark_zero_order(system:dict):
    """Determines the transition files needed by QOCT-RA for the transition between the ground state and each dark zero order state, expressed as a linear combination of eigenstates. Note: the handling of degenerated states must be turned off during execution of CONTROL LAUNCHER (the -dt / --degen_tresh command line argument must have a negative value).

    Parameters
    ----------
    system : dict
        Information extracted by the parsing function and derived from it.

    Returns
    -------
    transitions_list : list
        List of dictionaries containing six keys each: 

          - ``label`` is the label of the transition, which will be used for the name of the job directories.
          - ``init_file`` is the name of the initial state file, minus the number at the end.
          - ``init_content`` is the content of the initial state file.
          - ``target_file`` is the name of the target state file, minus the number at the end.
          - ``target_content`` is the content of the target state file.
          - ``momdip_key`` is the key of the transition dipole moments matrix used for this transition.
    """  

    # Initialize the list of dictionaries that will be returned by the function

    transitions_list = []

    # Define the initial density matrix file name and content

    gs_state = [state for state in system['zero_states_list'] if state['number'] == 0][0]

    init_file = gs_state['label'] + "_"

    init_content = np.zeros((len(system['zero_states_list']), len(system['zero_states_list'])),dtype=complex)  # Quick init of a zero-filled matrix
    init_content[gs_state['number']][gs_state['number']] = complex(1)

    # Iterate over the dark zero order states

    for state in system['zero_states_list']:

      if state['type'].lower() == "dark":

        # Build the target density matrix in the zero order states

        target_mtx = np.zeros((len(system['zero_states_list']), len(system['zero_states_list'])),dtype=complex)  # Quick init of a zero-filled matrix
        target_mtx[state['number']][state['number']] = complex(1)

        # Convert target density matrix to the eigenstates basis set and define the target density matrix file name 

        target_content = np.matmul(np.matmul(system['eigenvectors_inv'],target_mtx),system['eigenvectors'])

        target_file = state['label'] + "_"

        # Building the transition dictionary (one for each transition dipole moment matrix)

        for momdip_key in system['momdip_o_mtx']:

            # Building the transition dictionary

            transition = {
              "label" : momdip_key + "_" + gs_state['label'] + "-" + state['label'],
              "init_file" : init_file,
              "init_content" : init_content,
              "target_file" : target_file,
              "target_content" : target_content,
              "momdip_key" : momdip_key
              }

            # Pretty recap for the log file

            print("")
            print(''.center(70, '-'))
            print("{:<30} {:<40}".format("Label: ", transition["label"]))
            print(''.center(70, '-'))
            print(''.center(70, '-'))
            print("{:<30} {:<40}".format("Transition dipole matrix: ", momdip_key))
            print(''.center(70, '-'))
            print("Initial density matrix (%s) - Real parts only" % gs_state['label'])
            print(''.center(70, '-'))
            for line in init_content:
              for val in line:
                print(np.format_float_scientific(val.real,precision=3,unique=False,pad_left=2), end = " ")
              print("")
            print(''.center(70, '-'))
            print("Target density matrix (%s) - Real parts only" % state['label'])
            print(''.center(70, '-'))
            for line in target_content:
              for val in line:
                print(np.format_float_scientific(val.real,precision=3,unique=False,pad_left=2), end = " ")
              print("")
            print(''.center(70, '-'))

            transitions_list.append(transition)

    return transitions_list

