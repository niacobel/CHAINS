################################################################################################################################################
##                                                            Transition functions                                                            ##
##                                                                                                                                            ##
##                                     This script contains the transition functions for CONTROL LAUNCHER,                                    ##
##                                consult the documentation at https://chains-ulb.readthedocs.io/ for details                                 ##
################################################################################################################################################

import numpy as np
import scipy

import control_common

# =================================================================== #
# =================================================================== #
#                        Transition functions                         #
# =================================================================== #
# =================================================================== #

def brightests_to_darkests_and_reverse(system:dict):
    """Determines the transition files needed by QOCT-RA for the transition between each of the first B brightest states and each of the first D darkest states, and vice versa. The brightest states are identified by their transition dipole moment with the ground state while the darkest states are identified by their radiative lifetime. B and D are fixed parameters that can easily be changed at the beginning of the function definition.

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
    #! Keep in mind that you will end up with a total number of transitions corresponding to bmax * dmax * the number of transition dipole moments matrices!

    bmax = 1 # Maximum number of bright states to consider (e.g. if bmax = 2, then only the two brightest states will be considered)
    dmax = 1 # Maximum number of dark states to consider (e.g. if dmax = 2, then only the two darkest states will be considered)

    # ========================================================= #
    #                   Defining transitions                    #
    # ========================================================= #

    # We need to consider each transition dipole moments matrix separately
    for momdip_key in system['momdip_es_mtx']:

      # ========================================================= #
      #          Defining the bright-dark pair of states          #
      # ========================================================= #

      # Sorting the bright states
      # =========================

      # Consider the dipole moments for transitions involving the ground state (and make sure it's a list and not a NumPy array)
      gs_line = system['momdip_es_mtx'][momdip_key][0]
      if isinstance(gs_line, np.ndarray):
        gs_line = gs_line.tolist()

      # Sort the indices of the states by decreasing order of the absolute value of their transition dipole moments
      # e.g. if gs_line = [0.0, 0.0, 0.0, 0.0, 0.0, -1.2, 0.8, 1.0, -0.4] then bright_max_indices = [5, 7, 6, 8]
      abs_gs_line = [abs(mom) for mom in gs_line]
      bright_max_indices = [abs_gs_line.index(mom) for mom in sorted(abs_gs_line, reverse=True)]

      # Sorting the dark states
      # =======================

      # Sort the indices of the states by decreasing order of their radiative lifetime (not taking into account the lowest energy level which has an infinite lifetime)
      dark_max_indices = [system['eigenstates_list'].index(eigenstate) for eigenstate in sorted(system['eigenstates_list'], key = lambda i: i['lifetime'], reverse=True) if eigenstate['lifetime'] != float('inf')]

      # Iterate over the bright states
      # ==============================

      # Start iterating over the first Nth brightest states, where N is the bmax argument (see https://stackoverflow.com/questions/36106712/how-can-i-limit-iterations-of-a-loop-in-python for details)
      for iter_bright, bright_index in zip(range(bmax), bright_max_indices):

      # iter_bright is the number of the current iteration of this loop, e.g if bmax = 3, then iter_bright will be 0, then 1, then 2. (useful to label the transition)
      # bright_index is the number of the state currently considered, e.g. if bright_max_indices = [5, 7, 6, 8] and iter_bright = 1, then bright_index = 7.

        # Define the initial density matrix file name and content

        bright_label = system['eigenstates_list'][bright_index]["label"]

        init_file = bright_label + "_"

        init_content = np.zeros((len(system['eigenstates_list']), len(system['eigenstates_list'])),dtype=complex)  # Quick init of a zero-filled matrix
        init_content[bright_index][bright_index] = complex(1)

        # Iterate over the dark states
        # ============================

        # Start iterating over the first Nth darkest states, where N is the dmax argument (see https://stackoverflow.com/questions/36106712/how-can-i-limit-iterations-of-a-loop-in-python for details)
        for iter_dark, dark_index in zip(range(dmax), dark_max_indices):

        # iter_dark is the number of the current iteration of this loop, e.g if dmax = 2, then iter_dark will be 0, then 1. (useful to label the transition)
        # dark_index is the number of the state currently considered, e.g. if dark_max_indices = [4, 3, 2, 1] and iter_bright = 0, then dark_index = 4.

          # Exit if both states are the same (that transition will be skipped)
          if bright_index == dark_index:
            state_label = system['eigenstates_list'][bright_index]["label"]
            transition_label = momdip_key + "_" + str(iter_bright+1) + "B" + str(iter_dark+1) + "D"
            print("\nNOTICE: The transition %s will be skipped since both states are the same (%s)." % (transition_label,state_label))
            continue

          # Define the initial density matrix file name and content

          dark_label = system['eigenstates_list'][dark_index]["label"]

          target_file = dark_label + "_"

          target_content = np.zeros((len(system['eigenstates_list']), len(system['eigenstates_list'])),dtype=complex)  # Quick init of a zero-filled matrix
          target_content[dark_index][dark_index] = complex(1)

          # ========================================================= #
          #           Building the transition dictionary              #
          # ========================================================= #

          # Building the transition dictionary

          transition = {
            "label" : momdip_key + "_" + str(iter_bright+1) + "B" + str(iter_dark+1) + "D_" + bright_label + "-" + dark_label,
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

          transition = {
            "label" : "R_" + momdip_key + "_" + str(iter_dark+1) + "D" + str(iter_bright+1) + "B_" + dark_label + "-" +  bright_label,
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

    gs_state = [state for state in system['states_list'] if state['number'] == 0][0]

    init_file = gs_state['label'] + "_"

    init_content = np.zeros((len(system['states_list']), len(system['states_list'])),dtype=complex)  # Quick init of a zero-filled matrix
    init_content[gs_state['number']][gs_state['number']] = complex(1)

    # Iterate over the dark zero order states

    for state in system['states_list']:

      if state['type'].lower() == "dark":

        # Build the target density matrix in the zero order states

        target_mtx = np.zeros((len(system['states_list']), len(system['states_list'])),dtype=complex)  # Quick init of a zero-filled matrix
        target_mtx[state['number']][state['number']] = complex(1)

        # Convert target density matrix to the eigenstates basis set and define the target density matrix file name 

        transpose_inv = scipy.linalg.inv(system['transpose'])
        target_content = np.matmul(np.matmul(transpose_inv,target_mtx),system['transpose'])

        target_file = state['label'] + "_"

        # Building the transition dictionary (one for each transition dipole moment matrix)

        for momdip_key in system['momdip_mtx']:

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

