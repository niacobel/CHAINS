################################################################################################################################################
##                                                            Transition functions                                                            ##
##                                                                                                                                            ##
##                                     This script contains the transition functions for CONTROL LAUNCHER,                                    ##
##                                consult the documentation at https://chains-ulb.readthedocs.io/ for details                                 ##
################################################################################################################################################

import contextlib
import os

import numpy

import control_errors


def build_transition(init_state:int,target_state:int,init_label:str,target_label:str,momdip_key:str,system:dict):
    """Build the transition dictionary and the transition files for the transition going from the first given state to the second one.

    Parameters
    ----------
    init_state : int
        Number of the initial state.
    target_state : int
        Number of the target state.
    init_label : str
        Label of the initial state.
    target_label : str
        Label of the target state.
    momdip_key : str
        Key of the transition dipole moments matrix used for this transition.
    system : dict
        Information extracted by the parsing function and derived from it.

    Returns
    -------
    transition : dict
        Dictionary containing nine keys:

          - ``label`` is the label of the transition, which will be used for the name of the job directories.
          - ``init_state`` is the number of the initial state.
          - ``target_state`` is the number of the target state.
          - ``energy`` is the transition energy between the two states.
          - ``init_file`` is the name of the initial state file, minus the number at the end.
          - ``init_content`` is the content of the initial state file.
          - ``target_file`` is the name of the target state file, minus the number at the end.
          - ``target_content`` is the content of the target state file.
          - ``momdip_key`` is the key of the transition dipole moments matrix used for this transition.
    """

    # Defining the projector for the initial state

    init_proj = numpy.zeros((len(system['eigenvalues']), len(system['eigenvalues'])),dtype=complex)  # Quick init of a zero-filled matrix
    
    init_proj[init_state][init_state] = 1 + 0j

    # Creating the initial population file content by converting the projector to the eigenstates basis set

    init_file = init_label + "_"

    init_content = numpy.matmul(numpy.matmul(system['transpose'],init_proj),system['eigenvectors'])

    # Defining the projector for the target state

    target_proj = numpy.zeros((len(system['eigenvalues']), len(system['eigenvalues'])),dtype=complex)  # Quick init of a zero-filled matrix
    
    target_proj[target_state][target_state] = 1 + 0j

    # Creating the target population file content by converting the projector to the eigenstates basis set
    
    target_file = target_label + "_"

    target_content = numpy.matmul(numpy.matmul(system['transpose'],target_proj),system['eigenvectors'])

    # Calculate the transition energy (in Ha)

    energy = abs(system['mime'][init_state][init_state] - system['mime'][target_state][target_state]) #! Eigenstates or zero order states?

    # Building the transition dictionary

    transition = {
      "label" : momdip_key + "_" + init_label + "-" + target_label,
      "init_state" : init_state,
      "target_state" : target_state,
      "energy" : energy,
      "init_file" : init_file,
      "init_content" : init_content,
      "target_file" : target_file,
      "target_content" : target_content,
      "momdip_key" : momdip_key
      }
    
    return transition

# =================================================================== #
# =================================================================== #
#                        Transition functions                         #
# =================================================================== #
# =================================================================== #

def closest_bright_to_dark(system:dict):
    """Determines the transition files needed by QOCT-RA for the transition between the bright and dark states with the lowest transition energy between themselves, in the eigenstates basis set.

    Parameters
    ----------
    system : dict
        Information extracted by the parsing function and derived from it.

    Returns
    -------
    transitions_list : list
        List of dictionaries containing nine keys each: 

          - ``label`` is the label of the transition, which will be used for the name of the job directories.
          - ``init_state`` is the number of the initial state.
          - ``target_state`` is the number of the target state.
          - ``energy`` is the transition energy between the two states.
          - ``init_file`` is the name of the initial state file, minus the number at the end.
          - ``init_content`` is the content of the initial state file.
          - ``target_file`` is the name of the target state file, minus the number at the end.
          - ``target_content`` is the content of the target state file.
          - ``momdip_key`` is the key of the transition dipole moments matrix used for this transition.
    """

    # Initialize the dictionary that will be returned by the function

    transitions_list = []

    # ========================================================= #
    #            Identifying bright and dark states             #
    # ========================================================= #

    bright_list = []
    dark_list = []

    for state in system['states_list']:
      if state['number'] == 0:
        continue # exclude the ground state
      elif state['type'].lower() == "dark":
        dark_list.append(state['number'])
      elif state['type'].lower() == "bright":
        bright_list.append(state['number'])

    # ========================================================= #
    #      Transition between the two closest bright-dark       #
    # ========================================================= #    

    minimum = float('inf')

    for bright in bright_list:
      for dark in dark_list:
        energy = abs(system['eigenvalues'][bright] - system['eigenvalues'][dark]) # Using eigenvalues to better distinguish between degenerated zero order states
        if energy < minimum:
          minimum = energy
          bright_number = bright
          dark_number = dark

    bright_label = system['states_list'][bright_number]["label"]
    dark_label = system['states_list'][dark_number]["label"]

    # We need to consider each transition dipole moments matrix separately
    for momdip_key in system['momdip_es_mtx']:    

      transition = build_transition(bright_number,dark_number,bright_label,dark_label,momdip_key,system)

      print("")
      print(''.center(50, '-'))
      print("{:<20} {:<30}".format("Label: ", transition["label"]))
      print(''.center(50, '-'))
      print("{:<20} {:<30}".format("Initial state: ", "%s (%s)" % (bright_label,bright_number)))
      print("{:<20} {:<30}".format("Target state: ", "%s (%s)" % (dark_label,dark_number)))
      print("{:<20} {:<30}".format("Energy (Ha): ", "{:.4e}".format(transition["energy"])))
      print(''.center(50, '-'))

      transitions_list.append(transition)

    return transitions_list

def brightests_to_coupled_darks(system:dict):
    """Determines the transition files needed by QOCT-RA for transitions going from the brightest states of the molecule to their most coupled ('brightest') dark states.

    Parameters
    ----------
    system : dict
        Information extracted by the parsing function and derived from it.

    Returns
    -------
    transitions_list : list
        List of dictionaries containing eight keys each: 

          - ``label`` is the label of the transition, which will be used for the name of the job directories.
          - ``init_state`` is the number of the initial state.
          - ``target_state`` is the number of the target state.
          - ``energy`` is the transition energy between the two states.
          - ``init_file`` is the name of the initial state file, minus the number at the end.
          - ``init_content`` is the content of the initial state file.
          - ``target_file`` is the name of the target state file, minus the number at the end.
          - ``target_content`` is the content of the target state file.
    """

    bmax = 2 # Maximum number of bright states to consider (e.g. if bmax = 3, then only the three brightest states will be considered)
    dmax = 2 # Maximum number of coupled dark states to consider (e.g. if dmax = 3, then only the three most coupled dark states will be considered)

    # Initialize the dictionary that will be returned by the function

    transitions_list = []

    # ========================================================= #
    #            Identifying bright and dark states             #
    # ========================================================= #

    bright_list = []
    dark_list = []

    for state in system['states_list']:
      if state['number'] == 0:
        continue # exclude the ground state
      elif state['type'].lower() == "dark":
        dark_list.append(state['number'])
      elif state['type'].lower() == "bright":
        bright_list.append(state['number'])

    # Adjust bmax and dmax if needed

    if bmax > len(bright_list):
      bmax = len(bright_list)
    if dmax > len(dark_list):
      dmax = len(dark_list)

    # ========================================================= #
    #                        Subfunction                        #
    # ========================================================= #

    def brightest_to_most_coupled(bright_list:list,dark_list:list,momdip_mtx):
        """Determines the brightest state of the molecule and its most coupled dark state for a given transition dipole moments matrix.

        Parameters
        ----------
        bright_list : list
            List of the bright state numbers.
        dark_list : list
            List of the dark state numbers.
        momdip_mtx
            Either a list or a NumPy array representing the transition dipole moments matrix.

        Returns
        -------
        bright_number : int
          Number of the brightest state.
        dark_number : int
          Number of the dark state presenting the highest coupling with the brightest state.
        """

        # Ensure we are working with lists and not numpy arrays

        if isinstance(momdip_mtx, numpy.ndarray):
          matrix = momdip_mtx.tolist()
        if isinstance(momdip_mtx, list):
          matrix = momdip_mtx

        # Identifying the brightest state

        brightest_mom = 0

        for moment in matrix[0]:
          if (matrix[0].index(moment) in bright_list) and (abs(moment) > abs(brightest_mom)):
            brightest_mom = moment
            bright_number = matrix[0].index(moment)

        # Identifying the most coupled dark state

        best_dark = 0

        for moment in matrix[bright_number]:
          if (matrix[bright_number].index(moment) in dark_list) and (abs(moment) > abs(best_dark)):
            best_dark = moment
            dark_number = matrix[bright_number].index(moment)
        
        return bright_number,dark_number

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

      # Consider the dipole moments for transition involving the ground state (and make sure it's a list and not a NumPy array)
      gs_line = system['momdip_es_mtx'][momdip_key][0]
      if isinstance(gs_line, numpy.ndarray):
        gs_line = gs_line.tolist()

      # Sort the number of the bright states by decreasing order of the absolute value of their transition dipole moments
      # e.g. if gs_line = [0.0, 0.0, 0.0, 0.0, 0.0, -1.2, 0.8, 1.0, -0.4] and bright_list = [5, 6, 7, 8], then bright_max_numbers = [1, 3, 2, 4]
      abs_gs_line = [abs(mom) for mom in gs_line]
      bright_max_numbers = [abs_gs_line.index(mom) for mom in sorted(abs_gs_line, reverse=True) if abs_gs_line.index(mom) in bright_list]

      # Iterate over the bright states
      # ==============================

      # Start iterating over the first Nth brightest bright states, where N is the bmax argument (see https://stackoverflow.com/questions/36106712/how-can-i-limit-iterations-of-a-loop-in-python for details)
      for iter_bright, bright_number in zip(range(bmax), bright_max_numbers):

      # iter_bright is the index of the current iteration of this loop, e.g if bmax = 3, then iter_bright will be 0, then 1, then 2.
      # bright_number is the number of the bright state currently considered, e.g. if bright_max_numbers = [1, 3, 2, 4] and iter_bright = 1, then bright_number = 3.

        # Sorting the dark states
        # =======================

        # Consider the dipole moments for transition involving the current bright state (and make sure it's a list and not a NumPy array)
        bright_line = system['momdip_es_mtx'][momdip_key][bright_number]
        if isinstance(bright_line, numpy.ndarray):
          bright_line = bright_line.tolist()

        # Sort the number of the dark states by decreasing order of the absolute value of their transition dipole moments with the current bright state
        # e.g. if bright_line = [0.0, -0.2, -0.8, 1.1, -1.4, 0.0, 0.0, 0.0, 0.0, ] and dark_list = [1, 2, 3, 4], then dark_max_numbers = [4, 3, 2, 1]
        abs_bright_line = [abs(mom) for mom in bright_line]
        dark_max_numbers = [abs_bright_line.index(mom) for mom in sorted(abs_bright_line, reverse=True) if abs_bright_line.index(mom) in dark_list]

        # Iterate over the dark states
        # ============================

        # Start iterating over the first Nth most coupled dark states, where N is the dmax argument (see https://stackoverflow.com/questions/36106712/how-can-i-limit-iterations-of-a-loop-in-python for details)
        for iter_dark, dark_number in zip(range(dmax), dark_max_numbers):

        # iter_dark is the index of the current iteration of this loop, e.g if dmax = 2, then iter_dark will be 0, then 1.
        # dark_number is the number of the dark state currently considered, e.g. if dark_max_numbers = [4, 3, 2, 1] and iter_bright = 0, then dark_number = 4.

          # ========================================================= #
          #            Building the transition dictionary             #
          # ========================================================= #

          bright_label = system['states_list'][bright_number]["label"]
          dark_label = system['states_list'][dark_number]["label"]

          # Call the build_transition function using the current information for this transition

          transition = build_transition(bright_number,dark_number,bright_label,dark_label,momdip_key,system)

          # Update the label to something less generic and more informative
          # e.g. X_1B2D_S3-T2 indicates a transition based on the transition dipole moments matrix associated with the "X" key (X axis), involving the brightest bright state (1B) and its second most coupled dark state (2D), the labels of those states being S3 and T2, respectively.

          transition_label = momdip_key + "_" + str(iter_bright+1) + "B" + str(iter_dark+1) + "D_" + bright_label + "-" + dark_label 
          transition.update({ "label": transition_label}) 

          # Pretty recap for the log file

          print("")
          print(''.center(50, '-'))
          print("{:<20} {:<30}".format("Label: ", transition["label"]))
          print(''.center(50, '-'))
          print("{:<20} {:<30}".format("Initial state: ", "%s (%s)" % (bright_label,bright_number)))
          print("{:<20} {:<30}".format("Target state: ", "%s (%s)" % (dark_label,dark_number)))
          print("{:<20} {:<30}".format("Energy (Ha): ", "{:.4e}".format(transition["energy"])))
          print(''.center(50, '-'))

          # Add the transition to the transitions list, before proceeding with the next most coupled dark state

          transitions_list.append(transition)

      """  
      # ========================================================= #
      #      Brightest state to its most coupled dark state       #
      # ========================================================= #

      bright_number,dark_number = brightest_to_most_coupled(bright_list,dark_list,system['momdip_es_mtx'][key])

      bright_label = system['states_list'][bright_number]["label"]
      dark_label = system['states_list'][dark_number]["label"]

      transition = build_transition(bright_number,dark_number,bright_label,dark_label,system)
      transition.update({ "label": key + "_1B1D_" + bright_label + "-" + dark_label }) 

      print("")
      print(''.center(50, '-'))
      print("{:<20} {:<30}".format("Label: ", transition["label"]))
      print(''.center(50, '-'))
      print("{:<20} {:<30}".format("Initial state: ", "%s (%s)" % (bright_label,bright_number)))
      print("{:<20} {:<30}".format("Target state: ", "%s (%s)" % (dark_label,dark_number)))
      print("{:<20} {:<30}".format("Energy (Ha): ", "{:.4e}".format(transition["energy"])))
      print(''.center(50, '-'))

      transitions_list.append(transition)

      # ========================================================= #
      #     Brightest state to its 2nd most coupled dark state    #
      # ========================================================= #    

      copy_dark_list = dark_list[:]

      copy_dark_list.remove(dark_number)

      bright_number,dark_number = brightest_to_most_coupled(bright_list,copy_dark_list,system['momdip_es_mtx'][key])

      bright_label = system['states_list'][bright_number]["label"]
      dark_label = system['states_list'][dark_number]["label"]

      transition = build_transition(bright_number,dark_number,bright_label,dark_label,system)
      transition.update({ "label": key + "_1B2D_" + bright_label + "-" + dark_label }) 

      print("")
      print(''.center(50, '-'))
      print("{:<20} {:<30}".format("Label: ", transition["label"]))
      print(''.center(50, '-'))
      print("{:<20} {:<30}".format("Initial state: ", "%s (%s)" % (bright_label,bright_number)))
      print("{:<20} {:<30}".format("Target state: ", "%s (%s)" % (dark_label,dark_number)))
      print("{:<20} {:<30}".format("Energy (Ha): ", "{:.4e}".format(transition["energy"])))
      print(''.center(50, '-'))

      transitions_list.append(transition)

      # ========================================================= #
      #    2nd brightest state to its most coupled dark state     #
      # ========================================================= #

      copy_bright_list = bright_list[:]

      copy_bright_list.remove(bright_number)

      bright_number,dark_number = brightest_to_most_coupled(copy_bright_list,dark_list,system['momdip_es_mtx'][key])

      bright_label = system['states_list'][bright_number]["label"]
      dark_label = system['states_list'][dark_number]["label"]

      transition = build_transition(bright_number,dark_number,bright_label,dark_label,system)
      transition.update({ "label": key + "_2B1D_" + bright_label + "-" + dark_label }) 

      print("")
      print(''.center(50, '-'))
      print("{:<20} {:<30}".format("Label: ", transition["label"]))
      print(''.center(50, '-'))
      print("{:<20} {:<30}".format("Initial state: ", "%s (%s)" % (bright_label,bright_number)))
      print("{:<20} {:<30}".format("Target state: ", "%s (%s)" % (dark_label,dark_number)))
      print("{:<20} {:<30}".format("Energy (Ha): ", "{:.4e}".format(transition["energy"])))
      print(''.center(50, '-'))

      transitions_list.append(transition)

      # ========================================================= #
      #   2nd brightest state to its 2nd most coupled dark state  #
      # ========================================================= #    

      copy_dark_list = dark_list[:]

      copy_dark_list.remove(dark_number)

      bright_number,dark_number = brightest_to_most_coupled(copy_bright_list,copy_dark_list,system['momdip_es_mtx'][key])

      bright_label = system['states_list'][bright_number]["label"]
      dark_label = system['states_list'][dark_number]["label"]

      transition = build_transition(bright_number,dark_number,bright_label,dark_label,system)
      transition.update({ "label": key + "_2B2D_" + bright_label + "-" + dark_label }) 

      print("")
      print(''.center(50, '-'))
      print("{:<20} {:<30}".format("Label: ", transition["label"]))
      print(''.center(50, '-'))
      print("{:<20} {:<30}".format("Initial state: ", "%s (%s)" % (bright_label,bright_number)))
      print("{:<20} {:<30}".format("Target state: ", "%s (%s)" % (dark_label,dark_number)))
      print("{:<20} {:<30}".format("Energy (Ha): ", "{:.4e}".format(transition["energy"])))
      print(''.center(50, '-'))

      transitions_list.append(transition)

    # ========================================================= #
    #   Transition between the two closest bright-dark states   #
    # ========================================================= #    

    # Call the "closest_bright_to_dark" function, but without its standard output (https://stackoverflow.com/questions/2828953/silence-the-stdout-of-a-function-in-python-without-trashing-sys-stdout-and-resto)
    # Note that fifth_transition is a list containing only one dictionary (see closest_bright_to_dark function definition)
    with open(os.devnull, 'w') as devnull:
      with contextlib.redirect_stdout(devnull):
        fifth_transition = closest_bright_to_dark(system)

    if not any(transition["energy"] == fifth_transition[0]["energy"] for transition in transitions_list):

      transition = fifth_transition[0]
      transition.update({ "label": "LE_" + transition["label"] }) 

      bright_number = transition["init_state"]
      dark_number = transition["target_state"]

      bright_label = system['states_list'][bright_number]["label"]
      dark_label = system['states_list'][dark_number]["label"]

      print("")
      print(''.center(50, '-'))
      print("{:<20} {:<30}".format("Label: ", transition["label"]))
      print(''.center(50, '-'))
      print("{:<20} {:<30}".format("Initial state: ", "%s (%s)" % (bright_label,bright_number)))
      print("{:<20} {:<30}".format("Target state: ", "%s (%s)" % (dark_label,dark_number)))
      print("{:<20} {:<30}".format("Energy (Ha): ", "{:.4e}".format(transition["energy"])))
      print(''.center(50, '-'))

      transitions_list.append(transition)

    else:
      print("\nThe transition between the bright and dark pair with the lowest transition energy (%s) is already included so only four transitions will be considered." % fifth_transition[0]["label"])
    """

    return transitions_list
