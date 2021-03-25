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


def build_transition(init_state:int,target_state:int,init_label:str,target_label:str,system:dict):
    """Determines the transition dictionary and the transition files for the transition going from the first given state to the second one.

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
    system : dict
        Information extracted by the parsing function and derived from it.

    Returns
    -------
    transition : dict
        Dictionary containing eight keys:

          - ``label`` is the label of the transition, which will be used for the name of the job directories.
          - ``init_state`` is the number of the initial state.
          - ``target_state`` is the number of the target state.
          - ``energy`` is the transition energy between the two states.
          - ``init_file`` is the name of the initial state file, minus the number at the end.
          - ``init_content`` is the content of the initial state file.
          - ``target_file`` is the name of the target state file, minus the number at the end.
          - ``target_content`` is the content of the target state file.
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
      "label" : init_label + "-" + target_label,
      "init_state" : init_state,
      "target_state" : target_state,
      "energy" : energy,
      "init_file" : init_file,
      "init_content" : init_content,
      "target_file" : target_file,
      "target_content" : target_content
      }
    
    return transition

# =================================================================== #
# =================================================================== #
#                        Transition functions                         #
# =================================================================== #
# =================================================================== #

def closest_bright_to_dark(system:dict):
    """Determines the transition files needed by QOCT-RA for the transition between the bright and dark states with the lowest transition energy between themselves.

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
        
        Note that in this case, this list only contains one dictionary as only one transition is considered.
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
        continue
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

    transition = build_transition(bright_number,dark_number,bright_label,dark_label,system)

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

def brightests_to_coupled_darks_and_closest(system:dict):
    """Determines the transition files needed by QOCT-RA for transitions going from each of the two brightest states of the molecule to each of their two most coupled dark states. For good measure, the transition between the bright and dark states with the lowest transition energy is also added, if it does not correspond to one of the four transitions already considered.

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

    # Initialize the dictionary that will be returned by the function

    transitions_list = []

    # ========================================================= #
    #            Identifying bright and dark states             #
    # ========================================================= #

    bright_list = []
    dark_list = []

    for state in system['states_list']:
      if state['number'] == 0:
        continue
      elif state['type'].lower() == "dark":
        dark_list.append(state['number'])
      elif state['type'].lower() == "bright":
        bright_list.append(state['number'])

    # ========================================================= #
    #                        Subfunction                        #
    # ========================================================= #

    def brightest_to_most_coupled(bright_list:list,dark_list:list,system:dict):
        """Determines the brightest state of the molecule and its most coupled dark state.

        Parameters
        ----------
        bright_list : list
            List of the bright state numbers.
        dark_list : list
            List of the dark state numbers.
        system : dict
            Information extracted by the parsing function and derived from it.

        Returns
        -------
        bright_number : int
          Number of the brightest state.
        dark_number : int
          Number of the dark state presenting the highest coupling with the brightest state.
        """

        # Ensure we are working with lists and not numpy arrays

        if isinstance(system['momdip_mtx'], numpy.ndarray):
          momdip_mtx = system['momdip_mtx'].tolist()
        if isinstance(system['momdip_mtx'], list):
          momdip_mtx = system['momdip_mtx']

        if isinstance(system['mime'], numpy.ndarray):
          mime = system['mime'].tolist()
        if isinstance(system['mime'], list):
          mime = system['mime']

        # Identifying the brightest state

        highest_mom = -1

        for moment in momdip_mtx[0]:
          if (momdip_mtx[0].index(moment) in bright_list) and (moment > highest_mom):
            highest_mom = moment
            bright_number = momdip_mtx[0].index(moment)

        # Identifying the dark state with the strongest coupling

        strongest_cou = -1

        for coupling in mime[bright_number]:
          if (mime[bright_number].index(coupling) in dark_list) and (coupling > strongest_cou):
            strongest_cou = coupling
            dark_number = mime[bright_number].index(coupling)
        
        return bright_number,dark_number

    # ========================================================= #
    #      Brightest state to its most coupled dark state       #
    # ========================================================= #

    bright_number,dark_number = brightest_to_most_coupled(bright_list,dark_list,system)

    bright_label = system['states_list'][bright_number]["label"]
    dark_label = system['states_list'][dark_number]["label"]

    transition = build_transition(bright_number,dark_number,bright_label,dark_label,system)
    transition.update({ "label": "1B1C_" + bright_label + "-" + dark_label }) 

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

    bright_number,dark_number = brightest_to_most_coupled(bright_list,copy_dark_list,system)

    bright_label = system['states_list'][bright_number]["label"]
    dark_label = system['states_list'][dark_number]["label"]

    transition = build_transition(bright_number,dark_number,bright_label,dark_label,system)
    transition.update({ "label": "1B2C_" + bright_label + "-" + dark_label }) 

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

    bright_number,dark_number = brightest_to_most_coupled(copy_bright_list,dark_list,system)

    bright_label = system['states_list'][bright_number]["label"]
    dark_label = system['states_list'][dark_number]["label"]

    transition = build_transition(bright_number,dark_number,bright_label,dark_label,system)
    transition.update({ "label": "2B1C_" + bright_label + "-" + dark_label }) 

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

    bright_number,dark_number = brightest_to_most_coupled(copy_bright_list,copy_dark_list,system)

    bright_label = system['states_list'][bright_number]["label"]
    dark_label = system['states_list'][dark_number]["label"]

    transition = build_transition(bright_number,dark_number,bright_label,dark_label,system)
    transition.update({ "label": "2B2C_" + bright_label + "-" + dark_label }) 

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

    return transitions_list
