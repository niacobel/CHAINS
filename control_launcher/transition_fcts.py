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


def build_transition(state1:int,state2:int,label1:str,label2:str,system:dict):
    """Determines the transition dictionary and the transition files for the transition going from the first given state to the second one.

    Parameters
    ----------
    state1 : int
        Number of the initial state.
    state2 : int
        Number of the target state.
    label1 : str
        Label of the initial state.
    label2 : str
        Label of the target state.
    system : dict
        Information extracted by the parsing function and derived from it.

    Returns
    -------
    transition : dict
        Dictionary containing eight keys:

          - ``label`` is the label of the transition, which will be used for the name of the job directories.
          - ``state1`` is the number of the initial state.
          - ``state2`` is the number of the target state.
          - ``transition_energy`` is the transition energy between the two states.
          - ``init_file`` is the name of the initial state file, minus the number at the end.
          - ``init_content`` is the content of the initial state file.
          - ``target_file`` is the name of the target state file, minus the number at the end.
          - ``target_content`` is the content of the target state file.
    """

    # Defining the projector for the initial state

    init_proj = numpy.zeros((len(system['eigenvalues']), len(system['eigenvalues'])),dtype=complex)  # Quick init of a zero-filled matrix
    
    init_proj[state1][state1] = 1 + 0j

    # Creating the initial population file content by converting the projector to the eigenstates basis set

    init_file = label1 + "_"

    init_content = numpy.matmul(numpy.matmul(system['transpose'],init_proj),system['eigenvectors'])

    # Defining the projector for the target state

    target_proj = numpy.zeros((len(system['eigenvalues']), len(system['eigenvalues'])),dtype=complex)  # Quick init of a zero-filled matrix
    
    target_proj[state2][state2] = 1 + 0j

    # Creating the target population file content by converting the projector to the eigenstates basis set
    
    target_file = label2 + "_"

    target_content = numpy.matmul(numpy.matmul(system['transpose'],target_proj),system['eigenvectors'])

    # Calculate the transition energy (in Ha)

    energy = abs(system['mime'][state1][state1] - system['mime'][state2][state2])

    # Building the transition dictionary

    transition = {
      "label" : label1 + "-" + label2,
      "state1" : state1,
      "state2" : state2,
      "transition_energy" : energy,
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

def closest_singlet_to_triplet(system:dict):
    """Determines the transition files needed by QOCT-RA for the transition between the singlet and triplet states with the lowest transition energy between themselves.

    Parameters
    ----------
    system : dict
        Information extracted by the parsing function and derived from it.
        This dictionary needs an additionnal ``states_list`` key which is a list of dictionaries containing at least three keys each: ``number``, ``multiplicity`` and ``label`` where

          - ``number`` is the number of the state, starting at 0 (which is the ground state).
          - ``multiplicity`` is the multiplicity of the state (ex: Singlet, Triplet).
          - ``label`` is the label of the state, which will be used for the label of the transitions.

    Returns
    -------
    transitions_list : list
        List of dictionaries containing eight keys each: 

          - ``label`` is the label of the transition, which will be used for the name of the job directories.
          - ``state1`` is the number of the initial state.
          - ``state2`` is the number of the target state.
          - ``transition_energy`` is the transition energy between the two states.
          - ``init_file`` is the name of the initial state file, minus the number at the end.
          - ``init_content`` is the content of the initial state file.
          - ``target_file`` is the name of the target state file, minus the number at the end.
          - ``target_content`` is the content of the target state file.
        
        Note that in this case, this list only contains one dictionary as only one transition is considered.
    """

    # Initialize the dictionary that will be returned by the function

    transitions_list = []

    # ========================================================= #
    #             Identifying singlets and triplets             #
    # ========================================================= #

    singlet_list = []
    triplet_list = []

    for state in system['states_list']:
      if state['number'] == 0:
        continue
      elif state['multiplicity'].lower() == "triplet":
        triplet_list.append(state['number'])
      elif state['multiplicity'].lower() == "singlet":
        singlet_list.append(state['number'])

    # ========================================================= #
    #    Transition between the two closest singlet-triplet     #
    # ========================================================= #    

    minimum = float('inf')

    for singlet in singlet_list:
      for triplet in triplet_list:
        energy = abs(system['eigenvalues'][singlet] - system['eigenvalues'][triplet]) # Using eigenvalues to better distinguish between degenerated zero order states
        if energy < minimum:
          minimum = energy
          singlet_number = singlet
          triplet_number = triplet

    singlet_label = system['states_list'][singlet_number]["label"]
    triplet_label = system['states_list'][triplet_number]["label"]

    transition = build_transition(singlet_number,triplet_number,singlet_label,triplet_label,system)

    print("")
    print(''.center(50, '-'))
    print("{:<20} {:<30}".format("Label: ", transition["label"]))
    print(''.center(50, '-'))
    print("{:<20} {:<30}".format("Initial state: ", "%s (%s)" % (singlet_label,singlet_number)))
    print("{:<20} {:<30}".format("Target state: ", "%s (%s)" % (triplet_label,triplet_number)))
    print("{:<20} {:<30}".format("Energy (Ha): ", "{:.4e}".format(transition["transition_energy"])))
    print(''.center(50, '-'))

    transitions_list.append(transition)

    return transitions_list

def bright_singlets_to_coupled_triplets_and_closest(system:dict):
    """Determines the transition files needed by QOCT-RA for transitions going from each of the two brightest electronic singlet states of the molecule to each of their two most coupled triplet states. For good measure, the transition between the singlet and triplet states with the lowest transition energy is also added, if it does not correspond to one of the four transitions already considered.

    Parameters
    ----------
    system : dict
        Information extracted by the parsing function and derived from it.
        This dictionary needs an additionnal ``states_list`` key which is a list of dictionaries containing at least three keys each: ``number``, ``multiplicity`` and ``label`` where

          - ``number`` is the number of the state, starting at 0 (which is the ground state).
          - ``multiplicity`` is the multiplicity of the state (ex: Singlet, Triplet).
          - ``label`` is the label of the state, which will be used for the label of the transitions.

    Returns
    -------
    transitions_list : list
        List of dictionaries containing eight keys each: 

          - ``label`` is the label of the transition, which will be used for the name of the job directories.
          - ``state1`` is the number of the initial state.
          - ``state2`` is the number of the target state.
          - ``transition_energy`` is the transition energy between the two states.
          - ``init_file`` is the name of the initial state file, minus the number at the end.
          - ``init_content`` is the content of the initial state file.
          - ``target_file`` is the name of the target state file, minus the number at the end.
          - ``target_content`` is the content of the target state file.
    """

    # Initialize the dictionary that will be returned by the function

    transitions_list = []

    # ========================================================= #
    #             Identifying singlets and triplets             #
    # ========================================================= #

    singlet_list = []
    triplet_list = []

    for state in system['states_list']:
      if state['number'] == 0:
        continue
      elif state['multiplicity'].lower() == "triplet":
        triplet_list.append(state['number'])
      elif state['multiplicity'].lower() == "singlet":
        singlet_list.append(state['number'])

    # ========================================================= #
    #                        Subfunction                        #
    # ========================================================= #

    def brightest_to_most_coupled(singlet_list:list,triplet_list:list,system:dict):
        """Determines the brightest electronic singlet state of the molecule and its most coupled triplet state.

        Parameters
        ----------
        singlet_list : list
            List of the state numbers with the singlet multiplicity.
        triplet_list : list
            List of the state numbers with the triplet multiplicity.
        system : dict
            Information extracted by the parsing function and derived from it.

        Returns
        -------
        singlet_number : int
          Number of the brightest singlet state.
        triplet_number : int
          Number of the triplet state presenting the highest SOC with the brightest singlet state.
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

        # Identifying the brightest singlet

        highest_mom = -1

        for moment in momdip_mtx[0]:
          if (momdip_mtx[0].index(moment) in singlet_list) and (moment > highest_mom):
            highest_mom = moment
            singlet_number = momdip_mtx[0].index(moment)

        # Identifying the triplet with the strongest SOC

        strongest_soc = -1

        for soc in mime[singlet_number]:
          if (mime[singlet_number].index(soc) in triplet_list) and (soc > strongest_soc):
            strongest_soc = soc
            triplet_number = mime[singlet_number].index(soc)
        
        return singlet_number,triplet_number

    # ========================================================= #
    #       Brightest singlet to its most coupled triplet       #
    # ========================================================= #

    singlet_number,triplet_number = brightest_to_most_coupled(singlet_list,triplet_list,system)

    singlet_label = system['states_list'][singlet_number]["label"]
    triplet_label = system['states_list'][triplet_number]["label"]

    transition = build_transition(singlet_number,triplet_number,singlet_label,triplet_label,system)
    transition.update({ "label": "1B1C_" + singlet_label + "-" + triplet_label }) 

    print("")
    print(''.center(50, '-'))
    print("{:<20} {:<30}".format("Label: ", transition["label"]))
    print(''.center(50, '-'))
    print("{:<20} {:<30}".format("Initial state: ", "%s (%s)" % (singlet_label,singlet_number)))
    print("{:<20} {:<30}".format("Target state: ", "%s (%s)" % (triplet_label,triplet_number)))
    print("{:<20} {:<30}".format("Energy (Ha): ", "{:.4e}".format(transition["transition_energy"])))
    print(''.center(50, '-'))

    transitions_list.append(transition)

    # ========================================================= #
    #     Brightest singlet to its 2nd most coupled triplet     #
    # ========================================================= #    

    copy_triplet_list = triplet_list[:]

    copy_triplet_list.remove(triplet_number)

    singlet_number,triplet_number = brightest_to_most_coupled(singlet_list,copy_triplet_list,system)

    singlet_label = system['states_list'][singlet_number]["label"]
    triplet_label = system['states_list'][triplet_number]["label"]

    transition = build_transition(singlet_number,triplet_number,singlet_label,triplet_label,system)
    transition.update({ "label": "1B2C_" + singlet_label + "-" + triplet_label }) 

    print("")
    print(''.center(50, '-'))
    print("{:<20} {:<30}".format("Label: ", transition["label"]))
    print(''.center(50, '-'))
    print("{:<20} {:<30}".format("Initial state: ", "%s (%s)" % (singlet_label,singlet_number)))
    print("{:<20} {:<30}".format("Target state: ", "%s (%s)" % (triplet_label,triplet_number)))
    print("{:<20} {:<30}".format("Energy (Ha): ", "{:.4e}".format(transition["transition_energy"])))
    print(''.center(50, '-'))

    transitions_list.append(transition)

    # ========================================================= #
    #     2nd brightest singlet to its most coupled triplet     #
    # ========================================================= #

    copy_singlet_list = singlet_list[:]

    copy_singlet_list.remove(singlet_number)

    singlet_number,triplet_number = brightest_to_most_coupled(copy_singlet_list,triplet_list,system)

    singlet_label = system['states_list'][singlet_number]["label"]
    triplet_label = system['states_list'][triplet_number]["label"]

    transition = build_transition(singlet_number,triplet_number,singlet_label,triplet_label,system)
    transition.update({ "label": "2B1C_" + singlet_label + "-" + triplet_label }) 

    print("")
    print(''.center(50, '-'))
    print("{:<20} {:<30}".format("Label: ", transition["label"]))
    print(''.center(50, '-'))
    print("{:<20} {:<30}".format("Initial state: ", "%s (%s)" % (singlet_label,singlet_number)))
    print("{:<20} {:<30}".format("Target state: ", "%s (%s)" % (triplet_label,triplet_number)))
    print("{:<20} {:<30}".format("Energy (Ha): ", "{:.4e}".format(transition["transition_energy"])))
    print(''.center(50, '-'))

    transitions_list.append(transition)

    # ========================================================= #
    #   2nd brightest singlet to its 2nd most coupled triplet   #
    # ========================================================= #    

    copy_triplet_list = triplet_list[:]

    copy_triplet_list.remove(triplet_number)

    singlet_number,triplet_number = brightest_to_most_coupled(copy_singlet_list,copy_triplet_list,system)

    singlet_label = system['states_list'][singlet_number]["label"]
    triplet_label = system['states_list'][triplet_number]["label"]

    transition = build_transition(singlet_number,triplet_number,singlet_label,triplet_label,system)
    transition.update({ "label": "2B2C_" + singlet_label + "-" + triplet_label }) 

    print("")
    print(''.center(50, '-'))
    print("{:<20} {:<30}".format("Label: ", transition["label"]))
    print(''.center(50, '-'))
    print("{:<20} {:<30}".format("Initial state: ", "%s (%s)" % (singlet_label,singlet_number)))
    print("{:<20} {:<30}".format("Target state: ", "%s (%s)" % (triplet_label,triplet_number)))
    print("{:<20} {:<30}".format("Energy (Ha): ", "{:.4e}".format(transition["transition_energy"])))
    print(''.center(50, '-'))

    transitions_list.append(transition)

    # ========================================================= #
    #    Transition between the two closest singlet-triplet     #
    # ========================================================= #    

    # Call the "closest_singlet_to_triplet" function, but without its standard output (https://stackoverflow.com/questions/2828953/silence-the-stdout-of-a-function-in-python-without-trashing-sys-stdout-and-resto)
    # Note that fifth_transition is a list containing only one dictionary (see closest_singlet_to_triplet function definition)
    with open(os.devnull, 'w') as devnull:
      with contextlib.redirect_stdout(devnull):
        fifth_transition = closest_singlet_to_triplet(system)

    if not any(transition["transition_energy"] == fifth_transition[0]["transition_energy"] for transition in transitions_list):

      transition = fifth_transition[0]
      transition.update({ "label": "LE_" + transition["label"] }) 

      singlet_number = transition["state1"]
      triplet_number = transition["state2"]

      singlet_label = system['states_list'][singlet_number]["label"]
      triplet_label = system['states_list'][triplet_number]["label"]

      print("")
      print(''.center(50, '-'))
      print("{:<20} {:<30}".format("Label: ", transition["label"]))
      print(''.center(50, '-'))
      print("{:<20} {:<30}".format("Initial state: ", "%s (%s)" % (singlet_label,singlet_number)))
      print("{:<20} {:<30}".format("Target state: ", "%s (%s)" % (triplet_label,triplet_number)))
      print("{:<20} {:<30}".format("Energy (Ha): ", "{:.4e}".format(transition["transition_energy"])))
      print(''.center(50, '-'))

      transitions_list.append(transition)

    else:
      print("\nThe transition between the singlet and triplet pair with the lowest transition energy (%s) is already included so only four transitions will be considered." % fifth_transition[0]["label"])

    return transitions_list
