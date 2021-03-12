################################################################################################################################################
##                                                            Transition functions                                                            ##
##                                                                                                                                            ##
##                                     This script contains the transition functions for CONTROL LAUNCHER,                                    ##
##                                consult the documentation at https://chains-ulb.readthedocs.io/ for details                                 ##
################################################################################################################################################

import control_errors
import numpy as np
import os

def build_transition(state1:int,state2:int,label1:str,label2:str,system:dict,data_dir:str):
    """Creates the transition dictionary and the transition files for the transition going from the first given state to the second one.

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
    data_dir : str
        Path towards the data directory.

    Returns
    -------
    transition : dict
        Dictionary containing six keys: ``label``, ``state1``, ``state2``, ``transition_energy``, ``init_file`` and ``target_file`` where

          - ``label`` is the label of the transition, which will be used for the name of the job directories.
          - ``state1`` is the number of the initial state.
          - ``state2`` is the number of the target state.
          - ``transition_energy`` is the transition energy between the two states.
          - ``init_file`` is the name of the initial state file, minus the number at the end.
          - ``target_file`` is the name of the target state file, minus the number at the end.
    """

    # Creating the initial population file

    init_file = label1 + "_"

    init_pop = np.zeros((len(system['eigenvalues']), len(system['eigenvalues'])),dtype=complex)  # Quick init of a zero-filled matrix
    
    eigenvector = system['eigenvectors'][state1]

    for state in len(system['eigenvalues']):
      init_pop[state][state] = eigenvector[state] + 0j

    with open(os.path.join(data_dir, init_file + "1"), "w") as f:
      for line in init_pop:
        for val in line:
          print('( {0.real:.2f} , {0.imag:.2f} )'.format(val), end = " ", file = f)
        print('', file = f)

    # Creating the target population file
    
    target_file = label2 + "_"

    target_pop = np.zeros((len(system['eigenvalues']), len(system['eigenvalues'])),dtype=complex)  # Quick init of a zero-filled matrix

    eigenvector = system['eigenvectors'][state2]

    for state in len(system['eigenvalues']):
      target_pop[state][state] = eigenvector[state] + 0j

    with open(os.path.join(data_dir, target_file + "1"), "w") as f:
      for line in target_pop:
        for val in line:
          print('( {0.real:.2f} , {0.imag:.2f} )'.format(val), end = " ", file = f)
        print('', file = f)

    # Calculate the transition energy (in Ha)

    energy = abs(system['mime'][state1][state1] - system['mime'][state2][state2])

    # Building the transition dictionary

    transition = {
      "label" : label1 + "_" + label2,
      "state1" : state1,
      "state2" : state2,
      "transition_energy" : energy,
      "init_file" : init_file,
      "target_file" : target_file,
      }
    
    return transition

# =================================================================== #
# =================================================================== #
#                        Transition functions                         #
# =================================================================== #
# =================================================================== #

def closest_singlet_to_triplet(system:dict,data_dir:str):
    """Creates the transition files needed by QOCT-RA for the transition between the singlet and triplet states with the lowest transition energy between themselves.

    Parameters
    ----------
    system : dict
        Information extracted by the parsing function and derived from it.
        This dictionary needs an additionnal ``states_list`` key which is a list of dictionaries containing at least three keys each: ``number``, ``multiplicity`` and ``label`` where

          - ``number`` is the number of the state, starting at 0 (which is the ground state).
          - ``multiplicity`` is the multiplicity of the state (ex: Singlet, Triplet).
          - ``label`` is the label of the state, which will be used for the label of the transitions.

    data_dir : str
        Path towards the data directory.

    Returns
    -------
    transitions_list : list
        List of dictionaries containing six keys each: ``label``, ``state1``, ``state2``, ``transition_energy``, ``init_file`` and ``target_file`` where

          - ``label`` is the label of the transition, which will be used for the name of the job directories.
          - ``state1`` is the number of the initial state.
          - ``state2`` is the number of the target state.
          - ``transition_energy`` is the transition energy between the two states.
          - ``init_file`` is the name of the initial state file, minus the number at the end.
          - ``target_file`` is the name of the target state file, minus the number at the end.
        
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
        energy = abs(system['mime'][singlet][singlet] - system['mime'][triplet][triplet])
        if energy < minimum:
          minimum = energy
          singlet_number = singlet
          triplet_number = triplet

    singlet_label = system['states_list'][singlet_number]["label"]
    triplet_label = system['states_list'][triplet_number]["label"]

    transition = build_transition(singlet_number,triplet_number,singlet_label,triplet_label,system,data_dir)

    print("\nConsidered transition: from the %s singlet to the %s triplet with a transition energy of %s Ha" % (singlet_label,triplet_label,energy))
    print(transition)

    transitions_list.append(transition)

    return transitions_list

def bright_singlets_to_coupled_triplets_and_closest(system:dict,data_dir:str):
    """Creates the transition files needed by QOCT-RA for transitions going from each of the two brightest electronic singlet states of the molecule to each of their two most coupled triplet states. For good measure, the transition between the singlet and triplet states with the lowest transition energy is also added, if it does not correspond to one of the four transitions already considered.

    Parameters
    ----------
    system : dict
        Information extracted by the parsing function and derived from it.
        This dictionary needs an additionnal ``states_list`` key which is a list of dictionaries containing at least three keys each: ``number``, ``multiplicity`` and ``label`` where

          - ``number`` is the number of the state, starting at 0 (which is the ground state).
          - ``multiplicity`` is the multiplicity of the state (ex: Singlet, Triplet).
          - ``label`` is the label of the state, which will be used for the label of the transitions.

    data_dir : str
        Path towards the data directory.

    Returns
    -------
    transitions_list : list
        List of dictionaries containing six keys each: ``label``, ``state1``, ``state2``, ``transition_energy``, ``init_file`` and ``target_file`` where

          - ``label`` is the label of the transition, which will be used for the name of the job directories.
          - ``state1`` is the number of the initial state.
          - ``state2`` is the number of the target state.
          - ``transition_energy`` is the transition energy between the two states.
          - ``init_file`` is the name of the initial state file, minus the number at the end.
          - ``target_file`` is the name of the target state file, minus the number at the end.
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

        # Identifying the brightest singlet

        highest_mom = -1

        for moment in system['momdip_mtx'][0]:
          if (system['momdip_mtx'][0].index(moment) in singlet_list) and (moment > highest_mom):
            highest_mom = moment
            singlet_number = system['momdip_mtx'][0].index(moment)

        # Identifying the triplet with the strongest SOC

        strongest_soc = -1

        for soc in system['mime'][singlet_number]:
          if (system['mime'][singlet_number].index(soc) in triplet_list) and (soc > strongest_soc):
            strongest_soc = soc
            triplet_number = system['mime'][singlet_number].index(soc)
        
        return singlet_number,triplet_number

    # ========================================================= #
    #       Brightest singlet to its most coupled triplet       #
    # ========================================================= #

    singlet_number,triplet_number = brightest_to_most_coupled(singlet_list,triplet_list,system)

    singlet_label = system['states_list'][singlet_number]["label"]
    triplet_label = system['states_list'][triplet_number]["label"]

    transition = build_transition(singlet_number,triplet_number,singlet_label,triplet_label,system,data_dir)
    transition.update({ "label": "1B1C_" + singlet_label + "_" + triplet_label }) 

    print("\nFirst transition: brightest singlet (%s) to its most coupled triplet (%s) with a transition energy of %s Ha" % (singlet_label,triplet_label,transition["transition_energy"]))
    print(transition)

    transitions_list.append(transition)

    # ========================================================= #
    #     Brightest singlet to its 2nd most coupled triplet     #
    # ========================================================= #    

    copy_triplet_list = triplet_list[:]

    copy_triplet_list.remove(triplet_number)

    singlet_number,triplet_number = brightest_to_most_coupled(singlet_list,copy_triplet_list,system)

    singlet_label = system['states_list'][singlet_number]["label"]
    triplet_label = system['states_list'][triplet_number]["label"]

    transition = build_transition(singlet_number,triplet_number,singlet_label,triplet_label,system,data_dir)
    transition.update({ "label": "1B2C_" + singlet_label + "_" + triplet_label }) 

    print("\nSecond transition: brightest singlet (%s) to its second most coupled triplet (%s) with a transition energy of %s Ha" % (singlet_label,triplet_label,transition["transition_energy"]))
    print(transition)

    transitions_list.append(transition)

    # ========================================================= #
    #     2nd brightest singlet to its most coupled triplet     #
    # ========================================================= #

    copy_singlet_list = singlet_list[:]

    copy_singlet_list.remove(singlet_number)

    singlet_number,triplet_number = brightest_to_most_coupled(copy_singlet_list,triplet_list,system)

    singlet_label = system['states_list'][singlet_number]["label"]
    triplet_label = system['states_list'][triplet_number]["label"]

    transition = build_transition(singlet_number,triplet_number,singlet_label,triplet_label,system,data_dir)
    transition.update({ "label": "2B1C_" + singlet_label + "_" + triplet_label }) 

    print("\nThird transition: second brightest singlet (%s) to its most coupled triplet (%s) with a transition energy of %s Ha" % (singlet_label,triplet_label,transition["transition_energy"]))
    print(transition)

    transitions_list.append(transition)

    # ========================================================= #
    #   2nd brightest singlet to its 2nd most coupled triplet   #
    # ========================================================= #    

    copy_triplet_list = triplet_list[:]

    copy_triplet_list.remove(triplet_number)

    singlet_number,triplet_number = brightest_to_most_coupled(copy_singlet_list,copy_triplet_list,system)

    singlet_label = system['states_list'][singlet_number]["label"]
    triplet_label = system['states_list'][triplet_number]["label"]

    transition = build_transition(singlet_number,triplet_number,singlet_label,triplet_label,system,data_dir)
    transition.update({ "label": "2B2C_" + singlet_label + "_" + triplet_label }) 

    print("\nFourth transition: second brightest singlet (%s) to its second most coupled triplet (%s) with a transition energy of %s Ha" % (singlet_label,triplet_label,transition["transition_energy"]))
    print(transition)

    transitions_list.append(transition)

    # ========================================================= #
    #    Transition between the two closest singlet-triplet     #
    # ========================================================= #    

    fifth_transition = closest_singlet_to_triplet(system,data_dir)
    # fifth_transition is a list containing only one dictionary (see closest_singlet_to_triplet function definition)
    
    if not any(transition["transition_energy"] == fifth_transition[0]["transition_energy"] for transition in transitions_list):
      fifth_transition[0].update({ "label": "LE_" + singlet_label + "_" + triplet_label }) 
      print("\nFifth transition: the singlet and triplet pair %s with the lowest transition energy of %s Ha" % (fifth_transition[0]["label"],fifth_transition[0]["transition_energy"]))
      print(fifth_transition[0])
      transitions_list.append(fifth_transition[0])
    else:
      print("\nThe transition between the singlet and triplet pair with the lowest transition energy (%s) is already included so only four transitions will be considered." % fifth_transition[0]["label"])

    #! ========================================================= #
    #!                 Temporary states list file                #
    #! ========================================================= #

    print("{:<60}".format('\nCreating states.csv file ... '), end="")
    with open(os.path.join(data_dir, "states.csv"), "w") as f:
      print("Number;Multiplicity;Energy (cm-1);Label", file = f)
      for state in system['states_list']:
        # Print every item in state, separated by ";"
        print(";".join(map(str,state.values())), file = f)
    print('%12s' % "[ DONE ]")

    return transitions_list

def proj_ground_to_triplet(system:dict,data_dir:str):
    """Creates the transition files needed by QOCT-RA for transitions going from the ground state to each individual triplet state of the molecule, using projectors.

    Parameters
    ----------
    system : dict
        Information extracted by the parsing function and derived from it.
        This dictionary needs an additionnal ``states_list`` key which is a list of dictionaries containing at least three keys each: ``number``, ``multiplicity`` and ``label`` where

          - ``number`` is the number of the state, starting at 0 (which is the ground state).
          - ``multiplicity`` is the multiplicity of the state (ex: Singlet, Triplet).
          - ``label`` is the label of the state, which will be used for the name of the projectors.

    data_dir : str
        Path towards the data directory.

    Returns
    -------
    transitions_list : list
        List of dictionaries containing three keys each: ``label``, ``init_file`` and ``target_file`` where

          - ``label`` is the label of the transition, which will be used for the name of the job directories.
          - ``init_file`` is the name of the initial state file, minus the number at the end.
          - ``target_file`` is the name of the target state file, minus the number at the end.
    """

    # ========================================================= #
    #                     Initial state file                    #
    # ========================================================= #
    
    # Initial population in the ground state

    init_file = "ground_"

    print("{:<60}".format("\nCreating the initial population file ..."), end="") 

    init_pop = np.zeros((len(system['eigenvalues']), len(system['eigenvalues'])),dtype=complex)  # Quick init of a zero-filled matrix
    
    init_pop[0,0] = 1+0j

    with open(os.path.join(data_dir, init_file + "1"), "w") as f:
      for line in init_pop:
        for val in line:
          print('( {0.real:.2f} , {0.imag:.2f} )'.format(val), end = " ", file = f)
        print('', file = f)

    print('%12s' % "[ DONE ]")

    # ========================================================= #
    #                  Final state file (dummy)                 #
    # ========================================================= #

    # Dummy file but still needed by QOCT-RA

    final_file = "final_"

    print("{:<60}".format("\nCreating the dummy final population file ..."), end="") 

    final_pop = np.zeros((len(system['eigenvalues']), len(system['eigenvalues'])),dtype=complex)  # Quick init of a zero-filled matrix

    with open(os.path.join(data_dir, final_file + "1"), "w") as f:
      for line in final_pop:
        for val in line:
          print('( {0.real:.2f} , {0.imag:.2f} )'.format(val), end = " ", file = f)
        print('', file = f)

    print('%12s' % "[ DONE ]")

    # ========================================================= #
    #                      Projector files                      #
    # ========================================================= #

    # Each projector targets a different triplet

    transitions_list = [] # List of target state files (projectors)

    for state in system['states_list']:

      if state['multiplicity'].lower() == "triplet":

        proj_file = "projector" + state['label'] + "_"
        transitions_list.append({"label" : state['label'], "init_file" : init_file, "target_file" : proj_file})

        print("{:<60}".format("\nCreating the %s file ..." % (proj_file + "1")), end="")
        
        proj = np.zeros((len(system['states_list']),len(system['states_list'])),dtype=complex)
        
        proj[state['number'],state['number']] = 1+0j
        
        with open(os.path.join(data_dir, proj_file + "1"), "w") as f:
          for line in proj:
            for val in line:
              print('( {0.real:.2f} , {0.imag:.2f} )'.format(val), end = " ", file = f)
            print('', file = f)
        
        print('%12s' % "[ DONE ]")

    #! ========================================================= #
    #!                 Temporary states list file                #
    #! ========================================================= #

    print("{:<60}".format('\nCreating states.csv file ... '), end="")
    with open(os.path.join(data_dir, "states.csv"), "w") as f:
      print("Number;Multiplicity;Energy (cm-1);Label", file = f)
      #print(";".join(map(str,state.keys())), file = f)
      for state in system['states_list']:
        # Print every item in state, separated by ";"
        #print((str(state['number'])+";"+state['multiplicity']+";"+str(state['energy'])+";"+state['label']), file = f)
        print(";".join(map(str,state.values())), file = f)
    print('%12s' % "[ DONE ]")

    """
    print("{:<60}".format('\nCreating coupling_list.csv file ... '), end="")
    with open(os.path.join(data_dir, "coupling_list.csv"), "w") as f:
      print("State 1;State 2;Energy (cm-1)", file = f)
      for line in ori_coupling_list:
        # Print every item in line, separated by ";"
        print(";".join(map(str,line)), file = f)
    print('%12s' % "[ DONE ]")

    print("{:<60}".format('\nCreating momdip_list.csv file ... '), end="")
    with open(os.path.join(data_dir, "momdip_list.csv"), "w") as f:
      print("State 1;State 2;Dipole (a.u.)", file = f)
      for line in momdip_list:
        # Print every item in line, separated by ";"
        print(";".join(map(str,line)), file = f)
    print('%12s' % "[ DONE ]")
    """

    return transitions_list