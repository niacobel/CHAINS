################################################################################################################################################
##                                                            Transition functions                                                            ##
##                                                                                                                                            ##
##                                     This script contains the transition functions for CONTROL LAUNCHER,                                    ##
##                                consult the documentation at https://chains-ulb.readthedocs.io/ for details                                 ##
################################################################################################################################################

import control_errors
import numpy as np
import os

def singlet_to_triplet(system:dict,data_dir:str):
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

    def brightest_to_most_coupled(singlet_list:list,triplet_list:list,system:dict,data_dir:str):
        """Creates the transition files for the transition going from the brightest electronic singlet state of the molecule to its most coupled triplet state.

        Parameters
        ----------
        singlet_list : list
            List of the state numbers with the singlet multiplicity.
        triplet_list : list
            List of the state numbers with the triplet multiplicity.
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
        transition : dict
            Dictionary containing six keys: ``label``, ``state1``, ``state2``, ``transition_energy``, ``init_file`` and ``target_file`` where

              - ``label`` is the label of the transition, which will be used for the name of the job directories.
              - ``state1`` is the number of the initial state.
              - ``state2`` is the number of the target state.
              - ``transition_energy`` is the transition energy between the two states.
              - ``init_file`` is the name of the initial state file, minus the number at the end.
              - ``target_file`` is the name of the target state file, minus the number at the end.
        """

        # Identifying the brightest singlet

        highest_mom = 0

        for moment in system['momdip_mtx'][0]:
          if (system['momdip_mtx'][0].index(moment) in singlet_list) and (moment > highest_mom):
            singlet = system['momdip_mtx'][0].index(moment)
        
        singlet_label = system['states_list'][singlet]["label"]

        # Creating the initial population file

        #TODO

        # Identifying the triplet with the strongest SOC

        strongest_soc = 0

        for soc in system['mime'][singlet]:
          if (system['mime'][singlet].index(soc) in triplet_list) and (soc > strongest_soc):
            triplet = system['mime'][singlet].index(soc)
        
        triplet_label = system['states_list'][triplet]["label"]

        # Calculate the transition energy (in Ha)

        energy = abs(system['mime'][singlet][singlet] - system['mime'][triplet][triplet])

        # Creating the target population file

        #TODO

        # Adding the transition to the transitions list

        transition = {
          "label" : singlet_label + "_" + triplet_label,
          "state1" : singlet,
          "state2" : triplet,
          "transition_energy" : energy,
          "init_file" : ,
          "target_file" : ,
          }
        
        return transition

    # ========================================================= #
    #       Brightest singlet to its most coupled triplet       #
    # ========================================================= #

    transition = brightest_to_most_coupled(singlet_list,triplet_list,system,data_dir)

    print(transition)

    transitions_list.append(transition)

    # ========================================================= #
    #     Brightest singlet to its 2nd most coupled triplet     #
    # ========================================================= #    

    copy_triplet_list = triplet_list[:]

    copy_triplet_list.remove(transition['state2'])

    transition = brightest_to_most_coupled(singlet_list,copy_triplet_list,system,data_dir)

    print(transition)

    transitions_list.append(transition)

    # ========================================================= #
    #     2nd brightest singlet to its most coupled triplet     #
    # ========================================================= #

    copy_singlet_list = singlet_list[:]

    copy_singlet_list.remove(transition['state1'])

    transition = brightest_to_most_coupled(copy_singlet_list,triplet_list,system,data_dir)

    print(transition)

    transitions_list.append(transition)

    # ========================================================= #
    #   2nd brightest singlet to its 2nd most coupled triplet   #
    # ========================================================= #    

    copy_triplet_list = triplet_list[:]

    copy_triplet_list.remove(transition['state2'])

    transition = brightest_to_most_coupled(copy_singlet_list,copy_triplet_list,system,data_dir)

    print(transition)

    transitions_list.append(transition)

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