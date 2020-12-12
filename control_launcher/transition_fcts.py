################################################################################################################################################
##                                                            Transition functions                                                            ##
##                                                                                                                                            ##
##                                     This script contains the transition functions for CONTROL LAUNCHER,                                    ##
##                                consult the documentation at https://chains-ulb.readthedocs.io/ for details                                 ##
################################################################################################################################################

import control_errors
import numpy as np
import os

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

        print("{:<60}".format("\nCreating the %s file ..." % proj_file), end="")
        
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