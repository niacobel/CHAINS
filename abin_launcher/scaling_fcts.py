################################################################################################################################################
##                                                             Scaling functions                                                              ##
##                                                                                                                                            ##
##                                        This script contains the scaling functions for ABIN LAUNCHER,                                       ##
##                                consult the documentation at https://chains-ulb.readthedocs.io/ for details                                 ##
################################################################################################################################################

import abin_errors

def total_nb_elec(mendeleev:dict,file_data:dict):
    """Calculates the total number of electrons in a molecule.

    Parameters
    ----------
    mendeleev : dict
        Content of AlexGustafsson's Mendeleev Table YAML file, which can be found at https://github.com/AlexGustafsson/molecular-data.
    file_data : dict
        The extracted information of the geometry file.

    Returns
    -------
    total_elec : int
        Total number of electrons in the molecule.

    Raises
    ------
    AbinError
        If there is no atomic number defined in mendeleev for one of the constituting atoms of the molecule. *(This exception is raised in the* ``get_nb_elec_for_element`` *subfunction.)*
    """

    # Definition of the function that can fetch the number of electrons associated with each element in Gustafsson's table

    def get_nb_elec_for_element(symbol:str, mendeleev:dict):
        """Returns the number of electrons for a specific element.

        Parameters
        ----------
        symbol : str
            The atom symbol of the element.
        mendeleev : dict
            Content of AlexGustafsson's Mendeleev Table YAML file, which can be found at https://github.com/AlexGustafsson/molecular-data.

        Returns
        -------
        nb_elec : int
           Number of electrons.
        
        Raises
        ------
        AbinError
            If there is no atomic number defined in mendeleev for one of the constituting atoms of the molecule.
        """

        nb_elec = 0

        # Scan the mendeleev table and get the atomic number of our atom

        for element in mendeleev:
          if (element['symbol'] == symbol):
            nb_elec = element['number']
        
        if nb_elec == 0:
          raise abin_errors.AbinError ("ERROR: There is no atomic number defined for %s in AlexGustafsson's Mendeleev Table YAML file (mendeleev.yml)" % symbol)

        return nb_elec

    # Calculating the total number of electrons in the molecule

    total_elec = 0

    print("")
    print(''.center(69, '-'))
    print("{:<12} {:<16} {:<19} {:<22}".format('Atom Type','Atomic Number','Number of atoms','Number of electrons'))
    print(''.center(69, '-'))

    for atom,nb_atom in file_data['chemical_formula'].items():
      atomic_number = get_nb_elec_for_element(atom,mendeleev)
      subtotal_elec = nb_atom * atomic_number
      print("{:<12} {:<16} {:<19} {:<22}".format(atom, atomic_number, nb_atom, subtotal_elec))
      total_elec += subtotal_elec

    print(''.center(69, '-'))
    print("{:<29} {:<19} {:<22}".format('Total',sum(file_data['chemical_formula'].values()),total_elec))
    print(''.center(69, '-'))

    return total_elec

def total_nb_atoms(mendeleev:dict,file_data:dict):
    """Returns the total number of atoms in a molecule.

    Parameters
    ----------
    mendeleev : dict
        Content of AlexGustafsson's Mendeleev Table YAML file, which can be found at https://github.com/AlexGustafsson/molecular-data.
        Unused in this particular function.
    file_data : dict
        The extracted information of the geometry file.

    Returns
    -------
    total_atoms : int
        Total number of atoms in the molecule.
    """

    # Returns the total number of atoms in the molecule by summing all the values given in the chemical_formula dictionary

    total_atoms = sum(file_data['chemical_formula'].values())

    print("")
    print("Total number of atoms in the molecule: ",total_atoms)

    return total_atoms
