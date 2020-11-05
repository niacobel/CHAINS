################################################################################################################################################
##                                                             Geometry Files Scanner                                                         ##
##                                                                                                                                            ##
##                                        This script contains the scanning functions for ABIN LAUNCHER,                                      ##
##                                consult the documentation at https://chains-ulb.readthedocs.io/ for details                                 ##
################################################################################################################################################

import re
import abin_errors

def xyz_scan(mol_content:list):
    """Scans the content of an XYZ geometry file and extracts the chemical formula and atomic coordinates of the molecule.

    Parameters
    ----------
    mol_content : list
        Content of the XYZ geometry file. Each element of the list is a line of the file.

    Returns
    -------
    file_data : dict
        The extracted informations of the file, following the pattern { 'chemical_formula' : { }, 'atomic_coordinates' : [ ] }

    Raises
    ------
    AbinError
        If the number of atomic coordinates lines does not match the number of atoms mentioned in the first line of the .xyz file.
    """

    # Initialize the file_data dictionary that will be returned by the function

    file_data = {'chemical_formula':{}, 'atomic_coordinates':[]}
    
    # Determining the number of atoms (first line of the xyz file)
    
    nb_atoms = int(mol_content[0])

    # Initialize a variable will be used to check if the number of coordinate lines matches the number of atoms of the molecule

    checksum_nlines = 0 

    # Define the pattern of the atomic coordinates lines (They look like 'Si   -0.31438   1.89081   0.00000')
    # This is based on regular expressions (regex), consult https://docs.python.org/3/library/re.html for details
    # You can also paste everything inside the raw string (r'<here>') on https://regex101.com for an explanation of this particular regex (use your .xyz file as a test string on the site)

    pattern = re.compile(r'^\s*(?P<atomSymbol>[a-zA-Z]{1,3})(?:\s+-?\d+\.\d+){3}\s*$')

    # Scanning the content of the XYZ file to determine the chemical formula and atomic coordinates of the molecule
    # We only start at the 3rd line ([2:]) because the first two won't contain any coordinates
    
    for line in mol_content[2:]:                                        
      
      matching_line = pattern.match(line)

      # If the line matches our pattern

      if matching_line is not None:
        checksum_nlines += 1

        # Store the line in the 'atomic_coordinates' key to be rendered in the input file later on

        file_data['atomic_coordinates'].append(line)

        # Count the number of occurrences of the atom type

        atom_type = matching_line.group("atomSymbol")

        if atom_type not in file_data['chemical_formula']:
          file_data['chemical_formula'][atom_type] = 1
        else:
          file_data['chemical_formula'][atom_type] += 1

    # Check if the number of lines matches the number of atoms defined in the first line of the .xyz file
    
    if checksum_nlines != nb_atoms:
      raise abin_errors.AbinError("ERROR: Number of atomic coordinates lines (%s) doesn't match the number of atoms mentioned in the first line of the .xyz file (%s) !" % (checksum_nlines, nb_atoms))
  
    # Scanning complete, now return file_data

    return file_data