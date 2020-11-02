################################################################################################################################################
##                                                                Molecule Scanner                                                            ##
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
        Content of the XYZ geometry file.

    Returns
    -------
    file_data : dict
        The extracted informations of the file, following the pattern { 'chemical_formula' : { }, 'atomic_coordinates' : [ ] }

    Raises
    ------
    AbinError
        If the number of atomic coordinates lines does not match the number of atoms mentioned in the first line of the .xyz file.
    """

    file_data = {'chemical_formula':{}, 'atomic_coordinates':[]}
    
    # Determining the number of atoms (first line of the xyz file)
    
    nb_atoms = int(mol_content[0])

    # Scanning the content of the XYZ file to determine the chemical formula and atomic coordinates of the molecule
      
    lines_rx = {
        # Pattern for finding the atomic coordinates lines
        # They look like 'Si   -0.31438   1.89081   0.00000' 
        # This uses regex, for more information, see https://docs.python.org/3/library/re.html)
        'atomLine': re.compile(
            r'^\s{0,4}(?P<atomSymbol>[a-zA-Z]{0,3})\s+[-]?\d+\.\d+\s+[-]?\d+\.\d+\s+[-]?\d+\.\d+$')
    }
    
    checksum_nlines = 0                                                 # This variable will be used to check if the number of coordinate lines matches the number of atoms of the molecule 
    
    for line in mol_content[2:]:                                        # We only start at the 3rd line because the first two won't contain any coordinates
      m = lines_rx['atomLine'].match(line)
      if m is not None:                                                 # We only care if the line looks like an atom coordinates
        checksum_nlines += 1
        file_data['atomic_coordinates'].append(line)                    # All coordinates will be stored in this variable to be rendered in the input file later on
        if m.group("atomSymbol") not in file_data['chemical_formula']:
          file_data['chemical_formula'][m.group("atomSymbol")] = 1
        else:
          file_data['chemical_formula'][m.group("atomSymbol")] += 1
                
    # Check if the number of lines matches the number of atoms defined in the first line of the .xyz file
    
    if checksum_nlines != nb_atoms:
      raise abin_errors.AbinError("ERROR: Number of atomic coordinates lines (%s) doesn't match the number of atoms mentioned in the first line of the .xyz file (%s) !" % (checksum_nlines, nb_atoms))
  
    return file_data