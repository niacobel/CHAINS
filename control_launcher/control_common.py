#!/usr/bin/env python3

################################################################################################################################################
##                                   Common functions and exceptions of QOCT-GRAD Input Builder & Job Launcher                                ##
##                                                                                                                                            ##
##                                 This script contains all the custom exceptions and functions common to the                                 ##
##                                   QOCT-GRAD Input Builder & Job Launcher python script and its subscripts.                                 ##
################################################################################################################################################

import math
import os

import numpy as np

# =================================================================== #
# =================================================================== #
#                        EXCEPTIONS DEFINITIONS                       #
# =================================================================== #
# =================================================================== #

class Error(Exception):
    """Base class for exceptions in this module."""
    pass

#######################################################################

class ControlError(Error):
    """Exception raised for errors specific to certain instructions in CONTROL LAUNCHER and its subscripts.

    Attributes
    ----------
    message : str
        Proper error message explaining the error.
    """

    def __init__(self, message):
        self.message = message

# =================================================================== #
# =================================================================== #
#                        FUNCTIONS DEFINITIONS                        #
# =================================================================== #
# =================================================================== #

def check_abspath(path:str,context:str,type="either"):
    """Checks if a path towards a file or directory exists and is of the correct type. If it's a path towards a file, the function also checks that the file is not empty. If all goes well, the function then returns the absolute version of the path.

    Parameters
    ----------
    path : str
        The path towards the file or directory you want to test.
    context : str
        Message to show on screen to give more information in case of an exception (e.g. the role of the directory or file that was checked, where the checked path was given, etc.).
    type : str, optional
        The type of element for which you would like to test the path ('file', 'directory' or 'either').
        By default, checks if the path leads to either a file or a directory (type = 'either').
    
    Returns
    -------
    abspath : str
        Normalized absolute version of the path.

    Raises
    ------
    ValueError
        If the specified type when calling the function is not 'file', 'directory' or 'either'.
    ControlError
        If the type does not match what is given in the path, or if the path does not exist, or it's an empty file.
    """

    # Check "type" argument

    if type not in ["file","directory","either"]:
      raise ValueError ("The specified type for which the check_abspath function has been called is not one of 'file', 'directory' or 'either'")

    # Prepare to print a helpful error message in case of problem with the given path

    msg = "\nSomething went wrong when checking the path " + path + "\nContext: " + context + "\n"

    # Check path
    
    if not os.path.exists(path):
      raise ControlError (msg + "ERROR: %s does not seem to exist." % path)
    elif type == "file":
      if not os.path.isfile(path):
        raise ControlError (msg + "ERROR: %s is not a file" % path)
      elif os.stat(path).st_size == 0:
        raise ControlError (msg + "ERROR: %s is an empty file" % path)
    elif type == "directory":
      if not os.path.isdir(path):
        raise ControlError (msg + "ERROR: %s is not a directory" % path)
    elif type == "either":
      if not os.path.isdir(path) and not os.path.isfile(path):
        raise ControlError (msg + "ERROR: %s is neither a file nor a directory" % path)
      elif os.path.isfile(path) and os.stat(path).st_size == 0:
        raise ControlError (msg + "ERROR: %s is an empty file" % path)

    # If everything went well, get the normalized absolute version of the path
    
    abspath = os.path.abspath(path)

    return abspath

#######################################################################

def check_keys(keys:list,dicts,context:str):
    """Checks if a list of keys is present in a dictionary or list of dictionaries. In the latter case, this function also checks that each item of the list is indeed a dictionary.

    Parameters
    ----------
    keys : list
        The list of required keys.
    dicts : list or dict
        The dictionary or list of dictionaries that need to be checked.
    context : str
        Message describing the dictionary or dictionaries being checked (shown on screen in case of a missing key).

    Raises
    ------
    ControlError
        If one of the keys is missing or if one of the items of the 'dicts' list is not a dictionary.
    ValueError
        If dicts is neither a dictionary nor a list.
    """

    # Define a dictionary for correct English spelling during printing

    special_numbers = {1:"st", 2:"nd", 3:"rd"}

    # Case 1: dicts is a dictionary

    if isinstance(dicts, dict):
      for key in keys:
        if key not in dicts:
          raise ControlError ('\nContext: %s \nERROR: There is no defined "%s" key.' % (context, key))

    # Case 2: dicts is a list of dictionaries

    elif isinstance(dicts, list):
      for single_dict in dicts:
        idx = dicts.index(single_dict) + 1
        if not isinstance(single_dict, dict):
          raise ControlError ('\nContext: %s \nERROR: The %s%s item in the list is not a dictionary.' % (context, idx, ("th" if not idx in special_numbers else special_numbers[idx]))) 
        for key in keys:
          if key not in single_dict:
            raise ControlError ('\nContext: %s \nERROR: There is no defined "%s" key for the %s%s dictionary in the list.' % (context, key, idx, ("th" if not idx in special_numbers else special_numbers[idx])))

    # If dicts is neither a dictionary nor a list, raise an exception

    else:
      raise ValueError ("The type of the 'dicts' argument with which the check_keys function has been called is neither a dictionary nor a list.")

#######################################################################

def is_consecutive(numbers:list, context:str):
    """Check if a list of numbers contains consecutive values or not and if those values start at 0.

    Parameters
    ----------
    numbers : list
        List of numbers that need to be checked.
    context : str
        Message describing the list of numbers being checked (shown on screen in case of an exception).
    
    Raises
    ------
    ControlError
        If the values are not consecutive, or if the lowest value is not 0.

    """

    min_val = min(numbers)
    max_val = max(numbers)

    if min_val != 0:
      raise ControlError ('\nContext: %s \nList: %s \nERROR: The minimum value of this list (%s) is not 0.' % (context, numbers, min_val))

    if sorted(numbers) != list(range(min_val, max_val+1)):
      raise ControlError ('\nContext: %s \nSorted list: %s \nERROR: The values are not consecutive.' % (context, sorted(numbers)))      

#######################################################################

def energy_unit_conversion(value:float,init:str,target:str) -> float:
    """|  Converts an energy value from an initial unit to a target unit by using atomic units of energy (Hartree) as an intermediary.
    |  Currently supported units: Hartree, cm\ :sup:`-1`\ , eV, nm, Hz and Joules

    Parameters
    ----------
    value : float
        The energy value we need to convert.
    init : str
        The unit of the value we need to convert.
    target : str
        The unit we must convert the value to.
    
    Returns
    -------
    conv_value : float
        The converted energy value.

    Raises
    ------
    ControlError
        If either the initial or target unit are not supported.
    """

    # Define the dictionary of conversion factors, from atomic units (Hartree) to any unit you want. - Taken from the NIST website (https://physics.nist.gov/)

    conv_factors = {
      # 1 Hartree equals:
      "Ha" : 1,
      "cm-1" : 2.1947463136320e+05,
      "eV" : 27.211386245988,
      "nm" : 45.56337117,
      "Hz" : 6.579683920502e+15,
      "J" : 4.3597447222071e-18
      }

    # Put everything in lower cases, to make it case insensitive

    init_low = init.lower()
    target_low = target.lower()
    conv_factors_low = dict((key.lower(), value) for key, value in conv_factors.items())

    # Check if the desired units are supported

    if init_low not in conv_factors_low.keys():
      raise ControlError ("ERROR: The unit of the value you want to convert (%s) is currently not supported. Supported values include: %s" % (init, ', '.join(unit for unit in conv_factors.keys())))
    elif target_low not in conv_factors_low.keys():
      raise ControlError ("ERROR: The unit you want to convert the value to (%s) is currently not supported. Supported values include: %s" % (target, ', '.join(unit for unit in conv_factors.keys())))
    
    # Convert the value

    if target_low != 'nm' and init_low != 'nm':
      conv_value = (value / conv_factors_low[init_low]) * conv_factors_low[target_low] 
    elif target_low == 'nm' and value != 0:
      conv_value = conv_factors_low['nm'] / (value / conv_factors_low[init_low])
    elif init_low == 'nm' and value != 0:
      conv_value = (conv_factors_low['nm'] / value) * conv_factors_low[target_low]
    elif (target_low == 'nm' or init_low == 'nm') and value == 0:
      conv_value = 0

    return conv_value

#######################################################################

def is_indistinguishable(energy1:float, energy2:float):
    """Check if two states can be considered indistinguishable by comparing their energy differences (in cm-1) to the spectral resolution of the pulse shapers associated with their energy domain. The spectral resolutions are fetched from:
    - the specifications of the DAZZLER ultrafast pulse shapers provided by the FASTLITE company (https://fastlite.com/produits/dazzler-ultrafast-pulse-shaper/)
    - the Quickshaper IR provided by the PhaseTech company (http://phasetechspectroscopy.com/)
    - the compact optical parametric amplifier platform from Jakob et al. (https://doi.org/10.1364/oe.27.026979)

    Parameters
    ----------
    energy1 : float
        Energy of the first state in cm-1.
    energy2 : float
        Energy of the second state in cm-1.

    Returns
    -------
    bool
        "True" if the two states are indistinguishable, "False" otherwise.
    """  

    # Define the different energy domains and their associated spectral resolution

    resolutions = {
      40000: 15.99, # Dazzler Qz-250-400
      25000: 12.49, # Dazzler Qz-250-400
      20000:	8.00, # Dazzler WR-460-740
      14286:	6.12, # Dazzler WR-460-740
      11111:	6.17, # Dazzler WR-510-950
      9709:  	4.71, # Dazzler UHR-900-1700
      6452:	  4.16, # Dazzler UHR-900-1700
      4762:	  4.53, # Dazzler WB45-1450-3000
      3846:	  4.43, # Dazzler WB45-1450-3000
      2857:	  5.70, # Dazzler WB45-2000-3700
      1818:	  5.00, # Phasetech Quickshape IR
      926:    1.97  # Jakob et al.
    }

    # Define the domain of our values (see https://www.geeksforgeeks.org/python-find-closest-number-to-k-in-given-list/)

    domains = np.asarray(list(resolutions.keys()))

    domain1 = domains[(np.abs(domains - energy1)).argmin()]
    domain2 = domains[(np.abs(domains - energy2)).argmin()]

    # If the domains of the values differ, take the average resolution

    if domain1 != domain2:
      res = (resolutions[domain1] + resolutions[domain2]) / 2
    else:
      res = resolutions[domain1]

    # Check the energy difference of the states and compare it to the spectral resolution

    if math.isclose(energy1, energy2, abs_tol=res):
      return True
    else:
      return False
