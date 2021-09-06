#!/usr/bin/env python3

################################################################################################################################################
##                                      Errors and exceptions of Ab Initio Input Builder & Job Launcher                                       ##
##                                                                                                                                            ##
##                            This script contains all the custom exceptions and functions built to handle errors                             ##
##                               of the Ab Initio Input Builder & Job Launcher python script and its subscripts.                              ##
################################################################################################################################################

import os

# =================================================================== #
# =================================================================== #
#                        EXCEPTIONS DEFINITIONS                       #
# =================================================================== #
# =================================================================== #

class Error(Exception):
    """Base class for exceptions in this module."""
    pass

class AbinError(Error):
    """Exception raised for errors specific to certain instructions in ABIN LAUNCHER and its subscripts.

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
    """Checks if a path towards a file or directory exists and is of the correct type. If it is, the function returns its absolute version.

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
    AbinError
        If the type does not match what is given in the path, or if the path does not exist.
    """

    # Check "type" argument

    if type not in ["file","directory","either"]:
      raise ValueError ("The specified type for which the check_abspath function has been called is not one of 'file', 'directory' or 'either'")

    # Prepare to print a helpful error message in case of problem with the given path

    msg = "\nSomething went wrong when checking the path " + path + "\nContext: " + context + "\n"

    # Check path
    
    if not os.path.exists(path):
      raise AbinError (msg + "ERROR: %s does not seem to exist." % path)
    elif type == "file":
      if not os.path.isfile(path):
        raise AbinError (msg + "ERROR: %s is not a file" % path)
    elif type == "directory":
      if not os.path.isdir(path):
        raise AbinError (msg + "ERROR: %s is not a directory" % path)
    elif type == "either":
      if not os.path.isdir(path) and not os.path.isfile(path):
        raise AbinError (msg + "ERROR: %s is neither a file nor a directory" % path)

    # If everything went well, get the normalized absolute version of the path
    
    abspath = os.path.abspath(path)

    return abspath

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
    AbinError
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
          raise AbinError ('\nContext: %s \nERROR: There is no defined "%s" key.' % (context, key))

    # Case 2: dicts is a list of dictionaries

    elif isinstance(dicts, list):
      for single_dict in dicts:
        idx = dicts.index(single_dict) + 1
        if not isinstance(single_dict, dict):
          raise AbinError ('\nContext: %s \nERROR: The %s%s item in the list is not a dictionary.' % (context, idx, ("th" if not idx in special_numbers else special_numbers[idx]))) 
        for key in keys:
          if key not in single_dict:
            raise AbinError ('\nContext: %s \nERROR: There is no defined "%s" key for the %s%s dictionary in the list.' % (context, key, idx, ("th" if not idx in special_numbers else special_numbers[idx])))

    # If dicts is neither a dictionary nor a list, raise an exception

    else:
      raise ValueError ("The type of the 'dicts' argument with which the check_keys function has been called is neither a dictionary nor a list.")