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