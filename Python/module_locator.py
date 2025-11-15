"""
Module locator utility for finding module paths.

This module provides utilities to determine the path of the current module,
handling both frozen (e.g., py2exe) and unfrozen execution environments.
"""

import os
import sys


def we_are_frozen() -> bool:
    """
    Check if the application is running as a frozen executable.
    
    Returns:
        True if running as frozen executable (e.g., py2exe), False otherwise
    """
    # All of the modules are built-in to the interpreter, e.g., by py2exe
    return hasattr(sys, "frozen")


def module_path() -> str:
    """
    Get the path to the current module.
    
    Returns:
        Path to the module directory as a string
    """
    encoding = sys.getfilesystemencoding()
    if we_are_frozen():
        return os.path.dirname(str(sys.executable, encoding))
    return os.path.dirname(str(__file__))
