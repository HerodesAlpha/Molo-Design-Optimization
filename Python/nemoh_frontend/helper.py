# -*- coding: utf-8 -*-
"""
Helper module for validation and logging utilities.

This module provides helper methods for:
- Value validation (None checks, type checks, path checks)
- Logging utilities (entrance, exit, exception logging)

Updated since version 1.1:
    1. Added check_path_exist(), check_is_directory() and check_is_file().

Updated since version 1.2 (OpenWarp - Add Logging Functionality):
    Added support for logging
"""

__author__ = "Eivind Sonju"
__copyright__ = "Copyright (C) 2017-2019 Verbun AS. All rights reserved."
__version__ = "2.0"

import os
import traceback
from typing import Any, Dict, List, Optional


def check_not_none_nor_empty(val: Any, name: str) -> None:
    """
    Check if the given value is None or empty string.

    Args:
        val: The value to check
        name: Name of the value (for error messages)

    Raises:
        ValueError: If val is None or empty string
        TypeError: If val is not a string
    """
    if val is None:
        raise ValueError(f'Object {name} should not be None.')
    if not isinstance(val, str):
        raise TypeError(f'Object {name} should be a string.')
    if len(val.strip()) == 0:
        raise ValueError(f'Object {name} should not be empty string.')


def check_path_exist(val: str, name: str) -> None:
    """
    Check if the given value is a legal file or directory path.

    Args:
        val: The path to check
        name: Name of the value (for error messages)

    Raises:
        ValueError: If val is not a legal file or directory path
    """
    if not os.path.exists(val):
        raise ValueError(f'Path "{val}" does not exist.')


def check_is_directory(val: str, name: str) -> None:
    """
    Check if the given value is a legal directory path.

    Args:
        val: The path to check
        name: Name of the value (for error messages)

    Raises:
        ValueError: If val is not a legal directory path
    """
    check_path_exist(val, name)
    if not os.path.isdir(val):
        raise ValueError(f'Path "{val}" is not a legal directory.')


def check_is_file(val: str, name: str) -> None:
    """
    Check if the given value is a legal file path.

    Args:
        val: The path to check
        name: Name of the value (for error messages)

    Raises:
        ValueError: If val is not a legal file path
    """
    check_path_exist(val, name)
    if not os.path.isfile(val):
        raise ValueError(f'Path "{val}" is not a legal file.')


def check_type_value(val: Any, name: str, expected_type: type, allow_none: bool) -> None:
    """
    Check if the given value is of expected type and handle None values.

    Args:
        val: The value to check
        name: Name of the value (for error messages)
        expected_type: The expected type
        allow_none: Whether the val is allowed to be None

    Raises:
        ValueError: If val is None while not allow None
        TypeError: If val is not of expected type
    """
    if val is None and not allow_none:
        raise ValueError(f'Object {name} should not be None.')
    if not isinstance(val, expected_type):
        raise TypeError(f'Object {name} should be of {expected_type}.')


def log_entrance(logger: Any, signature: str, parasMap: Optional[Dict[str, Any]]) -> None:
    """
    Log entrance into public methods at DEBUG level.

    Args:
        logger: The logger object
        signature: The method signature
        parasMap: The passed parameters dictionary
    """
    logger.debug(f'[Entering method {signature}]')
    if parasMap is not None and len(list(parasMap.items())) > 0:
        paraStr = '[Input parameters['
        for (k, v) in list(parasMap.items()):
            paraStr += (str(k) + ':' + str(v) + ', '))
        paraStr += ']]'
        logger.debug(paraStr)


def log_exit(logger: Any, signature: str, parasList: Optional[List[Any]]) -> None:
    """
    Log exit from public methods at DEBUG level.

    Args:
        logger: The logger object
        signature: The method signature
        parasList: The objects to return
    """
    logger.debug(f'[Exiting method {signature}]')
    if parasList is not None and len(parasList) > 0:
        logger.debug(f'[Output parameter {parasList}]')


def log_exception(logger: Any, signature: str, e: Exception) -> Exception:
    """
    Log exception at ERROR level.

    Args:
        logger: The logger object
        signature: The method signature
        e: The error/exception

    Returns:
        The exception (for re-raising)
    """
    # This will log the traceback.
    logger.error(f'[Error in method {signature}: Details {e}]')
    logger.error(' Error stack:')
    logger.error(traceback.format_exc())
    return e
