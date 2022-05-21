#!/usr/bin/env python3
from __future__ import annotations

import argparse

#TODO Need to add argument for custom ouput path.
def args_parser(msg) -> argparse.Namespace:
    """Custom argument parser.

    Args:
        * `msg` (str): Description help message.

    Returns:
        `argparse.Namespace`: Namespace of input arguments.
    """

    parser = argparse.ArgumentParser(description = msg, formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-f", help = "Input .vcf file.")
    parser.add_argument("-o", help = "Optional argument: Name of output file. Must be a .txt file. Default is output.txt.")
    parser.add_argument("-dup", help = "Optional argument: Find only duplicate type cnv's. Default is False.")
    parser.add_argument("-del", help = "Optional argument: Find only deletion type cnv's. Default is False.")
    parser.add_argument("-inf", help = "Optional argument: Display information about a specific cnv. Needs to be a start or end position findings.")
    parser.add_argument("-csv", help = "Optional argument: Output findings as a .csv, default is True.")
    return parser.parse_args()

def bool_parser(var: any) -> bool:
    """Check if parameter is boolean, if not, convert it to boolean.
    Args:
        * `var` (Any): variable to check for boolean.

    Raises:
        TypeError: Unable to convert to boolean.

    Returns:
        boolean: True if var is boolean, False if not.
    """

    _true = ["true", "True", "1"]
    _false = ["false", "False", "0", None]
    if type(var) == bool:
        return var
    else:
        if var in _true:
            return True
        elif var in _false:
            return False
        else:
            raise TypeError(f"{var} must be true, True, 1, False, false, 0 or None.")