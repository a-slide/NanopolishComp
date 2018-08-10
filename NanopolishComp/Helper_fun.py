# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
import sys

#~~~~~~~~~~~~~~FUNCTIONS~~~~~~~~~~~~~~#

def to_string (*args, sep=" "):
    """make a string from a list of diverse element"""
    return sep.join (str(a) for a in args) + "\n"

def stderr_print (*args, sep=" "):
    """reproduce print with stderr.write
    """
    sys.stderr.write (to_string (*args, sep))

def stdout_print (*args, sep=" "):
    """reproduce print with stdout.write
    """
    sys.stdout.write (to_string (*args, sep))
