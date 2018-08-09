# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
import sys

#~~~~~~~~~~~~~~FUNCTIONS~~~~~~~~~~~~~~#

def stdout_print (*args, sep="\t"):
    """reproduce print with stdout.write
    """
    sys.stdout.write (sep.join (str(a) for a in args) + "\n")

def stderr_print (*args, sep="\t"):
    """reproduce print with stderr.write
    """
    sys.stderr.write (sep.join (str(a) for a in args) + "\n")
