# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
import sys
import os

#~~~~~~~~~~~~~~FUNCTIONS~~~~~~~~~~~~~~#

def stderr_print (*args):
    """reproduce print with stderr.write
    """
    sys.stderr.write(" ".join(str(a) for a in args))
    sys.stderr.flush()

def access_file (fn, **kwargs):
    """Check if the file is readable
    """
    return os.path.isfile (fn) and os.access (fn, os.R_OK)

#~~~~~~~~~~~~~~CUSTOM EXCEPTION AND WARN CLASSES~~~~~~~~~~~~~~#
class NanopolishCompError (Exception):
    """ Basic exception class for NanopolishComp package """
    pass
