# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
import sys
import os
from collections import *
import logging

#~~~~~~~~~~~~~~FUNCTIONS~~~~~~~~~~~~~~#

def stderr_print (*args):
    """reproduce print with stderr.write"""
    sys.stderr.write(" ".join(str(a) for a in args))
    sys.stderr.flush()

def file_readable (fn, **kwargs):
    """Check if the file is readable"""
    return os.path.isfile (fn) and os.access (fn, os.R_OK)

def dir_writable (fn, **kwargs):
    """Check if the file is readable"""
    if not os.path.isdir(fn):
        fn = os.path.dirname(fn)
    return os.path.dirname(fn) and os.access (fn, os.W_OK)

def mkdir (fn, exist_ok=False):
    """ Create directory recursivelly. Raise IO error if path exist or if error at creation """
    try:
        os.makedirs (fn, exist_ok=exist_ok)
    except:
        raise NanopolishCompError ("Error creating output folder `{}`".format(fn))

def numeric_cast_dict (d):
    """Cast str values to integer or float from a dict """
    for k, v in d.items():
        d[k] = numeric_cast(v)
    return d

def numeric_cast (v):
    """Try to cast values to int or to float"""
    if type(v)== str:
        try:
            v = int(v)
        except ValueError:
            try:
                v = float(v)
            except ValueError:
                pass
    return v

def dict_to_str (d, sep="\t", nsep=0, exclude_list=[]):
    """ Transform a multilevel dict to a tabulated str """
    m = ""

    if isinstance(d, Counter):
        for i, j in d.most_common():
            if not i in exclude_list:
                m += "{}{}: {:,}\n".format(sep*nsep, i, j)

    else:
        for i, j in d.items():
            if not i in exclude_list:
                if isinstance(j, dict):
                    j = dict_to_str(j, sep=sep, nsep=ntab+1)
                    m += "{}{}\n{}".format(sep*nsep, i, j)
                else:
                    m += "{}{}: {}\n".format(sep*nsep, i, j)
    return m[:-1]

def jhelp (f:"python function or method"):
    """
    Display a Markdown pretty help message for functions and class methods (default __init__ is a class is passed)
    jhelp also display default values and type annotations if available.
    Undocumented options are not displayed.
    The docstring synthax should follow the markdown formated convention below
    * f
        Function or method to display the help message for
    """
    # For some reason signature is not always importable. In these cases the build-in help is called instead
    try:
        from IPython.core.display import display, Markdown, HTML
        import inspect
    except (NameError, ImportError) as E:
        NanocomporeWarning ("jupyter notebook is required to use this function. Please verify your dependencies")
        help(f)
        return

    if inspect.isclass(f):
        f = f.__init__

    if inspect.isfunction(f) or inspect.ismethod(f):

        # Parse arguments default values and annotations
        sig_dict = OrderedDict()
        for name, p in inspect.signature(f).parameters.items():
            sig_dict[p.name] = []
            # Get Annotation
            if p.annotation != inspect._empty:
                sig_dict[p.name].append(": {}".format(p.annotation))
            # Get default value if available
            if p.default == inspect._empty:
                sig_dict[p.name].append("(required)")
            else:
                sig_dict[p.name].append("(default = {})".format(p.default))

        # Parse the docstring
        doc_dict = OrderedDict()
        descr = []
        lab=None
        for l in inspect.getdoc(f).split("\n"):
            l = l.strip()
            if l:
                if l.startswith("*"):
                    lab = l[1:].strip()
                    doc_dict[lab] = []
                elif lab:
                    doc_dict[lab].append(l)
                else:
                    descr.append(l)

        # Reformat collected information in Markdown synthax
        s = "---\n\n**{}.{}**\n\n{}\n\n---\n\n".format(f.__module__, f.__name__, " ".join(descr))
        for k, v in doc_dict.items():
            s+="* **{}** *{}*\n\n{}\n\n".format(k, " ".join(sig_dict[k]), " ".join(v))

        # Display in Jupyter
        display (Markdown(s))

def head (fp, n=10, sep="\t", comment=None):
    """
    Emulate linux head cmd. Handle gziped files and bam files
    * fp
        Path to the file to be parse.
    * n
        Number of lines to print starting from the begining of the file (Default 10)
    """
    line_list = []

    # Get lines
    try:
        with open(fp) as fh:
            line_num = 0
            while (line_num < n):
                l= next(fh).strip()
                if comment and l.startswith(comment):
                    continue
                if sep:
                    line_list.append (l.split(sep))
                else:
                    line_list.append (l)
                line_num+=1

    except StopIteration:
        pass

    # Add padding if sep given
    if sep:
        try:
            # Find longest elem per col
            col_len_list = [0 for _ in range (len(line_list[0]))]
            for ls in line_list:
                for i in range (len(ls)):
                    len_col = len(ls[i])
                    if len_col > col_len_list[i]:
                        col_len_list[i] = len_col

            # Add padding
            line_list_tab = []
            for ls in line_list:
                s = ""
                for i in range (len(ls)):
                    len_col = col_len_list[i]
                    len_cur_col = len(ls[i])
                    s += ls[i][0:len_col] + " "*(len_col-len_cur_col)+" "
                line_list_tab.append(s)
            line_list = line_list_tab

        # Fall back to non tabulated display
        except IndexError:
            return head (fp=fp, n=n, sep=None)

    for l in line_list:
        print (l)
    print()

def get_logger (name=None, verbose=False, quiet=False):
    """Set logger to appropriate log level"""

    logging.basicConfig(format='%(message)s')
    logger = logging.getLogger(name)

    # Define overall verbose level
    if verbose:
        logger.setLevel(logging.DEBUG)
    elif quiet:
        logger.setLevel(logging.WARNING)
    else:
        logger.setLevel(logging.INFO)

    return logger

class LineParser():
    """Simple line parser which returns namedtuples per line and does numeric
    type casting for file lines"""

    def __init__(self, header_line, sep="\t", cast_numeric_field=True):
        """"""
        self.cast_numeric_field = cast_numeric_field
        self.sep = sep
        self.header = header_line.strip().split(self.sep)
        self.line_tuple = namedtuple("Line", self.header)
        self.c = Counter()

    def __repr__ (self):
        m = "Header names: {}\n".format(self.header)
        m+="Line counter:{}".format(dict_to_str(self.c))
        return m

    def __call__(self, line):
        """"""
        line = line.strip().split(self.sep)

        # Try to cast numeric field to appropiate type
        if self.cast_numeric_field:
            for i in range(len(line)):
                line[i] = self._numeric_cast(line[i])

        if len(line) != len(self.header):
            self.c["Inconsistent field numbers"]+=1
            return None
        else:
            # Cast line_tuple object
            line = self.line_tuple(*line)
            self.c["Valid line parsed"]+=1
            return line

    def _numeric_cast(self, val):
        """Try to cast values to int or to float"""
        if type(val)== str:
            try:
                val = int(val)
            except ValueError:
                try:
                    val = float(val)
                except ValueError:
                    pass
        return val


#~~~~~~~~~~~~~~CUSTOM EXCEPTION AND WARN CLASSES~~~~~~~~~~~~~~#
class NanopolishCompError (Exception):
    """ Basic exception class for NanopolishComp package """
    pass
