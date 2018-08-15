  #!/usr/bin/env python3
# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
import argparse
import sys

# Local imports
from NanopolishComp import __version__ as package_version
from NanopolishComp import __name__ as package_name
from NanopolishComp.Eventalign_collapse import Eventalign_collapse
from NanopolishComp.Helper_fun import stderr_print

#~~~~~~~~~~~~~~TOP LEVEL ENTRY POINT~~~~~~~~~~~~~~#
def main ():
    # Simple triage function
    try:
        args = sys.argv
        if len(args) == 1:
            raise ValueError ("Error: Missing command\n")
        if args[1] == "Eventalign_collapse":
            Eventalign_collapse_argparse ()
        elif args[1] in ["-v", "--version"]:
            stderr_print ("{} v{}\n".format(package_name, package_version))
        elif args[1] in ["-h", "--help"]:
            raise ValueError ("NanopolishComp help\n")
        else:
            raise ValueError ("Error: Invalid command '{}'\n".format(args[1]))

    except ValueError as E:
        stderr_print (E)
        stderr_print ("Usage: NanopolishComp [command] [options]\n")
        stderr_print ("Valid command:\n\t-v/--version\n\tEventalign_collapse\n")
        stderr_print ("For help on given command, type NanopolishComp [command] -h\n")
        sys.exit()

#~~~~~~~~~~~~~~SUBPROGRAMS~~~~~~~~~~~~~~#
def Eventalign_collapse_argparse ():
    # Define parser object
    parser = argparse.ArgumentParser (description=
    "Collapse the nanopolish eventalign output by kmers rather that by events.\
    kmer level statistics (mean, median, std, var) are only computed if nanopolish is run with --samples option")
    parser.prog = "NanopolishComp Eventalign_collapse"
    # Define arguments
    parser.add_argument("subprogram")
    parser.add_argument("-i", "--input_fn", default=None, help="Path to a nanopolish eventalign tsv output file. If None, read from std input")
    parser.add_argument("-o", "--output_fn", default=None, help="Path the output tsv file. If None, write to std output")
    parser.add_argument("-s", "--write_samples", default=False, action='store_true', help="If given, will write the raw sample if eventalign is run with --samples option")
    parser.add_argument("-v", "--verbose", default=False, action='store_true', help="If given will be more chatty (default = False)")
    # Parse Arguments
    a = parser.parse_args()

    # Run command
    f = Eventalign_collapse(
        input_fn = a.input_fn,
        output_fn = a.output_fn,
        verbose = a.verbose)
