#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
import argparse
import sys

# Local imports
from NanopolishComp import __version__ as package_version
from NanopolishComp import __name__ as package_name
from NanopolishComp import __description__ as package_description
from NanopolishComp.Eventalign_collapse import Eventalign_collapse
from NanopolishComp.Freq_meth_calculate import Freq_meth_calculate

#~~~~~~~~~~~~~~TOP LEVEL ENTRY POINT~~~~~~~~~~~~~~#
def main(args=None):
    """
    Main entry point for NanoSnake command line interface
    """

    # Parser and subparsers for command
    parser = argparse.ArgumentParser (description=package_description)
    parser.add_argument("--version", action="version", version="{} v{}".format(package_name, package_version))
    subparsers = parser.add_subparsers (description="%(prog)s implements the following subcommands", dest="subcommands")
    subparsers.required = True

    # Eventalign_collapse subparser
    subparser_ec = subparsers.add_parser("Eventalign_collapse", description="Collapse the nanopolish eventalign output at kmers level and compute kmer level statistics")
    subparser_ec.set_defaults(func=Eventalign_collapse_main)
    subparser_ec_io = subparser_ec.add_argument_group("Input/Output options")
    subparser_ec_io.add_argument("-i", "--input_fn", default=0, help="Path to a nanopolish eventalign tsv output file. If '0' read from std input (default: %(default)s)")
    subparser_ec_io.add_argument("-o", "--outdir", type=str, default="./", help="Path to the output folder (will be created if it does exist yet) (default: %(default)s)")
    subparser_ec_io.add_argument("-p", "--outprefix", type=str, default="out", help="text outprefix for all the files generated (default: %(default)s)")
    subparser_ec_rp = subparser_ec.add_argument_group("Run parameters options")
    subparser_ec_rp.add_argument("-s", "--write_samples", default=False, action='store_true', help="If given, will write the raw sample if nanopolish eventalign was ran with --samples option (default: %(default)s)")
    subparser_ec_rp.add_argument("-r", "--max_reads", default=0 , type=int , help = "Maximum number of read to parse. 0 to deactivate (default: %(default)s)")
    subparser_ec_rp.add_argument("-f", "--stat_fields", default=["mean", "median", "num_signals"], type=str, nargs='+', help = "List of statistical fields to compute if nanopolish eventalign was ran with --sample option. Valid values = mean, std, median, mad, num_signals (default: %(default)s)")
    subparser_ec_other = subparser_ec.add_argument_group("Other options")
    subparser_ec_other.add_argument("-t", "--threads", default=4, type=int, help="Total number of threads. 1 thread is used for the reader and 1 for the writer (default: %(default)s)")

    # Freq_meth_calculate subparser
    subparser_fm = subparsers.add_parser("Freq_meth_calculate", description="Calculate methylation frequency at genomic CpG sites from the output of nanopolish call-methylation")
    subparser_fm.set_defaults(func=Freq_meth_calculate_main)
    subparser_fm_io = subparser_fm.add_argument_group("Input/Output options")
    subparser_fm_io.add_argument("-i", "--input_fn", default=0, help="Path to a nanopolish call_methylation tsv output file. If not specified read from std input")
    subparser_fm_io.add_argument("-b", "--output_bed_fn", type=str, default="", help="Path to write a summary result file in BED format (default: %(default)s)")
    subparser_fm_io.add_argument("-t", "--output_tsv_fn", type=str, default="", help="Path to write an more extensive result report in TSV format (default: %(default)s)")
    subparser_fm_fo = subparser_fm.add_argument_group("Filtering options")
    subparser_fm_fo.add_argument("-l", "--min_llr", type=float, default=2.5, help="Log likelihood ratio threshold (default: %(default)s)")
    subparser_fm_fo.add_argument("-d", "--min_depth", type=int, default=10, help="Minimal number of reads covering a site to be reported (default: %(default)s)")
    subparser_fm_fo.add_argument("-f", "--min_meth_freq", type=float, default=0.05, help="Minimal methylation frequency of a site to be reported (default: %(default)s)")

    # Add common group parsers
    for sp in [subparser_ec, subparser_fm]:
        sp_verbosity = sp.add_mutually_exclusive_group()
        sp_verbosity.add_argument("-v", "--verbose", action="store_true", default=False, help="Increase verbosity (default: %(default)s)")
        sp_verbosity.add_argument("-q", "--quiet", action="store_true", default=False, help="Reduce verbosity (default: %(default)s)")

    # Parse args and call subfunction
    args = parser.parse_args()
    args.func(args)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~SUBPARSERS FUNCTIONS~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

def Eventalign_collapse_main (args):
    """"""
    # Run corresponding class
    Eventalign_collapse (
        input_fn = args.input_fn,
        outdir = args.outdir,
        outprefix = args.outprefix,
        write_samples = args.write_samples,
        stat_fields= args.stat_fields,
        threads = args.threads,
        verbose = args.verbose,
        quiet = args.quiet)

def Freq_meth_calculate_main (args):
    """"""
    # Run corresponding class
    Freq_meth_calculate (
        input_fn = args.input_fn,
        output_bed_fn = args.output_bed_fn,
        output_tsv_fn = args.output_tsv_fn,
        min_llr = args.min_llr,
        min_depth = args.min_depth,
        min_meth_freq = args.min_meth_freq,
        verbose = args.verbose,
        quiet = args.quiet)
