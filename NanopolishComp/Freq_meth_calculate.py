# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Standard library imports
from collections import OrderedDict

# Third party imports


# Local imports
from NanopolishComp.common import file_readable, dir_writable, NanopolishCompError

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~LOGGING INFO~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
import logging
logging.basicConfig(level=logging.INFO, format="%(message)s")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~MAIN CLASS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

class Freq_meth_calculate():

    def __init__ (self,
        output_fn:"str",
        input_fn:"str",
        call_threshold:"float"=2.5,
        split_groups:"bool"=False,
        verbose:"bool"=False,
        quiet:"bool"=False):
        """
        Calculate methylation frequency at genomic CpG sites from the output of nanopolish call-methylation
        * output_fn
            Path the output reference position tsv file
        * input_fn
            Path to a nanopolish call_methylation tsv output file
        * call_threshold
            Log likelihood ratio threshold
        * split_groups
             split multi-cpg groups
        * verbose
            Increase verbosity
        * quiet
            Reduce verbosity
        """

        # Define overall verbose level
        self.log = logging.getLogger()
        if verbose:
            self.log.setLevel (logging.DEBUG)
            self.verbose_level = 2
        elif quiet:
            self.log.setLevel (logging.WARNING)
            self.verbose_level = 0
        else:
            self.log.setLevel (logging.INFO)
            self.verbose_level = 1

        # Verify parameters validity
        self.log.info ("Checking arguments")
        self.log.debug("Test input file readability")
        if input_fn != 0 and not file_readable (input_fn):
            raise IOError ("Cannot read input file")
        self.log.debug("\tTesting output file writability")
        if not dir_writable (output_fn):
            raise IOError ("Cannot write output file in indicated folder. Create the output folder if it does not exist yet")

        # Save args to self values
        self.output_fn = output_fn
        self.input_fn = input_fn
        self.call_threshold = call_threshold
        self.split_groups = split_groups


#
# #! /usr/bin/env python
#
# import sys
# import csv
# import argparse
#
#
# class SiteStats:
#     def __init__(self, g_size, g_seq):
#         self.num_reads = 0
#         self.called_sites = 0
#         self.called_sites_methylated = 0
#         self.group_size = g_size
#         self.sequence = g_seq
#
# def Freq_meth_calculate
#
#
# def update_call_stats(key, num_called_cpg_sites, is_methylated, sequence):
#     if key not in sites:
#         sites[key] = SiteStats(num_called_cpg_sites, sequence)
#
#     sites[key].num_reads += 1
#     sites[key].called_sites += num_called_cpg_sites
#     if is_methylated > 0:
#         sites[key].called_sites_methylated += num_called_cpg_sites
#
# parser = argparse.ArgumentParser( description='Calculate methylation frequency at genomic CpG sites')
# parser.add_argument('-c', '--call-threshold', type=float, required=False, default=2.5)
# parser.add_argument('-i', '--input', type=str, required=False)
# parser.add_argument('-s', '--split-groups', action='store_true')
# args = parser.parse_args()
# assert(args.call_threshold is not None)
#
# sites = dict()
#
# if args.input:
#     in_fh = open(args.input)
# else:
#     in_fh = sys.stdin
# csv_reader = csv.DictReader(in_fh, delimiter='\t')
#
# for record in csv_reader:
#
#     num_sites = int(record['num_motifs'])
#     llr = float(record['log_lik_ratio'])
#
#     # Skip ambiguous call
#     if abs(llr) < args.call_threshold:
#         continue
#     sequence = record['sequence']
#
#     is_methylated = llr > 0
#
#     # if this is a multi-cpg group and split_groups is set, break up these sites
#     if args.split_groups and num_sites > 1:
#         c = str(record['chromosome'])
#         s = int(record['start'])
#         e = int(record['end'])
#
#         # find the position of the first CG dinucleotide
#         sequence = record['sequence']
#         cg_pos = sequence.find("CG")
#         first_cg_pos = cg_pos
#         while cg_pos != -1:
#             key = (c, s + cg_pos - first_cg_pos, s + cg_pos - first_cg_pos)
#             update_call_stats(key, 1, is_methylated, "split-group")
#             cg_pos = sequence.find("CG", cg_pos + 1)
#     else:
#         key = (str(record['chromosome']), int(record['start']), int(record['end']))
#         update_call_stats(key, num_sites, is_methylated, sequence)
#
# # header
# print("\t".join(["chromosome", "start", "end", "num_motifs_in_group", "called_sites", "called_sites_methylated", "methylated_frequency", "group_sequence"]))
#
# sorted_keys = sorted(sites.keys(), key = lambda x: x)
#
# for key in sorted_keys:
#     if sites[key].called_sites > 0:
#         (c, s, e) = key
#         f = float(sites[key].called_sites_methylated) / sites[key].called_sites
#         print("%s\t%s\t%s\t%d\t%d\t%d\t%.3f\t%s" % (c, s, e, sites[key].group_size, sites[key].called_sites, sites[key].called_sites_methylated, f, sites[key].sequence))
