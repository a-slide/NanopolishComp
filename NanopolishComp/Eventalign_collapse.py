# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
import sys
from time import time

# Third party imports
import numpy as np

# Local imports
from NanopolishComp.Helper_fun import stdout_print, stderr_print, to_string

#~~~~~~~~~~~~~~FUNCTIONS~~~~~~~~~~~~~~#
def Eventalign_collapse (input_fn=None, output_fn=None, write_samples=False, skip_model_N=False, verbose=False):

    read_id_set = set()
    nkmers = nevents = 0

    if verbose:
        stderr_print ("Define input and output")
    input = open (input_fn, "r") if input_fn else sys.stdin
    output = open (output_fn, "w") if output_fn else sys.stdout

    if verbose:
        stderr_print ("Parse file")
    try:

        # Parse file header to detemine the option used by nanopolish
        ls = input.readline().rstrip().split("\t")
        header_list = ["ref_name", "ref_pos", "ref_kmer", "read_name", "n_events"]
        idx_agregate = False
        sample_agregate = False

        if "start_idx" in ls and "end_idx" in ls:
            if verbose:
                stderr_print ("\tFound signal index in nanopolish eventalign header")
            start_pos = ls.index("start_idx")
            end_pos = ls.index("end_idx")
            idx_agregate = True
            header_list.extend (["start_idx", "end_idx"])

        if "samples" in ls:
            if verbose:
                stderr_print ("\tFound samples in nanopolish eventalign header")
            sample_agregate = True
            sample_pos = ls.index("samples")
            header_list.extend (["mean", "median", "std", "var"])
            if write_samples:
                header_list.append ("samples")

        output.write (to_string (*header_list, sep="\t"))

        # First line exception
        n_events = 1
        ls = input.readline().rstrip().split("\t")
        ref_name, ref_pos, ref_kmer, read_name = ls[0], ls[1], ls[2], ls[3]

        if idx_agregate:
            start_idx, end_idx = ls[start_pos], ls[end_pos]
        if sample_agregate:
            raw_list = ls[sample_pos].split(",")

        # Iterate over all lines
        for line in input:

            # Extract important fields from the file
            ls = line.rstrip().split("\t")
            c_ref_name, c_ref_pos, c_ref_kmer, c_read_name = ls[0], ls[1], ls[2], ls[3]

            if idx_agregate:
                c_start_idx, c_end_idx = ls[start_pos], ls[end_pos]
            if sample_agregate:
                c_raw_list = ls[sample_pos].split(",")

            # Update values if same position
            if c_ref_name == ref_name and c_ref_pos == ref_pos:
                n_events += 1
                if idx_agregate:
                    start_idx = c_start_idx
                if sample_agregate:
                    raw_list.extend (c_raw_list) ######## Might have to reverse the order of the list but not very important for statistics

            # Write new kmer
            else:
                res_list = [ref_name, ref_pos, ref_kmer, read_name, n_events]
                if idx_agregate:
                    res_list.extend ([start_idx, end_idx])
                if sample_agregate:
                    mean, median, std, var = raw_stat (raw_list)
                    res_list.extend ([mean, median, std, var])
                    if write_samples:
                        res_list.append (";".join(raw_list))

                output.write (to_string (*res_list, sep="\t"))

                # Update Counters
                read_id_set.add (read_name)
                nkmers += 1

                # Initialise a new kmer
                n_events = 1
                ref_name, ref_pos, ref_kmer, read_name = c_ref_name, c_ref_pos, c_ref_kmer, c_read_name
                if idx_agregate:
                    start_idx, end_idx = c_start_idx, c_end_idx
                if sample_agregate:
                    raw_list = c_raw_list

        # Last line exception
        res_list = [ref_name, ref_pos, ref_kmer, read_name, n_events]
        if idx_agregate:
            res_list.extend ([start_idx, end_idx])
        if sample_agregate:
            mean, median, std, var = raw_stat (raw_list)
            res_list.extend ([mean, median, std, var])
            if write_samples:
                res_list.append (";".join(raw_list))

        output.write (to_string (*res_list, sep="\t"))

        # Update Counters
        read_id_set.add (read_name)
        nkmers += 1

    except (BrokenPipeError, IOError, KeyboardInterrupt) as E:
        print (E)
        pass

    # Print final counts
    stderr_print ("[NanopolishComp summary] Reads:{:,}\tKmers:{:,}".format (len(read_id_set), nkmers))

    # Close files
    if input_fn:
        input.close()
    if output_fn:
        output.close()

def raw_stat (raw_list):

    l = np.array (raw_list, dtype=np.float32)
    mean = round (np.mean (l), 3)
    median = round (np.median (l), 3)
    std = round (np.std (l), 3)
    var = round (np.var (l), 3)

    return mean, median, std, var
