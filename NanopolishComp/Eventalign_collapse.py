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
        header_list = ["ref_pos", "ref_kmer"]
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
            header_list.extend (["mean", "std"])
            if write_samples:
                header_list.append ("samples")

        output.write (to_string (*header_list, sep="\t"))

        # First line exception
        ls = input.readline().rstrip().split("\t")
        ref_name, ref_pos, ref_kmer, read_name = ls[0], ls[1], ls[2], ls[3]

        # Write first separator
        output.write ("#\t{}\t{}\n".format(read_name, ref_name))

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
            if c_ref_pos == ref_pos and c_ref_kmer == ref_kmer:
                if idx_agregate:
                    start_idx = c_start_idx
                if sample_agregate:
                    raw_list.extend (c_raw_list) ######## Might have to reverse the order of the list but not very important for statistics

            # Write new kmer
            else:
                res_list = [ref_pos, ref_kmer]
                if idx_agregate:
                    res_list.extend ([start_idx, end_idx])
                if sample_agregate:
                    l = np.array (raw_list, dtype=np.float32)
                    res_list.extend ([round (np.mean (l), 3), round (np.std (l), 3)])
                    if write_samples:
                        res_list.append (";".join(raw_list))

                output.write (to_string (*res_list, sep="\t"))

                # Update Counters
                read_id_set.add (read_name)
                nkmers += 1

                # Write new separator if needed
                if c_ref_name != ref_name or c_read_name != read_name:
                    output.write ("#\t{}\t{}\n".format(c_read_name, c_ref_name))

                # Initialise a new kmer
                ref_name, ref_pos, ref_kmer, read_name = c_ref_name, c_ref_pos, c_ref_kmer, c_read_name
                if idx_agregate:
                    start_idx, end_idx = c_start_idx, c_end_idx
                if sample_agregate:
                    raw_list = c_raw_list

        # Last line exception
        res_list = [ref_pos, ref_kmer]
        if idx_agregate:
            res_list.extend ([start_idx, end_idx])
        if sample_agregate:
            l = np.array (raw_list, dtype=np.float32)
            res_list.extend ([round (np.mean (l), 3), round (np.std (l), 3)])
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

    # Flag last line
    output.write ("#\n")

    # Close files
    if input_fn:
        input.close()
    if output_fn:
        output.close()
