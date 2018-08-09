# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
import sys

# Third party imports
import numpy as np

# Local imports
from NanopolishComp.Helper_fun import stdout_print, stderr_print

#~~~~~~~~~~~~~~FUNCTIONS~~~~~~~~~~~~~~#
def Eventalign_collapse ():

    read_id_set = set()
    nkmers = 0

    try:
        # Check header fields
        line = sys.stdin.readline()
        ls = line.rstrip().split("\t")

        if "start_idx" in ls and "end_idx" in ls:
            idx_agregate = True
            stdout_print ("ref_name", "ref_pos", "ref_kmer", "read_name", "kmer_mean", "kmer_std", "start_idx", "end_idx", sep="\t")
        else:
            idx_agregate = False
            stdout_print ("ref_name", "ref_pos", "ref_kmer", "read_name", "kmer_mean", "kmer_std", sep="\t")

        #First line exception
        ls = sys.stdin.readline().rstrip().split("\t")
        n_events = 1
        ref_name, ref_pos, ref_kmer, read_name = ls[0], ls[1], ls[2], ls[3]
        mean_list = [np.float32(ls[6])]
        std_list = [np.float32(ls[7])]
        len_list = [np.float32(ls[8])]
        if idx_agregate:
            start_idx, end_idx = ls[13], ls[14]

        for line in sys.stdin:

            # Extract important fields from the file
            ls = line.rstrip().split("\t")
            c_ref_name, c_ref_pos, c_ref_kmer, c_read_name = ls[0], ls[1], ls[2], ls[3]
            c_event_mean, c_event_std, c_event_length  = np.float32(ls[6]), np.float32(ls[7]), np.float32(ls[8])
            if idx_agregate:
                c_start_idx, c_end_idx = ls[13], ls[14]

            # Update start is still same poisition
            if c_ref_name == ref_name and c_ref_pos == ref_pos:
                n_events += 1
                mean_list.append (c_event_mean)
                std_list.append (c_event_std)
                len_list.append (c_event_length)
                if idx_agregate:
                    start_idx = c_start_idx

            # Write new kmer
            else:
                mean = mean_list[0] if n_events == 1 else round (np.average (mean_list, weights=len_list), 2)
                std = std_list[0] if n_events == 1 else round (np.average (std_list, weights=len_list), 3)
                if idx_agregate:
                    stdout_print (ref_name, ref_pos, ref_kmer, read_name, mean, std, start_idx, end_idx, sep="\t")
                else:
                    stdout_print (ref_name, ref_pos, ref_kmer, read_name, mean, std, sep="\t")

                # Update Counters
                read_id_set.add (read_name)
                nkmers += 1

                # Initialise a new kmer
                n_events = 1
                ref_name, ref_pos, ref_kmer, read_name = c_ref_name, c_ref_pos, c_ref_kmer, c_read_name
                mean_list = [c_event_mean]
                std_list = [c_event_std]
                len_list = [c_event_length]
                if idx_agregate:
                    start_idx, end_idx = ls[13], ls[14]

        # Last line exception
        mean = mean_list[0] if n_events == 1 else round (np.average (mean_list, weights=len_list), 2)
        std = std_list[0] if n_events == 1 else round (np.average (std_list, weights=len_list), 3)
        if idx_agregate:
            stdout_print (ref_name, ref_pos, ref_kmer, read_name, mean, std, start_idx, end_idx, sep="\t")
        else:
            stdout_print (ref_name, ref_pos, ref_kmer, read_name, mean, std, sep="\t")

        # Update Counters
        read_id_set.add (read_name)
        nkmers += 1

    except (BrokenPipeError, IOError, KeyboardInterrupt):
        pass

    stderr_print ("[NanopolishComp summary] total reads: {:,}, total kmer:{:,}".format (len (read_id_set), nkmers))
    sys.stderr.close()
    sys.stdout.close()
