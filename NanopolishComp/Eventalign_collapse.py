# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#

# Disable multithreading for MKL and openBlas
import os
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["MKL_THREADING_LAYER"] = "sequential"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ['OPENBLAS_NUM_THREADS'] = '1'

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
import multiprocessing as mp
from time import time
from collections import OrderedDict
import traceback

# Third party imports
import numpy as np
from tqdm import tqdm

# Local imports
from NanopolishComp.common import stderr_print, access_file, NanopolishCompError

#~~~~~~~~~~~~~~CLASS~~~~~~~~~~~~~~#
class Eventalign_collapse ():
    """
    Collapse the nanopolish eventalign output by kmers rather that by events.
    kmer level statistics (mean, median, std, mad) are only computed if nanopolish is run with --samples option
    """

    def __init__ (self,
        output_fn,
        input_fn=0,
        threads=4,
        max_reads=None,
        write_samples=False,
        stat_fields=["mean", "median", "num_signals"],
        verbose=False,):
        """
        * output_fn
            Path the output eventalign collapsed tsv file
        * input_fn
            Path to a nanopolish eventalign tsv output file. If '0' read from std input (default = 0)
        * threads
            Total number of threads. 1 thread is used for the reader and 1 for the writer (default = 4)
        * max_reads
            Maximum number of read to parse. 0 to deactivate (default = 0)
        * write_samples
            If given, will write the raw sample if eventalign is run with --samples option
        * stat_fields
            List  of statistical fields generated if nanopolish is ran with --sample option.
            Valid values = "mean", "std", "median", "mad", "num_signals"
        """

        if input_fn and not access_file (input_fn):
            raise IOError ("Cannot read input file")
        if threads < 3:
            raise ValueError ("At least 3 threads required")
        for field in stat_fields:
            if not field in ["mean", "std", "median", "mad", "num_signals"]:
                raise ValueError ("Invalid value in stat_field {}. Valid entries = mean, std, median, mad, num_signals".format(field))

        self.output_fn = output_fn
        self.input_fn = input_fn
        self.threads = threads-2 # Remove 2 threads for read and write
        self.max_reads = max_reads
        self.write_samples = write_samples
        self.stat_fields = stat_fields
        self.verbose = verbose

        if self.verbose: stderr_print ("Collapse file by read_id/ref_id\n")
        # Init Multiprocessing variables
        in_q = mp.Queue (maxsize = 1000)
        out_q = mp.Queue (maxsize = 1000)
        error_q = mp.Queue ()

        # Define processes
        ps_list = []
        ps_list.append (mp.Process (target=self._split_reads, args=(in_q, error_q)))
        for i in range (self.threads):
            ps_list.append (mp.Process (target=self._process_read, args=(in_q, out_q, error_q)))
        ps_list.append (mp.Process (target=self._write_output, args=(out_q, error_q)))

        try:
            # Start all processes
            for ps in ps_list:
                ps.start ()
            # Monitor error queue
            for E in iter (error_q.get, None):
                raise E
            # Join processes
            for ps in ps_list:
                ps.join ()

        # Kill processes if any error
        except (BrokenPipeError, KeyboardInterrupt, NanopolishCompError) as E:
            for ps in ps_list:
                ps.terminate ()
            stderr_print ("\nAn error occured. All processes were killed\n")
            raise E

    #~~~~~~~~~~~~~~PRIVATE METHODS~~~~~~~~~~~~~~#
    def _split_reads (self, in_q, error_q):
        """
        Mono-threaded reader
        """
        try:
            # Open input file or stdin if 0
            with open (self.input_fn) as fp:

                # Get header line and extract corresponding index
                input_header = fp.readline().rstrip().split("\t")
                idx = self._get_field_idx (input_header)
                n_reads = 0

                # First data line exception
                read_l = []
                event_l = fp.readline().rstrip().split("\t")
                cur_read_id = event_l[idx["read_id"]]
                cur_ref_id = event_l[idx["ref_id"]]
                event_d = self._event_list_to_dict (event_l, idx)
                read_l.append (event_d)

                for line in fp:
                    # Early ending if required
                    if self.max_reads and n_reads == self.max_reads:
                        break
                    # Get event line
                    event_l = line.rstrip().split("\t")
                    read_id = event_l[idx["read_id"]]
                    ref_id = event_l[idx["ref_id"]]
                    event_d = self._event_list_to_dict (event_l, idx)

                    # Line correspond to the same ids
                    if read_id != cur_read_id or ref_id != cur_ref_id:
                        in_q.put ((cur_read_id, cur_ref_id, read_l))
                        n_reads+=1
                        read_l = []
                        cur_read_id = read_id
                        cur_ref_id = ref_id

                    # In any case extend list corresponding to current read_id/ref_id
                    read_l.append(event_d)

                # Last data line exception
                in_q.put ((cur_read_id, cur_ref_id, read_l))
                n_reads+=1

        # Manage exceptions and deal poison pills
        except Exception:
            error_q.put (NanopolishCompError(traceback.format_exc()))
        finally:
            for i in range (self.threads):
                in_q.put(None)

    def _process_read (self, in_q, out_q, error_q):
        """
        Multi-threaded workers
        """
        try:
            # Collapse event at kmer level
            for read_id, ref_id, read_l in iter(in_q.get, None):

                # Write read header to str
                read_str = "#{}\t{}\n".format(read_id, ref_id)
                read_str+= "{}\n".format (self._make_ouput_header(event_d=read_l[0]))

                # Init values for first kmer
                kmer_d = self._init_kmer_dict(event_d=read_l[0])

                # Init read dictionary
                read_d = OrderedDict ()
                read_d["read_id"] = read_id
                read_d["ref_id"] = ref_id
                read_d["dwell_time"] = 0.0
                read_d["kmers"] = 0
                read_d["NNNNN_kmers"] = 0
                read_d["mismatch_kmers"] = 0
                read_d["missing_kmers"] = 0
                read_d["ref_start"] = kmer_d["ref_pos"]

                # Iterate over the rest of the lines
                for event_d in read_l [1:]:
                    pos_offset = event_d["ref_pos"]-kmer_d["ref_pos"]

                    # Same position = update current kmer
                    if pos_offset == 0:
                        kmer_d = self._update_kmer_dict(kmer_d=kmer_d, event_d=event_d)

                    # New position = write previous kmer and start new one
                    else:
                        # Update read counter
                        read_d["dwell_time"] += kmer_d["dwell_time"]
                        if kmer_d["NNNNN_dwell_time"]:
                            read_d["NNNNN_kmers"] += 1
                        if kmer_d ["mismatch_dwell_time"]:
                            read_d["mismatch_kmers"] += 1
                        if pos_offset >=2:
                            read_d["missing_kmers"] += (pos_offset-1)
                        read_d["kmers"] += 1
                        # Converts previous kmer to str and init new kmer
                        read_str += "{}\n".format(self._kmer_dict_to_str(kmer_d=kmer_d))
                        kmer_d = self._init_kmer_dict(event_d=event_d)

                # Last read_d update
                read_d["dwell_time"] += kmer_d["dwell_time"]
                if kmer_d ["NNNNN_dwell_time"]:
                    read_d["NNNNN_kmers"] += 1
                if kmer_d ["mismatch_dwell_time"]:
                    read_d["mismatch_kmers"] += 1
                if pos_offset >=2:
                    read_d["missing_kmers"] += (pos_offset-1)
                read_d["ref_end"] = kmer_d["ref_pos"]+1
                read_d["kmers"] += 1

                # Last kmer
                read_str += "{}\n".format(self._kmer_dict_to_str(kmer_d=kmer_d))

                # Add the current read details to queue
                out_q.put((read_d, read_str))

        # Manage exceptions and deal poison pills
        except Exception:
            error_q.put (NanopolishCompError(traceback.format_exc()))
        finally:
            out_q.put(None)

    def _write_output (self, out_q, error_q):
        """
        Mono-threaded Writer
        """
        byte_offset = n_reads = 0
        t = time()

        try:

            # Open output files
            with open (self.output_fn, "w") as output_fp,\
                 open (self.output_fn+".idx", "w") as idx_fp,\
                 tqdm (unit=" reads", mininterval=0.1, smoothing=0.1, disable= not self.verbose) as pbar:

                idx_fp.write ("ref_id\tref_start\tref_end\tread_id\tdwell_time\tkmers\tNNNNN_kmers\tmismatch_kmers\tmissing_kmers\tbyte_offset\tbyte_len\n")

                n_reads = 0
                for _ in range (self.threads):
                    for (read_d, read_str) in iter (out_q.get, None):
                        byte_len = len(read_str)

                        output_fp.write (read_str)
                        idx_fp.write ("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format (
                            read_d["ref_id"],
                            read_d["ref_start"],
                            read_d["ref_end"],
                            read_d["read_id"],
                            read_d["dwell_time"],
                            read_d["kmers"],
                            read_d["NNNNN_kmers"],
                            read_d["mismatch_kmers"],
                            read_d["missing_kmers"],
                            byte_offset,
                            byte_len-1))

                        byte_offset += byte_len
                        n_reads += 1
                        if self.verbose: pbar.update(1)

                # Flag last line
                output_fp.write ("#\n")

        # Manage exceptions and deal poison pills
        except Exception:
            error_q.put (NanopolishCompError(traceback.format_exc()))
        finally:
            error_q.put(None)

        stderr_print ("[Eventalign_collapse] total reads: {} [{} reads/s]\n".format(n_reads, round (n_reads/(time()-t), 2)))

    #~~~~~~~~~~~~~~HELPER PRIVATE METHODS~~~~~~~~~~~~~~#

    def _get_field_idx (self, input_header):
        """"""
        # Get index of fields to fetch
        idx = OrderedDict()
        idx["ref_id"] = input_header.index ("contig")
        if "read_name" in input_header:
            idx["read_id"] = input_header.index ("read_name")
        elif "read_index" in input_header:
            idx["read_id"] = input_header.index ("read_index")
        idx["ref_pos"] = input_header.index ("position")
        idx["ref_kmer"] = input_header.index ("reference_kmer")
        idx["mod_kmer"] = input_header.index ("model_kmer")
        idx["event_len"] = input_header.index ("event_length")
        # Facultative field start and end index
        if "start_idx" in input_header and "end_idx" in input_header:
            idx["start_idx"] = input_header.index ("start_idx")
            idx["end_idx"] = input_header.index ("end_idx")
        # Facultative field samples
        if "samples" in input_header:
            idx["samples"] = input_header.index ("samples")
        return idx

    def _event_list_to_dict (self, event_l, idx):
        """Get interesting fields from event list and cast in appropriate type"""
        event_d = OrderedDict()
        event_d["ref_pos"] = int(event_l[idx["ref_pos"]])
        event_d["ref_kmer"] = event_l[idx["ref_kmer"]]
        event_d["mod_kmer"] = event_l[idx["mod_kmer"]]
        event_d["event_len"] = float(event_l[idx["event_len"]])
        if "start_idx" in idx:
            event_d["start_idx"] = int(event_l[idx["start_idx"]])
            event_d["end_idx"] = int(event_l[idx["end_idx"]])
        if "samples" in idx:
            event_d["sample_list"] = event_l[idx["samples"]].split(",")
        return event_d

    def _init_kmer_dict (self, event_d):
        """Start a new kmer dict from first event values"""
        kmer_d = OrderedDict ()
        kmer_d["ref_pos"] = event_d["ref_pos"]
        kmer_d["ref_kmer"] = event_d["ref_kmer"]
        kmer_d["num_events"] = 1
        kmer_d["dwell_time"] = event_d["event_len"]
        kmer_d["NNNNN_dwell_time"] = 0.0
        kmer_d["mismatch_dwell_time"] = 0.0
        if event_d["mod_kmer"] == "NNNNN":
            kmer_d["NNNNN_dwell_time"] += event_d["event_len"]
        elif event_d["mod_kmer"] != event_d["ref_kmer"]:
            kmer_d["mismatch_dwell_time"] += event_d["event_len"]
        if "start_idx" in event_d:
            kmer_d["start_idx"] = event_d["start_idx"]
            kmer_d["end_idx"] = event_d["end_idx"]
        if "sample_list" in event_d:
            kmer_d["sample_list"] = event_d["sample_list"]
        return kmer_d

    def _update_kmer_dict (self, kmer_d, event_d):
        """Update kmer dict from subsequent event values"""
        kmer_d["num_events"] += 1
        kmer_d["dwell_time"] += event_d["event_len"]
        if event_d["mod_kmer"] == "NNNNN":
            kmer_d["NNNNN_dwell_time"] += event_d["event_len"]
        elif event_d["mod_kmer"] != event_d["ref_kmer"]:
            kmer_d["mismatch_dwell_time"] += event_d["event_len"]
        if "start_idx" in event_d:
            kmer_d["start_idx"] = event_d["start_idx"]
        if "sample_list" in event_d:
            kmer_d["sample_list"].extend(event_d["sample_list"])
        return kmer_d

    def _kmer_dict_to_str (self, kmer_d):
        """"""
        # Write base fields
        s = "{}\t{}\t{}\t{}\t{}\t{}".format(
            kmer_d["ref_pos"],
            kmer_d["ref_kmer"],
            kmer_d["num_events"],
            kmer_d["dwell_time"],
            kmer_d["NNNNN_dwell_time"],
            kmer_d["mismatch_dwell_time"])
        # Facultative index fields
        if "start_idx" in kmer_d:
            s += "\t{}\t{}".format(
                kmer_d["start_idx"],
                kmer_d["end_idx"])
        # Facultative samples fields
        if "sample_list" in kmer_d:
            sample_array = np.array (kmer_d["sample_list"], dtype=np.float32)
            if "mean" in self.stat_fields:
                s += "\t{}".format(np.mean (sample_array))
            if "std" in self.stat_fields:
                s += "\t{}".format(np.std (sample_array))
            if "median" in self.stat_fields:
                s += "\t{}".format(np.median (sample_array))
            if "mad" in self.stat_fields:
                s += "\t{}".format(np.median(np.abs(sample_array-np.median(sample_array))))
            if "num_signals" in self.stat_fields:
                s += "\t{}".format(len(sample_array))
            if self.write_samples:
                s += "\t{}".format(",".join(kmer_d["sample_list"]))
        return s

    def _make_ouput_header (self, event_d):
        """"""
        # Write base fields
        s = "ref_pos\tref_kmer\tnum_events\tdwell_time\tNNNNN_dwell_time\tmismatch_dwell_time"
        # Write extra fields
        if "start_idx" in event_d:
            s += "\tstart_idx\tend_idx"
        if "sample_list" in event_d:
            if "mean" in self.stat_fields:
                s += "\tmean"
            if "std" in self.stat_fields:
                s += "\tstd"
            if "median" in self.stat_fields:
                s += "\tmedian"
            if "mad" in self.stat_fields:
                s += "\tmad"
            if "num_signals" in self.stat_fields:
                s += "\tnum_signals"
            if self.write_samples:
                s += "\tsamples"
        return s
