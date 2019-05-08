# -*- coding: utf-8 -*-

# Script inspired from https://github.com/jts/nanopolish/blob/master/scripts/calculate_methylation_frequency.py

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Standard library imports
from collections import *
import csv
import datetime

# Third party imports
from tqdm import tqdm

# Local imports
from NanopolishComp.common import *
from NanopolishComp import __version__ as package_version
from NanopolishComp import __name__ as package_name

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~LOGGING INFO~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
import logging
logging.basicConfig(level=logging.INFO, format="%(message)s")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~MAIN CLASS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

class Freq_meth_calculate():

    def __init__ (self,
        input_fn:"str",
        outdir:"str"="./",
        outprefix:"str"="out",
        min_llr:"float"=2.5,
        min_depth:"int"=10,
        min_meth_freq:"float"=0.05,
        split_group:"bool"=False,
        motif:"{cpg,gpc,dam,dcm}"="cpg",
        verbose:"bool"=False,
        quiet:"bool"=False):
        """
        Calculate methylation frequency at genomic CpG sites from the output of nanopolish call-methylation
        * input_fn
            Path to a nanopolish call_methylation tsv output file
        * outdir
            Path to the output folder (will be created if it does exist yet)
        * outprefix
            text outprefix for all the files generated
        * min_llr
            Log likelihood ratio threshold
        * min_depth
            Minimal number of reads covering a site to be reported
        * min_meth_freq
            Minimal methylation frequency of a site to be reported
        * split_group
            If True, multi motif groups (sequence with close motifs) are split in individual site
        * motif
            Methylation motif type
        * verbose
            Increase verbosity
        * quiet
            Reduce verbosity
        """

        # Save init options in dict for later
        kwargs = locals()

        # Define overall verbose level
        self.log = logging.getLogger()
        if verbose:
            self.log.setLevel (logging.DEBUG)
        elif quiet:
            self.log.setLevel (logging.WARNING)
        else:
            self.log.setLevel (logging.INFO)

        # Collect args in dict for log report
        self.option_d = OrderedDict()
        self.option_d["package_name"] = package_name
        self.option_d["package_version"] = package_version
        self.option_d["timestamp"] = str(datetime.datetime.now())
        for i, j in kwargs.items():
            if i != "self":
                self.option_d[i]=j
        self.log.debug ("Options summary")
        self.log.debug (dict_to_str(self.option_d))

        # Verify parameters validity
        self.log.warning ("## Checking arguments ##")
        # Try to read input file if not a stream
        self.log.debug("\tTesting input file readability")
        if input_fn != 0 and not file_readable (input_fn):
            raise IOError ("Cannot read input file")
        # Try to create output folder
        self.log.debug("\tCreating output folder")
        mkdir(outdir, exist_ok=True)

        if motif == "cpg":
            motif_seq = "CG"
        elif motif == "gpc":
            motif_seq = "GC"
        elif motif == "dam":
            motif_seq = "GATC"
        elif motif == "dcm":
            motif_seq = "CCAGG"
        else:
            raise NanopolishCompError("Not a valid motif: {}".format(motif))

        # init object counters
        self.site_c = Counter()
        self.pos_c = Counter()

        self.log.warning ("## Parsing methylation_calls file ##")
        sites_d = self._parse_methylation_calls (
            fn = input_fn,
            min_llr = min_llr,
            min_depth = min_depth,
            min_meth_freq = min_meth_freq,
            split_group = split_group,
            motif_seq = motif_seq)

        if sites_d:
            self.log.warning ("## Write output files ##")
            self._write_output (sites_d=sites_d, outdir=outdir, outprefix=outprefix)
        else:
            self.log.warning("No valid site found in the imput file")

    def _parse_methylation_calls(self, fn, min_llr=2.5, min_depth=10, min_meth_freq=0.05, split_group=False, motif_seq="CG"):
        """"""
        # Create collection to store results
        self.log.info ("Starting to parse file methylation_calls file")
        sites_d = OrderedDict()
        with open (fn) as fp:
            for line_d in tqdm(csv.DictReader(fp, delimiter='\t'), unit=" lines", disable=self.log.level>=30):
                line_d = numeric_cast_dict (line_d)

                # Update counter if verbose only
                self.site_c["total"] +=1
                if line_d["log_lik_ratio"] >= min_llr:
                    self.site_c["methylated"] +=1
                elif line_d["log_lik_ratio"] <= min_llr:
                    self.site_c["unmethylated"] +=1
                else :
                    self.site_c["ambiguous"] +=1

                # In case multi motif group have to be split
                if line_d["num_motifs"] > 1 and split_group:

                    # Find all index corresponding to the motif
                    for pos in find_subseq_index (line_d["sequence"], motif_seq):
                        start = line_d["start"]+pos-5
                        seq = line_d["sequence"][pos-5:pos+6]

                        # Create entry in sites dict if it does not already exist
                        coord_key = "{}_{}_{}".format(line_d["chromosome"], start, line_d["strand"])
                        if not coord_key in sites_d:
                            sites_d[coord_key] = Site(seq, line_d["num_motifs"], line_d["chromosome"], start , start+1, line_d["strand"], min_llr=min_llr)
                        sites_d[coord_key].update(line_d["log_lik_ratio"])

                # Otherwise
                else:
                    # Create entry in sites dict if it does not already exist
                    coord_key = "{}_{}_{}".format(line_d["chromosome"], line_d["start"], line_d["strand"])
                    if not coord_key in sites_d:
                        sites_d[coord_key] = Site(line_d["sequence"], line_d["num_motifs"], line_d["chromosome"], line_d["start"] , line_d["end"]+1, line_d["strand"], min_llr=min_llr)
                    sites_d[coord_key].update(line_d["log_lik_ratio"])

        # Print read level counter summary
        self.log.debug ("Read sites summary")
        self.log.debug (dict_to_str(self.site_c))

        self.log.info ("Filtering out positions with low coverage or methylation frequency")
        filtered_sites_d = OrderedDict()
        for coord_key, site in tqdm(sites_d.items(), unit=" positions", disable=self.log.level>=30):

            # Update counter
            self.pos_c["total"] +=1
            if site.valid_call_reads < min_depth:
                self.pos_c["low_coverage"]+=1
            elif site.meth_freq < min_meth_freq:
                self.pos_c["low_meth_freq"]+=1
            else:
                self.pos_c["valid"] +=1

            if site.valid_call_reads >= min_depth and site.meth_freq >= min_meth_freq:
                filtered_sites_d[coord_key] = site

        # Print genomic positions level counter summary
        self.log.debug ("Genomic positions summary")
        self.log.debug (dict_to_str(self.pos_c))

        return filtered_sites_d

    def __repr__ (self):
        m = "General options:\n"
        m+=dict_to_str(self.option_d)
        m+="Read sites summary:\n"
        m+=dict_to_str(self.site_c)
        m+="Genomic positions summary:\n"
        m+=dict_to_str(self.pos_c)
        return m

    def _write_output (self, sites_d, outdir, outprefix):
        """"""
        # Bed file output
        self.log.info("Writing bed file")
        fn = os.path.join(outdir, outprefix+"_freq_meth_calculate.bed")
        with open (fn, "w") as fp:
            # Write header
            fp.write (Site.BED_header(outprefix)+"\n")
            # Sort by coordinates and write values
            for site in sorted(sites_d.values(), key=lambda s: (s.chrom, s.start, s.strand)):
                fp.write (site.to_bed()+"\n")

        # TSV file output
        self.log.info("Writing tsv file")
        fn = os.path.join(outdir, outprefix+"_freq_meth_calculate.tsv")
        with open (fn, "w") as fp:
            # Write header
            fp.write (Site.TSV_header()+"\n")
            # Sort by methylation frequency and write values
            for site in sorted(sites_d.values(), key=lambda s: s.meth_freq, reverse=True):
                fp.write (site.to_tsv()+"\n")

        # log file _write_output
        self.log.info("Writing log file")
        fn = os.path.join(outdir, outprefix+"_freq_meth_calculate.log")
        with open (fn, "w") as fp:
            fp.write (str(self))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~HELPER CLASS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

class Site ():
    """Structure like class to store site information"""

    # Class id dispatcher
    ID = 0
    @classmethod
    def next_id (cls):
        id = cls.ID
        cls.ID +=1
        return id

    @classmethod
    def BED_header (cls, track_name):
        return "track name='{}' description='Methylation frequency track generated with NanopolishComp' useScore=1".format(track_name)

    @classmethod
    def TSV_header (cls):
        return "\t".join(["chrom","start","end","strand","site_id","methylated_reads","unmethylated_reads","ambiguous_reads","sequence","num_motifs","meth_freq"])

    def __init__ (self, seq, num_motifs, chrom, start, end, strand, min_llr=2.5):
        """"""
        self.seq = seq
        self.num_motifs = num_motifs
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.methylated_reads = 0
        self.unmethylated_reads = 0
        self.ambiguous_reads = 0
        self.min_llr = min_llr
        self.id = self.next_id()

    def __repr__(self):
        return "{}:{}-{}({}) / id:{} / reads:{} / meth_freq:{:03}".format(
            self.chrom, self.start, self.end, self.strand, self.id, self.valid_call_reads, self.meth_freq)

    @property
    def valid_call_reads (self):
        return self.methylated_reads+self.unmethylated_reads

    @property
    def meth_freq (self):
        if self.valid_call_reads:
            return self.methylated_reads/self.valid_call_reads
        else:
            return 0

    def update (self, llr):
        """"""
        if llr >= self.min_llr:
            self.methylated_reads +=1
        elif llr <= -self.min_llr:
            self.unmethylated_reads +=1
        else:
            self.ambiguous_reads +=1

    def to_bed (self):
        """"""
        return "{}\t{}\t{}\t{}\t{:06}\t{}".format(
            self.chrom, self.start, self.end, self.id ,int(self.meth_freq*1000), self.strand)

    def to_tsv (self):
        """"""
        return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6f}".format(
            self.chrom, self.start, self.end, self.strand, self.id, self.methylated_reads, self.unmethylated_reads, self.ambiguous_reads, self.seq, self.num_motifs, self.meth_freq)
