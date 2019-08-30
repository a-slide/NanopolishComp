# -*- coding: utf-8 -*-

# Script inspired from https://github.com/jts/nanopolish/blob/master/scripts/calculate_methylation_frequency.py

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Standard library imports
from collections import *
import csv
import datetime
from random import randint

# Third party imports
from tqdm import tqdm
import numpy as np

# Local imports
from NanopolishComp.common import *
from NanopolishComp import __version__ as package_version
from NanopolishComp import __name__ as package_name

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~MAIN CLASS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

class Freq_meth_calculate():

    def __init__ (self,
        input_fn:"str",
        fasta_index:"str"="",
        output_bed_fn:"str"="",
        output_tsv_fn:"str"="",
        min_depth:"int"=10,
        sample_id:"str"="",
        verbose:"bool"=False,
        quiet:"bool"=False):
        """
        Calculate methylation frequency at genomic CpG sites from the output of nanopolish call-methylation
        * input_fn
            Path to a nanopolish call_methylation tsv output file
        * fasta_index
            fasta index file obtained with samtools faidx. Required for coordinate sorting
        * output_bed_fn
            Path to write a summary result file in BED format
        * output_tsv_fn
            Path to write an more extensive result report in TSV format
        * min_depth
            Minimal number of reads covering a site to be reported
        * sample_id
            Sample ID to be used for the bed track header
        * verbose
            Increase verbosity
        * quiet
            Reduce verbosity
        """

        # Save init options in dict for later
        kwargs = locals()

        # Define overall verbose level
        self.log = get_logger(name="Freq_meth_calculate", verbose=verbose, quiet=quiet)

        # Print option summary log
        self.log.debug ("## Options summary ##")
        self.log.debug ("\tpackage_name: {}".format(package_name))
        self.log.debug ("\tpackage_version: {}".format(package_version))
        self.log.debug ("\ttimestamp: {}".format(str(datetime.datetime.now())))
        self.log.debug (dict_to_str(kwargs, nsep=1, exclude_list=["self"]))

        # Verify parameters validity
        self.log.warning ("## Checking arguments ##")

        # Try to read input file if not a stream
        self.log.debug("\tTesting input file readability")
        if input_fn != 0 and not file_readable (input_fn):
            raise IOError ("Cannot read input file")

        # Verify that at least one output file is given:
        self.log.debug("\tCheck output file")
        if not output_bed_fn and not output_tsv_fn:
            raise NanopolishCompError("At least one output file should be given")
        if output_bed_fn:
            if os.path.dirname(output_bed_fn):
                mkdir (os.path.dirname(output_bed_fn), exist_ok=True)
            self.log.debug("\t\tOutput results in bed format")
        if output_tsv_fn:
            if os.path.dirname(output_tsv_fn):
                mkdir (os.path.dirname(output_tsv_fn), exist_ok=True)
            self.log.debug("\t\tOutput results in tsv format")

        # Create self variables
        self.counter = Counter()
        self.input_fn = input_fn
        self.fasta_index = fasta_index
        self.output_bed_fn = output_bed_fn
        self.output_tsv_fn = output_tsv_fn
        self.min_depth = min_depth
        self.sample_id = sample_id

        self.log.warning ("## Parsing methylation_calls file ##")
        self._parse_methylation_calls ()
        self.log.info ("## Results summary ##")
        self.log.info (dict_to_str(self.counter, nsep=1))

    def _parse_methylation_calls(self):
        """"""
        # Init SGC class with fasta_index
        if self.fasta_index:
            SGC.set_chrom_list(self.fasta_index)

        # Create collection to store results
        site_dict = defaultdict(list)

        try:
            input_fp = open (self.input_fn, "r")
            self.log.debug ("\tWrite output file header")
            if self.output_bed_fn:
                output_bed_fp = open (self.output_bed_fn, "w")
                output_bed_fp.write(Site.BED_header(self.sample_id)+"\n")
            if self.output_tsv_fn:
                output_tsv_fp = open (self.output_tsv_fn, "w")
                output_tsv_fp.write(Site.TSV_header()+"\n")

            self.log.info ("\tStarting to parse file Nanopolish methylation call file")
            header_line = input_fp.readline()
            byte_offset = len(header_line)
            lp = LineParser(header_line, sep="\t", cast_numeric_field=True)

            for line in tqdm(input_fp, desc="\t", unit=" lines", disable=self.log.level>=30):
                self.counter["Total read lines"]+=1
                byte_len = len(line)
                l = lp(line)

                if not l:
                    # Failsafe if line is malformed
                    self.counter["Invalid read line"]+=1
                else:
                    # Store byte offset corresponding to appropriate line
                    self.counter["Valid read lines"]+=1
                    site_dict[SGC(l.chromosome, l.start, l.end, l.strand)].append(byte_offset)
                    byte_offset += byte_len

            self.log.info ("\tFiltering out low coverage sites")
            filtered_site_dict = defaultdict(list)
            for k, offset_list in site_dict.items():
                self.counter["Total sites"]+=1

                # If low coverage unset list to release memory
                if len(offset_list) < self.min_depth:
                    self.counter["Low coverage sites"]+=1
                else:
                    self.counter["Valid sites"]+=1
                    filtered_site_dict[k]=offset_list
            del site_dict

            if self.fasta_index:
                self.log.info ("\tSorting by coordinates")

                filtered_site_dict = OrderedDict(sorted(filtered_site_dict.items(), key=lambda t: t[0]))

            self.log.info ("\tProcessing valid sites found")
            for k, offset_list in tqdm(filtered_site_dict.items(), desc="\t", unit=" sites", disable=self.log.level>=30):
                # Get all read lines corresponding to current site
                ll = []
                for offset in offset_list:
                    input_fp.seek(offset, 0)
                    ll.append (lp(input_fp.readline()))

                # Parse list with helper class Site
                site = Site (ll)
                self.counter["Valid sites"]+=1
                if self.output_bed_fn:
                    output_bed_fp.write(site.to_bed()+"\n")
                if self.output_tsv_fn:
                    output_tsv_fp.write(site.to_tsv()+"\n")
        finally:
            input_fp.close()
            if self.output_bed_fn:
                output_bed_fp.close()
            if self.output_tsv_fn:
                output_tsv_fp.close()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~HELPER CLASS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

class Site():
    """Structure like class to store site information"""

    # Class id dispatcher
    ID = 0
    @classmethod
    def next_id (cls):
        id = cls.ID
        cls.ID +=1
        return id

    @classmethod
    def BED_header (cls, sample_id=""):
        return "track name=Methylation_{} itemRgb=On".format(sample_id)

    @classmethod
    def TSV_header (cls):
        return "\t".join(["chromosome","start","end","strand","site_id","methylated_reads","unmethylated_reads","ambiguous_reads","sequence","num_motifs","llr_list"])

    def __init__ (self, ll, min_llr=2):
        """"""
        self.total = len(ll)
        self.methylated = 0
        self.unmethylated = 0
        self.ambiguous = 0
        self.id = self.next_id()
        self.sequence = ll[0].sequence
        self.num_motifs = ll[0].num_motifs
        self.chromosome = ll[0].chromosome
        self.start = ll[0].start
        self.end = ll[0].end+1
        self.strand = ll[0].strand

        llr_list = []
        for l in ll:
            llr_list.append(float(l.log_lik_ratio))
            # Count read methylation call per site
            if l.log_lik_ratio >= min_llr:
                self.methylated+=1
            elif l.log_lik_ratio <= -min_llr:
                self.unmethylated+=1
            else:
                self.ambiguous+=1

        self.med_llr = np.mean(llr_list)
        self.llr_list = llr_list
        if self.med_llr <= -min_llr:
            self.color = '8,121,207'
        elif self.med_llr < min_llr:
            self.color = '100,100,100'
        else:
            self.color = '235,5,79'

    def __repr__(self):
        return "{}:{}-{}({}) / id:{} / reads:{}".format(
            self.chromosome,
            self.start,
            self.end,
            self.strand,
            self.id,
            self.total)

    def to_bed (self):
        """"""
        return "{}\t{}\t{}\t{}\t{:.6f}\t{}\t{}\t{}\t'{}''".format(
            self.chromosome,
            self.start,
            self.end,
            self.id,
            self.med_llr,
            self.strand,
            self.start,
            self.end,
            self.color)

    def to_tsv (self):
        """"""
        return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
            self.chromosome,
            self.start,
            self.end,
            self.strand,
            self.id,
            self.methylated,
            self.unmethylated,
            self.ambiguous,
            self.sequence,
            self.num_motifs,
            ",".join([str(i) for i in self.llr_list]))

class SGC():
    """Sortable genomic coordinate object"""

    chr_list = OrderedDict()
    @classmethod
    def set_chrom_list (cls, index):
        if isinstance(index, dict):
            cls.chr_list = index
        else:
            with open(index) as fp:
                for i, line in enumerate(fp):
                    chrom = line.split()[0]
                    cls.chr_list[chrom]=i

    def __init__ (self, chrom, start=0, end=0, strand=""):
        # Verify and store chromosome name
        self.chrom = chrom
        self.start = int(start)
        if not end or int(end) < self.start:
            self.end = self.start
        else:
            self.end = int(end)
        self.strand = strand

    def __repr__ (self):
        return "{}:{}-{} ({})".format(self.chrom, self.start, self.end, self.strand)

    def __lt__(self, other):
        if self.chrom == other.chrom:
            if self.start == other.start:
                if self.end == other.end:
                    if self.strand == other.strand:
                        return randint(0,1)
                    else:
                        return self.strand < other.strand
                else:
                    return self.end < other.end
            else:
                return self.start < other.start
        else:
            return self.chr_list.get(self.chrom, 0) < self.chr_list.get(other.chrom, 0)

    def __hash__(self):
        return hash((self.chrom, self.start, self.end, self.strand))

    def __eq__(self, other):
        return (self.chrom, self.start, self.end, self.strand) == (other.chrom, other.start, other.end, other.strand)

    def __ne__(self, other):
        return not(self == other)
