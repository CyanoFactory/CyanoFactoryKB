"""
Copyright (c) 2014 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""
import re
from Bio import SeqIO

import cyano.models as cmodels
from cyano.helpers import slugify
from bioparser import BioParser
from django.db.transaction import commit_on_success

class FastaFeature(BioParser):
    """Parses a FASTA file format used to indicate Chromosome features.

    FASTA is a simple text format, the header starts with ``>`` and all
    non-``>`` lines is sequence data.

    This class expects files with header format:

    ``>ref|IDENTIFIER|:SeqBegin-SeqEnd|DESCRIPTION``

    Example:

    ``>ref|SNP_ID 0|:1-8|``

    Only start and end are used.
    """
    def __init__(self, wid, user, reason, chromosome, feature_type):
        super(FastaFeature, self).__init__(wid, user, reason)

        if chromosome is None:
            raise ValueError("chromosome argument is mandatory")
        if feature_type is None:
            raise ValueError("feature_type argument is mandatory")

        self.chromosome = cmodels.Genome.objects.for_species(self.species).for_wid(self.try_slugify("Chromosome", chromosome))
        self.feature_type = self.try_slugify("feature_type", feature_type)
        self.data = []

    @staticmethod
    def _parse_header(header):
        header = header.split("|")
        if len(header) < 3:
            raise ValueError("Invalid header {}".format(header))

        pos = re.match(r":([0-9]+)-([0-9]+)", header[2])
        if not pos:
            raise ValueError("Invalid location (was {})".format(header[2]))

        if len(header) > 3:
            description = header[3]
        else:
            description = ""

        start, end = pos.groups()

        return header[1], start, end, description

    def parse(self, handle):
        self.report_progress(current=0, total=1, message="Parsing FASTA file")

        for record in SeqIO.parse(handle, "fasta"):
            wid, start, end, description = FastaFeature._parse_header(record.description)
            wid = slugify(wid)

            self.data.append({
                "wid": wid,
                "start": int(start),
                "end": int(end),
                "description": description
            })

    @commit_on_success
    def apply(self):
        current = 0
        total = len(self.data)

        self.report_progress(current=current, total=total, message="Importing FASTA feature")

        self.detail.save()

        self.species.save(self.detail)

        typ = cmodels.Type.objects.for_wid(self.feature_type, create=True)
        typ.save(self.detail)
        typ.species.add(self.species)

        for i, item in enumerate(self.data):
            wid = item["wid"]
            start = item["start"]
            end = item["end"]

            direction = "f" if start < end else "r"
            length = abs(end - start)
            coordinate = start if direction == "f" else end

            cf = cmodels.ChromosomeFeature.objects.for_species(self.species).for_wid(wid, create=True)
            cf.save(self.detail)
            cf.type.add(typ)
            cf.species.add(self.species)

            cmodels.FeaturePosition.objects.get_or_create(chromosome_feature=cf,
                                                          chromosome=self.chromosome,
                                                          coordinate=coordinate,
                                                          length=length,
                                                          direction=direction)

            self.report_progress(current=i+1, total=total, message="Importing FASTA feature: {}".format(wid))
