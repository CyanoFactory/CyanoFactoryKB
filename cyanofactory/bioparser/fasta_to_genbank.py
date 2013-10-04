"""
Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

import re
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
from Bio import Seq
from Bio.SeqRecord import SeqRecord

class FastaToGenbank:
	"""Parses a FASTA file and converts it to a Genbank file.
	
	FASTA is a simple text format, the header starts with ``>`` and all
	non-``>`` lines is sequence data.
	
	This class expects files with header format
	
	``>ref|IDENTIFIER|:SeqBegin-SeqEnd|DESCRIPTION| [gene=NAME] [locus_tag=TAG]``
	
	Example:
	
	``>ref|NC_000913|:223771-225312|16S ribosomal RNA of rrnH operon| [gene=rrsH] [locus_tag=b0201]``
	
	The ``IDENTIFIER`` of the first FASTA item is used as ``id`` and ``name`` of the genbank file.
	
	For every gene a feature is created, containing the location, the locus_tag and the
	``DESCRIPTION`` as a note element. The type is always ``gene``.
	"""
	def __init__(self, filename):
		"""Parses the specified FASTA file"""
		self._file = filename
		self.errors = []

	def parse(self):
		with open(self._file) as handle:
			genbank = SeqRecord(Seq.UnknownSeq(0))
			header_pattern = re.compile(r"ref\|(?P<id>.*?)\|:(?P<start>[0-9]+)-(?P<end>[0-9]+)\|(?P<description>.*?)\|\s*\[gene=(?P<gene>\S+)\]\s*\[locus_tag=(?P<locus_tag>\S+)\]\s*")
			first = True			
			for record in SeqIO.parse(handle, "fasta"):
				header = record.description
				match = header_pattern.match(header)
				if not match:
					self.errors.append("Invalid header: >" + header)
					continue
				
				if first:
					first = False
					genbank.id = match.group("id")
					genbank.name = match.group("id")
				
				feature = SeqFeature(FeatureLocation(int(match.group("start")), int(match.group("end"))), type = "gene")
				feature.qualifiers = {"locus_tag": match.group("locus_tag"),
							"gene": match.group("gene"),
							"note": match.group("description"),
							"sequence": record.seq}
				genbank.features.append(feature)
			
			return genbank
		return None