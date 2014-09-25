"""
Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""
import cyano.models as cmodels
from cyano.helpers import slugify
from bioparser import BioParser
from django.db.transaction import atomic
from django.core.exceptions import ObjectDoesNotExist

import re
import xml.etree.ElementTree as ET
from collections import OrderedDict


class InterProScan(BioParser):
    def __init__(self, wid, user, reason):
        super(InterProScan, self).__init__(wid, user, reason)

        self.protein_monomers_cf = OrderedDict()
        self.xrefs = OrderedDict()
        self.feature_positions = OrderedDict()
        self.types = OrderedDict()
        self.type_cache = {}
        self.total = 1

    def parse(self, handle):
        if hasattr(self, "notify_progress"):
            self.notify_progress(current=0, total=1, message="Parsing InterProScan file...")
        
        xml = handle.read()
        # Remove xmlns namespace, makes working with ElementTree more complicated
        xml = re.sub(' xmlns="[^"]+"', '', xml, count=1)

        root = ET.fromstring(xml)

        all_proteins = root.findall("protein")
        all_proteins_len = len(all_proteins)

        for i, protein in enumerate(all_proteins):
            xref = protein.find("xref")
            wid = xref.get("id").split("|", 2)[0]

            self.total = all_proteins_len

            if hasattr(self, "notify_progress"):
                self.notify_progress(current=i+1, total=self.total, message="Parsing features of {} ({}/{})".
                                                                                        format(wid, i+1, self.total))

            try:
                protein_item = cmodels.ProteinMonomer.objects.for_species(self.species).for_wid(wid)
                self.protein_monomers_cf[protein_item] = []
            except ObjectDoesNotExist:
                # ToDo: Error reporting?
                continue

            for matches in protein.findall("matches"):
                for match in matches:
                    signature = match.find("signature")
                    wid = slugify(signature.get("ac"))
                    cf = cmodels.ChromosomeFeature(wid=wid)
                    cf.name = signature.get("name") or ""
                    cf.comments = signature.get("desc") or ""
                    self.xrefs[wid] = []
                    self.protein_monomers_cf[protein_item].append(cf)
                    self.types[wid] = match.tag.title()

                    for entry in signature.findall("entry"):
                        for xref in entry:
                            self.xrefs[wid].append([xref.get("id"), xref.get("db")])

                    self.feature_positions[wid] = []
                    locations = match.find("locations")
                    for location in locations:
                        start = int(location.get("start"))
                        end = int(location.get("end"))
                        direction = "f"
                        if start > end:
                            start, end = end, start
                            direction = "r"
                        length = end - start
                        self.feature_positions[wid].append({"chromosome": protein_item.gene.chromosome_id,
                                                            "coordinate": start + protein_item.gene.coordinate,
                                                            "length": length,
                                                            "direction": direction})

    @atomic
    def apply(self):
        self.detail.save()

        match_type = cmodels.Type.objects.for_wid("Match-Type", create=True)
        match_type.save(self.detail)
        match_type.species.add(self.species)

        for i, values in enumerate(self.protein_monomers_cf.items()):
            protein_monomer, chromosome_feature = values
            if hasattr(self, "notify_progress"):
                self.notify_progress(current=i+1, total=self.total, message="Importing features of {} ({}/{})".
                                                                        format(protein_monomer.wid, i+1, self.total))

            for cf in chromosome_feature:
                real_cf = cmodels.ChromosomeFeature.objects.for_species(self.species).for_wid(cf.wid, create=True)
                real_cf.name = cf.name
                real_cf.comments = cf.comments
                real_cf.save(self.detail)

                real_cf.species.add(self.species)

                if cf.wid in self.xrefs:
                    for item in self.xrefs[cf.wid]:
                        xid, source = item
                        x = cmodels.CrossReference.objects.get_or_create_with_revision(self.detail, xid=xid, source=source)
                        real_cf.cross_references.add(x)

                if cf.wid in self.feature_positions:
                    for fp in self.feature_positions[cf.wid]:
                        cmodels.FeaturePosition.objects.get_or_create_with_revision(self.detail,
                                                                                    chromosome_feature=real_cf,
                                                                                    chromosome_id=fp["chromosome"],
                                                                                    coordinate=fp["coordinate"],
                                                                                    length=fp["length"],
                                                                                    direction=fp["direction"])

                if cf.wid in self.types:
                    typ = self.types[cf.wid]
                    if typ in self.type_cache:
                        typ_obj = self.type_cache[typ]
                    else:
                        typ_obj = cmodels.Type.objects.for_wid(typ, create=True)
                        typ_obj.parent = match_type
                        typ_obj.save(self.detail)
                        self.type_cache[typ] = typ_obj
                    typ_obj.species.add(self.species)

                    real_cf.type.add(typ_obj)
