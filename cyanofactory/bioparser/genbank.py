"""
Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

from Bio import SeqIO

import cyano.models as cmodels
from cyano.helpers import slugify
from .bioparser import BioParser
from django.db.transaction import atomic


class Genbank(BioParser):    
    def __init__(self, wid, user, reason, chromosome, name):
        super(Genbank, self).__init__(wid, user, reason)
        
        if not chromosome:
            raise ValueError("chromosome argument is mandatory")
        if not name:
            raise ValueError("name argument is mandatory")
        
        self.chromosome = self.try_slugify("Chromosome", chromosome)
        self.name = name

    def parse(self, handle):
        if hasattr(self, "notify_progress"):
            self.notify_progress(current = 0, total = 1, message = "Parsing Genbank file...")
        
        f = SeqIO.parse(handle, "genbank")
        for record in f:
            if not record.annotations:
                raise ValueError("File lacks annotation block")

            self.annotation = record.annotations

            self.species.comments = ""

            self.is_chromosome = True

            if record.description is not None:
                if "plasmid" in record.description.lower():
                    self.is_chromosome = False
                self.species.comments = record.description
            
            if "comment" in self.annotation:
                self.species.comments += "\n" + self.annotation["comment"]

            self.species.genetic_code = ''

            self.record = record

            if len(record.features) > 0:
                if record.features[0].type == "source":
                    if "plasmid" in record.features[0].qualifiers:
                        self.is_chromosome = False
            
            # Genbank files only have one record
            break

    @atomic
    def apply(self):
        self.detail.save()
        
        self.species.save(self.detail)

        obj = cmodels.Chromosome if self.is_chromosome else cmodels.Plasmid
        chromosome = obj.objects.for_species(self.species).for_wid(self.chromosome, create=True)
        chromosome.name = self.name
        chromosome.sequence = str(self.record.seq) # Cast needed, otherwise revision-compare fails!
        chromosome.length = len(self.record.seq)
        chromosome.species = self.species
        chromosome.save(self.detail)

        if self.record.dbxrefs:
            for xref in self.record.dbxrefs:
                # BioPython doesnt always properly split the db xrefs
                xref = xref.split(" ")
                for x in xref:
                    if ":" in x:
                        source, xid = x.split(":")
                        x = cmodels.CrossReference.objects.get_or_create_with_revision(self.detail, source=source, xid=xid)
                        chromosome.cross_references.add(x)
    
        if "references" in self.annotation:
            for ref in self.annotation["references"]:
                # calculate the wid
                if ref.pubmed_id:
                    wid = "PUB_" + ref.pubmed_id
                    name = "Pubmed #" + ref.pubmed_id
                elif ref.medline_id:
                    wid = "MED_" + ref.medline_id
                    name = "Pubmed #" + ref.medline_id
                else:
                    publication = cmodels.PublicationReference.objects.filter(
                        authors__exact = ref.authors, title__exact = ref.title, publication__exact = ref.journal)
                    next_id = 0
                    if publication.exists():
                        wid = publication[0].wid
                        name = publication[0].name
                    else:
                        refs = cmodels.PublicationReference.objects.filter(wid__startswith = "REF_")
                        if refs.exists():
                            last = refs.reverse()[0]
                            next_id = int(last.wid[4:], 10) + 1
                            
                            wid = "REF_" + "%04d" % (next_id)
                            name = "Reference #%04d" % (next_id)
                        else:
                            wid = "REF_0001"
                            name = "Reference #0001"
            
                pubref = cmodels.PublicationReference.objects.for_wid(slugify(wid), create=True)
                pubref.name = name
                pubref.authors = ref.authors
                pubref.title = ref.title
                pubref.publication = ref.journal
                pubref.species = self.species
                pubref.save(self.detail)

                if ref.pubmed_id:
                    xref = cmodels.CrossReference.objects.get_or_create_with_revision(self.detail, source="PUBMED", xid=ref.pubmed_id)
                    pubref.cross_references.add(xref)

                if ref.medline_id:
                    xref = cmodels.CrossReference.objects.get_or_create_with_revision(self.detail, source="MEDLINE", xid=ref.medline_id)
                    pubref.cross_references.add(xref)

                chromosome.publication_references.add(pubref)
    
        if "gi" in self.annotation:
            xref = cmodels.CrossReference.objects.get_or_create_with_revision(self.detail, xid=self.annotation["gi"], source="GI")
            chromosome.cross_references.add(xref)
        
        features = self.record.features

        if len(features) > 0:
            if features[0].type == "source":
                if "db_xref" in features[0].qualifiers:
                    for xref in features[0].qualifiers["db_xref"]:
                        if ":" in xref:
                            source, xid = xref.split(":")
                            xref = cmodels.CrossReference.objects.get_or_create_with_revision(self.detail, source=source, xid=xid)
                            chromosome.cross_references.add(xref)
        
        gene_features = filter(lambda x : x.type == "gene", features)
        cds_features = filter(lambda x : x.type in ["CDS", "ncRNA", "rRNA", "tmRNA", "tRNA"], features)
        
        gene_map = {}
        for g in gene_features:
            if not "locus_tag" in g.qualifiers:
                self.stderr.write("WARN: " + str(g) + " without locus")
                continue
            loci = g.qualifiers["locus_tag"][0]
            if loci in gene_map:
                raise ValueError("locus_tag " + loci + " appeared twice")
            gene_map[loci] = g
        
        cds_map = {}
        for c in cds_features:
            if not "locus_tag" in c.qualifiers:
                self.stderr.write("WARN: " + str(c) + " without locus")
                continue
            loci = c.qualifiers["locus_tag"][0]
            if loci in cds_map:
                raise ValueError("locus_tag " + loci + " appeared twice")
            if loci in gene_map:
                cds_map[loci] = c
        
        sorted_cds_values = sorted(cds_map.values(), key = lambda x: x.qualifiers["locus_tag"])
        for i, v in enumerate(sorted_cds_values):

            qualifiers = v.qualifiers
            
            if not self.species.genetic_code:
                if "transl_table" in qualifiers:
                    self.species.genetic_code = qualifiers["transl_table"][0]
                    self.species.save(self.detail)
            
            g = cmodels.Gene.objects.for_species(self.species).for_wid(slugify(qualifiers["locus_tag"][0]), create = True)
            
            if hasattr(self, "notify_progress"):
                outstr = "Importing Gene %s (%d/%d)" % (g.wid, i + 1, len(cds_map.values()))
                self.notify_progress(current = i + 1, total = len(cds_map.values()), message = outstr)
    
            g.chromosome = chromosome
    
            if "gene" in qualifiers:
                g.name = qualifiers["gene"][0]
                g.symbol = qualifiers["gene"][0]

            g.direction = 'f' if v.location.strand == 1 else 'r'

            # __len__ because len() fails for numbers < 0
            # Joins output the wrong length
            if v.location.__len__() < 0:
                g.length = v.location.__len__() + len(self.record.seq)
            else:
                g.length = len(v.location)

            g.coordinate = v.location.start + 1 if 'f' else v.location.start

            if "note" in qualifiers:
                g.comments = "\n".join(qualifiers["note"])

            g.species = self.species
            g.save(self.detail)

            if "db_xref" in qualifiers:
                for xref in qualifiers["db_xref"]:
                    if ":" in xref:
                        source, xid = xref.split(":")
                        xref = cmodels.CrossReference.objects.get_or_create_with_revision(self.detail, source=source, xid=xid)
                        g.cross_references.add(xref)

            if "EC_number" in qualifiers:
                for ec in qualifiers["EC_number"]:
                    xref = cmodels.CrossReference.objects.get_or_create_with_revision(self.detail, source="EC", xid=ec)
                    g.cross_references.add(xref)

            if "gene_synonym" in qualifiers:
                for synonym in qualifiers["gene_synonym"]:
                    # Inconsistency: Multiple synonyms appear in one entry,
                    # why don't they split them like for all other items?
                    for syn in synonym.split(";"):
                        obj = cmodels.Synonym.objects.get_or_create_with_revision(self.detail, name=syn.strip())
                        g.synonyms.add(obj)

            if "protein_id" in qualifiers:
                protxref = qualifiers["protein_id"][0]
                wid = slugify(g.wid + "_Monomer")

                protein = cmodels.ProteinMonomer.objects.for_species(self.species).for_wid(wid, create=True)

                if "product" in qualifiers:
                    protein.name = qualifiers["product"][0]

                xref = cmodels.CrossReference.objects.get_or_create_with_revision(self.detail, source="RefSeq", xid=protxref)

                protein.gene = g
                protein.species = self.species
                protein.save(self.detail)


                protein.cross_references.add(xref)


            if v.type == "CDS":
                v.type = "mRNA"
            
            t = cmodels.Type.objects.for_wid(wid=slugify(v.type), create=True)
            t.name = v.type

            t.species = self.species
            t.save(self.detail)

            g.type.add(t)
            
        if hasattr(self, "notify_progress"):
            outstr = "Assigning KEGG pathways"
            self.notify_progress(current=len(cds_map.values()), total=len(cds_map.values()), message=outstr)
        
        cmodels.Pathway.add_kegg_pathway(self.species, self.detail)

