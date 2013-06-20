from django.core.management.base import BaseCommand
from Bio import SeqIO
import cyano.models as cmodels
from django.core.exceptions import ObjectDoesNotExist
from argparse import ArgumentError
import sys

class Command(BaseCommand):
    def handle(self, *args, **options):
        for arg in args:
            print arg
            handle = open(arg, "r")
            f = SeqIO.parse(handle, "genbank")
            for record in f:
                if not record.annotations:
                    return # TODO: FAIL
                
                anno = record.annotations
                
                revdetail = cmodels.RevisionDetail()
                revdetail.user = cmodels.UserProfile.objects.get(user__username__exact = "gabriel")
                revdetail.reason = "Update"

                try:
                    species = cmodels.Species.objects.get(wid = record.name)
                except ObjectDoesNotExist:
                    species = cmodels.Species(wid = record.name)

                if "organism" in anno:
                    species.name = anno["organism"]
                
                species.comments = ""
                if record.description != None:
                    species.comments = record.description
                if "comment" in anno:
                    species.comments += "\n" + anno["comment"]

                species.genetic_code = '11'
                species.save(revision_detail = revdetail)
                
                if record.dbxrefs:
                    for xref in record.dbxrefs:
                        # BioPython doesnt always properly split the db xrefs
                        xref = xref.split(" ")
                        for x in xref:
                            if ":" in x:
                                source, xid = x.split(":")
                                wid = source + ":" + xid
                                try:
                                    x = cmodels.CrossReference.objects.get(wid = wid)
                                except ObjectDoesNotExist:
                                    x = cmodels.CrossReference(wid = wid)

                                x.name = wid
                                x.xid = xid
                                x.source = source
                                x.save(revision_detail = revdetail)
                                x.species.add(species)
                                x.save(revision_detail = revdetail)
                                species.cross_references.add(x)
                    
                    species.save(revision_detail = revdetail)

                if "references" in anno:
                    for ref in anno["references"]:
                        # calculate the wid
                        if ref.pubmed_id:
                            wid = "PUB_" + ref.pubmed_id
                            name = "Pubmed #" + ref.pubmed_id
                        elif ref.medline_id:
                            wid = "MED_" + ref.medline_id
                            name = "Pubmed #" + ref.medline_id
                        else:
                            refs = cmodels.PublicationReference.objects.filter(wid__startswith = "REF_")
                            next_id = 0
                            if refs.exists():
                                last = refs.reverse()[0]
                                next_id = int(last.wid[4:], 10) + 1
                                wid = "REF_" + "%04d" % (next_id)
                                name = "Reference #%04d" % (next_id)
                            else:
                                wid = "REF_0001"
                                name = "Reference #0001"
                        
                        try:
                            pubref = cmodels.PublicationReference.objects.get(wid = wid)
                        except ObjectDoesNotExist:
                            pubref = cmodels.PublicationReference(wid = wid)
                        pubref.name = name
                        pubref.authors = ref.authors
                        pubref.title = ref.title
                        pubref.journal = ref.journal
                        pubref.save(revision_detail = revdetail)

                        if ref.pubmed_id:
                            wid = "PUBMED" + ":" + ref.pubmed_id
                            try:
                                xref = cmodels.CrossReference.objects.get(wid = wid)
                            except ObjectDoesNotExist:
                                xref = cmodels.CrossReference(wid = wid)
                            xref.name = wid
                            xref.xid = ref.pubmed_id
                            xref.source = "PUBMED"
                            xref.save(revision_detail = revdetail)
                            xref.species.add(species)
                            xref.save(revision_detail = revdetail)
                            pubref.cross_references.add(xref)

                        if ref.medline_id:
                            wid = "MEDLINE" + ":" + ref.medline_id
                            try:
                                xref = cmodels.CrossReference.objects.get(wid = wid)
                            except ObjectDoesNotExist:
                                xref = cmodels.CrossReference(wid = wid)
                            xref.name = wid
                            xref.xid = ref.medline_id
                            xref.source = "MEDLINE"
                            xref.save(revision_detail = revdetail)
                            xref.species.add(species)
                            xref.save(revision_detail = revdetail)
                            pubref.cross_references.add(xref)
                            
                        pubref.species.add(species)
                        pubref.save(revision_detail = revdetail)
                        species.publication_references.add(pubref)

                    species.save(revision_detail = revdetail)

                if "gi" in anno:
                    wid = "GI" + ":" + anno["gi"]
                    try:
                        xref = cmodels.CrossReference.objects.get(wid = wid)
                    except ObjectDoesNotExist:
                        xref = cmodels.CrossReference(wid = wid)
                    xref.name = wid
                    xref.xid = anno["gi"]
                    xref.source = "GI"
                    xref.save(revision_detail = revdetail)
                    xref.species.add(species)
                    xref.save(revision_detail = revdetail)
                    species.cross_references.add(xref)
                    species.save(revision_detail = revdetail)
                
                try:
                    chromosome = cmodels.Chromosome.objects.get(wid = "CHROMOSOME")
                except ObjectDoesNotExist:
                    chromosome = cmodels.Chromosome(wid = "CHROMOSOME")
                chromosome.name = "Chromosome"
                chromosome.sequence = record.seq
                chromosome.length = len(record.seq)
                chromosome.save(revision_detail = revdetail)
                chromosome.species.add(species)
                chromosome.save(revision_detail = revdetail)
                
                features = record.features
                
                gene_features = filter(lambda x : x.type == "gene", features)
                cds_features = filter(lambda x : x.type != "source" and x.type != "gene", features)
                
                gene_map = {}
                for g in gene_features:
                    if not "locus_tag" in g.qualifiers:
                        print "WARN: " + str(g) + " without locus"
                        continue
                    loci = g.qualifiers["locus_tag"][0]
                    if loci in gene_map:
                        raise ArgumentError("locus_tag " + loci + " appeared twice")
                    gene_map[loci] = g
                
                cds_map = {}
                for c in cds_features:
                    if not "locus_tag" in c.qualifiers:
                        print "WARN: " + str(c) + " without locus"
                        continue
                    loci = c.qualifiers["locus_tag"][0]
                    if loci in cds_map:
                        raise ArgumentError("locus_tag " + loci + " appeared twice")
                    if loci in gene_map:
                        cds_map[loci] = c
                
                for i, v in enumerate(cds_map.values()):
                    qualifiers = v.qualifiers
                    try:
                        g = cmodels.Gene.objects.get(wid = qualifiers["locus_tag"][0])
                    except ObjectDoesNotExist:
                        g = cmodels.Gene(wid = qualifiers["locus_tag"][0])

                    g.chromosome = chromosome

                    if "gene" in qualifiers:
                        g.name = qualifiers["gene"][0]
                        g.symbol = qualifiers["gene"][0]

                    g.direction = 'f' if v.location.strand == 1 else 'r'
                    
                    # __len__ because len() fails for numbers < 0
                    # Joins output the wrong length
                    if v.location.__len__() < 0:
                        g.length = v.location.__len__() + len(record.seq)
                    else:
                        g.length = len(v.location)
                    
                    g.coordinate = v.location.start if g.direction == 'f' else v.location.end

                    if "note" in qualifiers:
                        g.comments = "\n".join(qualifiers["note"])
                    
                    g.save(revision_detail = revdetail)
                    
                    if "db_xref" in qualifiers:
                        for xref in qualifiers["db_xref"]:
                            if ":" in xref:
                                source, xid = xref.split(":")
                                wid = source + ":" + xid
                                try:
                                    xref = cmodels.CrossReference.objects.get(wid = wid)
                                except ObjectDoesNotExist:
                                    xref = cmodels.CrossReference(wid = wid)
                                xref.name = wid
                                xref.xid = xid
                                xref.source = source
                                                                
                                xref.save(revision_detail = revdetail)
                                xref.species.add(species)
                                g.cross_references.add(xref)
                                xref.save(revision_detail = revdetail)
                    
                    if "EC_number" in qualifiers:
                        for ec in qualifiers["EC_number"]:
                            wid = "EC" + ":" + ec
                            try:
                                xref = cmodels.CrossReference.objects.get(wid = wid)
                            except ObjectDoesNotExist:
                                xref = cmodels.CrossReference(wid = wid)
                            xref.name = wid
                            xref.source = "EC"
                            xref.xid = ec
                            xref.save(revision_detail = revdetail)
                            xref.species.add(species)
                            g.cross_references.add(xref)
                            xref.save(revision_detail = revdetail)
                    
                    g.save(revision_detail = revdetail)
                    g.species.add(species)
                    
                    if v.type == "CDS":
                        v.type = "mRNA"
                    
                    try:
                        t = cmodels.Type.objects.get(wid = v.type)
                    except ObjectDoesNotExist:
                        t = cmodels.Type(wid = v.type, name = v.type)

                    t.save(revision_detail = revdetail)
                    t.species.add(species)
                    t.save(revision_detail = revdetail)

                    g.type.add(t)
                    
                    g.save(revision_detail = revdetail)
                    sys.stdout.write("Importing Gene %s (%d/%d)%10s\r" % (g.wid, i + 1, len(cds_map.values()), " "))

