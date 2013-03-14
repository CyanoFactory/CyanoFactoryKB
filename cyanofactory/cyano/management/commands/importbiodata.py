from django.core.management.base import BaseCommand
from Bio import SeqIO
from cyano.models import Entry, Species, RevisionDetail, UserProfile, Gene, Chromosome
from django.core.exceptions import ObjectDoesNotExist
from argparse import ArgumentError

class Command(BaseCommand):
    def handle(self, *args, **options):
        for arg in args:
            print arg
            handle = open(arg, "r")
            f = SeqIO.parse(handle, "genbank")
            for record in f:
                revdetail = RevisionDetail()
                revdetail.user = UserProfile.objects.get(user__username__exact = "gabriel")
                revdetail.reason = "Initial Load"

                try:
                    species = Species.objects.get(wid = record.name)
                except ObjectDoesNotExist:
                    species = Species()
                species.wid = record.name
                species.name = record.annotations["organism"]
                species.comments = record.description
                species.genetic_code = '2'
                species.save(revision_detail = revdetail)
                
                try:
                    chromosome = Chromosome.objects.get(wid = "CHR_1")
                except ObjectDoesNotExist:
                    chromosome = Chromosome()
                chromosome.wid = "CHR_1"
                chromosome.sequence = record.seq
                chromosome.length = len(record.seq)
                chromosome.save(revision_detail = revdetail)
                chromosome.species.append(species)
                chromosome.save(revision_detail = revdetail)
                
                features = record.features
                
                gene_features = filter(lambda x : x.type == "gene", features)
                cds_features = filter(lambda x : x.type == "CDS", features)
                
                gene_map = {}
                for g in gene_features:
                    loci = g.qualifiers["locus_tag"]
                    if loci in gene_map:
                        raise ArgumentError("locus_tag " + loci + " appeared twice")
                    gene_map[loci] = g
                
                cds_map = {}
                for c in cds_features:
                    loci = c.qualifiers["locus_tag"]
                    if loci in cds_map:
                        raise ArgumentError("locus_tag " + loci + " appeared twice")
                    if loci in gene_map:
                        cds_map[loci] = c
                
                for v in cds_map.values():
                    try:
                        g = Gene.objects.get(wid = g.qualifiers["locus_tag"][0])
                    except ObjectDoesNotExist:
                        g = Gene()

                    g.wid = g.qualifers["locus_tag"][0]
                    g.chromosome = chromosome
                    qualifiers = v.qualifiers
                    if "gene" in qualifiers:
                        g.symbol = qualifiers["gene"][0]
                    g.direction = 'f' if g.location.strand == 1 else 'r'
                    g.coordinate = g.location.start if g.direction == 'f' else g.location.end # FIXME Not for joins
                    g.save(revision_detail = revdetail)
                    g.species.append(g)
                    g.save(revision_detail = revdetail)
                    pass
                    
                # r.annotations
                # r.id = NC_000911.1
                # r.name = NC_000911
                # r.features
                
                # >>> a["comment"]
#'PROVISIONAL REFSEQ: This record has not yet been subject to final\nNCBI review. The reference sequence was derived from
# BA000022.\nCOMPLETENESS: full length.'
                # a["sequence_version"] = 1
                # a["source"] = 'Synechocystis sp. PCC 6803'
                # >>> a["taxonomy"] = ['Bacteria', 'Cyanobacteria', 'Chroococcales', 'Synechocystis']
                # a["references"] -> Liste -> authors, comment, title, journal
                # >>> a["accessions"] = ['NC_000911']
                # >>> a["date"] = '21-DEC-2012'
                # a["gi"] = 16329170
                
                #
                
                


