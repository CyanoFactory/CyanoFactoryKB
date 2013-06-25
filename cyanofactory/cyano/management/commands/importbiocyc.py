from django.core.management.base import BaseCommand
import cyano.models as cmodels
import biowarehouse.models as bmodels
from django.core.exceptions import ObjectDoesNotExist
from argparse import ArgumentError

class Command(BaseCommand):
    def handle(self, *args, **options):
        files = ["NC_000911.1_Chromosome.fasta", "NC_005229.1_Plasmid-1.fasta",
                 "NC_005232.1_Plasmid-2.fasta", "NC_005230.1_Plasmid-3.fasta",
                 "NC_005231.1_Plasmid-4.fasta"]
        chr_list = ["CHROMOSOME-1"] + list("PLASMID-" + str(x) for x in range(1,5))
        
        proteins = bmodels.Protein.objects.all()
        genes = bmodels.Gene.objects.all()
        tuc = bmodels.Transcriptionunitcomponent.objects.all()
        pathways = bmodels.Pathway.objects.all()
        
        gene_to_replicon = \
            lambda x: "CHROMOSOME-1" if x >= 1 and x <= 3230  \
                else "PLASMID-1" if x <= 3363 \
                else "PLASMID-2" if x <= 3374 \
                else "PLASMID-3" if x <= 3580 \
                else "PLASMID-4" if x <= 3630 \
                else list()[1] # raise IndexError
        
        revdetail = cmodels.RevisionDetail()
        revdetail.user = cmodels.UserProfile.objects.get(user__username__exact = "management")
        revdetail.reason = "Biocyc Pathways"
        
        try:
            species = cmodels.Species.objects.get(wid = "BioCyc_PPC6803")
        except ObjectDoesNotExist:
            species = cmodels.Species(wid = "BioCyc_PPC6803")
        species.name = "Synechosystis PPC6803 BioCyc"
        species.comments = ""
        species.genetic_code = '11'
        species.save(revision_detail = revdetail)
        
        #for protein in proteins:
        #    try:
        #        p = cmodels.Protein.objects.get(wid = protein.wid)
        #    except ObjectDoesNotExist:
        #        p = cmodels.Protein(wid = protein.wid)
        #    p.name = protein.name
        #    p.save(revision_detail = revdetail)
        #    p.species.add(species)
        #    p.save(revision_detail = revdetail)
        
        #proteins = cmodels.Protein.objects.filter(species = species)
        
        for fil, chr in zip(files, chr_list):
            f = open("../sequences/" + fil)
            
            # header
            f.readline()
            # sequence
            seq = ''.join(line.strip() for line in f)
            
            try:
                chromosome = cmodels.Chromosome.objects.get(wid = chr)
            except ObjectDoesNotExist:
                chromosome = cmodels.Chromosome(wid = chr)
 
            chromosome.name = chr
            chromosome.sequence = seq
            chromosome.length = len(seq)
            chromosome.save(revision_detail = revdetail)
            chromosome.species.add(species)
            chromosome.save(revision_detail = revdetail)
 
            for gene in genes:
                if gene_to_replicon(int(gene.genomeid.split("-")[1])) != chr:
                    continue
                
                try:
                    g = cmodels.Gene.objects.get(wid = gene.genomeid)
                except ObjectDoesNotExist:
                    g = cmodels.Gene(wid = gene.genomeid)
 
                g.name = gene.name
                g.chromosome = chromosome
                g.direction = gene.direction.lower()
                g.coordinate = gene.codingregionstart if g.direction == 'f' else gene.codingregionend
                g.length = abs(gene.codingregionstart - gene.codingregionend) # FIXME Not for joins
                g.save(revision_detail = revdetail)
                g.species.add(species)
 
                typ = gene.name.split("-")
                if len(typ) > 1:
                    if "RNA" in typ[0]:
                        typ = typ[0]
                    else:
                        typ = "mRNA"
                else:
                    typ = "mRNA"
                
                try:
                    t = cmodels.Type.objects.get(wid = typ)
                except ObjectDoesNotExist:
                    t = cmodels.Type(wid = typ, name = typ)
 
                t.save(revision_detail = revdetail)
                t.species.add(species)
                t.save(revision_detail = revdetail)
 
                g.type.add(t)
 
                g.save(revision_detail = revdetail)
 
            f.close()
 
        for t in tuc:
            if t.type != "gene":
                continue
            
            name = bmodels.Transcriptionunit.objects.get(Wid = t.transcriptionunitWid.Wid).name
            
            try:
                tu = cmodels.TranscriptionUnit.objects.get(wid = name)
            except ObjectDoesNotExist:
                tu = cmodels.TranscriptionUnit(wid = name)
 
            tu.name = name
            tu.save(revision_detail = revdetail)
            tu.species.add(species)
            
            gene = bmodels.Gene.objects.get(Wid = t.otherWid)
            cgene = cmodels.Gene.objects.get(wid = gene.genomeid)
            tu.genes.add(cgene)
            
            tu.save(revision_detail = revdetail)

        for pathway in pathways:
            wid = "P_" + str(pathway.Wid)
   
            try:
                p = cmodels.Pathway.objects.get(wid = wid)
            except ObjectDoesNotExist:
                p = cmodels.Pathway(wid = wid)
            
            p.name = pathway.name
            p.save(revision_detail = revdetail)
            p.species.add(species)
            p.save(revision_detail = revdetail)
            