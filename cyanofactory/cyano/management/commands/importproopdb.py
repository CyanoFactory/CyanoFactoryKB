# Imports Operon Prediction files from ProOpDB
# Reference:
# Taboada B., Ciria R., Martinez-Guerrer C.E., Merino E., 2012, ProOpDB:
# Prokaryotic Operon DataBase, Nucleic Acids Research, 40(D1), D627-D631

from django.core.exceptions import ObjectDoesNotExist

import cyano.models as cmodels
from cyano_command import CyanoCommand
    
class Command(CyanoCommand):
    args = '<file file ...>'
    help = 'Imports Operon Prediction Files retrieved from ProOpDB'

    def handle_command(self, species, revdetail, *args, **options):
        for arg in args:
            with file(arg, "r") as f:                
                count = 0
                for line in f:
                    if line.count("\t") == 7:
                        # 0: Operon ID
                        # 1: GeneName (unused)
                        # 2: Locus Name
                        # 3: GI (unused)
                        # 4: Strand (unused)
                        # 5: Left Pos (unused)
                        # 6: Right Pos (unused)
                        # 7: Gene Position in Operon (unused)
                        try:
                            line = line.split("\t")
                            operon = int(line[0])
                            count += 1
                        except ValueError:
                            continue
                
                f.seek(0)

                i = 0
            
                for line in f:
                    if line.count("\t") == 7:
                        try:
                            line = line.split("\t")
                            operon = int(line[0])
                            locus = line[2]
                            i += 1
                        except ValueError:
                            continue
                        
                        operon_str = "TU_%04d" % (operon)
                        
                        try:
                            tu = cmodels.TranscriptionUnit.objects.get(wid = operon_str)
                        except ObjectDoesNotExist:
                            tu = cmodels.TranscriptionUnit(wid = operon_str)
                        
                        try:
                            gene = cmodels.Gene.objects.get(wid = locus)
                        except ObjectDoesNotExist:
                            self.stderr.write("Gene with wid {} not found".format(locus))
                            continue
                        
                        tu.name = operon_str
                        tu.save(revdetail)
                        tu.species.add(species)
                        tu.genes.add(gene)
                        tu.save(revdetail)
                        self.stdout.write("Importing TU %s (%d/%d)" % (tu.wid, i, count))

