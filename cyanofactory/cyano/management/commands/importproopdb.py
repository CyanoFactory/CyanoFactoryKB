# Imports Operon Prediction files from ProOpDB
# Reference:
# Taboada B., Ciria R., Martinez-Guerrer C.E., Merino E., 2012, ProOpDB:
# Prokaryotic Operon DataBase, Nucleic Acids Research, 40(D1), D627-D631

from django.core.management.base import BaseCommand, CommandError
import cyano.models as cmodels
from django.core.exceptions import ObjectDoesNotExist
from optparse import make_option
    
class Command(BaseCommand):
    args = '<file file ...>'
    help = 'Imports Operon Prediction Files retrieved from ProOpDB'
    
    option_list = BaseCommand.option_list + (
        make_option('--wid', '-w',
            action='store',
            dest='wid',
            default=False,
            help='WID of the target species'),
        )
    
    def handle(self, *args, **options):
        if not options["wid"]:
            raise CommandError("wid argument is mandatory")
        
        wid = options["wid"]
        try:
            species = cmodels.Species.objects.get(wid = wid)
        except ObjectDoesNotExist:
            raise CommandError("No species with wid {} exists".format(wid))
        
        for arg in args:
            with file(arg, "r") as f:
                revdetail = cmodels.RevisionDetail()
                revdetail.user = cmodels.UserProfile.objects.get(user__username__exact = "management")
                revdetail.reason = "Import ProOpDB File " + arg
                
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
                        tu.save(revision_detail = revdetail)
                        tu.species.add(species)
                        tu.genes.add(gene)
                        tu.save(revision_detail = revdetail)
                        self.stdout.write("Importing TU %s (%d/%d)" % (tu.wid, i + 1, count))
