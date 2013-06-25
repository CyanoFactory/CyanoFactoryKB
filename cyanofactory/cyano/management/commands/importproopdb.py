# Imports Operon Prediction files from ProOpDB
# Reference:
# Taboada B., Ciria R., Martinez-Guerrer C.E., Merino E., 2012, ProOpDB:
# Prokaryotic Operon DataBase, Nucleic Acids Research, 40(D1), D627-D631

from django.core.management.base import BaseCommand, CommandError
import cyano.models as cmodels
import biowarehouse.models as bmodels
from django.core.exceptions import ObjectDoesNotExist
from argparse import ArgumentError
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
        
        self.stdout.write(str(options["wid"]))
        pass
        #for arg in args:
        #    print "Parsing %s" % (arg)
        #    with open(arg, "r") as handle:
        #        f = SeqIO.parse(handle, "genbank")