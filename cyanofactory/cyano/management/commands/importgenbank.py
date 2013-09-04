from optparse import make_option
from cyano.tasks import genbank
from cyano_command import CyanoCommand
from django.core.management.base import BaseCommand

class Command(BaseCommand):
    args = '<file file ...>'
    help = 'Imports NCBI GenBank Files'
    
    option_list = CyanoCommand.option_list + (
        make_option('--chromosome', '-c',
                    action='store',
                    dest='chromosome',
                    default=False,
                    help='Name of the Chromosome (or Plasmid) to assign the GenBank data to. Created if necessary.'),
        make_option('--name', '-n',
                    action='store',
                    dest='name',
                    default=False,
                    help='Human readable name of the species')
        )

    def handle(self, *args, **options):
        for arg in args:
            genbank.delay(filename = arg, wid = options["wid"], user = options["user"], reason = options["reason"],
                chromosome = options["chromosome"], name = options["name"])
