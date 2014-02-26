"""
Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

from optparse import make_option
from cyano.tasks import fastafeature
from cyano_command import CyanoCommand
from django.core.management.base import BaseCommand

class Command(BaseCommand):
    args = '<file file ...>'
    help = 'Imports FASTA files containing chromosome features'
    
    option_list = CyanoCommand.option_list + (
        make_option('--chromosome', '-c',
                    action='store',
                    dest='chromosome',
                    default=False,
                    help='Name of the Chromosome (or Plasmid) to assign the feature to.'),
        make_option('--feature-type', '-t',
                    action='store',
                    dest='feature-type',
                    default=False,
                    help='Type of the feature')
        )

    def handle(self, *args, **options):
        for arg in args:
            fastafeature.delay(filename = arg, wid = options["wid"], user = options["user"], reason = options["reason"],
                chromosome = options["chromosome"], feature_type = options["feature-type"])
