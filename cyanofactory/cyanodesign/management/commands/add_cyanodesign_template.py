"""
Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

from optparse import make_option
from django.core.exceptions import ObjectDoesNotExist
from django.core.management.base import BaseCommand
from cyanodesign.models import DesignTemplate
from cyanodesign.json_model import *
from PyNetMet2.metabolism import Metabolism
import os

class Command(BaseCommand):
    args = 'file'
    help = 'Add template model to CyanoDesign'

    option_list = BaseCommand.option_list + (
        make_option('--chromosome', '-c',
                    action='store',
                    dest='chromosome',
                    default=False,
                    help='Name of the Chromosome (or Plasmid) to assign the GenBank data to. Created if necessary.'),
        make_option('--name', '-n',
                    action='store',
                    dest='name',
                    default=False,
                    help='Human readable name of the species'),
        make_option('--description', '-d',
                    action='store',
                    dest='description',
                    default='',
                    help='Description of the model')
    )

    def handle(self, *args, **options):
        if len(args) > 0:
            m = Metabolism(args[0])
            j = JsonModel.from_model(m)
        
        try:
            DesignTemplate.objects.get(name=options["name"])
        except ObjectDoesNotExist:
            dt = DesignTemplate(name=options["name"])

        dt.description = options["description"]
        dt.filename = os.path.basename(args[0])
        dt.content = j.to_json()
        dt.save()
