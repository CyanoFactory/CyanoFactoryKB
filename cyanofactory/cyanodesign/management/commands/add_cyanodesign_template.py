"""
Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

from optparse import make_option
from django.core.exceptions import ObjectDoesNotExist
from django.core.management.base import BaseCommand
from cyanodesign.models import DesignTemplate
import os

import metabolic_model.sbml_parser as sbml_parser

class Command(BaseCommand):
    args = 'file'
    help = 'Add template model to CellDesign'

    def add_arguments(self, parser):
        parser.add_argument('--name', '-n',
                    action='store',
                    dest='name',
                    default=False,
                    help='Human readable name of the species'),
        parser.add_argument('--description', '-d',
                    action='store',
                    dest='description',
                    default='',
                    help='Description of the model')
        parser.add_argument('--file', '-f',
                    action='store',
                    dest='file',
                    default='')

    def handle(self, *args, **options):
        #if len(args) > 0:
        sbml_handler = sbml_parser.SbmlHandler()
        sbml_parser.push_handler(sbml_handler)
        sbml_parser.parser.parse(options["file"])
        #else:
        #    return
        
        try:
            dt = DesignTemplate.objects.get(name=options["name"])
            print("Updating existing template " + dt.name)
        except ObjectDoesNotExist:
            dt = DesignTemplate(name=options["name"])
            print("Creating new template " + dt.name)

        dt.description = options["description"]
        dt.filename = os.path.basename(options["file"])
        dt.content = sbml_handler.model.to_json()
        dt.save()
