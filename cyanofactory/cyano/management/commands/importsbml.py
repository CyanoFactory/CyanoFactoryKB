"""
Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

from cyano_command import CyanoCommand
from cyano.tasks import sbml
from django.core.management.base import BaseCommand

class Command(BaseCommand):  
    args = '<file file ...>'
    help = 'Imports SBML Files'
    
    option_list = CyanoCommand.option_list

    def handle(self, *args, **options):
        for arg in args:
            sbml.delay(filename = arg, wid = options["wid"], user = options["user"], reason = options["reason"])
