"""
Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

# Imports Operon Prediction files from ProOpDB
# Reference:
# Taboada B., Ciria R., Martinez-Guerrer C.E., Merino E., 2012, ProOpDB:
# Prokaryotic Operon DataBase, Nucleic Acids Research, 40(D1), D627-D631

from cyano_command import CyanoCommand
from cyano.tasks import proopdb
from django.core.management.base import BaseCommand
    
class Command(BaseCommand):
    args = '<file file ...>'
    help = 'Imports Operon Prediction Files retrieved from ProOpDB'
    
    option_list = CyanoCommand.option_list

    def handle(self, *args, **options):
        for arg in args:
            proopdb.delay(filename = arg, wid = options["wid"], user = options["user"], reason = options["reason"])
