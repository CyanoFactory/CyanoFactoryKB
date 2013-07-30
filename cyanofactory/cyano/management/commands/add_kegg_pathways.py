from django.core.exceptions import ObjectDoesNotExist

import kegg.models as kmodels

import cyano.models as cmodels
from cyano_command import CyanoCommand

class Command(CyanoCommand):
    help = 'Takes the data from the Kegg tables and registers it with the species'

    def handle_command(self, species, revdetail, *args, **options):        
        crs = cmodels.CrossReference.objects.filter(species = species, source = "EC").values_list('xid', flat=True)
        maps = kmodels.Map.objects.filter(ec_numbers__name__in = crs).distinct()
        
        for map_ in maps:
            try:
                x = cmodels.Pathway.objects.get(wid = map_.name)
            except ObjectDoesNotExist:
                x = cmodels.Pathway(wid = map_.name)
            
            x.name = map_.title
            x.save(revdetail)
            x.species.add(species)
            x.save(revdetail)
