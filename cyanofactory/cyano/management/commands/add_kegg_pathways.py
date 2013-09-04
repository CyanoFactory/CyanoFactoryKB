from cyano.models import Pathway
from cyano_command import CyanoCommand

class Command(CyanoCommand):
    help = 'Takes the data from the Kegg tables and registers it with the species'

    def handle_command(self, species, revdetail, *args, **options):        
        Pathway.add_kegg_pathway(species, revdetail)
