from cyano_command import CyanoCommand
from cyano.models import Pathway

class Command(CyanoCommand):
    help = 'Adds a reference to boehringer pathway to the species'

    def handle_command(self, species, revdetail, *args, **options):
        Pathway.add_boehringer_pathway(species, revdetail)

        self.stdout.write("Boehringer Pathway assigned")
