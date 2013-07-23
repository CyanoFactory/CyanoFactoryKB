from django.core.exceptions import ObjectDoesNotExist

import cyano.models as cmodels
from cyano_command import CyanoCommand

class Command(CyanoCommand):
    help = 'Adds a reference to boehringer pathway to the species'

    def handle_command(self, species, revdetail, *args, **options):
        wid = "Boehringer"
        try:
            x = cmodels.Pathway.objects.get(wid = wid, species = species)
            # No error: Already exists -> abort
            self.stderr.write("Boehringer Pathway already assigned")
            return
        except ObjectDoesNotExist:
            # Create new boehringer (if it was never created yet) or assign to one
            try:
                x = cmodels.Pathway.objects.get(wid = wid)
            except ObjectDoesNotExist:
                x = cmodels.Pathway(wid = wid)

        x.name = "Biochemical Pathways"
        x.save(revision_detail = revdetail)
        x.species.add(species)
        x.save(revision_detail = revdetail)

        self.stdout.write("Boehringer Pathway assigned")
