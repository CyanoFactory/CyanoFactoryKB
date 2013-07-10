from optparse import make_option
from django.core.management.base import BaseCommand, CommandError
from django.core.exceptions import ObjectDoesNotExist
import cyano.models as cmodels
from cyano.helpers import slugify

class Command(BaseCommand):
    help = 'Adds a reference to boehringer pathway to the species'
    
    option_list = BaseCommand.option_list + (
        make_option('--wid', '-w',
                    action='store',
                    dest='wid',
                    default=False,
                    help='WID of the target species'),
        make_option('--reason', '-r',
                    action='store',
                    dest='reason',
                    default=False,
                    help='Reason for this operation')
        )

    def handle(self, *args, **options):
        if not options["wid"]:
            raise CommandError("wid argument is mandatory")
        
        if not options["reason"]:
            raise CommandError("reason is mandatory")

        try:
            species_obj = cmodels.Species.objects.get(wid = slugify(options["wid"]))
        except:
            self.stderr.write("Error: Species not found")
            return
        
        revdetail = cmodels.RevisionDetail()
        revdetail.user = cmodels.UserProfile.objects.get(user__username__exact = "management")
        revdetail.reason = options["reason"]

        wid = "Boehringer"
        try:
            x = cmodels.Pathway.objects.get(wid = wid, species = species_obj)
            # No error: Already exists -> abort
            self.stderr.write("Error: Boehringer Pathway already assigned")
            return
        except ObjectDoesNotExist:
            # Create new boehringer (if it was never created yet) or assign to one
            try:
                x = cmodels.Pathway.objects.get(wid = wid)
            except ObjectDoesNotExist:
                x = cmodels.Pathway(wid = wid)

        x.name = "Biochemical Pathways"
        x.save(revision_detail = revdetail)
        x.species.add(species_obj)
        x.save(revision_detail = revdetail)

        self.stdout.write("Boehringer Pathway assigned")
