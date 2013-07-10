from optparse import make_option
from django.core.management.base import BaseCommand, CommandError
from django.core.exceptions import ObjectDoesNotExist
import kegg.models as kmodels
import cyano.models as cmodels
from cyano.helpers import slugify

class Command(BaseCommand):
    help = 'Takes the data from the Kegg tables and registers it with the species'
    
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
        
        revdetail = cmodels.RevisionDetail()
        revdetail.user = cmodels.UserProfile.objects.get(user__username__exact = "management")
        revdetail.reason = options["reason"]
        
        species_obj = cmodels.Species.objects.get(wid = slugify(options["wid"]))
        
        crs = cmodels.CrossReference.objects.filter(species = species_obj, source = "EC").values_list('xid', flat=True)
        maps = kmodels.Map.objects.filter(ec_numbers__name__in = crs).distinct()
        
        for map_ in maps:
            try:
                x = cmodels.Pathway.objects.get(wid = map_.name)
            except ObjectDoesNotExist:
                x = cmodels.Pathway(wid = map_.name)
            
            x.name = map_.title
            x.save(revision_detail = revdetail)
            x.species.add(species_obj)
            x.save(revision_detail = revdetail)

