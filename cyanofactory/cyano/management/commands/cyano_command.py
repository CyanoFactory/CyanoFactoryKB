from optparse import make_option
from django.core.management.base import BaseCommand, CommandError
from cyano.helpers import slugify
from cyano.models import Species, RevisionDetail, UserProfile

class CyanoCommand(BaseCommand):
    """
    Verifies arguments of all Cyanofactory management commands
    """
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
    
    verify_species_exists = True
    
    def handle(self, *args, **options):
        if not options["wid"]:
            raise CommandError("wid argument is mandatory")
        
        if not options["reason"]:
            raise CommandError("reason is mandatory")
        
        wid = slugify(options["wid"])
        reason = options["reason"]
        
        if options["wid"] != wid:
            raise CommandError("Wid {} contained invalid characters. Only letters, numbers and _ are allowed".format(wid))            
        
        try:
            species_obj = Species.objects.get(wid = wid)
        except:
            if self.verify_species_exists:
                raise CommandError("Species {} not found".format(wid))
            else:
                species_obj = Species(wid = wid)
        
        revdetail = RevisionDetail()
        revdetail.user = UserProfile.objects.get(user__username__exact = "management")
        revdetail.reason = reason
        
        self.handle_command(species_obj, revdetail, *args, **options)
        
    def handle_command(self, species, revision, *args, **options):
        raise CommandError("handle_command no implement")
