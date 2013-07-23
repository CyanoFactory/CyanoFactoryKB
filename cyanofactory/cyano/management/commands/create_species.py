from optparse import make_option

from django.core.management.base import CommandError

from cyano_command import CyanoCommand

class Command(CyanoCommand):
    help = 'Creates a new species'
    
    """
    Creates a new species (or alters the name when found)
    """
    option_list = CyanoCommand.option_list + (
        make_option('--name', '-n',
                    action='store',
                    dest='name',
                    default=False,
                    help='Human readable name of the species'),)
    #                Tuple syntax: Don't remove that comma ---^

    verify_species_exists = False
    
    def handle_command(self, species_obj, revdetail, *args, **options):
        if not options["name"]:
            raise CommandError("name argument is mandatory")
 
        species_obj.name = options["name"]
        species_obj.save(revision_detail = revdetail)
