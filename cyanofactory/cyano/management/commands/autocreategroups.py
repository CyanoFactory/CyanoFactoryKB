from django.core.management.base import BaseCommand
from django.contrib.auth.models import Group

from cyano.models import GroupProfile

class Command(BaseCommand):
    """Creates general groups
    """
    def handle(self, *args, **options):
        print "Creating groups"

        everybody, _ = Group.objects.get_or_create(name='Everybody')        
        registred, _ = Group.objects.get_or_create(name='Registred')        
        guest, _ = Group.objects.get_or_create(name='Guest')
        admin, _ = Group.objects.get_or_create(name='Administrator')
        cyano_leader, _ = Group.objects.get_or_create(name='Cyanofactory Leader')
        cyano_member, _ = Group.objects.get_or_create(name='Cyanofactory Member')
        non_cyano_member, _ = Group.objects.get_or_create(name='Non-Cyanofactory Member')        
        mw, _= Group.objects.get_or_create(name='Mittweida')
        
        GroupProfile.objects.get_or_create \
            (group=everybody, description="Any visitor (guest and registred)")
        
        GroupProfile.objects.get_or_create \
            (group=registred, description="Any user that is logged in")

        GroupProfile.objects.get_or_create \
            (group=guest, description="Any user that is not logged in")

        GroupProfile.objects.get_or_create \
            (group=admin, description="Super User. Can modify internal settings of the warehouse.")

        GroupProfile.objects.get_or_create \
            (group=cyano_leader, description="Workpackage Leader")
        
        GroupProfile.objects.get_or_create \
            (group=cyano_member, description="Works on the Cyanofactory project")
        
        GroupProfile.objects.get_or_create \
            (group=non_cyano_member, description="Not a member of the Cyanofactory project")
        
        GroupProfile.objects.get_or_create \
            (group=mw, description="From University of Applied Sciences Mittweida")
        
        
