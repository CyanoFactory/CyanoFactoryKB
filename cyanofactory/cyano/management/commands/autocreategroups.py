from django.core.management.base import BaseCommand
from cyano.models import UserProfile
from cyano.models import GroupProfile
from django.contrib.auth.models import User

from django.contrib.auth.models import Group

class Command(BaseCommand):
    def handle(self, *args, **options):
        everybody, created = Group.objects.get_or_create(name='Everybody')
        if created: everybody.save()
        
        registred, created = Group.objects.get_or_create(name='Registred')
        if created: registred.save()
        
        guest, created = Group.objects.get_or_create(name='Guest')
        if created: guest.save()
        
        admin, created = Group.objects.get_or_create(name='Administrator')
        if created: admin.save()
        
        cyano_leader, created = Group.objects.get_or_create(name='Cyanofactory Leader')
        if created: cyano_leader.save()
        
        cyano_member, created = Group.objects.get_or_create(name='Cyanofactory Member')
        if created: cyano_member.save()
        
        non_cyano_member, created = Group.objects.get_or_create(name='Non-Cyanofactory Member')
        if created: non_cyano_member.save()
        
        mw, created = Group.objects.get_or_create(name='Mittweida')
        if created: mw.save()
        
        everybody_profile, created = GroupProfile.objects.get_or_create \
            (group=everybody, description="Any visitor (guest and registred)")
        if created: everybody_profile.save()
        
        registred_profile, created = GroupProfile.objects.get_or_create \
            (group=registred, description="Any user that is logged in")
        if created: registred_profile.save()
        
        guest_profile, created = GroupProfile.objects.get_or_create \
            (group=guest, description="Any user that is not logged in")
        if created: guest_profile.save()
        
        admin_profile, created = GroupProfile.objects.get_or_create \
            (group=admin, description="Super User. Can modify internal settings of the warehouse.")
        if created: admin_profile.save()
        
        cyano_leader_profile, created = GroupProfile.objects.get_or_create \
            (group=cyano_leader, description="Workpackage Leader")
        if created: cyano_leader_profile.save()
        
        cyano_member_profile, created = GroupProfile.objects.get_or_create \
            (group=cyano_member, description="Works on the Cyanofactory project")
        if created: cyano_member_profile.save()
        
        non_cyano_member_profile, created = GroupProfile.objects.get_or_create \
            (group=non_cyano_member, description="Not a member of the Cyanofactory project")
        if created: non_cyano_member_profile.save()
        
        mw_profile, created = GroupProfile.objects.get_or_create \
            (group=mw, description="From University of Applied Sciences Mittweida")
        if created: mw_profile.save()

        gabriel = User.objects.get(username="gabriel")
        gabriel.groups = [cyano_member, mw]
        gabriel.save()
