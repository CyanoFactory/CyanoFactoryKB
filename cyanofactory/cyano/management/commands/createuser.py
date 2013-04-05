from django.core.management.base import BaseCommand
from cyano.models import UserProfile
from django.contrib.auth.models import User

class Command(BaseCommand):
    def handle(self, *args, **options):
        for name in args:
            user = UserProfile()
            print "Creating user " + name
            mail = raw_input("E-Mail: ")
            password = raw_input("Password: ")
            user.affiliation = raw_input(user._meta.get_field_by_name('affiliation')[0].verbose_name + ": ")
            user.website = raw_input(user._meta.get_field_by_name('website')[0].verbose_name + ": ")
            user.phone = raw_input(user._meta.get_field_by_name('phone')[0].verbose_name + ": ")
            user.address = raw_input(user._meta.get_field_by_name('address')[0].verbose_name + ": ")
            user.city = raw_input(user._meta.get_field_by_name('city')[0].verbose_name + ": ")
            user.state = raw_input(user._meta.get_field_by_name('state')[0].verbose_name + ": ")
            user.zip = raw_input(user._meta.get_field_by_name('zip')[0].verbose_name + ": ")
            user.country = raw_input(user._meta.get_field_by_name('country')[0].verbose_name + ": ")
            user.user = User.objects.create_user(username=name, email=mail,  password=password)
            user.save()
