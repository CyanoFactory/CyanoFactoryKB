"""
Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

from django.core.management.base import BaseCommand
from django.contrib.auth.models import User

from cyano.models import UserProfile

class Command(BaseCommand):
    def handle(self, *args, **options):
        for name in args:
            user = UserProfile()
            print("Creating user " + name)
            mail = input("E-Mail: ")
            password = input("Password: ")
            user.affiliation = input(user._meta.get_field_by_name('affiliation')[0].verbose_name + ": ")
            user.website = input(user._meta.get_field_by_name('website')[0].verbose_name + ": ")
            user.phone = input(user._meta.get_field_by_name('phone')[0].verbose_name + ": ")
            user.address = input(user._meta.get_field_by_name('address')[0].verbose_name + ": ")
            user.city = input(user._meta.get_field_by_name('city')[0].verbose_name + ": ")
            user.state = input(user._meta.get_field_by_name('state')[0].verbose_name + ": ")
            user.zip = input(user._meta.get_field_by_name('zip')[0].verbose_name + ": ")
            user.country = input(user._meta.get_field_by_name('country')[0].verbose_name + ": ")
            user.user = User.objects.create_user(username=name, email=mail,  password=password)
            user.save()
