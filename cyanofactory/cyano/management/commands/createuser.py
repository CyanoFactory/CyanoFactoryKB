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
            user.affiliation = input(user._meta.get_field('affiliation').verbose_name + ": ")
            user.website = input(user._meta.get_field('website').verbose_name + ": ")
            user.phone = input(user._meta.get_field('phone').verbose_name + ": ")
            user.address = input(user._meta.get_field('address').verbose_name + ": ")
            user.city = input(user._meta.get_field('city').verbose_name + ": ")
            user.state = input(user._meta.get_field('state').verbose_name + ": ")
            user.zip = input(user._meta.get_field('zip').verbose_name + ": ")
            user.country = input(user._meta.get_field('country').verbose_name + ": ")
            user.user = User.objects.create_user(username=name, email=mail,  password=password)
            user.save()
