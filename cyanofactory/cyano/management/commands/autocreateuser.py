from django.core.management.base import BaseCommand
from django.contrib.auth.models import User

class Command(BaseCommand):
    def handle(self, *args, **options):
        for name in args:
            print "Creating user " + name
            mail = name + "@example.com"
            password = "aaa"
            
            if User.objects.filter(username=name).count():
                print "User already exists. Aborting."
            else:
                new_user = User.objects.create_user(username=name, email=mail, password=password)
                new_user.save()
