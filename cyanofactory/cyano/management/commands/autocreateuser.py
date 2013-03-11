from django.core.management.base import BaseCommand
from cyano.models import UserProfile
from django.contrib.auth.models import User

class Command(BaseCommand):
    def handle(self, *args, **options):
        for name in args:
            user = UserProfile()
            print "Creating user " + name
            mail = name + "@example.com"
            password = "aaa"
            user.user = User.objects.create_user(username=name, email=mail,  password=password)
            user.save()
