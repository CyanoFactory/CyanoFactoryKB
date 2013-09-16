from django.core.management.base import BaseCommand
from django.contrib.auth.models import User

from django.contrib.sites.models import Site
from django.contrib.sessions.models import Session
from django.template import loader
from django.core.mail import send_mail

import string
import random
import settings

def all_sessions_for_user(user):
    user_sessions = []
    all_sessions = Session.objects.all()
    for session in all_sessions:
        session_data = session.get_decoded()
        if user.pk == session_data.get('_auth_user_id'):
            user_sessions.append(session.pk)
    return Session.objects.filter(pk__in=user_sessions)

# via http://stackoverflow.com/questions/2257441/
def id_generator(size=8, chars=string.ascii_uppercase + string.ascii_lowercase + string.digits):
    return ''.join(random.choice(chars) for x in range(size))

class Command(BaseCommand):
    help = 'Randomly changes the password of the provided users and mails it to them'
    
    def handle(self, *args, **options):
        site = Site.objects.get(pk = 1)
        site_name = site.name
        domain = site.domain
        
        for name in args:
            if not User.objects.filter(username=name).count():
                print "User {} not found".format(name)

            else:
                user = User.objects.get(username=name)                
                new_pass = id_generator()
                user.set_password(new_pass)
                user.profile.force_password_change = True
                
                from_email = getattr(settings, 'DEFAULT_FROM_EMAIL', "webmaster@localhost")
                
                c = {
                    'email': user.email,
                    'domain': domain,
                    'site_name': site_name,
                    'user': user,
                    'protocol': False and 'https' or 'http',
                    'password': new_pass
                }

                subject = loader.render_to_string("registration/password_created_subject.txt", c)
                # Email subject *must not* contain newlines
                subject = ''.join(subject.splitlines())
                email = loader.render_to_string("registration/password_created_email.html", c)
                
                # Invalidate active sessions (in case of an account hijack)
                all_sessions_for_user(user).delete()

                send_mail(subject, email, from_email, [user.email])
                
                user.save()
                user.profile.save()

                print "Sent new password for {} to {}".format(name, user.email)
