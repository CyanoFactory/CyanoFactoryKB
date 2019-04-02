import os
import sys
from django.core.wsgi import get_wsgi_application
import djcelery

paths = [
        os.path.dirname(os.path.realpath(__file__)) + '/../..',
        os.path.dirname(os.path.realpath(__file__)) + '/..',
        ]
for path in paths:
        if path not in sys.path:
                sys.path.append(path)

os.environ['DJANGO_SETTINGS_MODULE'] = 'settings'

djcelery.setup_loader()

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "settings")

application = get_wsgi_application()
