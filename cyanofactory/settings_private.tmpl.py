# Rename to settings_private.py

"""
Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

ROOT_URL = 'http://localhost:8000'

DEBUG = True
TEMPLATE_DEBUG = DEBUG

# URL prefix for static files.
# If loading static files fails try prepending ROOT_URL
STATIC_URL = '/static/'

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2',
          'NAME': 'name',
          'USER': 'user',
          'PASSWORD': 'passwort',
          'HOST': 'example.com',
          'PORT': '5432',
    }
}

# Generate a SECRET KEY for Django and place it here
SECRET_KEY = ''

# Host Header values allowed to access the site when DEBUG is False
ALLOWED_HOSTS = ['localhost:8000']

# Used by the redirect host middleware
REDIRECT_URL = "https://cyanofactory.hs-mittweida.de"

# Login url to a broker for the job system
BROKER_URL = 'amqp://user:password@cyanofactory.hs-mittweida.de:5672/cyanofactory-dev'
