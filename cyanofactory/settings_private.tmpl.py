# Rename to settings_private.py

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
