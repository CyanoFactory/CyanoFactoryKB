"""
Whole-cell project settings

Author: Jonathan Karr, jkarr@stanford.edu
Affiliation: Covert Lab, Department of Bioengineering, Stanford University
Last updated: 2012-07-17

Cyanofactory project settings

Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

import os
from sys import argv

import djcelery

from settings_private import *

UNIT_TEST_RUNNING = 'test' in argv

# Special case for unit testing (SQLite DB -> faster)

if UNIT_TEST_RUNNING:
    DATABASES = {'default': {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': 'cyanobase'
    }}

    TEST_RUNNER = 'django_nose.NoseTestSuiteRunner'

    DEBUG = False
else:
    # Workaround for transaction blocks (progress not visible)
    # Needs the database router
    DATABASES['djcelery'] = DATABASES['default']

ROOT_DIR = os.path.dirname(os.path.realpath(__file__))

MANAGERS = ADMINS

# Local time zone for this installation. Choices can be found here:
# http://en.wikipedia.org/wiki/List_of_tz_zones_by_name
# although not all choices may be available on all operating systems.
# On Unix systems, a value of None will cause Django to use the same
# timezone as the operating system.
# If running in a Windows environment this must be set to the same as your
# system time zone.
TIME_ZONE = 'Europe/London'

# Language code for this installation. All choices can be found here:
# http://www.i18nguy.com/unicode/language-identifiers.html
LANGUAGE_CODE = 'en-us'

SITE_ID = 1

# If you set this to False, Django will make some optimizations so as not
# to load the internationalization machinery.
USE_I18N = True

# If you set this to False, Django will not format dates, numbers and
# calendars according to the current locale
USE_L10N = True

# Absolute filesystem path to the directory that will hold user-uploaded files.
# Example: "/home/media/media.lawrence.com/media/"
MEDIA_ROOT = ''

# URL that handles the media served from MEDIA_ROOT. Make sure to use a
# trailing slash.
# Examples: "http://media.lawrence.com/media/", "http://example.com/media/"
MEDIA_URL = ''

# Absolute path to the directory static files should be collected to.
# Don't put anything in this directory yourself; store your static files
# in apps' "static/" subdirectories and in STATICFILES_DIRS.
# Example: "/home/media/media.lawrence.com/static/"
STATIC_ROOT = ''

# Additional locations of static files, use absolute paths
STATICFILES_DIRS = (
    ROOT_DIR + "/static",
)

# List of finder classes that know how to find static files in
# various locations.
STATICFILES_FINDERS = (
    'django.contrib.staticfiles.finders.FileSystemFinder',
    'django.contrib.staticfiles.finders.AppDirectoriesFinder',
    #'django.contrib.staticfiles.finders.DefaultStorageFinder',
)

# List of callables that know how to import templates from various sources.
TEMPLATE_LOADERS = (
    'django.template.loaders.filesystem.Loader',
    'django.template.loaders.app_directories.Loader',
    #'django.template.loaders.eggs.Loader',
)

MIDDLEWARE_CLASSES = (
    'cyano.middleware.RedirectAllowedHostMiddlware',
    #'johnny.middleware.LocalStoreClearMiddleware',
    #'johnny.middleware.QueryCacheMiddleware',
    #'django.middleware.cache.UpdateCacheMiddleware',
    'django.middleware.common.CommonMiddleware',
    #'django.middleware.cache.FetchFromCacheMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'cyano.middleware.PasswordChangeMiddleware',
    #'django.contrib.messages.middleware.MessageMiddleware',
)

#if DEBUG:
#    MIDDLEWARE_CLASSES += ('debug_toolbar.middleware.DebugToolbarMiddleware',)

ROOT_URLCONF = 'urls'

#use absolute directories
TEMPLATE_DIRS = (
    ROOT_DIR + "/templates",
)

HAYSTACK_CONNECTIONS = {
    'default': {
        'ENGINE': 'haystack.backends.simple_backend.SimpleEngine',
        'PATH': os.path.join(os.path.dirname(__file__), 'xapian_index'),
    },
}

djcelery.setup_loader()

SOUTH_TESTS_MIGRATE = False

INSTALLED_APPS = (
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.sites',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'django.contrib.admin',
    'django.contrib.admindocs',
    'django_extensions',
    'widget_tweaks',

    #apps
    'public',
    #'biosql',
    'cyano',
    'cyanodesign',
    'db_xref',
    #'biowarehouse',
    'bioparser',
    'boehringer',
    'kegg',
    'cyanointeraction',

    #helpers
    'django_dumpdb',
    'djcelery',
    'haystack',
    'endless_pagination',
    'rest_framework',
    'south',
    #'debug_toolbar'
)


if UNIT_TEST_RUNNING:
    INSTALLED_APPS += ('django_nose',)


# A sample logging configuration. The only tangible logging
# performed by this configuration is to send an email to
# the site admins on every HTTP 500 error.
# See http://docs.djangoproject.com/en/dev/topics/logging for
# more details on how to customize your logging configuration.

LOGGING = {
    'version': 1,
    'disable_existing_loggers': False,
    'filters': {
        'require_debug_false': {
            '()': 'django.utils.log.RequireDebugFalse'
        }
    },
    'formatters': {
        'simple': {
            'format': '%(levelname)s | %(message)s'
        },
    },
    'handlers': {
        'mail_admins': {
            'level': 'ERROR',
            'filters': ['require_debug_false'],
            'class': 'django.utils.log.AdminEmailHandler'
        },
        'console': {
            'class': 'logging.StreamHandler',
            'level': 'ERROR',
            'formatter': 'simple'
        },
    },
    'loggers': {
        'django.request': {
            'handlers': ['mail_admins'],
            'level': 'ERROR',
            'propagate': True,
        },
    }
}

CACHES = {
    'default': {
        'BACKEND': 'johnny.backends.locmem.LocMemCache',
        'LOCATION': 'cyano-cache',
        'TIMEOUT': 600,
        'OPTIONS': {
            'MAX_ENTRIES': 1000
        },
        'JOHNNY_CACHE': False,
    }
}

CACHE_MIDDLEWARE_ALIAS = "default"
CACHE_MIDDLEWARE_SECONDS = 600
CACHE_MIDDLEWARE_KEY_PREFIX = "cyano-pagecache"

INTERNAL_IPS = ('127.0.0.1',)

AUTH_PROFILE_MODULE = 'cyano.UserProfile'
LOGIN_URL = ROOT_URL + '/login/'
LOGIN_REDIRECT_URL = ROOT_URL + '/'

DATABASE_ROUTERS = ['public.router.WarehouseRouter']

ABSOLUTE_URL_OVERRIDES = {
    'auth.user': lambda o: ROOT_URL + '/user/' + o.username,
}

from django.conf.global_settings import TEMPLATE_CONTEXT_PROCESSORS

TEMPLATE_CONTEXT_PROCESSORS += (
    'django.core.context_processors.request',
    'cyano.context_processor.process',
)

GOOGLE_SEARCH_ENABLED = False

CELERY_ACCEPT_CONTENT = ['pickle', 'json', 'msgpack', 'yaml']
CELERY_RESULT_BACKEND = "djcelery.backends.database:DatabaseBackend"

JOHNNY_TABLE_BLACKLIST = ["djcelery_crontabschedule", "djcelery_intervalschedule", "djcelery_periodictask", "djcelery_periodictasks", "djcelery_taskstate", "djcelery_workerstate", "celery_taskmeta", "celery_tasksetmeta"]

REST_FRAMEWORK = {
      'DEFAULT_RENDERER_CLASSES': (
        'rest_framework.renderers.JSONRenderer',
        'rest_framework.renderers.XMLRenderer',
        'rest_framework.renderers.BrowsableAPIRenderer',
      ),
    'DEFAULT_FILTER_BACKENDS': ('rest_framework.filters.DjangoFilterBackend',)
}
