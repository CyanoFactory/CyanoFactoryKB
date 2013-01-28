import django.core.handlers.wsgi
import os
import sys

paths = [
	os.path.dirname(os.path.realpath(__file__)) + '/../..',
	os.path.dirname(os.path.realpath(__file__)) + '/..',
	]
for path in paths:
	if path not in sys.path:
		sys.path.append(path)

os.environ['DJANGO_SETTINGS_MODULE'] = 'knowledgebase_2.settings'

application = django.core.handlers.wsgi.WSGIHandler()

import monitor
monitor.start(interval = 1.0)
