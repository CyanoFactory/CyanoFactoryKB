from django.conf.urls import patterns, url

urlpatterns = patterns('db_xref.views',
    url(r'^$', 'index'),
    url(r'^(?P<database>[a-zA-Z0-9_-]+)\s*:\s*(?P<organism>[a-zA-Z0-9_-]+)$', 'dbxref'),
)
