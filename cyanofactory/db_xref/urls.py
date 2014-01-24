from django.conf.urls import patterns, url

urlpatterns = patterns('db_xref.views',
    url(r'^$', 'index'),
    url(r'^(?P<source>[a-zA-Z0-9_\-]+)\s*:\s*(?P<xid>[a-zA-Z0-9_.\-:.+]+)$', 'dbxref'),
)
