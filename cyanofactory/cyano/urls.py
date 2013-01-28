from django.conf.urls import patterns, include, url
from django.views.generic import DetailView, ListView

_accession = r'(?P<organism>[A-Z]{2}_?\d{1,10})'

urlpatterns = patterns('cyano.views',
    url(r'^$', 'index'),
    url(r'^' + _accession + r'/$', 'organism'),
    url(r'^' + _accession + r'/genes$', 'genes'),
    url(r'^' + _accession + r'/genes/(?P<gene>.+)', 'gene'),
    ##url(r'^' + _accession + r'/genes', 'genes'),
#    url(r'^(?P<poll_id>\d+)/$', 'detail'),
#    url(r'^(?P<poll_id>\d+)/results/$', 'results'),
#    url(r'^(?P<poll_id>\d+)/vote/$', 'vote'),
)
