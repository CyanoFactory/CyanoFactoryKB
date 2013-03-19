from django.conf.urls import patterns, url

_accession = r'(?P<species_wid>[A-Z]{2}_?\d{1,10})'

urlpatterns = patterns('cyano.views',
    url(r'^$', 'index'),
    
    url(r'^login$', 'login'),
    url(r'^' + _accession + r'/login$', 'login'),
    url(r'^logout$', 'logout'),
    url(r'^' + _accession + r'/logout', 'logout'),

    url(r'^' + _accession + r'/permission$', 'permission'),
    url(r'^' + _accession + r'/?$', 'species'),
    
    url(r'^' + _accession + r'/(?P<model_type>\w+)/permission$', 'permission'),
    url(r'^' + _accession + r'/(?P<model_type>\w+)/?$', 'listing'),

    url(r'^' + _accession + r'/(?P<model_type>\w+)/(?P<wid>.*?)/permission$', 'permission'),
    url(r'^' + _accession + r'/(?P<model_type>\w+)/(?P<wid>.*)', 'detail')

    #url(r'^' + _accession + r'/genes$', 'genes'),
    #url(r'^' + _accession + r'/genes/(?P<gene>.+)', 'gene'),
    ##url(r'^' + _accession + r'/genes', 'genes'),
#    url(r'^(?P<poll_id>\d+)/$', 'detail'),
#    url(r'^(?P<poll_id>\d+)/results/$', 'results'),
#    url(r'^(?P<poll_id>\d+)/vote/$', 'vote'),
)
