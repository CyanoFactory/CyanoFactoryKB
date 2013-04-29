from django.conf.urls import patterns, url

_species_wid = r'(?P<species_wid>[a-zA-Z0-9_\-]+)'
_wid = r'(?P<wid>[a-zA-Z0-9_\-]+)'

urlpatterns = patterns('cyano.views',         
    #url(r'^user/?$', 'user'),
    #url(r'^')
                                    
    url(r'^login$', 'login'),
    url(_species_wid + r'/login$', 'login'),

    url(r'^logout$', 'logout'),
    url(_species_wid + r'/logout', 'logout'),
	
	url(r'^about/*$', 'about'),
	url('^about/' + _species_wid + '/*$', 'about'),
	
	url(r'^tutorial/*$', 'tutorial'),
	url(r'^tutorial/' + _species_wid + '/*$', 'tutorial'),
	
	url(r'^users/*$', 'users'),
	url(r'^users/' + _species_wid + '/*$', 'users'),
	
	url(r'^user/(?P<username>[\w\d]+)/*$', 'user'),
	url(r'^user/(?P<username>[\w\d]+)/' + _species_wid + '/*$', 'user'),
	
	url(r'^search/*', 'search'),
	url(r'^search/' + _species_wid + '/*', 'search'),
	url(r'^edit/' + _species_wid + '/' + _wid + '/*$', 'edit'),
	url(r'^delete/' + _species_wid + '/' + _wid + '/*$', 'delete'),
	url(r'^add/' + _species_wid + '/(?P<model_type>[a-zA-Z0-9_]+)/*$', 'add'),
	url(r'^add/(?P<model_type>[a-zA-Z0-9_]+)/*$', 'add'),
	
	url(r'^export/*$', 'exportData'),
	url(r'^export/' + _species_wid + '/*$', 'exportData'),
	
	url(r'^import/*$', 'importData'),
	url(_species_wid + '/import(/.*)?$', 'importData'),
	
    url(r'^import_submit/*$', 'importSubmitData'),
    url(_species_wid + '/import_submit(/.*)?$', 'importSubmitData'),
    
	url(r'^validate/' + _species_wid + '/*$', 'validate'),
	
	url(r'^sitemap.xml$', 'sitemap'),
	url(r'^sitemap_toplevel.xml$', 'sitemap_toplevel'),
	url(r'^sitemap_' + _species_wid + '.xml$', 'sitemap_species'),

	url(r'^' + _species_wid + r'/?$', 'species'),
	
    url(r'^' + _species_wid + r'/history/?', 'history'),
    url(r'^' + _species_wid + r'/(?P<model_type>\w+)/' + r'history/?', 'history'),
    url(r'^' + _species_wid + r'/(?P<model_type>\w+)/' + _wid + '/history/?', 'history'),
    
	url(r'^' + _species_wid + r'/permission$', 'permission'),

    url(r'^' + _species_wid + r'/(?P<model_type>\w+)/?$', 'list'),
    url(r'^' + _species_wid + r'/(?P<model_type>\w+)/' + _wid + '/?$', 'detail'),

	url(r'^$', 'index'),
)
