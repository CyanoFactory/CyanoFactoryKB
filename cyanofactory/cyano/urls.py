from django.conf.urls import patterns, url

_species_wid = r'(?P<species_wid>[a-zA-Z0-9_\-]+)'
_wid = r'(?P<wid>[a-zA-Z0-9_\-]+)'

_species_wid_model_type = _species_wid + r'/(?P<model_type>[a-zA-Z0-9_\-]+)'
_species_wid_model_type_wid = _species_wid + r'/(?P<model_type>[a-zA-Z0-9_\-]+)/' + _wid

urlpatterns = patterns('cyano.views',         
                     
    url(r'^login/?$', 'login'),
    url(r'^' + _species_wid + r'/login/?$', 'login'),

    url(r'^logout/?$', 'logout'),
    url(r'^' + _species_wid + r'/logout/?', 'logout'),
	
	url(r'^about/?$', 'about'),
	url(r'^' + _species_wid + '/about/?$', 'about'),
	
	url(r'^tutorial/?$', 'tutorial'),
	url(r'^' + _species_wid + '/tutorial/?$', 'tutorial'),

	url(r'^users/?$', 'users'),
	url(r'^' + _species_wid + '/users/?$', 'users'),
	
	url(r'^user/(?P<username>[\w\d]+)/?$', 'user'),
	url(r'^' + _species_wid + '/user/(?P<username>[\w\d]+)/?$', 'user'),
	
	url(r'^search/?$', 'search'),
	url(r'^' + _species_wid + '/search/?$', 'search'),	

    url(r'^' + _species_wid + r'/delete/?$', 'delete'),
	url(r'^' + _species_wid_model_type_wid + '/delete/?$', 'delete'),

    url(r'^(?P<model_type>Species)/?', 'add'), # Special case to create new species
    url(r'^' + _species_wid + r'/add/?$', 'add'),
	url(r'^' + _species_wid_model_type + r'/add/?$', 'add'),

    url(r'^' + _species_wid + r'/edit/?$', 'edit'),
    url(r'^' + _species_wid_model_type_wid + '/edit/?$', 'edit'),
	
	url(r'^search/?$', 'search'),
	url(r'^' + _species_wid + '/search/?$', 'search'),

	url(r'^export/?$', 'exportData'),
	url(r'^' + _species_wid + '/export/?$', 'exportData'),

	url(r'^import/?$', 'importData'),
	url(r'^' + _species_wid + '/import/?$', 'importData'),
	
    url(r'^import_submit/*$', 'importSubmitData'),
    url(_species_wid + '/import_submit/?$', 'importSubmitData'),
    
	url(r'^' + _species_wid + '/validate/?$', 'validate'),
	
	url(r'^sitemap.xml$', 'sitemap'),
	url(r'^sitemap_toplevel.xml$', 'sitemap_toplevel'),
	url(r'^sitemap_' + _species_wid + '.xml$', 'sitemap_species'),

	url(r'^' + _species_wid + r'/?$', 'species'),

    url(r'^' + _species_wid + r'/history/?$', 'history'),
    url(r'^' + _species_wid_model_type + r'/history/?$', 'history'),
    url(r'^' + _species_wid_model_type_wid + r'/history/?$', 'history'),

    #url(r'^' + _species_wid + r'/history/(?P<detail_id>[0-9]+)/?$', 'history_detail'),
    #url(r'^' + _species_wid_model_type + r'/history/(?P<detail_id>[0-9]+)/?$', 'history_detail'),
    url(r'^' + _species_wid_model_type_wid + r'/history/(?P<detail_id>[0-9]+)/?$', 'history_detail'),

    url(r'^' + _species_wid + r'/permission/?$', 'permission'),
    #url(r'^' + _species_wid_model_type + r'/permission/?$', 'permission'),
    url(r'^' + _species_wid_model_type_wid + r'/permission/?$', 'permission'),

    url(r'^' + _species_wid + r'/permission/edit/?$', 'permission_edit'),
    #url(r'^' + _species_wid_model_type + r'/permission/edit/?$', 'permission_edit'),
    url(r'^' + _species_wid_model_type_wid + r'/permission/edit/?$', 'permission_edit'),

    url(r'^' + _species_wid_model_type + r'/?$', 'listing'),
    url(r'^' + _species_wid_model_type_wid + r'/?$', 'detail'),

	url(r'^$', 'index'),
)
