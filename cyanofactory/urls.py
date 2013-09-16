'''
URL patterns
'''

from django.conf.urls import patterns, include, url

# enable the admin:
#from django.contrib import admin
#admin.autodiscover()

# database crossreference resolver
urlpatterns = patterns('db_xref.views',
    url(r'^dbxref/', include("db_xref.urls")),
)

# admin interface
#urlpatterns = patterns('',
#	url(r'^admin/doc/', include('django.contrib.admindocs.urls')),
#	url(r'^admin/', include(admin.site.urls)),
#)

# boehringer map
urlpatterns += patterns('boehringer.views',
    url(r'^boehringer/', include("boehringer.urls")),
)

# authentication
urlpatterns += patterns('django.contrib.auth.views',
    url(r'^login/password_reset_done/$', 'password_reset_done'),
    url(r'^login/password_reset_confirm/(?P<uidb36>[0-9A-Za-z]+)-(?P<token>.+)/$', 'password_reset_confirm', name='password_reset_confirm'),
    url(r'^login/password_change/$', 'password_change',  name = "password_change"),
    url(r'^login/password_change_done/$', 'password_change_done', name = "password_change_done"),
    url(r'^login/password_reset/$', 'password_reset', name = "password_reset"),
    url(r'^login/password_reset_confirm/$', 'password_reset_confirm', name = "password_reset_confirm"),
    url(r'^login/password_reset_complete/$', 'password_reset_complete', name = "password_reset_complete"),
)

# cyanofactory project
urlpatterns += patterns('cyano.views',
    url(r'^', include("cyano.urls")),
)

# authentication
#urlpatterns += patterns('public.views',
#	url(r'^login/*$', 'login'),
#	url(r'^login/(?P<species_wid>[a-zA-Z0-9_\-]+)/*$', 'login'),
#	
#	url(r'^logout/*$', 'logout'),
#	url(r'^logout/(?P<species_wid>[a-zA-Z0-9_\-]+)/*$', 'logout'),
#)

# public interface
# urlpatterns += patterns('public.views',	
#	url(r'^about/*$', 'about'),
#	url(r'^about/(?P<species_wid>[a-zA-Z0-9_\-]+)/*$', 'about'),
#	
#	url(r'^tutorial/*$', 'tutorial'),
#	url(r'^tutorial/(?P<species_wid>[a-zA-Z0-9_\-]+)/*$', 'tutorial'),
#	
#	url(r'^contributors/*$', 'contributors'),
#	url(r'^contributors/(?P<species_wid>[a-zA-Z0-9_\-]+)/*$', 'contributors'),
#	
#	url(r'^contributor/(?P<username>[\w\d]+)/*$', 'contributor'),
#	url(r'^contributor/(?P<username>[\w\d]+)/(?P<species_wid>[a-zA-Z0-9_\-]+)/*$', 'contributor'),
#	
#	url(r'^search/*', 'search'),
#	url(r'^search/(?P<species_wid>[a-zA-Z0-9_\-]+)/*', 'search'),
#	url(r'^list/(?P<species_wid>[a-zA-Z0-9_\-]+)/(?P<model_type>\w+)/*$', 'list'),	
#	url(r'^detail/(?P<species_wid>[a-zA-Z0-9_\-]+)/(?P<wid>[a-zA-Z0-9_\-]+)/*$', 'detail'),
#	url(r'^edit/(?P<species_wid>[a-zA-Z0-9_\-]+)/(?P<wid>[a-zA-Z0-9_\-]+)/*$', 'edit'),
#	url(r'^delete/(?P<species_wid>[a-zA-Z0-9_\-]+)/(?P<wid>[a-zA-Z0-9_\-]+)/*$', 'delete'),
#	url(r'^add/(?P<species_wid>[a-zA-Z0-9_\-]+)/(?P<model_type>[a-zA-Z0-9_]+)/*$', 'add'),
#	url(r'^add/(?P<model_type>[a-zA-Z0-9_]+)/*$', 'add'),
#	
#	url(r'^export/*$', 'exportData'),
#	url(r'^export/(?P<species_wid>[a-zA-Z0-9_\-]+)/*$', 'exportData'),
#	
#	url(r'^import/*$', 'importData'),
#	url(r'^import/(?P<species_wid>[a-zA-Z0-9_\-]+)/*$', 'importData'),
#	
#	url(r'^validate/(?P<species_wid>[a-zA-Z0-9_\-]+)/*$', 'validate'),
#	
#	url(r'^sitemap.xml$', 'sitemap'),
#	url(r'^sitemap_toplevel.xml$', 'sitemap_toplevel'),
#	url(r'^sitemap_(?P<species_wid>[a-zA-Z0-9_\-]+).xml$', 'sitemap_species'),
# 
#	url(r'^$', 'index'),
#	url(r'^(?P<species_wid>[a-zA-Z0-9_\-]+)/*$', 'index'),
# )
