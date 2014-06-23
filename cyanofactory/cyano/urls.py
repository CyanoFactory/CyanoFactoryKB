"""
Author: Jonathan Karr, jkarr@stanford.edu
Affiliation: Covert Lab, Department of Bioengineering, Stanford University
Last updated: 2012-07-17

Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

from django.conf.urls import patterns, url

_species_wid = r'(?P<species_wid>[a-zA-Z0-9_\-]+)'
_wid = r'(?P<wid>[a-zA-Z0-9_\-]+)'

_species_wid_model_type = _species_wid + r'/(?P<model_type>[a-zA-Z0-9_\-]+)'
_species_wid_model_type_wid = _species_wid + r'/(?P<model_type>[a-zA-Z0-9_\-]+)/' + _wid

urlpatterns = patterns('cyano.views',
    url(r'^login/$', 'login', name="login"),
    url(r'^' + _species_wid + r'/login/$', 'login', name="login"),

    url(r'^logout/$', 'logout', name="logout"),
    url(r'^' + _species_wid + r'/logout/$', 'logout', name="logout"),

    url(r'^sbgn/$', 'sbgn'),

    url(r'^account/password_change_required/$', 'password_change_required', name="password_change_required"),

    url(r'^license/$', 'licensing'),
    url(r'^' + _species_wid + '/license/$', 'licensing'),
    
    url(r'^' + _species_wid + '/basket/$', 'basket'),
    url(r'^basket/$', 'basket'),
    
    url(r'^' + _species_wid_model_type + r'/basket/op/$', 'basket_op', name="cyano-basket-op"),

    #url(r'^about/$', 'about'),
    #url(r'^' + _species_wid + '/about/$', 'about'),

    #url(r'^tutorial/$', 'tutorial'),
    #url(r'^' + _species_wid + '/tutorial/$', 'tutorial'),

    url(r'^users/$', 'users'),
    url(r'^' + _species_wid + '/users/$', 'users'),

    url(r'^user/(?P<username>[\w\d]+)/$', 'user'),
    url(r'^' + _species_wid + '/user/(?P<username>[\w\d]+)/$', 'user'),

    url(r'^jobs/$', 'jobs'),
    url(r'^' + _species_wid + '/jobs/$', 'jobs'),

    url(r'^search/$', 'search'),
    url(r'^' + _species_wid + '/search/$', 'search'),

    url(r'^' + _species_wid + r'/delete/$', 'delete'),
    url(r'^' + _species_wid_model_type_wid + '/delete/$', 'delete'),

    url(r'^(?P<model_type>Species)/add/', 'add'),  # Special case to create new species
    url(r'^' + _species_wid_model_type + r'/add/$', 'add'),

    url(r'^' + _species_wid + r'/edit/$', 'edit'),
    url(r'^' + _species_wid_model_type_wid + '/edit/$', 'edit'),

    url(r'^' + _species_wid + '/export/$', 'exportData'),

    url(r'^import/data/$', 'importData'),
    url(r'^' + _species_wid + '/import/data/$', 'importData'),

    url(r'^import/species/$', 'importSpeciesData'),
    url(r'^' + _species_wid + '/import/species/$', 'importSpeciesData'),

    #url(r'^' + _species_wid + '/validate/$', 'validate'),

    url(r'^sitemap\.xml$', 'sitemap'),
    url(r'^sitemap_toplevel\.xml$', 'sitemap_toplevel'),
    url(r'^sitemap/' + _species_wid + '\.xml$', 'sitemap_species'),

    url(r'^' + _species_wid + r'/$', 'species'),

    url(r'^' + _species_wid + r'/history/$', 'history'),
    url(r'^' + _species_wid_model_type + r'/history/$', 'history'),
    url(r'^' + _species_wid_model_type_wid + r'/history/$', 'history'),

    #url(r'^' + _species_wid + r'/history/(?P<detail_id>[0-9]+)/$', 'history_detail'),
    #url(r'^' + _species_wid_model_type + r'/history/(?P<detail_id>[0-9]+)/$', 'history_detail'),
    url(r'^' + _species_wid_model_type_wid + r'/history/(?P<detail_id>[0-9]+)/$', 'history_detail'),

    url(r'^' + _species_wid + r'/permission/$', 'permission'),
    #url(r'^' + _species_wid_model_type + r'/permission/$', 'permission'),
    url(r'^' + _species_wid_model_type_wid + r'/permission/$', 'permission'),

    url(r'^' + _species_wid + r'/permission/edit/$', 'permission_edit'),
    #url(r'^' + _species_wid_model_type + r'/permission/edit/$', 'permission_edit'),
    url(r'^' + _species_wid_model_type_wid + r'/permission/edit/$', 'permission_edit'),

    url(r'^' + _species_wid_model_type + r'/$', 'listing'),
    url(r'^' + _species_wid_model_type_wid + r'/$', 'detail'),
    url(r'^' + _species_wid_model_type_wid + r'/field/$', 'detail_field'),

    url(r'^$', 'index'),
)
