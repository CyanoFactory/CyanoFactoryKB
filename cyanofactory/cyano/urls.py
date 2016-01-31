"""
Author: Jonathan Karr, jkarr@stanford.edu
Affiliation: Covert Lab, Department of Bioengineering, Stanford University
Last updated: 2012-07-17

Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

from django.conf.urls import patterns, url
from rest_framework.urlpatterns import format_suffix_patterns
from cyano import views

_species_wid = r'(?P<species_wid>[a-zA-Z0-9_\-]+)'
_wid = r'(?P<wid>[a-zA-Z0-9_\-]+)'

_species_wid_model_type = _species_wid + r'/(?P<model_type>[a-zA-Z0-9_\-]+)'
_species_wid_model_type_wid = _species_wid + r'/(?P<model_type>[a-zA-Z0-9_\-]+)/' + _wid

def cyano_pattern(uri, view, root=False, species=False, model_type=False, item=False, not_found=False, name=None):
    if uri[-1] != "/":
        raise ValueError("URL must end with /")
    elif len(uri) == 1:
        uri = ""

    if name is None:
        name = "cyano.views." + view.__name__

    urls = []
    if item:
        last_url = "{}/{}".format(_species_wid_model_type_wid, uri)
        urls.append(url("^{}$".format(last_url), view, name=name))
        if not_found:
             urls.append(url("^{}.*".format(last_url), views.not_found))
    if model_type:
        last_url = "{}/{}".format(_species_wid_model_type, uri)
        urls.append(url("^{}$".format(last_url), view, name=name))
        if not_found:
             urls.append(url("^{}.*".format(last_url), views.not_found))
    if species:
        last_url = "{}/{}".format(_species_wid, uri)
        urls.append(url("^{}$".format(last_url), view, name=name))
        if not_found:
             urls.append(url("^{}.*".format(last_url), views.not_found))
    if root:
        last_url = uri
        urls.append(url("^{}$".format(last_url), view, name=name))
        if not_found:
             urls.append(url("^{}.*".format(last_url), views.not_found))

    if len(urls) == 0:
        raise ValueError("url list empty")

    return urls

urlpatterns = patterns('',
    url(r'^api/basket/$', views.BasketList.as_view()),
    url(r'^api/basket/(?P<basket_id>[0-9]+)/$', views.BasketDetail.as_view()),
    url(r'^api/$', views.ApiIndex.as_view(), name="cyano-api-index"),
    url(r'^api/' + _species_wid + '/$', views.ApiSpecies.as_view(), name="cyano-api-species"),
    url(r'^api/' + _species_wid_model_type + '/$', views.ApiEntryList.as_view(), name="cyano-api-list"),
    url(r'^api/' + _species_wid_model_type_wid + '/$', views.ApiEntryDetail.as_view(), name="cyano-api-detail")
)

urlpatterns = format_suffix_patterns(urlpatterns, allowed=['json', 'xml'])

urlpatterns += patterns('cyano.views',
    url(r'^sitemap\.xml$', 'sitemap'),
    url(r'^sitemap_toplevel\.xml$', 'sitemap_toplevel'),
    url(r'^sitemap/' + _species_wid + '\.xml$', 'sitemap_species'),
)

urlpatterns += cyano_pattern("login/", views.login, name="login", root=True, species=True, not_found=True)
urlpatterns += cyano_pattern("logout/", views.logout, name="logout", root=True, species=True, not_found=True)
urlpatterns += cyano_pattern("sbgn/", views.sbgn, root=True, not_found=True)
urlpatterns += cyano_pattern("account/password_change_required/", views.password_change_required, name="password_change_required", root=True, not_found=True)
urlpatterns += cyano_pattern("license/", views.licensing, name="license", root=True, species=True, not_found=True)
urlpatterns += cyano_pattern("basket/", views.basket, name="cyano-basket", root=True, species=True, not_found=True)
urlpatterns += cyano_pattern(r'basket/(?P<basket_id>[0-9]+)/', views.basket, name="cyano-basket", root=True, species=True, not_found=True)
urlpatterns += cyano_pattern('basket/create/', views.basket_create, name="cyano-basket-create-save-form", root=True, not_found=True)
urlpatterns += cyano_pattern(r'basket/rename/(?P<basket_id>[0-9]+)/', views.basket_rename, name="cyano-basket-rename-save-form", root=True, not_found=True)
urlpatterns += cyano_pattern('basket/op/', views.basket_op, name="cyano-basket-op", species=True, not_found=True)
#urlpatterns += cyano_pattern("about/", views.about, root=True, species=True, not_found=True)
#urlpatterns += cyano_pattern("tutorial/", views.tutorial, root=True, species=True, not_found=True)
urlpatterns += cyano_pattern(r"user/(?P<username>[\w\d]+)/", views.user, root=True, species=True, not_found=True)
urlpatterns += cyano_pattern("user/", views.users, root=True, species=True, not_found=True)
urlpatterns += cyano_pattern("tutorial/", views.tutorial, root=True, species=True, not_found=True)
urlpatterns += cyano_pattern("jobs/", views.jobs, root=True, species=True, not_found=True)
urlpatterns += cyano_pattern("search/", views.search, root=True, species=True, not_found=True)
urlpatterns += cyano_pattern("add/", views.add, model_type=True, not_found=True)
urlpatterns += cyano_pattern("delete/", views.delete, species=True, item=True, not_found=True)
urlpatterns += cyano_pattern("edit/", views.edit, species=True, item=True, not_found=True)
urlpatterns += cyano_pattern("export/", views.exportData, species=True, not_found=True)
urlpatterns += cyano_pattern("import/data/", views.importData, species=True, not_found=True)
urlpatterns += cyano_pattern("import/species/", views.importSpeciesData, root=True, not_found=True)
urlpatterns += cyano_pattern("permission/", views.global_permission, root=True, not_found=True)
urlpatterns += cyano_pattern("history/(?P<detail_id>[0-9]+)/", views.history_detail, item=True, not_found=True)
urlpatterns += cyano_pattern("history/", views.history, species=True, model_type=True, item=True, not_found=True)
urlpatterns += cyano_pattern("permission/edit/", views.permission_edit, species=True, item=True, not_found=True)
urlpatterns += cyano_pattern("permission/", views.permission, species=True, item=True, not_found=True)
urlpatterns += cyano_pattern("field/", views.detail_field, item=True)
urlpatterns += cyano_pattern("/", views.detail, item=True)
urlpatterns += cyano_pattern("/", views.listing, model_type=True)
urlpatterns += cyano_pattern("/", views.species, species=True)

urlpatterns += patterns('cyano.views',
    url(r'^$', 'index'),
)
