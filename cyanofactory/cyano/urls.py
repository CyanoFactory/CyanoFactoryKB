"""
Author: Jonathan Karr, jkarr@stanford.edu
Affiliation: Covert Lab, Department of Bioengineering, Stanford University
Last updated: 2012-07-17

Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

from django.conf.urls import url
from rest_framework.urlpatterns import format_suffix_patterns
from . import views

app_name = "cyano"

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

urlpatterns = [
    #url(r'^api/basket/$', views.BasketList.as_view()),
    #url(r'^api/basket/(?P<basket_id>[0-9]+)/$', views.BasketDetail.as_view()),
    url(r'^api/$', views.ApiIndex.as_view(), name="api-index"),
    url(r'^api/' + _species_wid + '/$', views.ApiSpecies.as_view(), name="api-species"),
    url(r'^api/' + _species_wid_model_type + '/$', views.ApiEntryList.as_view(), name="api-list"),
    url(r'^api/' + _species_wid_model_type_wid + '/$', views.ApiEntryDetail.as_view(), name="api-detail")
]

urlpatterns += format_suffix_patterns(urlpatterns, allowed=['json', 'xml'])

urlpatterns += [
    url(r'^sitemap\.xml$', views.sitemap),
    url(r'^sitemap_toplevel\.xml$', views.sitemap_toplevel),
    url(r'^sitemap/' + _species_wid + '\.xml$', views.sitemap_species),
]

urlpatterns += cyano_pattern("login/", views.login, name="login", root=True, species=True)
urlpatterns += cyano_pattern("logout/", views.logout, name="logout", root=True, species=True)
urlpatterns += cyano_pattern("register/", views.register, name="register", root=True)
urlpatterns += cyano_pattern("sbgn/", views.sbgn, name="sbgn", root=True)
urlpatterns += cyano_pattern("account/edit/", views.account_change, name="account_change", root=True, species=True)
urlpatterns += cyano_pattern("account/password_change_required/", views.password_change_required, name="password_change_required", root=True)
urlpatterns += cyano_pattern("license/", views.licensing, name="license", root=True, species=True)
urlpatterns += cyano_pattern("basket/", views.basket, name="basket", root=True, species=True)
urlpatterns += cyano_pattern(r'basket/(?P<basket_id>[0-9]+)/', views.basket, name="basket", root=True, species=True)
urlpatterns += cyano_pattern('basket/create/', views.basket_create, name="basket-create", root=True)
urlpatterns += cyano_pattern(r'basket/rename/(?P<basket_id>[0-9]+)/', views.basket_rename, name="basket-rename", root=True)
urlpatterns += cyano_pattern('basket/op/', views.basket_op, name="basket-op", root=True, species=True)
#urlpatterns += cyano_pattern("about/", views.about, root=True, species=True)
#urlpatterns += cyano_pattern("tutorial/", views.tutorial, root=True, species=True)

# Ab hier
urlpatterns += cyano_pattern(r"user/(?P<username>[\w\d]+)/", views.user, name="user", root=True, species=True)
urlpatterns += cyano_pattern("user/", views.users, name="users", root=True, species=True)
urlpatterns += cyano_pattern("tutorial/", views.tutorial, name="tutorial", root=True, species=True)
urlpatterns += cyano_pattern("jobs/", views.jobs, name="jobs", root=True, species=True)
urlpatterns += cyano_pattern("search/", views.search, name="search", root=True, species=True)
urlpatterns += cyano_pattern("add/", views.add, name="add", model_type=True)
urlpatterns += cyano_pattern("delete/", views.delete, name="delete", species=True, item=True)
urlpatterns += cyano_pattern("edit/", views.edit, name="edit", species=True, item=True)
urlpatterns += cyano_pattern("export/", views.exportData, name="exportData", species=True)
urlpatterns += cyano_pattern("import/data/", views.importData, name="importData", species=True)
urlpatterns += cyano_pattern("import/species/", views.importSpeciesData, name="importSpeciesData", root=True)
urlpatterns += cyano_pattern("permission/", views.global_permission, name="global_permission", root=True)
urlpatterns += cyano_pattern("history/(?P<detail_id>[0-9]+)/", views.history_detail, name="history_detail", item=True)
urlpatterns += cyano_pattern("history/", views.history, name="history", species=True, model_type=True, item=True)
urlpatterns += cyano_pattern("permission/edit/", views.permission_edit, name="permission_edit", species=True, item=True)
urlpatterns += cyano_pattern("permission/", views.permission, name="permission", species=True, item=True)
urlpatterns += cyano_pattern("field/", views.detail_field, name="detail_field", item=True)
urlpatterns += cyano_pattern("/", views.detail, name="detail", item=True)
urlpatterns += cyano_pattern("/", views.listing, name="listing", model_type=True)
urlpatterns += cyano_pattern("/", views.species, name="species", species=True)

urlpatterns += [
    url(r'^$', views.index, name="index"),
]
