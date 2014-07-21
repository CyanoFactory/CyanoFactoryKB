"""
Copyright (c) 2014 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

from django.conf.urls import patterns, url

pk_str = r"(?P<pk>[0-9]+)/$"

urlpatterns = patterns('cyanodesign.views',
    url(r'^$', 'index', name="cyano-design-index"),
    url(r'^upload/$', 'upload', name="cyano-design-upload"),
    url(r'^edit/' + pk_str, 'design', name="cyano-design-design"),
    url(r'^get_reactions/' + pk_str, 'get_reactions', name="cyano-design-get-reactions"),
    url(r'^simulate/' + pk_str, 'simulate', name="cyano-design-simulate"),
    url(r'^export/$', 'export'),
    url(r'^delete/$', 'delete', name="cyano-design-delete")
)
