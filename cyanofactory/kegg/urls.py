"""
Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""
from django.conf.urls import patterns, url

urlpatterns = patterns('kegg.views',
    url(r'^$', 'index'),
    url(r'^(?P<map_id>map[0-9]{5})/$', 'map_view')
)
