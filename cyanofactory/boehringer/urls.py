"""
Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""
from django.conf.urls import patterns, url

urlpatterns = patterns('boehringer.views',
    url(r'^$', 'index'),
    url(r'^ajax/$', 'index_ajax'),
    url(r'^legacy/$', 'legacy'),
)
