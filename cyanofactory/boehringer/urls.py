"""
Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""
from django.conf.urls import url
from . import views

app_name = "boehringer"

urlpatterns = [
    url(r'^$', views.index, name="index"),
    url(r'^ajax/$', views.index_ajax, name="index_ajax"),
    url(r'^legacy/$', views.legacy, name="legacy"),
]
