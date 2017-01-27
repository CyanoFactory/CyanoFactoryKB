"""
Copyright (c) 2014 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

from django.conf.urls import url
from . import views

app_name = "cyanodesign"

pk_str = r"(?P<pk>[0-9]+)/$"

urlpatterns = [
    url(r'^$', views.index, name="index"),
    url(r'^upload/' + pk_str, views.upload, name="upload"),
    url(r'^edit/' + pk_str, views.design, name="design"),
    url(r'^history/' + pk_str, views.history, name="history"),
    url(r'^get_reactions/' + pk_str, views.get_reactions, name="get_reactions"),
    url(r'^simulate/' + pk_str, views.simulate, name="simulate"),
    url(r'^export/' + pk_str, views.export, name="export"),
    url(r'^save/' + pk_str, views.save, name="save"),
    url(r'^save_as/' + pk_str, views.save_as, name="save_as"),
    url(r'^delete/$', views.delete, name="delete")
]
