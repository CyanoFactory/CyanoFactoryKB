"""
Copyright (c) 2014 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

from django.conf.urls import url
from . import views

pk_str = r"(?P<pk>[0-9]+)/$"

urlpatterns = [
    url(r'^$', views.index, name="cyano-design-index"),
    url(r'^upload/' + pk_str, views.upload, name="cyano-design-save-upload-form"),
    url(r'^edit/' + pk_str, views.design, name="cyano-design-design"),
    url(r'^history/' + pk_str, views.history, name="cyano-design-history"),
    url(r'^get_reactions/' + pk_str, views.get_reactions, name="cyano-design-get-reactions"),
    url(r'^simulate/' + pk_str, views.simulate, name="cyano-design-simulate"),
    url(r'^export/' + pk_str, views.export, name="cyano-design-export"),
    url(r'^save/' + pk_str, views.save, name="cyano-design-save"),
    url(r'^save_as/' + pk_str, views.save_as, name="cyano-design-saveas"),
    url(r'^delete/$', views.delete, name="cyano-design-delete")
]
