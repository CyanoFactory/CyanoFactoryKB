from django.conf.urls import url
from . import views

urlpatterns = [
    url(r'^$', views.index),
    url(r'^(?P<source>[a-zA-Z0-9_\-]+)\s*:\s*(?P<xid>[a-zA-Z0-9_.\-:.+]+)$', views.dbxref),
]
