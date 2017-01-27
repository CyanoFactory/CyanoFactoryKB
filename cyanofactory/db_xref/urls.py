from django.conf.urls import url
from . import views

app_name = "db_xref"

urlpatterns = [
    url(r'^$', views.index, name="index"),
    url(r'^(?P<source>[a-zA-Z0-9_\-]+)\s*:\s*(?P<xid>[a-zA-Z0-9_.\-:.+]+)$', views.dbxref, name="index"),
]
