from django.conf.urls import patterns, url
from stringdb import views

urlpatterns = patterns('',
    #url(r'^$', proteingraph, name='graph'),
    url(r'^graph$', views.proteingraph2json),
    url(r'^index/$', views.index),
    url(r'^interaction/(?P<protID>\d+)/$', views.checkInteraction),
    url(r'^interaction/(?P<protID>\d+)_(?P<limit>\d+)/$', views.checkInteraction),
    url(r'^interaction/ajax/$', views.onlygraph),
    url(r'^whole', views.getwholenetwork)
)