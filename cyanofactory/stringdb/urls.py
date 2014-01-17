from django.conf.urls import patterns, url
from stringdb.views import *

urlpatterns = patterns('',
    url(r'^$', graph, name='index'),
    url(r'^/graph$', graph, name='test')
)