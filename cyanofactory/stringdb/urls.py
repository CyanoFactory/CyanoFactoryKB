from django.conf.urls import patterns, url
from stringdb import views

urlpatterns = patterns('',
    #url(r'^$', proteingraph, name='graph'),
    url(r'^calcTime', 'stringdb.views.calcTimes'),
    url(r'^index/$', 'stringdb.views.calcTimes'),
    url(r'^interaction/(?P<protID>\w+\d+)/$', 'stringdb.views.checkInteraction'),
    url(r'^interaction/(?P<protID>\w+\d+)_(?P<limit>\d+)/$', 'stringdb.views.checkInteraction'),
    url(r'^interaction/ajax/prot$', 'stringdb.views.onlyProtGraph'),
    url(r'^interaction/ajax/chem$', 'stringdb.views.onlyChemGraph'),
    url(r'^interaction/request/$', 'stringdb.views.checkRequest'),
    url(r'^whole', 'stringdb.views.getwholenetwork'),
    url(r'^interaction/pathway_start/$', 'stringdb.views.startpathway'),
    url(r'^interaction/ajax/$', 'stringdb.views.findChems'),
    url(r'^interaction/pathway/$', 'stringdb.views.getpath')
)