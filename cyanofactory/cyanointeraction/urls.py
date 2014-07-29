from django.conf.urls import patterns, url
from cyanointeraction import views

urlpatterns = patterns('',
    #url(r'^$', proteingraph, name='graph'),
    url(r'^calcTime', 'cyanointeraction.views.calcTimes'),
    url(r'^index/$', 'cyanointeraction.views.calcTimes'),
    url(r'^interaction/(?P<protID>\w+\d+)/$', 'cyanointeraction.views.checkInteraction'),
    url(r'^interaction/(?P<protID>\w+\d+)_(?P<limit>\d+)/$', 'cyanointeraction.views.checkInteraction'),
    url(r'^interaction/ajax/prot$', 'cyanointeraction.views.onlyProtGraph'),
    url(r'^interaction/ajax/chem$', 'cyanointeraction.views.onlyChemGraph'),
    url(r'^interaction/request/$', 'cyanointeraction.views.checkRequest'),
    url(r'^whole', 'cyanointeraction.views.getwholenetwork'),
    url(r'^interaction/pathway_start/$', 'cyanointeraction.views.startpathway'),
    url(r'^interaction/ajax/$', 'cyanointeraction.views.findChems'),
    url(r'^interaction/pathway/$', 'cyanointeraction.views.getpath')
)