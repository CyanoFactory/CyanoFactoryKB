from django.conf.urls import url
from . import views

urlpatterns = [
    #url(r'^$', proteingraph, name='graph),
    url(r'^calcTime', views.calcTimes),
    url(r'^index/$', views.calcTimes),
    url(r'^interaction/(?P<protID>\w+\d+)/$', views.checkInteraction),
    url(r'^interaction/(?P<protID>\w+\d+)_(?P<limit>\d+)/$', views.checkInteraction),
    url(r'^interaction/ajax/prot$', views.onlyProtGraph),
    url(r'^interaction/ajax/chem$', views.onlyChemGraph),
    url(r'^interaction/request/$', views.checkRequest),
    url(r'^whole', views.getwholenetwork),
    url(r'^interaction/pathway_start/$', views.startpathway),
    url(r'^interaction/ajax/$', views.findChems),
    url(r'^interaction/pathway/$', views.getpath),
    url(r'^interaction/test/$', views.testSTITCH)
]
