from django.conf.urls import url
from . import views

app_name="cyanointeraction"

urlpatterns = [
    #url(r'^$', proteingraph, name='graph),
    url(r'^calcTime', views.calcTimes, name="calcTimes"),
    url(r'^index/$', views.calcTimes, name="calcTimes"),
    url(r'^interaction/(?P<protID>\w+\d+)/$', views.checkInteraction, name="checkInteraction"),
    url(r'^interaction/(?P<protID>\w+\d+)_(?P<limit>\d+)/$', views.checkInteraction, name="checkInteraction"),
    url(r'^interaction/ajax/prot$', views.onlyProtGraph, name="onlyProtGraph"),
    url(r'^interaction/ajax/chem$', views.onlyChemGraph, name="onlyChemGraph"),
    url(r'^interaction/request/$', views.checkRequest, name="checkRequest"),
    url(r'^whole', views.getwholenetwork, name="getwholenetwork"),
    url(r'^interaction/pathway_start/$', views.startpathway, name="startpathway"),
    url(r'^interaction/ajax/$', views.findChems, name="findChems"),
    url(r'^interaction/pathway/$', views.getpath, name="getpath"),
    url(r'^interaction/test/$', views.testSTITCH, name="testSTITCH")
]
