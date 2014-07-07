from django.conf.urls import patterns, url

urlpatterns = patterns('cyanodesign.views',
    url(r'^$', 'index'),
    url(r'^design/$', 'design'),
    url(r'^get_reactions/$', 'get_reactions'),
    url(r'^simulate/$', 'simulate'),
    url(r'^export/$', 'export')
)
