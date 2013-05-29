from django.conf.urls import patterns, url

urlpatterns = patterns('boehringer.views',
    url(r'^$', 'index'),
)
