"""
Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

from django.conf import settings

#import sys
import re

from django.shortcuts import redirect
from django.http.response import HttpResponseRedirect
from django.core.urlresolvers import reverse

from .helpers import render_queryset_to_response_error


class RedirectAllowedHostMiddlware(object):
    """
    Tests if the provided HOST header is the first one in ALLOWED_HOSTS.
    If not tests if HOST is any of the others. If yes redirect to REDIRECT_URL.
    """
    def process_request(self, request):
        redirect_url = self.get_redirect_url_options_value()
        if redirect_url is None:
            return None

        # Use X_FORWARDED or HOST depending on settings
        use_x = self.get_forwarded_host_options_value()
        host = None
        if use_x:
            host = request.META.get('HTTP_X_FORWARDED_FOR', None)
        if host is None:
            host = request.META.get('HTTP_HOST', None)
        if host is None:
            return None
        allowed_hosts = self.get_allowed_hosts_options_value()

        #print >>sys.stderr, "X_FOR: " + str(request.META.get('HTTP_X_FORWARDED_FOR', None))
        #print >>sys.stderr, "HOST : " + str(request.META.get('HTTP_HOST', None))

        if len(allowed_hosts) == 0:
            return None

        if host != allowed_hosts[0]:
            for ahost in allowed_hosts[1:]:
                if ahost == host:
                    #print >>sys.stderr, host
                    #print >>sys.stderr, "Redi: " + redirect_url + request.path
                    return redirect(redirect_url + request.path)

        # Host suspicious, this is handled by django
        return None

    def get_allowed_hosts_options_value(self):
        return getattr(settings, 'ALLOWED_HOSTS', [])

    def get_forwarded_host_options_value(self):
        return getattr(settings, 'USE_X_FORWARDED_HOST', False)

    def get_redirect_url_options_value(self):
        return getattr(settings, 'REDIRECT_URL', None)

# via http://stackoverflow.com/questions/2093593/


class PasswordChangeMiddleware(object):
    def process_request(self, request):
        if request.user.is_authenticated() and \
            not re.match(r'^' + reverse("cyano:password_change_required"), request.path) and \
            not re.match(r'^' + reverse("cyano:logout"), request.path):
            #re.match(r'^/admin/?', request.path) and

            profile = request.user.profile
            if profile.force_password_change:
                return HttpResponseRedirect(reverse("cyano:password_change_required"))


class ProtectedUserMiddleware(object):
    def process_request(self, request):
        if request.user.username == "cyanofactorist" and \
                ("/user/" in request.path_info or "/account/" in request.path_info):
            return render_queryset_to_response_error(request, msg="Not allowed for this user")
