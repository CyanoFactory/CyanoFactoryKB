"""
Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""
from django.contrib.auth.models import User

import cyano.models as cmodels
from django.contrib.auth.models import Permission

def process(request):
    if request.is_ajax():
        return {}

    data = {
        'request': request,
        'email': "wuenschi@hs-mittweida.de"
    }

    if not request.user.is_authenticated():
        # Special handling for guests
        user = User.objects.get(username="guest")
    else:
        user = request.user

    data['species_list'] = cmodels.Species.objects.for_permission(cmodels.PermissionEnum.READ_NORMAL, user.profile)

    class GlobalPermHelper(list):
        def __init__(self, *args, **kwargs):
            super(GlobalPermHelper, self).__init__(*args, **kwargs)

        def any_of_cyanomaps(self):
            return any(x in ["access_boehringer", "access_kegg", "access_sbgn"] for x in self)

        def any_of_tools(self):
            return self.any_of_cyanomaps() or "access_cyanodesign" in self


    data['global_permissions'] = GlobalPermHelper(user.profile.get_global_perms().values_list("codename", flat=True))

    return data
