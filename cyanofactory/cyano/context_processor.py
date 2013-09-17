"""
Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

import cyano.models as cmodels

def process(request):
    data = {}

    data['request'] = request
    data['email'] = "wuenschi@hs-mittweida.de"
    data['species_list'] = cmodels.Species.objects.all()
    
    #species = request.resolver_match.kwargs.get("species_wid")

    return data
