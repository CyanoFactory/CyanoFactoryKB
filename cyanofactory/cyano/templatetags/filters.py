"""
Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

import re
from django import template

register = template.Library()

@register.filter
def index(value, arg):
    """Access an array (or hash) with the given value"""
    try:
        return value[arg]
    except KeyError:
        from django.conf import settings

        value = settings.TEMPLATE_STRING_IF_INVALID

#via http://stackoverflow.com/questions/844746/performing-a-getattr-style-lookup-in-a-django-template
numeric_test = re.compile("^\d+$")

@register.filter(name='getattr')
def getattribute(value, arg):
    """Gets an attribute of an object dynamically from a string name"""

    if hasattr(value, str(arg)):
        return getattr(value, arg)
    elif hasattr(value, 'has_key') and value.has_key(arg):
        return value[arg]
    elif numeric_test.match(str(arg)) and len(value) > int(arg):
        return value[int(arg)]
    else:
        from django.conf import settings
        return settings.TEMPLATE_STRING_IF_INVALID

