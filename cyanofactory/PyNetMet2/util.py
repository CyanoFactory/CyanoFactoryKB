import re
import sys

PY3 = sys.version_info[0] == 3

def python_2_unicode_compatible(cls):
    """
    The implementation comes from django.utils.encoding.
    """
    if not PY3:
        cls.__unicode__ = cls.__str__
        cls.__str__ = lambda self: self.__unicode__().encode('utf-8')
    return cls

def create_sid(name):
    name = re.sub('[^0-9a-zA-Z]+', '_', name)
    if len(name) > 0 and name[0].isdigit():
        name = "_" + name
    return name
