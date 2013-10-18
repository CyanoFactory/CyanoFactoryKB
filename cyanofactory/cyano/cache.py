"""
Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

from django.core.cache import cache as djcache


class Cache(object):
    @staticmethod
    def get(name):
        return djcache.get(name)
 
    @staticmethod
    def set(name, value, timeout=-1):
        # No idea how to specify the default timeout
        if timeout == -1:
            return djcache.set(name, value)
        return djcache.set(name, value, timeout)
    
    @staticmethod
    def try_get(name, function, timeout=-1):
        """
        Tries to get the item "name" out of the cache. If this fails
        function is executed and the result added to the cache and returned 
        """
        item = Cache.get(name)
        if item is None:
            item = function()
            Cache.set(name, item, timeout)

        return item
