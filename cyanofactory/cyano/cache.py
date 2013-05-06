import django.core.cache
import time

class Cache(object):
    @staticmethod
    def get(name):
        return django.core.cache.cache.get(name)
 
    @staticmethod
    def set(name, value):
        django.core.cache.cache.set(name, value)
    
    @staticmethod
    def try_get(name, function):
        """
        Tries to get the item "name" out of the cache. If this fails
        function is executed and the result added to the cache and returned 
        """
        start = time.time()
        item = Cache.get(name)
        if item == None:
            #print "Cache miss for " + name
            item = function()
            Cache.set(name, item)
        #else:
            #print "Cache hit for " + name
        end = time.time()
        elapsed = end - start
        #print "Operation took: ", elapsed, "seconds."
        return item
