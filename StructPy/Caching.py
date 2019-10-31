"""
Some functions are 'expensive' (they require quite a bit of compute time), and if the structure hasn't been modified, then these types of functions can be cached and the cache can be cleared if the cached property becomes invalid
Changing a property to be invalid only requires you to set `self.__cache__{method_name} = None` to clear the cached variable.
"""

def cached_property(expensive_function):
    """
    A decorator function like @property that will cache return a cached quantity if it is available
    """
    @property
    def caching_function(self):
        cacheName = f"__cache__{expensive_function.__name__}"
        
		
        try: # check if the cache has been initialized
            cacheExists = True
            cache = getattr(self, cacheName)
        except AttributeError:
            cacheExists = False
            cache = None
        
		# Check if the cache is valid (not None), caching is requested, and that it exists
        if ( cache is not None ) and ( self.withCaching == True ) and (cacheExists == True):
            return cache
        else:
			#worst case, now we have to compute the quantity
            computed = expensive_function(self)
            setattr(self, cacheName, computed)
            return computed
        
    return caching_function