from cyano.models import RevisionDetail, UserProfile, Species
from django.core.exceptions import ObjectDoesNotExist
from cyano.helpers import slugify

class BioParser(object):
    """
    Base class for all biological file parsers
    """

    def __init__(self, wid = None, user = None, reason = ""):
        if wid is None:
            raise ValueError("wid argument is mandatory")
        
        if user is None:
            raise ValueError("user argument is mandatory")
        
        if reason is None or len(reason) == 0:
            raise ValueError("reason is mandatory")
        
        self.wid = self.try_slugify("Wid", wid)
       
        try:
            self.species = Species.objects.get(wid = wid)
        except ObjectDoesNotExist:
            raise ValueError("Species {} not found".format(wid))

        if isinstance(user, UserProfile):
            self.user = user
        else:
            try:
                self.user = UserProfile.objects.get(user__username = user)
            except:
                try:           
                    self.user = UserProfile.objects.get(user__pk = int(user))
                except:
                    raise ValueError("Invalid username " + str(user))

        self.detail = RevisionDetail(user = self.user, reason = reason)
    
    def try_slugify(self, name, not_slug):
        slug = slugify(not_slug)
        
        if slug != not_slug:
            raise ValueError("{} {} contained invalid characters. Only letters, numbers and _ are allowed".format(name, not_slug))
        
        return slug
        
    def parse(self, handle):
        """
        Parses the provided handle
        """
        raise NotImplementedError("Subclasses must implement me")
    
    def apply(self):
        """
        Applies changes of the parsed data (in parse) do the database.
        This operation should be in a transaction.
        """
        raise NotImplementedError("Subclasses must implement me")