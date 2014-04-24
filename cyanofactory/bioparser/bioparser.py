"""
Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

from cyano.models import RevisionDetail, UserProfile, Species
from cyano.helpers import slugify

class BioParser(object):
    """
    Base class for all biological file parsers
    """

    def __init__(self, wid = None, user = None, reason = ""):
        if not wid:
            raise ValueError("wid argument is mandatory")
        
        if not user:
            raise ValueError("user argument is mandatory")
        
        if not reason or len(reason) == 0:
            raise ValueError("reason is mandatory")
        
        self.wid = self.try_slugify("Wid", wid)
        self.species = Species.objects.for_wid(wid, create=True)

        try:
            self.user = UserProfile.objects.get(user__username = user)
        except:
            try:
                self.user = UserProfile.objects.get(user__pk = int(user))
            except:
                raise ValueError("Invalid username " + str(user))

        self.detail = RevisionDetail(user = self.user, reason = reason)

    def report_progress(self, current, total, message):
        if hasattr(self, "notify_progress"):
            self.notify_progress(current=current, total=total, message=message)

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


def notify_progress_print(current, total, message):
    percent = current * 100 / total
    remaining = 100 - percent
    print("[{}{}] - {}".format("*"*percent, " "*remaining, message))
