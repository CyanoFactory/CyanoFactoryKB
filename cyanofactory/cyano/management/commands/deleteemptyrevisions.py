"""
Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

from django.core.management.base import BaseCommand
from django.db.models.query_utils import Q

from cyano.models import Revision, RevisionDetail

class Command(BaseCommand):
    """Deletes empty revisions
    """
    def handle(self, *args, **options):
        revisions = Revision.objects.filter(Q(new_data="{}") & ~Q(action="D"))
        print("Deleting {} revisions".format(revisions.count()))
        revisions.delete()
        rev_details = RevisionDetail.objects.filter(revisions__isnull=True)
        print("Deleting {} revision details".format(rev_details.count()))
        rev_details.delete()
