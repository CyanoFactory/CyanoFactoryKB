from django.core.management.base import BaseCommand
from Bio import SeqIO
from cyano.models import Entry
from cyano.models import Species
from cyano.models import RevisionDetail
from cyano.models import UserProfile
from django.core.exceptions import ObjectDoesNotExist

class Command(BaseCommand):
    def handle(self, *args, **options):
        for arg in args:
            print arg
            handle = open(arg, "r")
            f = SeqIO.parse(handle, "genbank")
            for record in f:
                revdetail = RevisionDetail()
                revdetail.user = UserProfile.objects.get(user__username__exact = "gabriel")
                revdetail.reason = "Initial Load"
                try:
                    species = Species.objects.get(wid = record.name)
                except ObjectDoesNotExist:
                    species = Species()

                species.wid = record.name
                species.name = record.annotations["source"]
                species.comments = None
                species.genetic_code = '1'
                species.save(revision_detail = revdetail)
