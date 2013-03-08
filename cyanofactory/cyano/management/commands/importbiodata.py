from django.core.management.base import BaseCommand
from Bio import SeqIO
from cyano.models import Entry
from cyano.models import Species
from cyano.models import RevisionDetail
from cyano.models import UserProfile

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
                species = Species()
                species.wid = record.name
                species.name = record.annotations["source"]
                species.detail = revdetail
                species.save()
