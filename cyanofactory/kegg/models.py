from django.db import models
from django.db.models.fields.related import ManyToManyField

class EcNumber(models.Model):
    name = models.CharField(max_length=255)

    def __unicode__(self):
        return self.name

class Map(models.Model):
    name = models.CharField(max_length=255)
    title = models.CharField(max_length=255)
    ec_numbers = ManyToManyField(EcNumber, verbose_name = 'Ec numbers belonging to this map')

    def __unicode__(self):
        return self.name
