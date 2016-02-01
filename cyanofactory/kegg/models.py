"""
Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

from django.db import models
from django.db.models.fields.related import ManyToManyField
from cyano.models import UserProfile


class GlobalPermission(models.Model):
    class Meta:
        permissions = (
            ("access_kegg", "Can access Kegg maps"),
        )


class EcNumber(models.Model):
    name = models.CharField(max_length=255, blank=False, null=False)

    def __str__(self):
        return self.name


class Map(models.Model):
    name = models.CharField(max_length=255, blank=False, null=False)
    title = models.CharField(max_length=255, blank=True, null=False)
    overview = models.BooleanField(default=False, null=False, verbose_name='Is an overview pathway map')
    ec_numbers = ManyToManyField(EcNumber, verbose_name='Ec numbers belonging to this map')

    def __str__(self):
        return self.name


class Query(models.Model):
    user = models.ForeignKey(UserProfile, verbose_name="Saved by", related_name='+', editable=False)
    name = models.CharField(max_length=255, null=False, blank=False, verbose_name="Name of the query")
    query = models.TextField(verbose_name="Query text", null=False, blank=False)
