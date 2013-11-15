"""
Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

from django.db.models import Model, CharField, IntegerField, ForeignKey, TextField
from cyano.models import UserProfile


class Color(Model):
    name = CharField(max_length=255, blank=True, default='', verbose_name='color')


class BioMolecule(Model):
    x = IntegerField()
    y = IntegerField()
    w = IntegerField()
    h = IntegerField()
    color = ForeignKey(Color)
    title = CharField(max_length=255, blank=True, default='', verbose_name='title')


class Enzyme(BioMolecule):
    ec = CharField(max_length=255, blank=True, default='', verbose_name='ec-number')


class Metabolite(BioMolecule):
    pass


class Query(Model):
    user = ForeignKey(UserProfile, verbose_name="Saved by", related_name='+', editable=False)
    name = CharField(max_length=255, null=False, blank=False, verbose_name="Name of the query")
    query = TextField(verbose_name="Query text", null=False, blank=False)
