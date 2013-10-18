"""
Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

from django.db.models import F, Model, OneToOneField, CharField, IntegerField, URLField, PositiveIntegerField, FloatField, ForeignKey, BooleanField, SlugField, ManyToManyField, TextField, DateTimeField, options, permalink, SET_NULL, Min, Max

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
