"""
Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

# Model used to dump Kegg Data from Microbes online

from django.db import models

class KEGGCompound(models.Model):
    compound = models.CharField(primary_key=True, max_length=7, db_column='compound')
    name = models.TextField(db_column = "name", blank = True)
    formula = models.TextField(db_column = "formula", blank = True)
    mass = models.FloatField(db_column = "mass", null = True)
    
    class Meta:
        db_table = 'KEGGCompound'

class KEGGConf(models.Model):
    map_id = models.CharField(primary_key=True, max_length=20, db_column='mapId', blank = True)
    object = models.CharField(max_length=20, db_column='object')
    type = models.IntegerField(db_column = "type", null = True)
    url = models.CharField(max_length=255, db_column='url')
    coord = models.CharField(max_length=100, db_column='coord', blank = True)
    
    class Meta:
        db_table = 'KEGGConf'

class KEGGInfo(models.Model):
    hit_id = models.CharField(primary_key=True, max_length=100, db_column='hitId')
    s_length = models.IntegerField(db_column = "sLength")
    description = models.TextField(db_column = "description")
    
    class Meta:
        db_table = 'KEGGInfo'

class KEGGMap(models.Model):
    map_id = models.CharField(primary_key=True, max_length=20, db_column='mapId', blank = True)
    title = models.CharField(max_length=255, db_column='title')
    
    class Meta:
        db_table = 'KEGGMap'
