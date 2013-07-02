from django.db import models

class Compound(models.Model):
    compound = models.CharField(primary_key=True, max_length=7)
    name = models.TextField(blank = True)
    formula = models.TextField(blank = True)
    mass = models.FloatField(null = True)

class Conf(models.Model):
    # pk key id used
    map_id = models.CharField(max_length=20, blank = True)
    object = models.CharField(max_length=20)
    type = models.IntegerField(null = True)
    url = models.CharField(max_length=255)
    coord = models.CharField(max_length=100, blank = True)

class Info(models.Model):
    hit_id = models.CharField(primary_key=True, max_length=100)
    length = models.IntegerField()
    description = models.TextField()

class Map(models.Model):
    map_id = models.CharField(primary_key=True, max_length=20,blank = True)
    title = models.CharField(max_length=255)
