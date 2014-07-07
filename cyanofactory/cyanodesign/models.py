"""
Copyright (c) 2014 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

from django.db import models
from cyano.models import UserProfile


class DesignModel(models.Model):
    user = models.ForeignKey(UserProfile, verbose_name="Saved by", related_name='+', editable=False)
    name = models.CharField(max_length=255, null=False, blank=False, verbose_name="Name")
    comments = models.TextField(blank=True, default='', verbose_name='Comments')
    filename = models.CharField(max_length=255, null=False, blank=False, verbose_name="Filename")
    content = models.TextField(verbose_name="BioOpt file content", null=False, blank=False)
