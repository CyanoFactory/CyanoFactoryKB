"""
Copyright (c) 2014 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

from datetime import datetime

from django.db import models
from django_extensions.db.fields.json import JSONField
from cyano.models import UserProfile


class GlobalPermission(models.Model):
    class Meta:
        permissions = (
            ("access_cyanodesign", "Can access CyanoDesign"),
        )


class DesignModel(models.Model):
    user = models.ForeignKey(UserProfile, verbose_name="Saved by", related_name='+', editable=False)
    name = models.CharField(max_length=255, null=False, blank=False, verbose_name="Name")
    filename = models.CharField(max_length=255, null=False, blank=False, verbose_name="Filename")
    content = models.TextField(verbose_name="BioOpt file content", null=False, blank=False)

    def get_latest_revision(self):
        return self.revisions.order_by("date").last()

class Revision(models.Model):
    model = models.ForeignKey(DesignModel, related_name='revisions', verbose_name='Model')
    content = models.TextField(verbose_name="BioOpt file content", null=False, blank=False)
    date = models.DateTimeField(default=datetime.now, verbose_name = "Modification date")
    changes = JSONField(verbose_name='Summary of changes')
    reason = models.TextField(blank=True, default='', verbose_name='Description of changes')

    class Meta:
        ordering = ["-date"]
