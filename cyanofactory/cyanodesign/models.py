"""
Copyright (c) 2014 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

import base64

from django.db import models
from cyano.models import UserProfile


class GlobalPermission(models.Model):
    class Meta:
        permissions = (
            ("access_cyanodesign", "Can access CyanoDesign"),
        )


class DesignModel(models.Model):
    user = models.ForeignKey(UserProfile, verbose_name="Saved by", related_name='+', editable=False)
    name = models.CharField(max_length=255, null=False, blank=False, verbose_name="Name")
    comments = models.TextField(blank=True, default='', verbose_name='Comments')
    filename = models.CharField(max_length=255, null=False, blank=False, verbose_name="Filename")
    _content = models.TextField(verbose_name="BioOpt file content", null=False, blank=False, db_column="content")

    def set_content(self, content):
        self._content = base64.encodestring(content)

    def get_content(self):
        return base64.decodestring(self._content)

    content = property(get_content, set_content)
