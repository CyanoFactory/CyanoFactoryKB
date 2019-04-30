"""
Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

from django.core.management.base import BaseCommand
from django.apps import apps
from cyano.models import TableMeta


class Command(BaseCommand):
    def handle(self, *args, **options):
        app = apps.get_app_config('cyano')
        for model in app.get_models():
            obj = model()

            TableMeta.objects.get_or_create(table_name=obj._meta.db_table,
                                            model_name=obj._meta.object_name)
