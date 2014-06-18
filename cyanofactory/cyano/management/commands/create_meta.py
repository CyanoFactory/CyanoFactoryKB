"""
Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

from django.core.management.base import BaseCommand
from django.db.models import get_app, get_models
from cyano.models import TableMeta


class Command(BaseCommand):
    def handle(self, *args, **options):
        app = get_app('cyano')
        for model in get_models(app):
            obj = model()

            TableMeta.objects.get_or_create(table_name=obj._meta.db_table,
                                            model_name=obj._meta.object_name)
