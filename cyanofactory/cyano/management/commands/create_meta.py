from django.core.management.base import BaseCommand
from django.db.models import get_app, get_models
from cyano.models import TableMeta, TableMetaColumn, TableMetaManyToMany
from django.db import connection

class Command(BaseCommand):
    def handle(self, *args, **options):
        app = get_app('cyano')
        for model in get_models(app):
            obj = model()
            
            table_name = obj._meta.db_table
            
            table, _ = TableMeta.objects.get_or_create(table_name = table_name, model_name = obj._meta.object_name)

            for field in obj._meta.local_fields:                
                column_name = field.column
                cursor = connection.cursor()
                cursor.execute("SELECT ordinal_position FROM information_schema.columns WHERE table_name = %s AND column_name = %s", [table_name, column_name])
                column_id = cursor.fetchone()[0]
                TableMetaColumn.objects.get_or_create(table = table, column_name = column_name, column_id = column_id)
            
            for m2m in obj._meta.local_many_to_many:
                TableMetaManyToMany.objects.get_or_create(table = table, m2m_name = m2m.name)
