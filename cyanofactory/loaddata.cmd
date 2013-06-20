python manage.py syncdb
python manage.py loaddata cyano/fixtures/metadata.json --traceback
python manage.py restoredb < d:\db.dump
python manage.py autocreateinitial
python manage.py importbiodata ../sequences/NC_000911.gb --traceback

