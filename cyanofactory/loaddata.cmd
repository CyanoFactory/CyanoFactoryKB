python manage.py syncdb
python manage.py loaddata cyano/fixtures/metadata.json --traceback
python manage.py restoredb < d:\db.dump
python manage.py autocreateuser gabriel
python manage.py importbiodata ../sequences/NC_000911.gb
