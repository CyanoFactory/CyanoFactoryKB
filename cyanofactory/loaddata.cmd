d:\programme\python\python manage.py syncdb
d:\programme\python\python manage.py loaddata cyano/fixtures/metadata.json --traceback
:: d:\programme\python\python manage.py restoredb < d:\db.dump
d:\programme\python\python manage.py autocreateuser gabriel
d:\programme\python\python manage.py importbiodata d:/programmierung/cyanofactory/sequences/NC_000911.gb