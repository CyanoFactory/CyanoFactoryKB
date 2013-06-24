python manage.py syncdb
python manage.py loaddata cyano/fixtures/metadata.json --traceback
python manage.py restoredb < boehringer/fixtures/data.dump
:: python manage.py restoredb < d:\db.dump
python manage.py autocreateinitial
python manage.py importbiodata ../sequences/NC_000911.gb --traceback
python manage.py importbiocyc ../sample_data/iSyn811_v2-2_sbml_fixed.xml --traceback
python manage.py importbiodata ../sequences/NC_000913.gb --traceback
