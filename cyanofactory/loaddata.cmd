python manage.py syncdb
python manage.py loaddata cyano/fixtures/metadata.json --traceback
python manage.py restoredb < boehringer/fixtures/data.dump
:: python manage.py restoredb < d:\db.dump
python manage.py autocreateinitial
python manage.py importgenbank ../sequences/NC_000911.gb --traceback
python manage.py importsbml ../sample_data/iSyn811_v2-2_sbml_fixed.xml --wid=NC_000911 --traceback
python manage.py importproopdb ../sample_data/SynOperonPrediction.txt --wid=NC_000911 --traceback
python manage.py importgenbank ../sequences/NC_000913.gb --traceback
