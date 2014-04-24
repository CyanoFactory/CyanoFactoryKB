@echo off
echo Cyanofactory KB Initial Import script
@echo on

python manage.py syncdb
python manage.py migrate
python manage.py loaddata cyano/fixtures/metadata.json
python manage.py restoredb < boehringer/fixtures/data.dump
python manage.py autocreateinitial
python manage.py create_meta

python manage.py import_kegg_ec

python manage.py create_species --wid=Synechocystis --name="Synechocystis sp. PCC 6803" --reason="Create Synechocystis"
python manage.py create_species --wid=Ecoli --name="Escherichia coli str. K-12 substr. MG1655" --reason="Create Ecoli"

python manage.py importgenbank ../sequences/NC_000911_Syn.gb --wid=Synechocystis --reason="Import Chromosome of Synechocystis" --chromosome=Chromosome-Syn --name="Chromosome" --user=management
python manage.py importgenbank ../sequences/NC_005229_Syn_Plasmid_pSYSM.gb --wid=Synechocystis --reason="Import Plasmid pSYSM of Synechocystis" --chromosome=pSYSM --name="Plasmid pSYSM" --user=management
python manage.py importgenbank ../sequences/NC_005230_Syn_Plasmid_pSYSA.gb --wid=Synechocystis --reason="Import Plasmid pSYSA of Synechocystis" --chromosome=pSYSA --name="Plasmid pSYSA" --user=management
python manage.py importgenbank ../sequences/NC_005231_Syn_Plasmid_pSYSG.gb --wid=Synechocystis --reason="Import Plasmid pSYSG of Synechocystis" --chromosome=pSYSG --name="Plasmid pSYSG" --user=management
python manage.py importgenbank ../sequences/NC_005232_Syn_Plasmid_pSYSX.gb --wid=Synechocystis --reason="Import Plasmid pSYSX of Synechocystis" --chromosome=pSYSX --name="Plasmid pSYSX" --user=management

python manage.py importgenbank ../sequences/NC_000913_Ecoli.gb --wid=Ecoli --reason="Import Chromosome of Ecoli" --chromosome=Chromosome-Ecoli --name="Chromosome" --user=management

python manage.py importsbml ../sample_data/iSyn811_v2-2_sbml_fixed.xml --wid=Synechocystis --reason="Import SBML iSyn from Valencia Dataset" --user=management

python manage.py importproopdb ../sample_data/SynOperonPrediction.txt --wid=Synechocystis --reason="Import Operon Prediction Files" --user=management
