DROP SCHEMA public CASCADE;
CREATE SCHEMA "public";

CREATE TABLE "Warehouse"
(
  "Version"  numeric(6,3)  NOT NULL,
  "LoadDate"  timestamp  NOT NULL,
  "MaxSpecialWID"  bigint  NOT NULL,
  "MaxReservedWID"  bigint  NOT NULL,
  "Description"  text
  --
);



INSERT INTO "Warehouse" ("Version", "LoadDate", "MaxSpecialWID", "MaxReservedWID")
        VALUES (4.5, now(), 999, 999999);
    

CREATE TABLE "WIDTable"
(
  "PreviousWID"  bigserial  NOT NULL ,
  --
  CONSTRAINT "PK_WIDTable" PRIMARY KEY ("PreviousWID")
);



INSERT INTO "WIDTable" ("PreviousWID") values (999999);

CREATE TABLE "SpecialWIDTable"
(
  "PreviousWID"  bigserial  NOT NULL ,
  --
  CONSTRAINT "PK_SpecialWIDTable" PRIMARY KEY ("PreviousWID")
);



INSERT INTO "SpecialWIDTable" ("PreviousWID") values (1);

CREATE TABLE "Enumeration"
(
  "TableName"  character varying(50)  NOT NULL,
  "ColumnName"  character varying(50)  NOT NULL,
  "Value"  character varying(50)  NOT NULL,
  "Meaning"  text
  --
);



CREATE TABLE "DataSet"
(
  "WID"  bigint  NOT NULL,
  "Name"  character varying(255)  NOT NULL,
  "Version"  character varying(50),
  "ReleaseDate"  timestamp,
  "LoadDate"  timestamp  NOT NULL,
  "ChangeDate"  timestamp,
  "HomeURL"  character varying(255),
  "QueryURL"  character varying(255),
  "LoadedBy"  character varying(255),
  "Application"  character varying(255),
  "ApplicationVersion"  character varying(255),
  --
  CONSTRAINT "PK_DataSet1" PRIMARY KEY ("WID")
);



CREATE TABLE "DataSetHierarchy"
(
  "SuperWID"  bigint  NOT NULL,
  "SubWID"  bigint  NOT NULL
  --
);



CREATE TABLE "Entry"
(
  "OtherWID"  bigint  NOT NULL,
  "InsertDate"  timestamp  NOT NULL,
  "CreationDate"  timestamp,
  "ModifiedDate"  timestamp,
  "LoadError"  character(1)  NOT NULL,
  "LineNumber"  bigint,
  "ErrorMessage"  text,
  "DataSetWID"  bigint  NOT NULL,
  --
  CONSTRAINT "PK_Entry" PRIMARY KEY ("OtherWID")
);



CREATE TABLE "GeneExpressionData"
(
  "B"  smallint  NOT NULL,
  "D"  smallint  NOT NULL,
  "Q"  smallint  NOT NULL,
  "Value"  character varying(100)  NOT NULL,
  "BioAssayValuesWID"  bigint  NOT NULL
  --
);



CREATE TABLE "Element"
(
  "WID"  bigint  NOT NULL,
  "Name"  character varying(15)  NOT NULL,
  "ElementSymbol"  character varying(2)  NOT NULL,
  "AtomicWeight"  numeric  NOT NULL,
  "AtomicNumber"  smallint  NOT NULL,
  --
  CONSTRAINT "PK_Element" PRIMARY KEY ("WID")
);



CREATE TABLE "Valence"
(
  "OtherWID"  bigint  NOT NULL,
  "Valence"  smallint  NOT NULL
  --
);



CREATE TABLE "LightSource"
(
  "WID"  bigint  NOT NULL,
  "Wavelength"  numeric,
  "Type"  character varying(100),
  "InstrumentWID"  bigint,
  "DataSetWID"  bigint  NOT NULL,
  --
  CONSTRAINT "PK_LightSource1" PRIMARY KEY ("WID")
);

CREATE INDEX "LightSource_DWID" ON "LightSource"("DataSetWID");

INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('LightSource', 'Type', 'arc-lamp', 'A type of light source');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('LightSource', 'Type', 'argon-laser', 'An ion laser that uses argon gas to produce light');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('LightSource', 'Type', 'krypton-laser', 'An ion laser that uses krypton gas to produce light');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('LightSource', 'Type', 'argon-krypton-laser', 'An ion laser that uses a combination of krypton and argon gases to produce light');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('LightSource', 'Type', 'helium-neon-laser', 'A laser that uses helium and neon gases to produce light');

CREATE TABLE "FlowCytometrySample"
(
  "WID"  bigint  NOT NULL,
  "BioSourceWID"  bigint,
  "FlowCytometryProbeWID"  bigint,
  "MeasurementWID"  bigint,
  "ManufacturerWID"  bigint,
  "DataSetWID"  bigint  NOT NULL,
  --
  CONSTRAINT "PK_FlowCytometrySample1" PRIMARY KEY ("WID")
);

CREATE INDEX "FlowCytometrySample_DWID" ON "FlowCytometrySample"("DataSetWID");


CREATE TABLE "FlowCytometryProbe"
(
  "WID"  bigint  NOT NULL,
  "Type"  character varying(100)  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  --
  CONSTRAINT "PK_FlowCytometryProbe1" PRIMARY KEY ("WID")
);

CREATE INDEX "FlowCytometryProbe_DWID" ON "FlowCytometryProbe"("DataSetWID");

INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('FlowCytometryProbe', 'Type', 'acid-dye', 'An acid dye probe');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('FlowCytometryProbe', 'Type', 'antibody', 'An antibody probe');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('FlowCytometryProbe', 'Type', 'reporter', 'A protein probe');

CREATE TABLE "TranscriptionUnit"
(
  "WID"  bigint  NOT NULL,
  "Name"  character varying(255)  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  --
  CONSTRAINT "PK_TranscriptionUnit1" PRIMARY KEY ("WID")
);

CREATE INDEX "TranscriptionUnit_DWID" ON "TranscriptionUnit"("DataSetWID");


CREATE TABLE "TranscriptionUnitComponent"
(
  "Type"  character varying(100)  NOT NULL,
  "TranscriptionUnitWID"  bigint  NOT NULL,
  "OtherWID"  bigint  NOT NULL
  --
);


INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('TranscriptionUnitComponent', 'Type', 'gene', '');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('TranscriptionUnitComponent', 'Type', 'binding site', '');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('TranscriptionUnitComponent', 'Type', 'promoter', '');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('TranscriptionUnitComponent', 'Type', 'terminator', '');

CREATE TABLE "Chemical"
(
  "WID"  bigint  NOT NULL,
  "Name"  character varying(255)  NOT NULL,
  "Class"  character(1),
  "BeilsteinName"  character varying(50),
  "SystematicName"  character varying(255),
  "CAS"  character varying(50),
  "Charge"  smallint,
  "EmpiricalFormula"  character varying(50),
  "MolecularWeightCalc"  numeric,
  "MolecularWeightExp"  numeric,
  "OctH2OPartitionCoeff"  character varying(50),
  "PKA1"  numeric,
  "PKA2"  numeric,
  "PKA3"  numeric,
  "WaterSolubility"  character(1),
  "Smiles"  character varying(255),
  "DataSetWID"  bigint  NOT NULL,
  --
  CONSTRAINT "PK_Chemical1" PRIMARY KEY ("WID")
);

CREATE INDEX "CHEMICAL_DWID_NAME" ON "Chemical"("DataSetWID", "Name");
CREATE INDEX "CHEMICAL_NAME" ON "Chemical"("Name");


CREATE TABLE "Reaction"
(
  "WID"  bigint  NOT NULL,
  "Name"  character varying(250),
  "DeltaG"  character varying(50),
  "ECNumber"  character varying(50),
  "ECNumberProposed"  character varying(50),
  "Spontaneous"  character(1),
  "DataSetWID"  bigint  NOT NULL,
  --
  CONSTRAINT "PK_Reaction" PRIMARY KEY ("WID")
);

CREATE INDEX "REACTION_DWID" ON "Reaction"("DataSetWID");


CREATE TABLE "Interaction"
(
  "WID"  bigint  NOT NULL,
  "Type"  character varying(100),
  "Name"  character varying(120),
  "DataSetWID"  bigint  NOT NULL,
  --
  CONSTRAINT "PK_Interaction" PRIMARY KEY ("WID")
);

CREATE INDEX "INTERACTION_DWID" ON "Interaction"("DataSetWID");


CREATE TABLE "Protein"
(
  "WID"  bigint  NOT NULL,
  "Name"  text,
  "AASequence"  text,
  "Length"  bigint,
  "LengthApproximate"  character varying(10),
  "Charge"  smallint,
  "Fragment"  character(1),
  "MolecularWeightCalc"  numeric,
  "MolecularWeightExp"  numeric,
  "PICalc"  character varying(50),
  "PIExp"  character varying(50),
  "DataSetWID"  bigint  NOT NULL,
  --
  CONSTRAINT "PK_Protein" PRIMARY KEY ("WID")
);

CREATE INDEX "PROTEIN_DWID" ON "Protein"("DataSetWID");


CREATE TABLE "Feature"
(
  "WID"  bigint  NOT NULL,
  "Description"  character varying(1300),
  "Type"  character varying(50),
  "Class"  character varying(50),
  "SequenceType"  character(1)  NOT NULL,
  "SequenceWID"  bigint,
  "Variant"  text,
  "RegionOrPoint"  character varying(10),
  "PointType"  character varying(10),
  "StartPosition"  bigint,
  "EndPosition"  bigint,
  "StartPositionApproximate"  character varying(10),
  "EndPositionApproximate"  character varying(10),
  "ExperimentalSupport"  character(1),
  "ComputationalSupport"  character(1),
  "DataSetWID"  bigint  NOT NULL,
  --
  CONSTRAINT "PK_Feature" PRIMARY KEY ("WID")
);


INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Feature', 'Class', 'binding site', 'Identifies the presence of a DNA binding site');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Feature', 'Class', 'promoter', 'Identifies the presence of a promoter');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Feature', 'Class', 'terminator', 'Identifies the presence of a terminator');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Feature', 'Class', 'pseudogene', 'Identifies a pseudogene, whether non-transcribed or processed');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Feature', 'Class', 'ORF', 'Identifies a truly unknown open reading frame according to warehouse definition (no strong evidence that a product is produced)');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Feature', 'Class', 'partial', 'Qualifier: States that feature is not complete');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Feature', 'Class', 'unknown product', 'Identifies an unspecified product is produced from this genomic location, as stated in dataset');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Feature', 'Class', 'notable', 'Qualifier: Characterizes the feature value as notable (same as ''Exceptional'' in GB dataset)');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Feature', 'Class', 'by similarity', 'Qualifier: states that feature was derived by similarity analysis');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Feature', 'Class', 'potential', 'Qualifier: states that feature may be incorrect');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Feature', 'Class', 'probably', 'Qualifier: states that feature is probably correct');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Feature', 'SequenceType', 'P', 'Feature resides on a protein. Implies SequenceWID (if nonNULL) references a Protein');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Feature', 'SequenceType', 'S', 'Feature resides on a nucleic acid. Implies SequenceWID is nonNULL and references a Subsequence');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Feature', 'SequenceType', 'N', 'Feature resides on a nucleic acid. Implies SequenceWID (if nonNULL) references a NucleicAcid');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Feature', 'RegionOrPoint', 'region', 'Feature is specified by a start point and an end point on the sequence');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Feature', 'RegionOrPoint', 'point', 'Feature is specified by a single point on the sequence');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Feature', 'PointType', 'center', 'Feature is centered at location.');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Feature', 'PointType', 'left', 'Feature extends to the left (decreasing position) of location.');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Feature', 'PointType', 'right', 'Feature extends to the right (increasing position) of location.');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Feature', 'StartPositionApproximate', 'gt', 'The start position of the feature is greater than the actual position specified.');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Feature', 'StartPositionApproximate', 'lt', 'The start position of the feature is less than the actual position specified.');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Feature', 'StartPositionApproximate', 'ne', 'The start position of the feature is less than or greater than the actual position. All we know is that its not the exact position.');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Feature', 'EndPositionApproximate', 'gt', 'The end position of the feature is greater than the actual position specified.');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Feature', 'EndPositionApproximate', 'lt', 'The end position of the feature is less than the actual position specified.');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Feature', 'EndPositionApproximate', 'ne', 'The end position of the feature is less than or greater than the actual position. All we know is that its not the exact position.');

CREATE TABLE "Function"
(
  "WID"  bigint  NOT NULL,
  "Name"  character varying(255),
  "DataSetWID"  bigint  NOT NULL,
  --
  CONSTRAINT "PK_Function" PRIMARY KEY ("WID")
);



CREATE TABLE "EnzymaticReaction"
(
  "WID"  bigint  NOT NULL,
  "ReactionWID"  bigint  NOT NULL,
  "ProteinWID"  bigint  NOT NULL,
  "ComplexWID"  bigint,
  "ReactionDirection"  character varying(30),
  "DataSetWID"  bigint  NOT NULL,
  --
  CONSTRAINT "PK_EnzymaticReaction" PRIMARY KEY ("WID")
);

CREATE INDEX "ER_DATASETWID" ON "EnzymaticReaction"("DataSetWID");

INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('EnzymaticReaction', 'ReactionDirection', 'reversible', 'Reaction occurs in both directions in physiological settings.');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('EnzymaticReaction', 'ReactionDirection', 'physiol-left-to-right', 'Reaction occurs in the left-to-right direction in physiological settings.');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('EnzymaticReaction', 'ReactionDirection', 'physiol-right-to-left', 'Reaction occurs in the right-to-left direction in physiological settings.');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('EnzymaticReaction', 'ReactionDirection', 'irreversible-left-to-right', 'For all practical purposes, the reaction occurs only in the left-to-right direction in physiological settings, because of chemical properties of the reaction.');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('EnzymaticReaction', 'ReactionDirection', 'irreversible-right-to-left', 'For all practical purposes, the reaction occurs only in the right-to-left direction in physiological settings, because of chemical properties of the reaction.');

CREATE TABLE "GeneticCode"
(
  "WID"  bigint  NOT NULL,
  "NCBIID"  character varying(2),
  "Name"  character varying(100),
  "TranslationTable"  character varying(64),
  "StartCodon"  character varying(64),
  "DataSetWID"  bigint  NOT NULL,
  --
  CONSTRAINT "PK_GeneticCode" PRIMARY KEY ("WID")
);



CREATE TABLE "Division"
(
  "WID"  bigint  NOT NULL,
  "Code"  character varying(10),
  "Name"  character varying(100),
  "DataSetWID"  bigint  NOT NULL,
  --
  CONSTRAINT "PK_Division" PRIMARY KEY ("WID")
);



CREATE TABLE "Taxon"
(
  "WID"  bigint  NOT NULL,
  "ParentWID"  bigint,
  "Name"  character varying(100),
  "Rank"  character varying(100),
  "DivisionWID"  bigint,
  "InheritedDivision"  character(1),
  "GencodeWID"  bigint,
  "InheritedGencode"  character(1),
  "MCGencodeWID"  bigint,
  "InheritedMCGencode"  character(1),
  "DataSetWID"  bigint  NOT NULL,
  --
  CONSTRAINT "PK_Taxon" PRIMARY KEY ("WID")
);


INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Taxon', 'Rank', 'no rank', 'Origin:NCBI- Taxonomy databaseOrigin:NCBI- Taxonomy database. Parent: none');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Taxon', 'Rank', 'forma', 'Origin:NCBI- Taxonomy database.');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Taxon', 'Rank', 'varietas', 'Origin:NCBI- Taxonomy database');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Taxon', 'Rank', 'subspecies', 'Origin:NCBI- Taxonomy database Parent: species');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Taxon', 'Rank', 'species', 'Origin:NCBI- Taxonomy database. Parent: species subgroup');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Taxon', 'Rank', 'species subgroup', 'Origin:NCBI- Taxonomy database. Parent: species group');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Taxon', 'Rank', 'species group', 'Origin:NCBI- Taxonomy database. Parent: subgenus');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Taxon', 'Rank', 'subgenus', 'Origin:NCBI- Taxonomy database. Parent: genus');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Taxon', 'Rank', 'genus', 'Origin:NCBI- Taxonomy database. Parent: subtribe');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Taxon', 'Rank', 'subtribe', 'Origin:NCBI- Taxonomy database. Parent: tribe');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Taxon', 'Rank', 'tribe', 'Origin:NCBI- Taxonomy database. Parent: subfamily');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Taxon', 'Rank', 'subfamily', 'Origin:NCBI- Taxonomy database. Parent: family');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Taxon', 'Rank', 'family', 'Origin:NCBI- Taxonomy database Parent:superfamily');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Taxon', 'Rank', 'superfamily', 'Origin:NCBI- Taxonomy database. Parent: infraorder');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Taxon', 'Rank', 'infraorder', 'Origin:NCBI- Taxonomy database. Parent: parvorder');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Taxon', 'Rank', 'parvorder', 'Origin:NCBI- Taxonomy database. Parent: suborder');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Taxon', 'Rank', 'suborder', 'Origin:NCBI- Taxonomy database. Parent:order');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Taxon', 'Rank', 'order', 'Origin:NCBI- Taxonomy database Parent:superorder');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Taxon', 'Rank', 'superorder', 'Origin:NCBI- Taxonomy database. Parent: infraclass');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Taxon', 'Rank', 'infraclass', 'Origin:NCBI- Taxonomy database. Parent: subclass');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Taxon', 'Rank', 'subclass', 'Origin:NCBI- Taxonomy database.Parent: class');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Taxon', 'Rank', 'class', 'Origin:NCBI- Taxonomy database. Parent:superclass');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Taxon', 'Rank', 'superclass', 'Origin:NCBI- Taxonomy database. Parent: subphylum');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Taxon', 'Rank', 'subphylum', 'Origin:NCBI- Taxonomy database. Parent: phylum.');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Taxon', 'Rank', 'phylum', 'Origin:NCBI- Taxonomy database.Parent:superphylum');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Taxon', 'Rank', 'superphylum', 'Origin:NCBI- Taxonomy database. Parent:kingdom');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Taxon', 'Rank', 'subkingdom', 'Origin:NCBI- Taxonomy database. Parent: kingdom');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Taxon', 'Rank', 'kingdom', 'Origin:NCBI- Taxonomy database. Parent:superkingdom');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Taxon', 'Rank', 'superkingdom', 'Origin:NCBI- Taxonomy database. Parent: root');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Taxon', 'Rank', 'root', 'Origin:NCBI- Taxonomy database');

CREATE TABLE "BioSource"
(
  "WID"  bigint  NOT NULL,
  "MAGEClass"  character varying(100),
  "TaxonWID"  bigint,
  "Name"  character varying(200),
  "Strain"  character varying(220),
  "Organ"  character varying(50),
  "Organelle"  character varying(50),
  "Tissue"  character varying(100),
  "CellType"  character varying(50),
  "CellLine"  character varying(50),
  "ATCCId"  character varying(50),
  "Diseased"  character(1),
  "Disease"  character varying(250),
  "DevelopmentStage"  character varying(50),
  "Sex"  character varying(15),
  "DataSetWID"  bigint  NOT NULL,
  --
  CONSTRAINT "PK_BioSource" PRIMARY KEY ("WID")
);


INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('BioSource', 'Sex', 'Male', 'Sex of organism is male');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('BioSource', 'Sex', 'Female', 'Sex of organism is female');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('BioSource', 'Sex', 'Hermaphrodite', 'Sex of organism is hermaphrodite');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('BioSource', 'Sex', 'DoesNotApply', 'The notion of a sex does not apply to this organism');

CREATE TABLE "BioSubtype"
(
  "WID"  bigint  NOT NULL,
  "Type"  character varying(25),
  "Value"  character varying(50)  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  --
  CONSTRAINT "PK_BioSubtype" PRIMARY KEY ("WID")
);



CREATE TABLE "NucleicAcid"
(
  "WID"  bigint  NOT NULL,
  "Name"  character varying(200),
  "Type"  character varying(30)  NOT NULL,
  "Class"  character varying(30),
  "Topology"  character varying(30),
  "Strandedness"  character varying(30),
  "SequenceDerivation"  character varying(30),
  "Fragment"  character(1),
  "FullySequenced"  character(1),
  "MoleculeLength"  bigint,
  "MoleculeLengthApproximate"  character varying(10),
  "CumulativeLength"  bigint,
  "CumulativeLengthApproximate"  character varying(10),
  "GeneticCodeWID"  bigint,
  "BioSourceWID"  bigint,
  "DataSetWID"  bigint  NOT NULL,
  --
  CONSTRAINT "PK_NucleicAcid" PRIMARY KEY ("WID")
);


INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('NucleicAcid', 'Type', 'DNA', 'The molecule is composed of DNA');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('NucleicAcid', 'Type', 'RNA', 'The molecule is composed of RNA');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('NucleicAcid', 'Type', 'NA', 'The molecule is specified as a nucleic acid but whether of type DNA or RNA is not known');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('NucleicAcid', 'Class', 'pre-RNA', 'As it exists within the organism, the molecule is a pre-RNA molecule. NCBI: used when there is no evidence that mature RNA is produced');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('NucleicAcid', 'Class', 'mRNA', 'As it exists within the organism, the molecule is a mature mRNA. NCBI: used when there is evidence that mature mRNA is produced');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('NucleicAcid', 'Class', 'rRNA', 'As it exists within the organism, the molecule is an rRNA; NCBI used when there is evidence that mature rRNA is produced');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('NucleicAcid', 'Class', 'tRNA', 'As it exists within the organism, the molecule is a tRNA; NCBI: used when there is evidence that mature tRNA is produced');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('NucleicAcid', 'Class', 'snRNA', 'As it exists within the organism, the molecule is a small nuclear RNA. NCBI: used when there is evidence that snRNA is produced');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('NucleicAcid', 'Class', 'scRNA', 'As it exists within the organism, the molecule is an RNA which encodes small cytoplasmic ribonucleic proteins. NCBI: used when there is evidence that gene codes of small cytoplasmic (sc) ribonucleoproteins (RNP)s');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('NucleicAcid', 'Class', 'snoRNA', 'As it exists within the organism, the molecule is a small nucleolar RNA. NCBI: used when there is evidence that transcript is a small nucleolar RNA');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('NucleicAcid', 'Class', 'other', 'Check notes as to how we map this...');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('NucleicAcid', 'Class', 'chromosome', 'As it exists within the organism, the molecule is a chromosome');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('NucleicAcid', 'Class', 'plasmid', 'As it exists within the organism, the molecule is a plasmid');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('NucleicAcid', 'Class', 'organelle-chromosome', 'As it exists within the organism, the molecule is the chromosome of an organelle');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('NucleicAcid', 'Class', 'transposon', 'As it exists within the organism, the molecule is a transposon');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('NucleicAcid', 'Class', 'virus', 'As it exists within the organism, the molecule is a virus');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('NucleicAcid', 'Class', 'unknown', 'Used when the class is unknown');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('NucleicAcid', 'Topology', 'linear', 'The topology of the molecule is Linear');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('NucleicAcid', 'Topology', 'circular', 'The topology of the molecule is Circular');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('NucleicAcid', 'Topology', 'other', 'The topology of the molecule is neither circular nor linear but is known');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('NucleicAcid', 'Strandedness', 'ss', 'The molecule is single stranded');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('NucleicAcid', 'Strandedness', 'ds', 'The molecule is double stranded');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('NucleicAcid', 'Strandedness', 'mixed', 'The molecule is composed of single and double stranded regions');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('NucleicAcid', 'SequenceDerivation', 'virtual', 'NCBI: The sequence of the molecule is NOT known. NOTE: this should map to null.');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('NucleicAcid', 'SequenceDerivation', 'raw', 'POORLY DEFINED; seems to be used to indicate that the sequence of the molecule was actually generated from one single continuous molecule, as opposed to assembled together from sequences from different molecules (e.g., different clones)');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('NucleicAcid', 'SequenceDerivation', 'seg', 'NCBI: The sequence of the molecule is made up of collection of segments arranged according to specified coordinates, e.g., sequence was derived from assembling sequences from different clones');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('NucleicAcid', 'SequenceDerivation', 'reference', 'NCBI: The sequence of the molecule is constructed from existing Bioseqs;It behaves exactly like a segmented Bioseq in taking it''s data and character from the Bioseq to which it points.');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('NucleicAcid', 'SequenceDerivation', 'constructed', 'NCBI: The sequence of the molecule is constructed by assembling other Bioseqs');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('NucleicAcid', 'SequenceDerivation', 'consensus', 'NCBI: The sequence of the molecule represents a pattern typical of a sequence region or family of sequences;'' it summarizes attributes of an aligned collection of real sequences. Note that this is NOT A REAL OBJECT');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('NucleicAcid', 'SequenceDerivation', 'map', 'NCBI: The molecule does not have a sequence describing it, but rather a set of coordinates (restriction fragment order, genetic markers, physical map, etc)');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('NucleicAcid', 'MoleculeLengthApproximate', 'gt', 'The length of the molecule''s sequence is greater than the actual length specified');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('NucleicAcid', 'MoleculeLengthApproximate', 'lt', 'The length of the molecule''s sequence is less than the actual length specified');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('NucleicAcid', 'MoleculeLengthApproximate', 'ne', 'The length of the molecule''s sequence is less than or greater than the actual length. All we know is that its not the exact length');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('NucleicAcid', 'CumulativeLengthApproximate', 'gt', 'The total length of the molecule''s sequence is greater than the actual length');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('NucleicAcid', 'CumulativeLengthApproximate', 'lt', 'The total length of the molecule''s sequence is less than the actual length');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('NucleicAcid', 'CumulativeLengthApproximate', 'ne', 'The total length of the molecule''s sequence is less than or greater than the actual length. All we know is that its not the exact length');

CREATE TABLE "Subsequence"
(
  "WID"  bigint  NOT NULL,
  "NucleicAcidWID"  bigint  NOT NULL,
  "FullSequence"  character(1),
  "Sequence"  text,
  "Length"  bigint,
  "LengthApproximate"  character varying(10),
  "PercentGC"  numeric,
  "Version"  character varying(30),
  "DataSetWID"  bigint  NOT NULL,
  --
  CONSTRAINT "PK_Subsequence" PRIMARY KEY ("WID")
);



CREATE TABLE "Gene"
(
  "WID"  bigint  NOT NULL,
  "Name"  character varying(255),
  "NucleicAcidWID"  bigint,
  "SubsequenceWID"  bigint,
  "Type"  character varying(100),
  "GenomeID"  character varying(35),
  "CodingRegionStart"  bigint,
  "CodingRegionEnd"  bigint,
  "CodingRegionStartApproximate"  character varying(10),
  "CodingRegionEndApproximate"  character varying(10),
  "Direction"  character varying(25),
  "Interrupted"  character(1),
  "DataSetWID"  bigint  NOT NULL,
  --
  CONSTRAINT "PK_Gene" PRIMARY KEY ("WID")
);

CREATE INDEX "GENE_DATASETWID" ON "Gene"("DataSetWID");

INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Gene', 'Type', 'unknown', 'used when transcriptional status of gene is unknown');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Gene', 'Type', 'pre-RNA', 'used when there is no evidence that a mature RNA is ultimately produced by the gene');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Gene', 'Type', 'mRNA', 'used when there is evidence that a mature mRNA is ultimately produced by the gene');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Gene', 'Type', 'rRNA', 'used when there is evidence that a mature rRNA is ultimately produced by the gene');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Gene', 'Type', 'tRNA', 'used when there is evidence that a mature tRNA is ultimately produced by the gene');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Gene', 'Type', 'snRNA', 'used when there is evidence that an snRNA is ultimately produced by the gene');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Gene', 'Type', 'scRNA', '????');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Gene', 'Type', 'polypeptide', 'used when there is evidence that the ultimate product of this gene is proteinaceous');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Gene', 'Type', 'snoRNA', '???');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Gene', 'Type', 'other', 'catchall ?');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Gene', 'Direction', 'unknown', 'unknown the strand being transcribed is unknown');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Gene', 'Direction', 'forward', 'the plus strand is being transcribed');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Gene', 'Direction', 'reverse', 'the minus strand is being transcribed');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Gene', 'Direction', 'forward_and_reverse', 'both strands are being transcribed');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Gene', 'Direction', 'undefined_value', 'UNDEFINED; the NCBI documentation does not define this value');

CREATE TABLE "Pathway"
(
  "WID"  bigint  NOT NULL,
  "Name"  character varying(255)  NOT NULL,
  "Type"  character(1)  NOT NULL,
  "BioSourceWID"  bigint,
  "DataSetWID"  bigint  NOT NULL,
  --
  CONSTRAINT "PK_Pathway" PRIMARY KEY ("WID")
);

CREATE INDEX "PATHWAY_BSWID_WID_DWID" ON "Pathway"("BioSourceWID", "WID", "DataSetWID");
CREATE INDEX "PATHWAY_TYPE_WID_DWID" ON "Pathway"("Type", "WID", "DataSetWID");
CREATE INDEX "PATHWAY_DWID" ON "Pathway"("DataSetWID");

INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Pathway', 'Name', 'unknown', 'Name assigned when it is unknown or missing');

CREATE TABLE "Term"
(
  "WID"  bigint  NOT NULL,
  "Name"  character varying(255)  NOT NULL,
  "Definition"  text,
  "Hierarchical"  character(1),
  "Root"  character(1),
  "Obsolete"  character(1),
  "DataSetWID"  bigint  NOT NULL,
  --
  CONSTRAINT "PK_Term" PRIMARY KEY ("WID")
);



CREATE TABLE "Computation"
(
  "WID"  bigint  NOT NULL,
  "Name"  character varying(50)  NOT NULL,
  "Description"  text,
  "DataSetWID"  bigint  NOT NULL,
  --
  CONSTRAINT "PK_Computation" PRIMARY KEY ("WID")
);



CREATE TABLE "Citation"
(
  "WID"  bigint  NOT NULL,
  "Citation"  text,
  "PMID"  numeric,
  "Title"  character varying(255),
  "Authors"  character varying(255),
  "Publication"  character varying(255),
  "Publisher"  character varying(255),
  "Editor"  character varying(255),
  "Year"  character varying(255),
  "Volume"  character varying(255),
  "Issue"  character varying(255),
  "Pages"  character varying(255),
  "URI"  character varying(255),
  "DataSetWID"  bigint  NOT NULL,
  --
  CONSTRAINT "PK_Citation" PRIMARY KEY ("WID")
);

CREATE INDEX "CITATION_PMID" ON "Citation"("PMID");
--CREATE INDEX "CITATION_CITATION" ON "Citation"("Citation"(20));


CREATE TABLE "Archive"
(
  "WID"  bigint  NOT NULL,
  "OtherWID"  bigint  NOT NULL,
  "Format"  character varying(10)  NOT NULL,
  "Contents"  bytea,
  "URL"  text,
  "ToolName"  character varying(50),
  "DataSetWID"  bigint  NOT NULL,
  --
  CONSTRAINT "PK_ARCHIVE" PRIMARY KEY ("WID")
);



CREATE TABLE "Experiment"
(
  "WID"  bigint  NOT NULL,
  "Type"  character varying(50)  NOT NULL,
  "ContactWID"  bigint,
  "ArchiveWID"  bigint,
  "StartDate"  timestamp,
  "EndDate"  timestamp,
  "Description"  text,
  "GroupWID"  bigint,
  "GroupType"  character varying(50),
  "GroupSize"  int  NOT NULL,
  "GroupIndex"  bigint,
  "TimePoint"  bigint,
  "TimeUnit"  character varying(20),
  "DataSetWID"  bigint  NOT NULL,
  "BioSourceWID"  bigint,
  --
  CONSTRAINT "PK_Experiment" PRIMARY KEY ("WID")
);


INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Experiment', 'GroupType', 'replicate', 'All subexperiments attempt to replicate identical experimental conditions and parameters');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Experiment', 'GroupType', 'variant', 'Subexperiments are variations upon an experimental procedure, technique, conditions, and/or parameters');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Experiment', 'GroupType', 'time-series', 'Subexperiments consist of observations according to a specific schedule');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Experiment', 'GroupType', 'step', 'Subexperiments comprise an experimental procedure and consist of an ordered sequence');

CREATE TABLE "ExperimentData"
(
  "WID"  bigint  NOT NULL,
  "ExperimentWID"  bigint  NOT NULL,
  "Data"  text,
  "MageData"  bigint,
  "Role"  character varying(50)  NOT NULL,
  "Kind"  character(1)  NOT NULL,
  "DateProduced"  timestamp,
  "OtherWID"  bigint,
  "DataSetWID"  bigint  NOT NULL,
  --
  CONSTRAINT "PK_ExpData" PRIMARY KEY ("WID")
);


INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('ExperimentData', 'Kind', 'O', 'Data is an observation');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('ExperimentData', 'Kind', 'C', 'Data is computed from an observation');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('ExperimentData', 'Kind', 'P', 'Data is a parameter to a procedure or compuatation');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('ExperimentData', 'Kind', 'M', 'Data is metadata describing other data');

CREATE TABLE "Support"
(
  "WID"  bigint  NOT NULL,
  "OtherWID"  bigint  NOT NULL,
  "Type"  character varying(100),
  "EvidenceType"  character varying(100),
  "Confidence"  numeric,
  "DataSetWID"  bigint  NOT NULL,
  --
  CONSTRAINT "CK_Support" CHECK ("Confidence" > 0 AND "Confidence" <= 1),
  CONSTRAINT "PK_Support" PRIMARY KEY ("WID")
);


INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Support', 'Type', 'protein function', 'The evidence supports protein function');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Support', 'Type', 'protein existence', 'The evidence supports protein existence');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Support', 'EvidenceType', 'computational', 'Valid when type is ''protein function''. Protein function is supported by a computational prediction (example: existence of a gene could be supported by a gene-finding program).');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Support', 'EvidenceType', 'experimental', 'Valid when type is ''protein function''. Protein function is supported by data from a wet-lab experiment (example: existence of a gene could be supported observation of the gene''s mRNA product in a microarray experiment).');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Support', 'EvidenceType', 'Evidence at protein level', 'Valid when type is ''protein existence''. Protein existence is supported by the evidences at protein level (e.g. clear identification by mass spectrometry).');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Support', 'EvidenceType', 'Evidence at transcript level', 'Valid when type is ''protein existence''. Protein existence is supported by evidences at transcript level (e.g. Northern blot).');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Support', 'EvidenceType', 'Inferred from homology', 'Valid when type is ''protein existence''. Protein existence is inferred by homology (strong sequence similarity to known proteins in related species).');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Support', 'EvidenceType', 'Predicted', 'Valid when type is ''protein existence''. Protein existence is predicted');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('Support', 'EvidenceType', 'Uncertain', 'Valid when type is ''protein existence''. Protein existence is uncertain (e.g. dubious sequences that could be the erroneous translation of a pseudogene).');

CREATE TABLE "ChemicalAtom"
(
  "ChemicalWID"  bigint  NOT NULL,
  "AtomIndex"  smallint  NOT NULL,
  "Atom"  character varying(2)  NOT NULL,
  "Charge"  smallint  NOT NULL,
  "X"  numeric,
  "Y"  numeric,
  "Z"  numeric,
  "StereoParity"  numeric,
  --
  CONSTRAINT "UN_ChemicalAtom" UNIQUE ("ChemicalWID", "AtomIndex")
);


INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('ChemicalAtom', 'StereoParity', '0', 'Not stereo.');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('ChemicalAtom', 'StereoParity', '1', 'Odd parity.');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('ChemicalAtom', 'StereoParity', '2', 'Even parity.');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('ChemicalAtom', 'StereoParity', '3', 'Either or unmarked stereo center.');

CREATE TABLE "ChemicalBond"
(
  "ChemicalWID"  bigint  NOT NULL,
  "Atom1Index"  smallint  NOT NULL,
  "Atom2Index"  smallint  NOT NULL,
  "BondType"  smallint  NOT NULL,
  "BondStereo"  numeric
  --
);


INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('ChemicalBond', 'BondType', '1', 'Single bond.');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('ChemicalBond', 'BondType', '2', 'Double bond.');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('ChemicalBond', 'BondType', '3', 'Triple bond.');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('ChemicalBond', 'BondType', '4', 'Aromatic bond.');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('ChemicalBond', 'BondType', '5', 'Single or double bond.');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('ChemicalBond', 'BondType', '6', 'Single or aromatic bond.');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('ChemicalBond', 'BondType', '7', 'Double or aromatic bond.');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('ChemicalBond', 'BondType', '8', 'Any bond.');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('ChemicalBond', 'BondStereo', '0', 'For single bonds, 0 = not stereo.   Four double bonds, 0 = use X,Y,Z coords to determine cis or trans.');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('ChemicalBond', 'BondStereo', '1', 'For single bonds, 1 = up.');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('ChemicalBond', 'BondStereo', '3', 'For double bonds, 3 = cis or trans (either) (presumably meaning unspecified or a mixture).');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('ChemicalBond', 'BondStereo', '4', 'For single bonds, 4 = either.');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('ChemicalBond', 'BondStereo', '6', 'For single bonds, 6 = down.');

CREATE TABLE "EnzReactionCofactor"
(
  "EnzymaticReactionWID"  bigint  NOT NULL,
  "ChemicalWID"  bigint  NOT NULL,
  "Prosthetic"  character(1)
  --
);



CREATE TABLE "EnzReactionAltCompound"
(
  "EnzymaticReactionWID"  bigint  NOT NULL,
  "PrimaryWID"  bigint  NOT NULL,
  "AlternativeWID"  bigint  NOT NULL,
  "Cofactor"  character(1)
  --
);



CREATE TABLE "EnzReactionInhibitorActivator"
(
  "EnzymaticReactionWID"  bigint  NOT NULL,
  "CompoundWID"  bigint  NOT NULL,
  "InhibitOrActivate"  character(1),
  "Mechanism"  character(1),
  "PhysioRelevant"  character(1)
  --
);


INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('EnzReactionInhibitorActivator', 'InhibitOrActivate', 'I', 'Specifies a compound that inhibits an enzyme.');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('EnzReactionInhibitorActivator', 'InhibitOrActivate', 'A', 'Specifies a compound that activates an enzyme.');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('EnzReactionInhibitorActivator', 'Mechanism', 'A', 'The mechanism of inhibition or activation is allosteric.');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('EnzReactionInhibitorActivator', 'Mechanism', 'I', 'The mechanism of inhibition or activation is irreversible.');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('EnzReactionInhibitorActivator', 'Mechanism', 'C', 'The mechanism of inhibition or activation is competitive.');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('EnzReactionInhibitorActivator', 'Mechanism', 'N', 'The mechanism of inhibition or activation is neither allosteric nor competitive.');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('EnzReactionInhibitorActivator', 'Mechanism', 'O', 'The mechanism of inhibition or activation is known but not in this controlled vocabulary.');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('EnzReactionInhibitorActivator', 'Mechanism', 'U', 'The mechanism of inhibition or activation is unknown.');

CREATE TABLE "Location"
(
  "ProteinWID"  bigint  NOT NULL,
  "Location"  character varying(100)  NOT NULL
  --
);



CREATE TABLE "PathwayLink"
(
  "Pathway1WID"  bigint  NOT NULL,
  "Pathway2WID"  bigint  NOT NULL,
  "ChemicalWID"  bigint  NOT NULL
  --
);



CREATE TABLE "PathwayReaction"
(
  "PathwayWID"  bigint  NOT NULL,
  "ReactionWID"  bigint  NOT NULL,
  "PriorReactionWID"  bigint,
  "Hypothetical"  character(1)  NOT NULL
  --
);

CREATE INDEX "PR_PATHWID_REACTIONWID" ON "PathwayReaction"("PathwayWID", "ReactionWID");


CREATE TABLE "Product"
(
  "ReactionWID"  bigint  NOT NULL,
  "OtherWID"  bigint  NOT NULL,
  "Coefficient"  smallint  NOT NULL
  --
);



CREATE TABLE "Reactant"
(
  "ReactionWID"  bigint  NOT NULL,
  "OtherWID"  bigint  NOT NULL,
  "Coefficient"  smallint  NOT NULL
  --
);



CREATE TABLE "InteractionParticipant"
(
  "InteractionWID"  bigint  NOT NULL,
  "OtherWID"  bigint  NOT NULL,
  "Coefficient"  smallint
  --
);

CREATE INDEX "PR_INTERACTIONWID_OTHERWID" ON "InteractionParticipant"("InteractionWID", "OtherWID");


CREATE TABLE "SequenceMatch"
(
  "QueryWID"  bigint  NOT NULL,
  "MatchWID"  bigint  NOT NULL,
  "ComputationWID"  bigint  NOT NULL,
  "EValue"  numeric,
  "PValue"  numeric,
  "PercentIdentical"  numeric,
  "PercentSimilar"  numeric,
  "Rank"  smallint,
  "Length"  bigint,
  "QueryStart"  bigint,
  "QueryEnd"  bigint,
  "MatchStart"  bigint,
  "MatchEnd"  int
  --
);



CREATE TABLE "Subunit"
(
  "ComplexWID"  bigint  NOT NULL,
  "SubunitWID"  bigint  NOT NULL,
  "Coefficient"  smallint
  --
);



CREATE TABLE "SuperPathway"
(
  "SubPathwayWID"  bigint  NOT NULL,
  "SuperPathwayWID"  bigint  NOT NULL
  --
);



CREATE TABLE "TermRelationship"
(
  "TermWID"  bigint  NOT NULL,
  "RelatedTermWID"  bigint  NOT NULL,
  "Relationship"  character varying(10)  NOT NULL
  --
);



CREATE TABLE "RelatedTerm"
(
  "TermWID"  bigint  NOT NULL,
  "OtherWID"  bigint  NOT NULL,
  "Relationship"  character varying(50)
  --
);


INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('RelatedTerm', 'Relationship', 'keyword', 'The term is a keyword that characterizes the object.');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('RelatedTerm', 'Relationship', 'superclass', 'The term names a class, and the object is an instance or a subclass of that class.');

CREATE TABLE "CitationWIDOtherWID"
(
  "OtherWID"  bigint  NOT NULL,
  "CitationWID"  bigint  NOT NULL
  --
);



CREATE TABLE "CommentTable"
(
  "OtherWID"  bigint  NOT NULL,
  "Comm"  text
  --
);



CREATE TABLE "CrossReference"
(
  "OtherWID"  bigint  NOT NULL,
  "CrossWID"  bigint,
  "XID"  character varying(50),
  "Type"  character varying(20),
  "Version"  character varying(10),
  "Relationship"  character varying(50),
  "DataBaseName"  character varying(255)
  --
);


INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('CrossReference', 'Type', 'Accession', 'Type of XID is an accession number (not guaranteed to be unique across datasets)');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('CrossReference', 'Type', 'GUID', 'Type of XID is a global unique identifier (guaranteed to be unique across datasets for a given database provider)');

CREATE TABLE "Description"
(
  "OtherWID"  bigint  NOT NULL,
  "TableName"  character varying(30)  NOT NULL,
  "Comm"  text
  --
);



CREATE TABLE "DBID"
(
  "OtherWID"  bigint  NOT NULL,
  "XID"  character varying(150)  NOT NULL,
  "Type"  character varying(20),
  "Version"  character varying(10)
  --
);

CREATE INDEX "DBID_XID_OTHERWID" ON "DBID"("XID", "OtherWID");
CREATE INDEX "DBID_OTHERWID" ON "DBID"("OtherWID");

INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('DBID', 'Type', 'Accession', 'Type of XID is an accession number (not guaranteed to be unique across datasets)');
INSERT INTO "Enumeration" ("TableName", "ColumnName", "Value", "Meaning") VALUES ('DBID', 'Type', 'GUID', 'Type of XID is a global unique identifier (guaranteed to be unique across datasets for a given database provider)');

CREATE TABLE "SynonymTable"
(
  "OtherWID"  bigint  NOT NULL,
  "Syn"  character varying(255)  NOT NULL
  --
);

CREATE INDEX "SYNONYM_OTHERWID_SYN" ON "SynonymTable"("OtherWID", "Syn");


CREATE TABLE "ToolAdvice"
(
  "OtherWID"  bigint  NOT NULL,
  "ToolName"  character varying(50)  NOT NULL,
  "Advice"  text
  --
);



CREATE TABLE "BioSourceWIDBioSubtypeWID"
(
  "BioSourceWID"  bigint  NOT NULL,
  "BioSubtypeWID"  bigint  NOT NULL
  --
);



CREATE TABLE "BioSourceWIDGeneWID"
(
  "BioSourceWID"  bigint  NOT NULL,
  "GeneWID"  bigint  NOT NULL
  --
);



CREATE TABLE "BioSourceWIDProteinWID"
(
  "BioSourceWID"  bigint  NOT NULL,
  "ProteinWID"  bigint  NOT NULL
  --
);



CREATE TABLE "GeneWIDNucleicAcidWID"
(
  "GeneWID"  bigint  NOT NULL,
  "NucleicAcidWID"  bigint  NOT NULL
  --
);



CREATE TABLE "GeneWIDProteinWID"
(
  "GeneWID"  bigint  NOT NULL,
  "ProteinWID"  bigint  NOT NULL
  --
);



CREATE TABLE "ProteinWIDFunctionWID"
(
  "ProteinWID"  bigint  NOT NULL,
  "FunctionWID"  bigint  NOT NULL
  --
);



CREATE TABLE "ExperimentRelationship"
(
  "ExperimentWID"  bigint  NOT NULL,
  "RelatedExperimentWID"  bigint  NOT NULL
  --
);



CREATE TABLE "GelLocation"
(
  "WID"  bigint  NOT NULL,
  "SpotWID"  bigint  NOT NULL,
  "Xcoord"  numeric,
  "Ycoord"  numeric,
  "refGel"  character varying(1),
  "ExperimentWID"  bigint  NOT NULL,
  "DatasetWID"  bigint  NOT NULL,
  --
  CONSTRAINT "PK_GelLocation" PRIMARY KEY ("WID")
);



CREATE TABLE "ProteinWIDSpotWID"
(
  "ProteinWID"  bigint  NOT NULL,
  "SpotWID"  bigint  NOT NULL
  --
);



CREATE TABLE "Spot"
(
  "WID"  bigint  NOT NULL,
  "SpotId"  character varying(25),
  "MolecularWeightEst"  numeric,
  "PIEst"  character varying(50),
  "DatasetWID"  bigint  NOT NULL,
  --
  CONSTRAINT "PK_Spot" PRIMARY KEY ("WID")
);



CREATE TABLE "SpotIdMethod"
(
  "WID"  bigint  NOT NULL,
  "MethodName"  character varying(100)  NOT NULL,
  "MethodDesc"  character varying(500),
  "MethodAbbrev"  character varying(10),
  "DatasetWID"  bigint  NOT NULL,
  --
  CONSTRAINT "PK_SpotIdMethod" PRIMARY KEY ("WID")
);



CREATE TABLE "SpotWIDSpotIdMethodWID"
(
  "SpotWID"  bigint  NOT NULL,
  "SpotIdMethodWID"  bigint  NOT NULL
  --
);



CREATE TABLE "NameValueType"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "Name"  character varying(255),
  "Value"  character varying(255),
  "Type_"  character varying(255),
  "NameValueType_PropertySets"  bigint,
  "OtherWID"  bigint,
  --
  CONSTRAINT "PK_NameValueType" PRIMARY KEY ("WID")
);



CREATE TABLE "Contact"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "MAGEClass"  character varying(100)  NOT NULL,
  "Identifier"  character varying(255),
  "Name"  character varying(255),
  "URI"  character varying(255),
  "Address"  character varying(255),
  "Phone"  character varying(255),
  "TollFreePhone"  character varying(255),
  "Email"  character varying(255),
  "Fax"  character varying(255),
  "Parent"  bigint,
  "LastName"  character varying(255),
  "FirstName"  character varying(255),
  "MidInitials"  character varying(255),
  "Affiliation"  bigint,
  --
  CONSTRAINT "PK_Contact" PRIMARY KEY ("WID")
);



CREATE TABLE "ArrayDesign"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "MAGEClass"  character varying(100)  NOT NULL,
  "Identifier"  character varying(255),
  "Name"  character varying(255),
  "Version"  character varying(255),
  "NumberOfFeatures"  smallint,
  "SurfaceType"  bigint,
  --
  CONSTRAINT "PK_ArrayDesign" PRIMARY KEY ("WID")
);



CREATE TABLE "DesignElementGroup"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "MAGEClass"  character varying(100)  NOT NULL,
  "Identifier"  character varying(255),
  "Name"  character varying(255),
  "ArrayDesign_FeatureGroups"  bigint,
  "DesignElementGroup_Species"  bigint,
  "FeatureWidth"  numeric,
  "FeatureLength"  numeric,
  "FeatureHeight"  numeric,
  "FeatureGroup_TechnologyType"  bigint,
  "FeatureGroup_FeatureShape"  bigint,
  "FeatureGroup_DistanceUnit"  bigint,
  --
  CONSTRAINT "PK_DesignElementGroup" PRIMARY KEY ("WID")
);



CREATE TABLE "Zone"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "Identifier"  character varying(255),
  "Name"  character varying(255),
  "Row_"  smallint,
  "Column_"  smallint,
  "UpperLeftX"  numeric,
  "UpperLeftY"  numeric,
  "LowerRightX"  numeric,
  "LowerRightY"  numeric,
  "Zone_DistanceUnit"  bigint,
  "ZoneGroup_ZoneLocations"  bigint,
  --
  CONSTRAINT "PK_Zone" PRIMARY KEY ("WID")
);



CREATE TABLE "ZoneGroup"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "PhysicalArrayDesign_ZoneGroups"  bigint,
  "SpacingsBetweenZonesX"  numeric,
  "SpacingsBetweenZonesY"  numeric,
  "ZonesPerX"  smallint,
  "ZonesPerY"  smallint,
  "ZoneGroup_DistanceUnit"  bigint,
  "ZoneGroup_ZoneLayout"  bigint,
  --
  CONSTRAINT "PK_ZoneGroup" PRIMARY KEY ("WID")
);



CREATE TABLE "ZoneLayout"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "NumFeaturesPerRow"  smallint,
  "NumFeaturesPerCol"  smallint,
  "SpacingBetweenRows"  numeric,
  "SpacingBetweenCols"  numeric,
  "ZoneLayout_DistanceUnit"  bigint,
  --
  CONSTRAINT "PK_ZoneLayout" PRIMARY KEY ("WID")
);



CREATE TABLE "ExperimentDesign"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "Experiment_ExperimentDesigns"  bigint,
  "QualityControlDescription"  bigint,
  "NormalizationDescription"  bigint,
  "ReplicateDescription"  bigint,
  --
  CONSTRAINT "PK_ExperimentDesign" PRIMARY KEY ("WID")
);



CREATE TABLE "ExperimentalFactor"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "Identifier"  character varying(255),
  "Name"  character varying(255),
  "ExperimentDesign"  bigint,
  "ExperimentalFactor_Category"  bigint,
  --
  CONSTRAINT "PK_ExperimentalFactor" PRIMARY KEY ("WID")
);



CREATE TABLE "FactorValue"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "Identifier"  character varying(255),
  "Name"  character varying(255),
  "ExperimentalFactor"  bigint,
  "ExperimentalFactor2"  bigint,
  "FactorValue_Measurement"  bigint,
  "FactorValue_Value"  bigint,
  --
  CONSTRAINT "PK_FactorValue" PRIMARY KEY ("WID")
);



CREATE TABLE "QuantitationType"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "MAGEClass"  character varying(100)  NOT NULL,
  "Identifier"  character varying(255),
  "Name"  character varying(255),
  "IsBackground"  character(1),
  "Channel"  bigint,
  "QuantitationType_Scale"  bigint,
  "QuantitationType_DataType"  bigint,
  "TargetQuantitationType"  bigint,
  --
  CONSTRAINT "PK_QuantitationType" PRIMARY KEY ("WID")
);



CREATE TABLE "Database_"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "Identifier"  character varying(255),
  "Name"  character varying(255),
  "Version"  character varying(255),
  "URI"  character varying(255),
  --
  CONSTRAINT "PK_Database_" PRIMARY KEY ("WID")
);



CREATE TABLE "Array_"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "Identifier"  character varying(255),
  "Name"  character varying(255),
  "ArrayIdentifier"  character varying(255),
  "ArrayXOrigin"  numeric,
  "ArrayYOrigin"  numeric,
  "OriginRelativeTo"  character varying(255),
  "ArrayDesign"  bigint,
  "Information"  bigint,
  "ArrayGroup"  bigint,
  --
  CONSTRAINT "PK_Array_" PRIMARY KEY ("WID")
);



CREATE TABLE "ArrayGroup"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "Identifier"  character varying(255),
  "Name"  character varying(255),
  "Barcode"  character varying(255),
  "ArraySpacingX"  numeric,
  "ArraySpacingY"  numeric,
  "NumArrays"  smallint,
  "OrientationMark"  character varying(255),
  "OrientationMarkPosition"  character varying(25),
  "Width"  numeric,
  "Length"  numeric,
  "ArrayGroup_SubstrateType"  bigint,
  "ArrayGroup_DistanceUnit"  bigint,
  --
  CONSTRAINT "PK_ArrayGroup" PRIMARY KEY ("WID")
);



CREATE TABLE "ArrayManufacture"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "Identifier"  character varying(255),
  "Name"  character varying(255),
  "ManufacturingDate"  character varying(255),
  "Tolerance"  numeric,
  --
  CONSTRAINT "PK_ArrayManufacture" PRIMARY KEY ("WID")
);



CREATE TABLE "ArrayManufactureDeviation"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "Array_"  bigint,
  --
  CONSTRAINT "PK_ArrayManufactureDeviation" PRIMARY KEY ("WID")
);



CREATE TABLE "FeatureDefect"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "ArrayManufactureDeviation"  bigint,
  "FeatureDefect_DefectType"  bigint,
  "FeatureDefect_PositionDelta"  bigint,
  "Feature"  bigint,
  --
  CONSTRAINT "PK_FeatureDefect" PRIMARY KEY ("WID")
);



CREATE TABLE "Fiducial"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "ArrayGroup_Fiducials"  bigint,
  "Fiducial_FiducialType"  bigint,
  "Fiducial_DistanceUnit"  bigint,
  "Fiducial_Position"  bigint,
  --
  CONSTRAINT "PK_Fiducial" PRIMARY KEY ("WID")
);



CREATE TABLE "ManufactureLIMS"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "MAGEClass"  character varying(100)  NOT NULL,
  "ArrayManufacture_FeatureLIMSs"  bigint,
  "Quality"  character varying(255),
  "Feature"  bigint,
  "BioMaterial"  bigint,
  "ManufactureLIMS_IdentifierLIMS"  bigint,
  "BioMaterialPlateIdentifier"  character varying(255),
  "BioMaterialPlateRow"  character varying(255),
  "BioMaterialPlateCol"  character varying(255),
  --
  CONSTRAINT "PK_ManufactureLIMS" PRIMARY KEY ("WID")
);



CREATE TABLE "PositionDelta"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "DeltaX"  numeric,
  "DeltaY"  numeric,
  "PositionDelta_DistanceUnit"  bigint,
  --
  CONSTRAINT "PK_PositionDelta" PRIMARY KEY ("WID")
);



CREATE TABLE "ZoneDefect"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "ArrayManufactureDeviation"  bigint,
  "ZoneDefect_DefectType"  bigint,
  "ZoneDefect_PositionDelta"  bigint,
  "Zone"  bigint,
  --
  CONSTRAINT "PK_ZoneDefect" PRIMARY KEY ("WID")
);



CREATE TABLE "SeqFeatureLocation"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "SeqFeature_Regions"  bigint,
  "StrandType"  character varying(255),
  "SeqFeatureLocation_Subregions"  bigint,
  "SeqFeatureLocation_Coordinate"  bigint,
  --
  CONSTRAINT "PK_SeqFeatureLocation" PRIMARY KEY ("WID")
);



CREATE TABLE "SequencePosition"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "MAGEClass"  character varying(100)  NOT NULL,
  "Start_"  smallint,
  "End"  smallint,
  "CompositeCompositeMap"  bigint,
  "Composite"  bigint,
  "ReporterCompositeMap"  bigint,
  "Reporter"  bigint,
  --
  CONSTRAINT "PK_SequencePosition" PRIMARY KEY ("WID")
);



CREATE TABLE "DesignElement"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "MAGEClass"  character varying(100)  NOT NULL,
  "Identifier"  character varying(255),
  "Name"  character varying(255),
  "FeatureGroup_Features"  bigint,
  "DesignElement_ControlType"  bigint,
  "Feature_Position"  bigint,
  "Zone"  bigint,
  "Feature_FeatureLocation"  bigint,
  "FeatureGroup"  bigint,
  "Reporter_WarningType"  bigint,
  --
  CONSTRAINT "PK_DesignElement" PRIMARY KEY ("WID")
);



CREATE TABLE "FeatureInformation"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "Feature"  bigint,
  "FeatureReporterMap"  bigint,
  --
  CONSTRAINT "PK_FeatureInformation" PRIMARY KEY ("WID")
);



CREATE TABLE "FeatureLocation"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "Row_"  smallint,
  "Column_"  smallint,
  --
  CONSTRAINT "PK_FeatureLocation" PRIMARY KEY ("WID")
);



CREATE TABLE "MismatchInformation"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "CompositePosition"  bigint,
  "FeatureInformation"  bigint,
  "StartCoord"  smallint,
  "NewSequence"  character varying(255),
  "ReplacedLength"  smallint,
  "ReporterPosition"  bigint,
  --
  CONSTRAINT "PK_MismatchInformation" PRIMARY KEY ("WID")
);



CREATE TABLE "Position_"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "X"  numeric,
  "Y"  numeric,
  "Position_DistanceUnit"  bigint,
  --
  CONSTRAINT "PK_Position_" PRIMARY KEY ("WID")
);



CREATE TABLE "BioEvent"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "MAGEClass"  character varying(100)  NOT NULL,
  "Identifier"  character varying(255),
  "Name"  character varying(255),
  "CompositeSequence"  bigint,
  "Reporter"  bigint,
  "CompositeSequence2"  bigint,
  "BioAssayMapTarget"  bigint,
  "TargetQuantitationType"  bigint,
  "DerivedBioAssayDataTarget"  bigint,
  "QuantitationTypeMapping"  bigint,
  "DesignElementMapping"  bigint,
  "Transformation_BioAssayMapping"  bigint,
  "BioMaterial_Treatments"  bigint,
  "Order_"  smallint,
  "Treatment_Action"  bigint,
  "Treatment_ActionMeasurement"  bigint,
  "Array_"  bigint,
  "PhysicalBioAssayTarget"  bigint,
  "PhysicalBioAssay"  bigint,
  "Target"  bigint,
  "PhysicalBioAssaySource"  bigint,
  "MeasuredBioAssayTarget"  bigint,
  "PhysicalBioAssay2"  bigint,
  --
  CONSTRAINT "PK_BioEvent" PRIMARY KEY ("WID")
);



CREATE TABLE "BioAssayData"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "MAGEClass"  character varying(100)  NOT NULL,
  "Identifier"  character varying(255),
  "Name"  character varying(255),
  "BioAssayDimension"  bigint,
  "DesignElementDimension"  bigint,
  "QuantitationTypeDimension"  bigint,
  "BioAssayData_BioDataValues"  bigint,
  "ProducerTransformation"  bigint,
  --
  CONSTRAINT "PK_BioAssayData" PRIMARY KEY ("WID")
);



CREATE TABLE "BioAssayDimension"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "Identifier"  character varying(255),
  "Name"  character varying(255),
  --
  CONSTRAINT "PK_BioAssayDimension" PRIMARY KEY ("WID")
);



CREATE TABLE "BioAssayMapping"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  --
  CONSTRAINT "PK_BioAssayMapping" PRIMARY KEY ("WID")
);



CREATE TABLE "BioAssayTuple"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "BioAssay"  bigint,
  "BioDataTuples_BioAssayTuples"  bigint,
  --
  CONSTRAINT "PK_BioAssayTuple" PRIMARY KEY ("WID")
);



CREATE TABLE "BioDataValues"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "MAGEClass"  character varying(100)  NOT NULL,
  "Order_"  character varying(25),
  "BioDataCube_DataInternal"  bigint,
  "BioDataCube_DataExternal"  bigint,
  --
  CONSTRAINT "PK_BioDataValues" PRIMARY KEY ("WID")
);



CREATE TABLE "DataExternal"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "DataFormat"  character varying(255),
  "DataFormatInfoURI"  character varying(255),
  "FilenameURI"  character varying(255),
  --
  CONSTRAINT "PK_DataExternal" PRIMARY KEY ("WID")
);



CREATE TABLE "DataInternal"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  --
  CONSTRAINT "PK_DataInternal" PRIMARY KEY ("WID")
);



CREATE TABLE "Datum"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "Value"  character varying(255),
  --
  CONSTRAINT "PK_Datum" PRIMARY KEY ("WID")
);



CREATE TABLE "DesignElementDimension"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "MAGEClass"  character varying(100)  NOT NULL,
  "Identifier"  character varying(255),
  "Name"  character varying(255),
  --
  CONSTRAINT "PK_DesignElementDimension" PRIMARY KEY ("WID")
);



CREATE TABLE "DesignElementMapping"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  --
  CONSTRAINT "PK_DesignElementMapping" PRIMARY KEY ("WID")
);



CREATE TABLE "DesignElementTuple"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "BioAssayTuple"  bigint,
  "DesignElement"  bigint,
  --
  CONSTRAINT "PK_DesignElementTuple" PRIMARY KEY ("WID")
);



CREATE TABLE "QuantitationTypeDimension"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "Identifier"  character varying(255),
  "Name"  character varying(255),
  --
  CONSTRAINT "PK_QuantitationTypeDimension" PRIMARY KEY ("WID")
);



CREATE TABLE "QuantitationTypeMapping"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  --
  CONSTRAINT "PK_QuantitationTypeMapping" PRIMARY KEY ("WID")
);



CREATE TABLE "QuantitationTypeTuple"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "DesignElementTuple"  bigint,
  "QuantitationType"  bigint,
  "QuantitationTypeTuple_Datum"  bigint,
  --
  CONSTRAINT "PK_QuantitationTypeTuple" PRIMARY KEY ("WID")
);



CREATE TABLE "BioMaterialMeasurement"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "BioMaterial"  bigint,
  "Measurement"  bigint,
  "Treatment"  bigint,
  "BioAssayCreation"  bigint,
  --
  CONSTRAINT "PK_BioMaterialMeasurement" PRIMARY KEY ("WID")
);



CREATE TABLE "CompoundMeasurement"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "Compound_ComponentCompounds"  bigint,
  "Compound"  bigint,
  "Measurement"  bigint,
  "Treatment_CompoundMeasurements"  bigint,
  --
  CONSTRAINT "PK_CompoundMeasurement" PRIMARY KEY ("WID")
);



CREATE TABLE "BioAssay"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "MAGEClass"  character varying(100)  NOT NULL,
  "Identifier"  character varying(255),
  "Name"  character varying(255),
  "DerivedBioAssay_Type"  bigint,
  "FeatureExtraction"  bigint,
  "BioAssayCreation"  bigint,
  --
  CONSTRAINT "PK_BioAssay" PRIMARY KEY ("WID")
);



CREATE TABLE "Channel"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "Identifier"  character varying(255),
  "Name"  character varying(255),
  --
  CONSTRAINT "PK_Channel" PRIMARY KEY ("WID")
);



CREATE TABLE "Image"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "Identifier"  character varying(255),
  "Name"  character varying(255),
  "URI"  character varying(255),
  "Image_Format"  bigint,
  "PhysicalBioAssay"  bigint,
  --
  CONSTRAINT "PK_Image" PRIMARY KEY ("WID")
);



CREATE TABLE "BioAssayDataCluster"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "Identifier"  character varying(255),
  "Name"  character varying(255),
  "ClusterBioAssayData"  bigint,
  --
  CONSTRAINT "PK_BioAssayDataCluster" PRIMARY KEY ("WID")
);



CREATE TABLE "Node"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "BioAssayDataCluster_Nodes"  bigint,
  "Node_Nodes"  bigint,
  --
  CONSTRAINT "PK_Node" PRIMARY KEY ("WID")
);



CREATE TABLE "NodeContents"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "Node_NodeContents"  bigint,
  "BioAssayDimension"  bigint,
  "DesignElementDimension"  bigint,
  "QuantitationDimension"  bigint,
  --
  CONSTRAINT "PK_NodeContents" PRIMARY KEY ("WID")
);



CREATE TABLE "NodeValue"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "Node_NodeValue"  bigint,
  "Name"  character varying(255),
  "Value"  character varying(255),
  "NodeValue_Type"  bigint,
  "NodeValue_Scale"  bigint,
  "NodeValue_DataType"  bigint,
  --
  CONSTRAINT "PK_NodeValue" PRIMARY KEY ("WID")
);



CREATE TABLE "Measurement"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "Type_"  character varying(25),
  "Value"  character varying(255),
  "KindCV"  character varying(25),
  "OtherKind"  character varying(255),
  "Measurement_Unit"  bigint,
  --
  CONSTRAINT "PK_Measurement" PRIMARY KEY ("WID")
);



CREATE TABLE "Unit"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "MAGEClass"  character varying(100)  NOT NULL,
  "UnitName"  character varying(255),
  "UnitNameCV"  character varying(25),
  "UnitNameCV2"  character varying(25),
  "UnitNameCV3"  character varying(25),
  "UnitNameCV4"  character varying(25),
  "UnitNameCV5"  character varying(25),
  "UnitNameCV6"  character varying(25),
  "UnitNameCV7"  character varying(25),
  --
  CONSTRAINT "PK_Unit" PRIMARY KEY ("WID")
);



CREATE TABLE "Parameter"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "Identifier"  character varying(255),
  "Name"  character varying(255),
  "Parameter_DefaultValue"  bigint,
  "Parameter_DataType"  bigint,
  "Parameterizable_ParameterTypes"  bigint,
  --
  CONSTRAINT "PK_Parameter" PRIMARY KEY ("WID")
);



CREATE TABLE "ParameterValue"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "Value"  character varying(255),
  "ParameterType"  bigint,
  "ParameterizableApplication"  bigint,
  --
  CONSTRAINT "PK_ParameterValue" PRIMARY KEY ("WID")
);



CREATE TABLE "Parameterizable"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "MAGEClass"  character varying(100)  NOT NULL,
  "Identifier"  character varying(255),
  "Name"  character varying(255),
  "URI"  character varying(255),
  "Model"  character varying(255),
  "Make"  character varying(255),
  "Hardware_Type"  bigint,
  "Text"  character varying(1000),
  "Title"  character varying(255),
  "Protocol_Type"  bigint,
  "Software_Type"  bigint,
  "Hardware"  bigint,
  --
  CONSTRAINT "PK_Parameterizable" PRIMARY KEY ("WID")
);



CREATE TABLE "ParameterizableApplication"
(
  "WID"  bigint  NOT NULL,
  "DataSetWID"  bigint  NOT NULL,
  "MAGEClass"  character varying(100)  NOT NULL,
  "ArrayDesign"  bigint,
  "ArrayManufacture"  bigint,
  "BioEvent_ProtocolApplications"  bigint,
  "SerialNumber"  character varying(255),
  "Hardware"  bigint,
  "ActivityDate"  character varying(255),
  "ProtocolApplication"  bigint,
  "ProtocolApplication2"  bigint,
  "Protocol"  bigint,
  "Version"  character varying(255),
  "ReleaseDate"  timestamp,
  "Software"  bigint,
  --
  CONSTRAINT "PK_ParameterizableApplication" PRIMARY KEY ("WID")
);



CREATE TABLE "ArrayDesignWIDReporterGroupWID"
(
  "ArrayDesignWID"  bigint  NOT NULL,
  "ReporterGroupWID"  bigint  NOT NULL
  --
);



CREATE TABLE "ArrayDesignWIDCompositeGrpWID"
(
  "ArrayDesignWID"  bigint  NOT NULL,
  "CompositeGroupWID"  bigint  NOT NULL
  --
);



CREATE TABLE "ArrayDesignWIDContactWID"
(
  "ArrayDesignWID"  bigint  NOT NULL,
  "ContactWID"  bigint  NOT NULL
  --
);



CREATE TABLE "ComposGrpWIDComposSequenceWID"
(
  "CompositeGroupWID"  bigint  NOT NULL,
  "CompositeSequenceWID"  bigint  NOT NULL
  --
);



CREATE TABLE "ReporterGroupWIDReporterWID"
(
  "ReporterGroupWID"  bigint  NOT NULL,
  "ReporterWID"  bigint  NOT NULL
  --
);



CREATE TABLE "ExperimentWIDContactWID"
(
  "ExperimentWID"  bigint  NOT NULL,
  "ContactWID"  bigint  NOT NULL
  --
);



CREATE TABLE "ExperimWIDBioAssayDataClustWID"
(
  "ExperimentWID"  bigint  NOT NULL,
  "BioAssayDataClusterWID"  bigint  NOT NULL
  --
);



CREATE TABLE "ExperimentWIDBioAssayDataWID"
(
  "ExperimentWID"  bigint  NOT NULL,
  "BioAssayDataWID"  bigint  NOT NULL
  --
);



CREATE TABLE "ExperimentWIDBioAssayWID"
(
  "ExperimentWID"  bigint  NOT NULL,
  "BioAssayWID"  bigint  NOT NULL
  --
);



CREATE TABLE "ExperimentDesignWIDBioAssayWID"
(
  "ExperimentDesignWID"  bigint  NOT NULL,
  "BioAssayWID"  bigint  NOT NULL
  --
);



CREATE TABLE "QuantTypeWIDConfidenceIndWID"
(
  "QuantitationTypeWID"  bigint  NOT NULL,
  "ConfidenceIndicatorWID"  bigint  NOT NULL
  --
);



CREATE TABLE "QuantTypeWIDQuantTypeMapWID"
(
  "QuantitationTypeWID"  bigint  NOT NULL,
  "QuantitationTypeMapWID"  bigint  NOT NULL
  --
);



CREATE TABLE "DatabaseWIDContactWID"
(
  "DatabaseWID"  bigint  NOT NULL,
  "ContactWID"  bigint  NOT NULL
  --
);



CREATE TABLE "ArrayGroupWIDArrayWID"
(
  "ArrayGroupWID"  bigint  NOT NULL,
  "ArrayWID"  bigint  NOT NULL
  --
);



CREATE TABLE "ArrayManufactureWIDArrayWID"
(
  "ArrayManufactureWID"  bigint  NOT NULL,
  "ArrayWID"  bigint  NOT NULL
  --
);



CREATE TABLE "ArrayManufactureWIDContactWID"
(
  "ArrayManufactureWID"  bigint  NOT NULL,
  "ContactWID"  bigint  NOT NULL
  --
);



CREATE TABLE "CompositeSeqWIDBioSeqWID"
(
  "CompositeSequenceWID"  bigint  NOT NULL,
  "BioSequenceWID"  bigint  NOT NULL
  --
);



CREATE TABLE "ComposSeqWIDRepoComposMapWID"
(
  "CompositeSequenceWID"  bigint  NOT NULL,
  "ReporterCompositeMapWID"  bigint  NOT NULL
  --
);



CREATE TABLE "ComposSeqWIDComposComposMapWID"
(
  "CompositeSequenceWID"  bigint  NOT NULL,
  "CompositeCompositeMapWID"  bigint  NOT NULL
  --
);



CREATE TABLE "FeatureWIDFeatureWID"
(
  "FeatureWID1"  bigint  NOT NULL,
  "FeatureWID2"  bigint  NOT NULL
  --
);



CREATE TABLE "FeatureWIDFeatureWID2"
(
  "FeatureWID1"  bigint  NOT NULL,
  "FeatureWID2"  bigint  NOT NULL
  --
);



CREATE TABLE "ReporterWIDBioSequenceWID"
(
  "ReporterWID"  bigint  NOT NULL,
  "BioSequenceWID"  bigint  NOT NULL
  --
);



CREATE TABLE "ReporterWIDFeatureReporMapWID"
(
  "ReporterWID"  bigint  NOT NULL,
  "FeatureReporterMapWID"  bigint  NOT NULL
  --
);



CREATE TABLE "BioAssayDimensioWIDBioAssayWID"
(
  "BioAssayDimensionWID"  bigint  NOT NULL,
  "BioAssayWID"  bigint  NOT NULL
  --
);



CREATE TABLE "BioAssayMapWIDBioAssayWID"
(
  "BioAssayMapWID"  bigint  NOT NULL,
  "BioAssayWID"  bigint  NOT NULL
  --
);



CREATE TABLE "BAssayMappingWIDBAssayMapWID"
(
  "BioAssayMappingWID"  bigint  NOT NULL,
  "BioAssayMapWID"  bigint  NOT NULL
  --
);



CREATE TABLE "ComposSeqDimensWIDComposSeqWID"
(
  "CompositeSequenceDimensionWID"  bigint  NOT NULL,
  "CompositeSequenceWID"  bigint  NOT NULL
  --
);



CREATE TABLE "DesnElMappingWIDDesnElMapWID"
(
  "DesignElementMappingWID"  bigint  NOT NULL,
  "DesignElementMapWID"  bigint  NOT NULL
  --
);



CREATE TABLE "FeatureDimensionWIDFeatureWID"
(
  "FeatureDimensionWID"  bigint  NOT NULL,
  "FeatureWID"  bigint  NOT NULL
  --
);



CREATE TABLE "QuantTypeDimensWIDQuantTypeWID"
(
  "QuantitationTypeDimensionWID"  bigint  NOT NULL,
  "QuantitationTypeWID"  bigint  NOT NULL
  --
);



CREATE TABLE "QuantTypeMapWIDQuantTypeWID"
(
  "QuantitationTypeMapWID"  bigint  NOT NULL,
  "QuantitationTypeWID"  bigint  NOT NULL
  --
);



CREATE TABLE "QuantTyMapWIDQuantTyMapWI"
(
  "QuantitationTypeMappingWID"  bigint  NOT NULL,
  "QuantitationTypeMapWID"  bigint  NOT NULL
  --
);



CREATE TABLE "ReporterDimensWIDReporterWID"
(
  "ReporterDimensionWID"  bigint  NOT NULL,
  "ReporterWID"  bigint  NOT NULL
  --
);



CREATE TABLE "TransformWIDBioAssayDataWID"
(
  "TransformationWID"  bigint  NOT NULL,
  "BioAssayDataWID"  bigint  NOT NULL
  --
);



CREATE TABLE "BioSourceWIDContactWID"
(
  "BioSourceWID"  bigint  NOT NULL,
  "ContactWID"  bigint  NOT NULL
  --
);



CREATE TABLE "LabeledExtractWIDCompoundWID"
(
  "LabeledExtractWID"  bigint  NOT NULL,
  "CompoundWID"  bigint  NOT NULL
  --
);



CREATE TABLE "BioAssayWIDChannelWID"
(
  "BioAssayWID"  bigint  NOT NULL,
  "ChannelWID"  bigint  NOT NULL
  --
);



CREATE TABLE "BioAssayWIDFactorValueWID"
(
  "BioAssayWID"  bigint  NOT NULL,
  "FactorValueWID"  bigint  NOT NULL
  --
);



CREATE TABLE "ChannelWIDCompoundWID"
(
  "ChannelWID"  bigint  NOT NULL,
  "CompoundWID"  bigint  NOT NULL
  --
);



CREATE TABLE "DerivBioAWIDDerivBioADataWID"
(
  "DerivedBioAssayWID"  bigint  NOT NULL,
  "DerivedBioAssayDataWID"  bigint  NOT NULL
  --
);



CREATE TABLE "DerivBioAssayWIDBioAssayMapWID"
(
  "DerivedBioAssayWID"  bigint  NOT NULL,
  "BioAssayMapWID"  bigint  NOT NULL
  --
);



CREATE TABLE "ImageWIDChannelWID"
(
  "ImageWID"  bigint  NOT NULL,
  "ChannelWID"  bigint  NOT NULL
  --
);



CREATE TABLE "ImageAcquisitionWIDImageWID"
(
  "ImageAcquisitionWID"  bigint  NOT NULL,
  "ImageWID"  bigint  NOT NULL
  --
);



CREATE TABLE "MeasBAssayWIDMeasBAssayDataWID"
(
  "MeasuredBioAssayWID"  bigint  NOT NULL,
  "MeasuredBioAssayDataWID"  bigint  NOT NULL
  --
);



CREATE TABLE "HardwareWIDSoftwareWID"
(
  "HardwareWID"  bigint  NOT NULL,
  "SoftwareWID"  bigint  NOT NULL
  --
);



CREATE TABLE "HardwareWIDContactWID"
(
  "HardwareWID"  bigint  NOT NULL,
  "ContactWID"  bigint  NOT NULL
  --
);



CREATE TABLE "ProtocolWIDHardwareWID"
(
  "ProtocolWID"  bigint  NOT NULL,
  "HardwareWID"  bigint  NOT NULL
  --
);



CREATE TABLE "ProtocolWIDSoftwareWID"
(
  "ProtocolWID"  bigint  NOT NULL,
  "SoftwareWID"  bigint  NOT NULL
  --
);



CREATE TABLE "ProtocolApplWIDPersonWID"
(
  "ProtocolApplicationWID"  bigint  NOT NULL,
  "PersonWID"  bigint  NOT NULL
  --
);



CREATE TABLE "SoftwareWIDSoftwareWID"
(
  "SoftwareWID1"  bigint  NOT NULL,
  "SoftwareWID2"  bigint  NOT NULL
  --
);



CREATE TABLE "SoftwareWIDContactWID"
(
  "SoftwareWID"  bigint  NOT NULL,
  "ContactWID"  bigint  NOT NULL
  --
);



ALTER TABLE "DataSetHierarchy"
ADD CONSTRAINT "FK_DataSetHierarchy1"
FOREIGN KEY  ("SuperWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "DataSetHierarchy"
ADD CONSTRAINT "FK_DataSetHierarchy2"
FOREIGN KEY  ("SubWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Entry"
ADD CONSTRAINT "FK_Entry"
FOREIGN KEY  ("DatasetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "GeneExpressionData"
ADD CONSTRAINT "FK_GEDBAV"
FOREIGN KEY  ("BioAssayValuesWID")
REFERENCES  "BioDataValues" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Valence"
ADD CONSTRAINT "FK_Valence"
FOREIGN KEY  ("OtherWID")
REFERENCES  "Element" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "LightSource"
ADD CONSTRAINT "FK_LightSource1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "FlowCytometrySample"
ADD CONSTRAINT "FK_FlowCytometrySample1"
FOREIGN KEY  ("BioSourceWID")
REFERENCES  "BioSource" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "FlowCytometrySample"
ADD CONSTRAINT "FK_FlowCytometrySample2"
FOREIGN KEY  ("FlowCytometryProbeWID")
REFERENCES  "FlowCytometryProbe" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "FlowCytometrySample"
ADD CONSTRAINT "FK_FlowCytometrySample3"
FOREIGN KEY  ("MeasurementWID")
REFERENCES  "Measurement" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "FlowCytometrySample"
ADD CONSTRAINT "FK_FlowCytometrySample4"
FOREIGN KEY  ("ManufacturerWID")
REFERENCES  "Contact" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "FlowCytometrySample"
ADD CONSTRAINT "FK_FlowCytometrySampleDS"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "FlowCytometryProbe"
ADD CONSTRAINT "FK_Probe1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "TranscriptionUnit"
ADD CONSTRAINT "FK_TranscriptionUnit1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "TranscriptionUnitComponent"
ADD CONSTRAINT "FK_TranscriptionUnitComponent1"
FOREIGN KEY  ("TranscriptionUnitWID")
REFERENCES  "TranscriptionUnit" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Chemical"
ADD CONSTRAINT "FK_Chemical1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Reaction"
ADD CONSTRAINT "FK_Reaction"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Interaction"
ADD CONSTRAINT "FK_Interaction1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Protein"
ADD CONSTRAINT "FK_Protein"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Feature"
ADD CONSTRAINT "FK_Feature"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Function"
ADD CONSTRAINT "FK_Function"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "EnzymaticReaction"
ADD CONSTRAINT "FK_EnzymaticReaction1"
FOREIGN KEY  ("ReactionWID")
REFERENCES  "Reaction" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "EnzymaticReaction"
ADD CONSTRAINT "FK_EnzymaticReaction2"
FOREIGN KEY  ("ProteinWID")
REFERENCES  "Protein" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "EnzymaticReaction"
ADD CONSTRAINT "FK_EnzymaticReaction3"
FOREIGN KEY  ("ComplexWID")
REFERENCES  "Protein" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "EnzymaticReaction"
ADD CONSTRAINT "FK_EnzymaticReaction4"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "GeneticCode"
ADD CONSTRAINT "FK_GeneticCode"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Division"
ADD CONSTRAINT "FK_Division"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Taxon"
ADD CONSTRAINT "FK_Taxon_Division"
FOREIGN KEY  ("DivisionWID")
REFERENCES  "Division" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Taxon"
ADD CONSTRAINT "FK_Taxon_GeneticCode"
FOREIGN KEY  ("GencodeWID")
REFERENCES  "GeneticCode" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Taxon"
ADD CONSTRAINT "FK_Taxon"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioSource"
ADD CONSTRAINT "FK_BioSource1"
FOREIGN KEY  ("TaxonWID")
REFERENCES  "Taxon" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioSource"
ADD CONSTRAINT "FK_BioSource2"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioSubtype"
ADD CONSTRAINT "FK_BioSubtype2"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "NucleicAcid"
ADD CONSTRAINT "FK_NucleicAcid1"
FOREIGN KEY  ("GeneticCodeWID")
REFERENCES  "GeneticCode" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "NucleicAcid"
ADD CONSTRAINT "FK_NucleicAcid2"
FOREIGN KEY  ("BioSourceWID")
REFERENCES  "BioSource" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "NucleicAcid"
ADD CONSTRAINT "FK_NucleicAcid3"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Subsequence"
ADD CONSTRAINT "FK_Subsequence1"
FOREIGN KEY  ("NucleicAcidWID")
REFERENCES  "NucleicAcid" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Subsequence"
ADD CONSTRAINT "FK_Subsequence2"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Gene"
ADD CONSTRAINT "FK_Gene1"
FOREIGN KEY  ("NucleicAcidWID")
REFERENCES  "NucleicAcid" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Gene"
ADD CONSTRAINT "FK_Gene2"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Pathway"
ADD CONSTRAINT "FK_Pathway1"
FOREIGN KEY  ("BioSourceWID")
REFERENCES  "BioSource" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Pathway"
ADD CONSTRAINT "FK_Pathway2"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Computation"
ADD CONSTRAINT "FK_Computation"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Citation"
ADD CONSTRAINT "FK_Citation"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Experiment"
ADD CONSTRAINT "FK_Experiment3"
FOREIGN KEY  ("ContactWID")
REFERENCES  "Contact" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Experiment"
ADD CONSTRAINT "FK_Experiment4"
FOREIGN KEY  ("ArchiveWID")
REFERENCES  "Archive" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Experiment"
ADD CONSTRAINT "FK_Experiment2"
FOREIGN KEY  ("GroupWID")
REFERENCES  "Experiment" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Experiment"
ADD CONSTRAINT "FK_Experiment5"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Experiment"
ADD CONSTRAINT "FK_Experiment6"
FOREIGN KEY  ("BioSourceWID")
REFERENCES  "BioSource" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ExperimentData"
ADD CONSTRAINT "FK_ExpData1"
FOREIGN KEY  ("ExperimentWID")
REFERENCES  "Experiment" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ExperimentData"
ADD CONSTRAINT "FK_ExpDataMD"
FOREIGN KEY  ("MageData")
REFERENCES  "ParameterValue" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ExperimentData"
ADD CONSTRAINT "FK_ExpData2"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ChemicalAtom"
ADD CONSTRAINT "FK_ChemicalAtom"
FOREIGN KEY  ("ChemicalWID")
REFERENCES  "Chemical" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ChemicalBond"
ADD CONSTRAINT "FK_ChemicalBond"
FOREIGN KEY  ("ChemicalWID")
REFERENCES  "Chemical" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "EnzReactionCofactor"
ADD CONSTRAINT "FK_EnzReactionCofactor1"
FOREIGN KEY  ("EnzymaticReactionWID")
REFERENCES  "EnzymaticReaction" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "EnzReactionCofactor"
ADD CONSTRAINT "FK_EnzReactionCofactor2"
FOREIGN KEY  ("ChemicalWID")
REFERENCES  "Chemical" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "EnzReactionAltCompound"
ADD CONSTRAINT "FK_ERAC1"
FOREIGN KEY  ("EnzymaticReactionWID")
REFERENCES  "EnzymaticReaction" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "EnzReactionAltCompound"
ADD CONSTRAINT "FK_ERAC2"
FOREIGN KEY  ("PrimaryWID")
REFERENCES  "Chemical" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "EnzReactionAltCompound"
ADD CONSTRAINT "FK_ERAC3"
FOREIGN KEY  ("AlternativeWID")
REFERENCES  "Chemical" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "EnzReactionInhibitorActivator"
ADD CONSTRAINT "FK_EnzReactionIA1"
FOREIGN KEY  ("EnzymaticReactionWID")
REFERENCES  "EnzymaticReaction" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Location"
ADD CONSTRAINT "FK_Location"
FOREIGN KEY  ("ProteinWID")
REFERENCES  "Protein" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "PathwayLink"
ADD CONSTRAINT "FK_PathwayLink1"
FOREIGN KEY  ("Pathway1WID")
REFERENCES  "Pathway" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "PathwayLink"
ADD CONSTRAINT "FK_PathwayLink2"
FOREIGN KEY  ("Pathway2WID")
REFERENCES  "Pathway" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "PathwayLink"
ADD CONSTRAINT "FK_PathwayLink3"
FOREIGN KEY  ("ChemicalWID")
REFERENCES  "Chemical" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "PathwayReaction"
ADD CONSTRAINT "FK_PathwayReaction1"
FOREIGN KEY  ("PathwayWID")
REFERENCES  "Pathway" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "PathwayReaction"
ADD CONSTRAINT "FK_PathwayReaction3"
FOREIGN KEY  ("PriorReactionWID")
REFERENCES  "Reaction" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Product"
ADD CONSTRAINT "FK_Product"
FOREIGN KEY  ("ReactionWID")
REFERENCES  "Reaction" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Reactant"
ADD CONSTRAINT "FK_Reactant"
FOREIGN KEY  ("ReactionWID")
REFERENCES  "Reaction" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "InteractionParticipant"
ADD CONSTRAINT "FK_InteractionParticipant1"
FOREIGN KEY  ("InteractionWID")
REFERENCES  "Interaction" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "SequenceMatch"
ADD CONSTRAINT "FK_SequenceMatch"
FOREIGN KEY  ("ComputationWID")
REFERENCES  "Computation" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Subunit"
ADD CONSTRAINT "FK_Subunit1"
FOREIGN KEY  ("ComplexWID")
REFERENCES  "Protein" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Subunit"
ADD CONSTRAINT "FK_Subunit2"
FOREIGN KEY  ("SubunitWID")
REFERENCES  "Protein" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "SuperPathway"
ADD CONSTRAINT "FK_SuperPathway1"
FOREIGN KEY  ("SubPathwayWID")
REFERENCES  "Pathway" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "SuperPathway"
ADD CONSTRAINT "FK_SuperPathway2"
FOREIGN KEY  ("SuperPathwayWID")
REFERENCES  "Pathway" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "TermRelationship"
ADD CONSTRAINT "FK_TermRelationship1"
FOREIGN KEY  ("TermWID")
REFERENCES  "Term" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "TermRelationship"
ADD CONSTRAINT "FK_TermRelationship2"
FOREIGN KEY  ("RelatedTermWID")
REFERENCES  "Term" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "RelatedTerm"
ADD CONSTRAINT "FK_RelatedTerm1"
FOREIGN KEY  ("TermWID")
REFERENCES  "Term" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "CitationWIDOtherWID"
ADD CONSTRAINT "FK_CitationWIDOtherWID"
FOREIGN KEY  ("CitationWID")
REFERENCES  "Citation" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioSourceWIDBioSubtypeWID"
ADD CONSTRAINT "FK_BioSourceWIDBioSubtypeWID1"
FOREIGN KEY  ("BioSourceWID")
REFERENCES  "BioSource" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioSourceWIDBioSubtypeWID"
ADD CONSTRAINT "FK_BioSourceWIDBioSubtypeWID2"
FOREIGN KEY  ("BioSubtypeWID")
REFERENCES  "BioSubtype" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioSourceWIDGeneWID"
ADD CONSTRAINT "FK_BioSourceWIDGeneWID1"
FOREIGN KEY  ("BioSourceWID")
REFERENCES  "BioSource" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioSourceWIDGeneWID"
ADD CONSTRAINT "FK_BioSourceWIDGeneWID2"
FOREIGN KEY  ("GeneWID")
REFERENCES  "Gene" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioSourceWIDProteinWID"
ADD CONSTRAINT "FK_BioSourceWIDProteinWID1"
FOREIGN KEY  ("BioSourceWID")
REFERENCES  "BioSource" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioSourceWIDProteinWID"
ADD CONSTRAINT "FK_BioSourceWIDProteinWID2"
FOREIGN KEY  ("ProteinWID")
REFERENCES  "Protein" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "GeneWIDNucleicAcidWID"
ADD CONSTRAINT "FK_GeneWIDNucleicAcidWID1"
FOREIGN KEY  ("GeneWID")
REFERENCES  "Gene" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "GeneWIDNucleicAcidWID"
ADD CONSTRAINT "FK_GeneWIDNucleicAcidWID2"
FOREIGN KEY  ("NucleicAcidWID")
REFERENCES  "NucleicAcid" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "GeneWIDProteinWID"
ADD CONSTRAINT "FK_GeneWIDProteinWID1"
FOREIGN KEY  ("GeneWID")
REFERENCES  "Gene" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "GeneWIDProteinWID"
ADD CONSTRAINT "FK_GeneWIDProteinWID2"
FOREIGN KEY  ("ProteinWID")
REFERENCES  "Protein" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ProteinWIDFunctionWID"
ADD CONSTRAINT "FK_ProteinWIDFunctionWID2"
FOREIGN KEY  ("ProteinWID")
REFERENCES  "Protein" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ProteinWIDFunctionWID"
ADD CONSTRAINT "FK_ProteinWIDFunctionWID3"
FOREIGN KEY  ("FunctionWID")
REFERENCES  "Function" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ExperimentRelationship"
ADD CONSTRAINT "FK_ExpRelationship1"
FOREIGN KEY  ("ExperimentWID")
REFERENCES  "Experiment" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ExperimentRelationship"
ADD CONSTRAINT "FK_ExpRelationship2"
FOREIGN KEY  ("RelatedExperimentWID")
REFERENCES  "Experiment" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "GelLocation"
ADD CONSTRAINT "FK_GelLocSpotWid"
FOREIGN KEY  ("SpotWID")
REFERENCES  "Spot" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "GelLocation"
ADD CONSTRAINT "FK_GelLocExp"
FOREIGN KEY  ("ExperimentWID")
REFERENCES  "Experiment" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "GelLocation"
ADD CONSTRAINT "FK_GelLocDataset"
FOREIGN KEY  ("DatasetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ProteinWIDSpotWID"
ADD CONSTRAINT "FK_ProteinWIDSpotWID1"
FOREIGN KEY  ("ProteinWID")
REFERENCES  "Protein" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ProteinWIDSpotWID"
ADD CONSTRAINT "FK_ProteinWIDSpotWID2"
FOREIGN KEY  ("SpotWID")
REFERENCES  "Spot" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Spot"
ADD CONSTRAINT "FK_Spot"
FOREIGN KEY  ("DatasetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "SpotIdMethod"
ADD CONSTRAINT "FK_SpotIdMethDataset"
FOREIGN KEY  ("DatasetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "SpotWIDSpotIdMethodWID"
ADD CONSTRAINT "FK_SpotWIDMethWID1"
FOREIGN KEY  ("SpotWID")
REFERENCES  "Spot" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "SpotWIDSpotIdMethodWID"
ADD CONSTRAINT "FK_SpotWIDMethWID2"
FOREIGN KEY  ("SpotIdMethodWID")
REFERENCES  "SpotIdMethod" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "NameValueType"
ADD CONSTRAINT "FK_NameValueType1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "NameValueType"
ADD CONSTRAINT "FK_NameValueType66"
FOREIGN KEY  ("NameValueType_PropertySets")
REFERENCES  "NameValueType" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Contact"
ADD CONSTRAINT "FK_Contact1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Contact"
ADD CONSTRAINT "FK_Contact3"
FOREIGN KEY  ("Parent")
REFERENCES  "Contact" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Contact"
ADD CONSTRAINT "FK_Contact4"
FOREIGN KEY  ("Affiliation")
REFERENCES  "Contact" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ArrayDesign"
ADD CONSTRAINT "FK_ArrayDesign1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ArrayDesign"
ADD CONSTRAINT "FK_ArrayDesign3"
FOREIGN KEY  ("SurfaceType")
REFERENCES  "Term" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "DesignElementGroup"
ADD CONSTRAINT "FK_DesignElementGroup1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "DesignElementGroup"
ADD CONSTRAINT "FK_DesignElementGroup3"
FOREIGN KEY  ("ArrayDesign_FeatureGroups")
REFERENCES  "ArrayDesign" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "DesignElementGroup"
ADD CONSTRAINT "FK_DesignElementGroup4"
FOREIGN KEY  ("DesignElementGroup_Species")
REFERENCES  "Term" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "DesignElementGroup"
ADD CONSTRAINT "FK_DesignElementGroup5"
FOREIGN KEY  ("FeatureGroup_TechnologyType")
REFERENCES  "Term" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "DesignElementGroup"
ADD CONSTRAINT "FK_DesignElementGroup6"
FOREIGN KEY  ("FeatureGroup_FeatureShape")
REFERENCES  "Term" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "DesignElementGroup"
ADD CONSTRAINT "FK_DesignElementGroup7"
FOREIGN KEY  ("FeatureGroup_DistanceUnit")
REFERENCES  "Unit" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Zone"
ADD CONSTRAINT "FK_Zone1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Zone"
ADD CONSTRAINT "FK_Zone3"
FOREIGN KEY  ("Zone_DistanceUnit")
REFERENCES  "Unit" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Zone"
ADD CONSTRAINT "FK_Zone4"
FOREIGN KEY  ("ZoneGroup_ZoneLocations")
REFERENCES  "ZoneGroup" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ZoneGroup"
ADD CONSTRAINT "FK_ZoneGroup1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ZoneGroup"
ADD CONSTRAINT "FK_ZoneGroup2"
FOREIGN KEY  ("PhysicalArrayDesign_ZoneGroups")
REFERENCES  "ArrayDesign" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ZoneGroup"
ADD CONSTRAINT "FK_ZoneGroup3"
FOREIGN KEY  ("ZoneGroup_DistanceUnit")
REFERENCES  "Unit" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ZoneGroup"
ADD CONSTRAINT "FK_ZoneGroup4"
FOREIGN KEY  ("ZoneGroup_ZoneLayout")
REFERENCES  "ZoneLayout" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ZoneLayout"
ADD CONSTRAINT "FK_ZoneLayout1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ZoneLayout"
ADD CONSTRAINT "FK_ZoneLayout2"
FOREIGN KEY  ("ZoneLayout_DistanceUnit")
REFERENCES  "Unit" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ExperimentDesign"
ADD CONSTRAINT "FK_ExperimentDesign1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ExperimentDesign"
ADD CONSTRAINT "FK_ExperimentDesign3"
FOREIGN KEY  ("Experiment_ExperimentDesigns")
REFERENCES  "Experiment" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ExperimentalFactor"
ADD CONSTRAINT "FK_ExperimentalFactor1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ExperimentalFactor"
ADD CONSTRAINT "FK_ExperimentalFactor3"
FOREIGN KEY  ("ExperimentDesign")
REFERENCES  "ExperimentDesign" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ExperimentalFactor"
ADD CONSTRAINT "FK_ExperimentalFactor4"
FOREIGN KEY  ("ExperimentalFactor_Category")
REFERENCES  "Term" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "FactorValue"
ADD CONSTRAINT "FK_FactorValue1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "FactorValue"
ADD CONSTRAINT "FK_FactorValue3"
FOREIGN KEY  ("ExperimentalFactor")
REFERENCES  "ExperimentalFactor" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "FactorValue"
ADD CONSTRAINT "FK_FactorValue4"
FOREIGN KEY  ("ExperimentalFactor2")
REFERENCES  "ExperimentalFactor" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "FactorValue"
ADD CONSTRAINT "FK_FactorValue5"
FOREIGN KEY  ("FactorValue_Measurement")
REFERENCES  "Measurement" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "FactorValue"
ADD CONSTRAINT "FK_FactorValue6"
FOREIGN KEY  ("FactorValue_Value")
REFERENCES  "Term" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "QuantitationType"
ADD CONSTRAINT "FK_QuantitationType1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "QuantitationType"
ADD CONSTRAINT "FK_QuantitationType3"
FOREIGN KEY  ("Channel")
REFERENCES  "Channel" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "QuantitationType"
ADD CONSTRAINT "FK_QuantitationType4"
FOREIGN KEY  ("QuantitationType_Scale")
REFERENCES  "Term" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "QuantitationType"
ADD CONSTRAINT "FK_QuantitationType5"
FOREIGN KEY  ("QuantitationType_DataType")
REFERENCES  "Term" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "QuantitationType"
ADD CONSTRAINT "FK_QuantitationType6"
FOREIGN KEY  ("TargetQuantitationType")
REFERENCES  "QuantitationType" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Database_"
ADD CONSTRAINT "FK_Database_1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Array_"
ADD CONSTRAINT "FK_Array_1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Array_"
ADD CONSTRAINT "FK_Array_3"
FOREIGN KEY  ("ArrayDesign")
REFERENCES  "ArrayDesign" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Array_"
ADD CONSTRAINT "FK_Array_4"
FOREIGN KEY  ("Information")
REFERENCES  "ArrayManufacture" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Array_"
ADD CONSTRAINT "FK_Array_5"
FOREIGN KEY  ("ArrayGroup")
REFERENCES  "ArrayGroup" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ArrayGroup"
ADD CONSTRAINT "FK_ArrayGroup1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ArrayGroup"
ADD CONSTRAINT "FK_ArrayGroup3"
FOREIGN KEY  ("ArrayGroup_SubstrateType")
REFERENCES  "Term" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ArrayGroup"
ADD CONSTRAINT "FK_ArrayGroup4"
FOREIGN KEY  ("ArrayGroup_DistanceUnit")
REFERENCES  "Unit" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ArrayManufacture"
ADD CONSTRAINT "FK_ArrayManufacture1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ArrayManufactureDeviation"
ADD CONSTRAINT "FK_ArrayManufactureDeviation1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ArrayManufactureDeviation"
ADD CONSTRAINT "FK_ArrayManufactureDeviation2"
FOREIGN KEY  ("Array_")
REFERENCES  "Array_" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "FeatureDefect"
ADD CONSTRAINT "FK_FeatureDefect1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "FeatureDefect"
ADD CONSTRAINT "FK_FeatureDefect2"
FOREIGN KEY  ("ArrayManufactureDeviation")
REFERENCES  "ArrayManufactureDeviation" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "FeatureDefect"
ADD CONSTRAINT "FK_FeatureDefect3"
FOREIGN KEY  ("FeatureDefect_DefectType")
REFERENCES  "Term" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "FeatureDefect"
ADD CONSTRAINT "FK_FeatureDefect4"
FOREIGN KEY  ("FeatureDefect_PositionDelta")
REFERENCES  "PositionDelta" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "FeatureDefect"
ADD CONSTRAINT "FK_FeatureDefect5"
FOREIGN KEY  ("Feature")
REFERENCES  "DesignElement" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Fiducial"
ADD CONSTRAINT "FK_Fiducial1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Fiducial"
ADD CONSTRAINT "FK_Fiducial3"
FOREIGN KEY  ("ArrayGroup_Fiducials")
REFERENCES  "ArrayGroup" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Fiducial"
ADD CONSTRAINT "FK_Fiducial4"
FOREIGN KEY  ("Fiducial_FiducialType")
REFERENCES  "Term" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Fiducial"
ADD CONSTRAINT "FK_Fiducial5"
FOREIGN KEY  ("Fiducial_DistanceUnit")
REFERENCES  "Unit" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Fiducial"
ADD CONSTRAINT "FK_Fiducial6"
FOREIGN KEY  ("Fiducial_Position")
REFERENCES  "Position_" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ManufactureLIMS"
ADD CONSTRAINT "FK_ManufactureLIMS1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ManufactureLIMS"
ADD CONSTRAINT "FK_ManufactureLIMS3"
FOREIGN KEY  ("ArrayManufacture_FeatureLIMSs")
REFERENCES  "ArrayManufacture" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ManufactureLIMS"
ADD CONSTRAINT "FK_ManufactureLIMS4"
FOREIGN KEY  ("Feature")
REFERENCES  "DesignElement" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ManufactureLIMS"
ADD CONSTRAINT "FK_ManufactureLIMS5"
FOREIGN KEY  ("BioMaterial")
REFERENCES  "BioSource" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "PositionDelta"
ADD CONSTRAINT "FK_PositionDelta1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "PositionDelta"
ADD CONSTRAINT "FK_PositionDelta2"
FOREIGN KEY  ("PositionDelta_DistanceUnit")
REFERENCES  "Unit" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ZoneDefect"
ADD CONSTRAINT "FK_ZoneDefect1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ZoneDefect"
ADD CONSTRAINT "FK_ZoneDefect2"
FOREIGN KEY  ("ArrayManufactureDeviation")
REFERENCES  "ArrayManufactureDeviation" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ZoneDefect"
ADD CONSTRAINT "FK_ZoneDefect3"
FOREIGN KEY  ("ZoneDefect_DefectType")
REFERENCES  "Term" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ZoneDefect"
ADD CONSTRAINT "FK_ZoneDefect4"
FOREIGN KEY  ("ZoneDefect_PositionDelta")
REFERENCES  "PositionDelta" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ZoneDefect"
ADD CONSTRAINT "FK_ZoneDefect5"
FOREIGN KEY  ("Zone")
REFERENCES  "Zone" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "SeqFeatureLocation"
ADD CONSTRAINT "FK_SeqFeatureLocation1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "SeqFeatureLocation"
ADD CONSTRAINT "FK_SeqFeatureLocation2"
FOREIGN KEY  ("SeqFeature_Regions")
REFERENCES  "Feature" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "SeqFeatureLocation"
ADD CONSTRAINT "FK_SeqFeatureLocation3"
FOREIGN KEY  ("SeqFeatureLocation_Subregions")
REFERENCES  "SeqFeatureLocation" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "SeqFeatureLocation"
ADD CONSTRAINT "FK_SeqFeatureLocation4"
FOREIGN KEY  ("SeqFeatureLocation_Coordinate")
REFERENCES  "SequencePosition" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "SequencePosition"
ADD CONSTRAINT "FK_SequencePosition1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "SequencePosition"
ADD CONSTRAINT "FK_SequencePosition2"
FOREIGN KEY  ("CompositeCompositeMap")
REFERENCES  "BioEvent" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "SequencePosition"
ADD CONSTRAINT "FK_SequencePosition3"
FOREIGN KEY  ("Composite")
REFERENCES  "DesignElement" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "SequencePosition"
ADD CONSTRAINT "FK_SequencePosition4"
FOREIGN KEY  ("ReporterCompositeMap")
REFERENCES  "BioEvent" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "SequencePosition"
ADD CONSTRAINT "FK_SequencePosition5"
FOREIGN KEY  ("Reporter")
REFERENCES  "DesignElement" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "DesignElement"
ADD CONSTRAINT "FK_DesignElement1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "DesignElement"
ADD CONSTRAINT "FK_DesignElement3"
FOREIGN KEY  ("FeatureGroup_Features")
REFERENCES  "DesignElementGroup" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "DesignElement"
ADD CONSTRAINT "FK_DesignElement4"
FOREIGN KEY  ("DesignElement_ControlType")
REFERENCES  "Term" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "DesignElement"
ADD CONSTRAINT "FK_DesignElement5"
FOREIGN KEY  ("Feature_Position")
REFERENCES  "Position_" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "DesignElement"
ADD CONSTRAINT "FK_DesignElement6"
FOREIGN KEY  ("Zone")
REFERENCES  "Zone" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "DesignElement"
ADD CONSTRAINT "FK_DesignElement7"
FOREIGN KEY  ("Feature_FeatureLocation")
REFERENCES  "FeatureLocation" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "DesignElement"
ADD CONSTRAINT "FK_DesignElement8"
FOREIGN KEY  ("FeatureGroup")
REFERENCES  "DesignElementGroup" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "DesignElement"
ADD CONSTRAINT "FK_DesignElement9"
FOREIGN KEY  ("Reporter_WarningType")
REFERENCES  "Term" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "FeatureInformation"
ADD CONSTRAINT "FK_FeatureInformation1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "FeatureInformation"
ADD CONSTRAINT "FK_FeatureInformation2"
FOREIGN KEY  ("Feature")
REFERENCES  "DesignElement" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "FeatureInformation"
ADD CONSTRAINT "FK_FeatureInformation3"
FOREIGN KEY  ("FeatureReporterMap")
REFERENCES  "BioEvent" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "FeatureLocation"
ADD CONSTRAINT "FK_FeatureLocation1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "MismatchInformation"
ADD CONSTRAINT "FK_MismatchInformation1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "MismatchInformation"
ADD CONSTRAINT "FK_MismatchInformation2"
FOREIGN KEY  ("CompositePosition")
REFERENCES  "SequencePosition" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "MismatchInformation"
ADD CONSTRAINT "FK_MismatchInformation3"
FOREIGN KEY  ("FeatureInformation")
REFERENCES  "FeatureInformation" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "MismatchInformation"
ADD CONSTRAINT "FK_MismatchInformation4"
FOREIGN KEY  ("ReporterPosition")
REFERENCES  "SequencePosition" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Position_"
ADD CONSTRAINT "FK_Position_1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Position_"
ADD CONSTRAINT "FK_Position_2"
FOREIGN KEY  ("Position_DistanceUnit")
REFERENCES  "Unit" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioEvent"
ADD CONSTRAINT "FK_BioEvent1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioEvent"
ADD CONSTRAINT "FK_BioEvent3"
FOREIGN KEY  ("CompositeSequence")
REFERENCES  "DesignElement" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioEvent"
ADD CONSTRAINT "FK_BioEvent4"
FOREIGN KEY  ("Reporter")
REFERENCES  "DesignElement" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioEvent"
ADD CONSTRAINT "FK_BioEvent5"
FOREIGN KEY  ("CompositeSequence2")
REFERENCES  "DesignElement" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioEvent"
ADD CONSTRAINT "FK_BioEvent6"
FOREIGN KEY  ("BioAssayMapTarget")
REFERENCES  "BioAssay" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioEvent"
ADD CONSTRAINT "FK_BioEvent7"
FOREIGN KEY  ("TargetQuantitationType")
REFERENCES  "QuantitationType" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioEvent"
ADD CONSTRAINT "FK_BioEvent8"
FOREIGN KEY  ("DerivedBioAssayDataTarget")
REFERENCES  "BioAssayData" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioEvent"
ADD CONSTRAINT "FK_BioEvent9"
FOREIGN KEY  ("QuantitationTypeMapping")
REFERENCES  "QuantitationTypeMapping" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioEvent"
ADD CONSTRAINT "FK_BioEvent10"
FOREIGN KEY  ("DesignElementMapping")
REFERENCES  "DesignElementMapping" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioEvent"
ADD CONSTRAINT "FK_BioEvent11"
FOREIGN KEY  ("Transformation_BioAssayMapping")
REFERENCES  "BioAssayMapping" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioEvent"
ADD CONSTRAINT "FK_BioEvent12"
FOREIGN KEY  ("BioMaterial_Treatments")
REFERENCES  "BioSource" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioEvent"
ADD CONSTRAINT "FK_BioEvent13"
FOREIGN KEY  ("Treatment_Action")
REFERENCES  "Term" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioEvent"
ADD CONSTRAINT "FK_BioEvent14"
FOREIGN KEY  ("Treatment_ActionMeasurement")
REFERENCES  "Measurement" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioEvent"
ADD CONSTRAINT "FK_BioEvent15"
FOREIGN KEY  ("Array_")
REFERENCES  "Array_" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioEvent"
ADD CONSTRAINT "FK_BioEvent16"
FOREIGN KEY  ("PhysicalBioAssayTarget")
REFERENCES  "BioAssay" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioEvent"
ADD CONSTRAINT "FK_BioEvent17"
FOREIGN KEY  ("PhysicalBioAssay")
REFERENCES  "BioAssay" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioEvent"
ADD CONSTRAINT "FK_BioEvent18"
FOREIGN KEY  ("Target")
REFERENCES  "BioAssay" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioEvent"
ADD CONSTRAINT "FK_BioEvent19"
FOREIGN KEY  ("PhysicalBioAssaySource")
REFERENCES  "BioAssay" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioEvent"
ADD CONSTRAINT "FK_BioEvent20"
FOREIGN KEY  ("MeasuredBioAssayTarget")
REFERENCES  "BioAssay" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioEvent"
ADD CONSTRAINT "FK_BioEvent21"
FOREIGN KEY  ("PhysicalBioAssay2")
REFERENCES  "BioAssay" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioAssayData"
ADD CONSTRAINT "FK_BioAssayData1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioAssayData"
ADD CONSTRAINT "FK_BioAssayData3"
FOREIGN KEY  ("BioAssayDimension")
REFERENCES  "BioAssayDimension" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioAssayData"
ADD CONSTRAINT "FK_BioAssayData4"
FOREIGN KEY  ("DesignElementDimension")
REFERENCES  "DesignElementDimension" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioAssayData"
ADD CONSTRAINT "FK_BioAssayData5"
FOREIGN KEY  ("QuantitationTypeDimension")
REFERENCES  "QuantitationTypeDimension" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioAssayData"
ADD CONSTRAINT "FK_BioAssayData6"
FOREIGN KEY  ("BioAssayData_BioDataValues")
REFERENCES  "BioDataValues" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioAssayData"
ADD CONSTRAINT "FK_BioAssayData7"
FOREIGN KEY  ("ProducerTransformation")
REFERENCES  "BioEvent" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioAssayDimension"
ADD CONSTRAINT "FK_BioAssayDimension1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioAssayMapping"
ADD CONSTRAINT "FK_BioAssayMapping1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioAssayTuple"
ADD CONSTRAINT "FK_BioAssayTuple1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioAssayTuple"
ADD CONSTRAINT "FK_BioAssayTuple2"
FOREIGN KEY  ("BioAssay")
REFERENCES  "BioAssay" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioAssayTuple"
ADD CONSTRAINT "FK_BioAssayTuple3"
FOREIGN KEY  ("BioDataTuples_BioAssayTuples")
REFERENCES  "BioDataValues" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioDataValues"
ADD CONSTRAINT "FK_BioDataValues1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioDataValues"
ADD CONSTRAINT "FK_BioDataValues2"
FOREIGN KEY  ("BioDataCube_DataInternal")
REFERENCES  "DataInternal" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioDataValues"
ADD CONSTRAINT "FK_BioDataValues3"
FOREIGN KEY  ("BioDataCube_DataExternal")
REFERENCES  "DataExternal" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "DataExternal"
ADD CONSTRAINT "FK_DataExternal1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "DataInternal"
ADD CONSTRAINT "FK_DataInternal1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Datum"
ADD CONSTRAINT "FK_Datum1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "DesignElementDimension"
ADD CONSTRAINT "FK_DesignElementDimension1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "DesignElementMapping"
ADD CONSTRAINT "FK_DesignElementMapping1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "DesignElementTuple"
ADD CONSTRAINT "FK_DesignElementTuple1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "DesignElementTuple"
ADD CONSTRAINT "FK_DesignElementTuple2"
FOREIGN KEY  ("BioAssayTuple")
REFERENCES  "BioAssayTuple" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "DesignElementTuple"
ADD CONSTRAINT "FK_DesignElementTuple3"
FOREIGN KEY  ("DesignElement")
REFERENCES  "DesignElement" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "QuantitationTypeDimension"
ADD CONSTRAINT "FK_QuantitationTypeDimension1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "QuantitationTypeMapping"
ADD CONSTRAINT "FK_QuantitationTypeMapping1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "QuantitationTypeTuple"
ADD CONSTRAINT "FK_QuantitationTypeTuple1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "QuantitationTypeTuple"
ADD CONSTRAINT "FK_QuantitationTypeTuple2"
FOREIGN KEY  ("DesignElementTuple")
REFERENCES  "DesignElementTuple" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "QuantitationTypeTuple"
ADD CONSTRAINT "FK_QuantitationTypeTuple3"
FOREIGN KEY  ("QuantitationType")
REFERENCES  "QuantitationType" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "QuantitationTypeTuple"
ADD CONSTRAINT "FK_QuantitationTypeTuple4"
FOREIGN KEY  ("QuantitationTypeTuple_Datum")
REFERENCES  "Datum" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioMaterialMeasurement"
ADD CONSTRAINT "FK_BioMaterialMeasurement1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioMaterialMeasurement"
ADD CONSTRAINT "FK_BioMaterialMeasurement2"
FOREIGN KEY  ("BioMaterial")
REFERENCES  "BioSource" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioMaterialMeasurement"
ADD CONSTRAINT "FK_BioMaterialMeasurement3"
FOREIGN KEY  ("Measurement")
REFERENCES  "Measurement" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioMaterialMeasurement"
ADD CONSTRAINT "FK_BioMaterialMeasurement4"
FOREIGN KEY  ("Treatment")
REFERENCES  "BioEvent" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioMaterialMeasurement"
ADD CONSTRAINT "FK_BioMaterialMeasurement5"
FOREIGN KEY  ("BioAssayCreation")
REFERENCES  "BioEvent" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "CompoundMeasurement"
ADD CONSTRAINT "FK_CompoundMeasurement1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "CompoundMeasurement"
ADD CONSTRAINT "FK_CompoundMeasurement2"
FOREIGN KEY  ("Compound_ComponentCompounds")
REFERENCES  "Chemical" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "CompoundMeasurement"
ADD CONSTRAINT "FK_CompoundMeasurement3"
FOREIGN KEY  ("Compound")
REFERENCES  "Chemical" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "CompoundMeasurement"
ADD CONSTRAINT "FK_CompoundMeasurement4"
FOREIGN KEY  ("Measurement")
REFERENCES  "Measurement" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "CompoundMeasurement"
ADD CONSTRAINT "FK_CompoundMeasurement5"
FOREIGN KEY  ("Treatment_CompoundMeasurements")
REFERENCES  "BioEvent" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioAssay"
ADD CONSTRAINT "FK_BioAssay1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioAssay"
ADD CONSTRAINT "FK_BioAssay3"
FOREIGN KEY  ("DerivedBioAssay_Type")
REFERENCES  "Term" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioAssay"
ADD CONSTRAINT "FK_BioAssay4"
FOREIGN KEY  ("FeatureExtraction")
REFERENCES  "BioEvent" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioAssay"
ADD CONSTRAINT "FK_BioAssay5"
FOREIGN KEY  ("BioAssayCreation")
REFERENCES  "BioEvent" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Channel"
ADD CONSTRAINT "FK_Channel1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Image"
ADD CONSTRAINT "FK_Image1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Image"
ADD CONSTRAINT "FK_Image3"
FOREIGN KEY  ("Image_Format")
REFERENCES  "Term" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Image"
ADD CONSTRAINT "FK_Image4"
FOREIGN KEY  ("PhysicalBioAssay")
REFERENCES  "BioAssay" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioAssayDataCluster"
ADD CONSTRAINT "FK_BioAssayDataCluster1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioAssayDataCluster"
ADD CONSTRAINT "FK_BioAssayDataCluster3"
FOREIGN KEY  ("ClusterBioAssayData")
REFERENCES  "BioAssayData" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Node"
ADD CONSTRAINT "FK_Node1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Node"
ADD CONSTRAINT "FK_Node3"
FOREIGN KEY  ("BioAssayDataCluster_Nodes")
REFERENCES  "BioAssayDataCluster" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Node"
ADD CONSTRAINT "FK_Node4"
FOREIGN KEY  ("Node_Nodes")
REFERENCES  "Node" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "NodeContents"
ADD CONSTRAINT "FK_NodeContents1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "NodeContents"
ADD CONSTRAINT "FK_NodeContents3"
FOREIGN KEY  ("Node_NodeContents")
REFERENCES  "Node" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "NodeContents"
ADD CONSTRAINT "FK_NodeContents4"
FOREIGN KEY  ("BioAssayDimension")
REFERENCES  "BioAssayDimension" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "NodeContents"
ADD CONSTRAINT "FK_NodeContents5"
FOREIGN KEY  ("DesignElementDimension")
REFERENCES  "DesignElementDimension" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "NodeContents"
ADD CONSTRAINT "FK_NodeContents6"
FOREIGN KEY  ("QuantitationDimension")
REFERENCES  "QuantitationTypeDimension" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "NodeValue"
ADD CONSTRAINT "FK_NodeValue1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "NodeValue"
ADD CONSTRAINT "FK_NodeValue2"
FOREIGN KEY  ("Node_NodeValue")
REFERENCES  "Node" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "NodeValue"
ADD CONSTRAINT "FK_NodeValue3"
FOREIGN KEY  ("NodeValue_Type")
REFERENCES  "Term" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "NodeValue"
ADD CONSTRAINT "FK_NodeValue4"
FOREIGN KEY  ("NodeValue_Scale")
REFERENCES  "Term" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "NodeValue"
ADD CONSTRAINT "FK_NodeValue5"
FOREIGN KEY  ("NodeValue_DataType")
REFERENCES  "Term" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Measurement"
ADD CONSTRAINT "FK_Measurement1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Measurement"
ADD CONSTRAINT "FK_Measurement2"
FOREIGN KEY  ("Measurement_Unit")
REFERENCES  "Unit" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Unit"
ADD CONSTRAINT "FK_Unit1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Parameter"
ADD CONSTRAINT "FK_Parameter1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Parameter"
ADD CONSTRAINT "FK_Parameter3"
FOREIGN KEY  ("Parameter_DefaultValue")
REFERENCES  "Measurement" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Parameter"
ADD CONSTRAINT "FK_Parameter4"
FOREIGN KEY  ("Parameter_DataType")
REFERENCES  "Term" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Parameter"
ADD CONSTRAINT "FK_Parameter5"
FOREIGN KEY  ("Parameterizable_ParameterTypes")
REFERENCES  "Parameterizable" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ParameterValue"
ADD CONSTRAINT "FK_ParameterValue1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ParameterValue"
ADD CONSTRAINT "FK_ParameterValue2"
FOREIGN KEY  ("ParameterType")
REFERENCES  "Parameter" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ParameterValue"
ADD CONSTRAINT "FK_ParameterValue3"
FOREIGN KEY  ("ParameterizableApplication")
REFERENCES  "ParameterizableApplication" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Parameterizable"
ADD CONSTRAINT "FK_Parameterizable1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Parameterizable"
ADD CONSTRAINT "FK_Parameterizable3"
FOREIGN KEY  ("Hardware_Type")
REFERENCES  "Term" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Parameterizable"
ADD CONSTRAINT "FK_Parameterizable4"
FOREIGN KEY  ("Protocol_Type")
REFERENCES  "Term" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Parameterizable"
ADD CONSTRAINT "FK_Parameterizable5"
FOREIGN KEY  ("Software_Type")
REFERENCES  "Term" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "Parameterizable"
ADD CONSTRAINT "FK_Parameterizable6"
FOREIGN KEY  ("Hardware")
REFERENCES  "Parameterizable" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ParameterizableApplication"
ADD CONSTRAINT "FK_ParameterizableApplicatio1"
FOREIGN KEY  ("DataSetWID")
REFERENCES  "DataSet" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ParameterizableApplication"
ADD CONSTRAINT "FK_ParameterizableApplicatio3"
FOREIGN KEY  ("ArrayDesign")
REFERENCES  "ArrayDesign" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ParameterizableApplication"
ADD CONSTRAINT "FK_ParameterizableApplicatio4"
FOREIGN KEY  ("ArrayManufacture")
REFERENCES  "ArrayManufacture" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ParameterizableApplication"
ADD CONSTRAINT "FK_ParameterizableApplicatio5"
FOREIGN KEY  ("BioEvent_ProtocolApplications")
REFERENCES  "BioEvent" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ParameterizableApplication"
ADD CONSTRAINT "FK_ParameterizableApplicatio6"
FOREIGN KEY  ("Hardware")
REFERENCES  "Parameterizable" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ParameterizableApplication"
ADD CONSTRAINT "FK_ParameterizableApplicatio7"
FOREIGN KEY  ("ProtocolApplication")
REFERENCES  "ParameterizableApplication" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ParameterizableApplication"
ADD CONSTRAINT "FK_ParameterizableApplicatio8"
FOREIGN KEY  ("ProtocolApplication2")
REFERENCES  "ParameterizableApplication" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ParameterizableApplication"
ADD CONSTRAINT "FK_ParameterizableApplicatio9"
FOREIGN KEY  ("Protocol")
REFERENCES  "Parameterizable" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ParameterizableApplication"
ADD CONSTRAINT "FK_ParameterizableApplicatio10"
FOREIGN KEY  ("Software")
REFERENCES  "Parameterizable" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ArrayDesignWIDReporterGroupWID"
ADD CONSTRAINT "FK_ArrayDesignWIDReporterGro1"
FOREIGN KEY  ("ArrayDesignWID")
REFERENCES  "ArrayDesign" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ArrayDesignWIDReporterGroupWID"
ADD CONSTRAINT "FK_ArrayDesignWIDReporterGro2"
FOREIGN KEY  ("ReporterGroupWID")
REFERENCES  "DesignElementGroup" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ArrayDesignWIDCompositeGrpWID"
ADD CONSTRAINT "FK_ArrayDesignWIDCompositeGr1"
FOREIGN KEY  ("ArrayDesignWID")
REFERENCES  "ArrayDesign" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ArrayDesignWIDCompositeGrpWID"
ADD CONSTRAINT "FK_ArrayDesignWIDCompositeGr2"
FOREIGN KEY  ("CompositeGroupWID")
REFERENCES  "DesignElementGroup" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ArrayDesignWIDContactWID"
ADD CONSTRAINT "FK_ArrayDesignWIDContactWID1"
FOREIGN KEY  ("ArrayDesignWID")
REFERENCES  "ArrayDesign" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ArrayDesignWIDContactWID"
ADD CONSTRAINT "FK_ArrayDesignWIDContactWID2"
FOREIGN KEY  ("ContactWID")
REFERENCES  "Contact" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ComposGrpWIDComposSequenceWID"
ADD CONSTRAINT "FK_ComposGrpWIDComposSequenc1"
FOREIGN KEY  ("CompositeGroupWID")
REFERENCES  "DesignElementGroup" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ComposGrpWIDComposSequenceWID"
ADD CONSTRAINT "FK_ComposGrpWIDComposSequenc2"
FOREIGN KEY  ("CompositeSequenceWID")
REFERENCES  "DesignElement" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ReporterGroupWIDReporterWID"
ADD CONSTRAINT "FK_ReporterGroupWIDReporterW1"
FOREIGN KEY  ("ReporterGroupWID")
REFERENCES  "DesignElementGroup" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ReporterGroupWIDReporterWID"
ADD CONSTRAINT "FK_ReporterGroupWIDReporterW2"
FOREIGN KEY  ("ReporterWID")
REFERENCES  "DesignElement" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ExperimentWIDContactWID"
ADD CONSTRAINT "FK_ExperimentWIDContactWID1"
FOREIGN KEY  ("ExperimentWID")
REFERENCES  "Experiment" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ExperimentWIDContactWID"
ADD CONSTRAINT "FK_ExperimentWIDContactWID2"
FOREIGN KEY  ("ContactWID")
REFERENCES  "Contact" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ExperimWIDBioAssayDataClustWID"
ADD CONSTRAINT "FK_ExperimWIDBioAssayDataClu1"
FOREIGN KEY  ("ExperimentWID")
REFERENCES  "Experiment" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ExperimWIDBioAssayDataClustWID"
ADD CONSTRAINT "FK_ExperimWIDBioAssayDataClu2"
FOREIGN KEY  ("BioAssayDataClusterWID")
REFERENCES  "BioAssayDataCluster" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ExperimentWIDBioAssayDataWID"
ADD CONSTRAINT "FK_ExperimentWIDBioAssayData1"
FOREIGN KEY  ("ExperimentWID")
REFERENCES  "Experiment" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ExperimentWIDBioAssayDataWID"
ADD CONSTRAINT "FK_ExperimentWIDBioAssayData2"
FOREIGN KEY  ("BioAssayDataWID")
REFERENCES  "BioAssayData" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ExperimentWIDBioAssayWID"
ADD CONSTRAINT "FK_ExperimentWIDBioAssayWID1"
FOREIGN KEY  ("ExperimentWID")
REFERENCES  "Experiment" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ExperimentWIDBioAssayWID"
ADD CONSTRAINT "FK_ExperimentWIDBioAssayWID2"
FOREIGN KEY  ("BioAssayWID")
REFERENCES  "BioAssay" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ExperimentDesignWIDBioAssayWID"
ADD CONSTRAINT "FK_ExperimentDesignWIDBioAss1"
FOREIGN KEY  ("ExperimentDesignWID")
REFERENCES  "ExperimentDesign" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ExperimentDesignWIDBioAssayWID"
ADD CONSTRAINT "FK_ExperimentDesignWIDBioAss2"
FOREIGN KEY  ("BioAssayWID")
REFERENCES  "BioAssay" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "QuantTypeWIDConfidenceIndWID"
ADD CONSTRAINT "FK_QuantTypeWIDConfidenceInd1"
FOREIGN KEY  ("QuantitationTypeWID")
REFERENCES  "QuantitationType" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "QuantTypeWIDConfidenceIndWID"
ADD CONSTRAINT "FK_QuantTypeWIDConfidenceInd2"
FOREIGN KEY  ("ConfidenceIndicatorWID")
REFERENCES  "QuantitationType" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "QuantTypeWIDQuantTypeMapWID"
ADD CONSTRAINT "FK_QuantTypeWIDQuantTypeMapW1"
FOREIGN KEY  ("QuantitationTypeWID")
REFERENCES  "QuantitationType" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "QuantTypeWIDQuantTypeMapWID"
ADD CONSTRAINT "FK_QuantTypeWIDQuantTypeMapW2"
FOREIGN KEY  ("QuantitationTypeMapWID")
REFERENCES  "BioEvent" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "DatabaseWIDContactWID"
ADD CONSTRAINT "FK_DatabaseWIDContactWID1"
FOREIGN KEY  ("DatabaseWID")
REFERENCES  "Database_" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "DatabaseWIDContactWID"
ADD CONSTRAINT "FK_DatabaseWIDContactWID2"
FOREIGN KEY  ("ContactWID")
REFERENCES  "Contact" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ArrayGroupWIDArrayWID"
ADD CONSTRAINT "FK_ArrayGroupWIDArrayWID1"
FOREIGN KEY  ("ArrayGroupWID")
REFERENCES  "ArrayGroup" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ArrayGroupWIDArrayWID"
ADD CONSTRAINT "FK_ArrayGroupWIDArrayWID2"
FOREIGN KEY  ("ArrayWID")
REFERENCES  "Array_" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ArrayManufactureWIDArrayWID"
ADD CONSTRAINT "FK_ArrayManufactureWIDArrayW1"
FOREIGN KEY  ("ArrayManufactureWID")
REFERENCES  "ArrayManufacture" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ArrayManufactureWIDArrayWID"
ADD CONSTRAINT "FK_ArrayManufactureWIDArrayW2"
FOREIGN KEY  ("ArrayWID")
REFERENCES  "Array_" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ArrayManufactureWIDContactWID"
ADD CONSTRAINT "FK_ArrayManufactureWIDContac1"
FOREIGN KEY  ("ArrayManufactureWID")
REFERENCES  "ArrayManufacture" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ArrayManufactureWIDContactWID"
ADD CONSTRAINT "FK_ArrayManufactureWIDContac2"
FOREIGN KEY  ("ContactWID")
REFERENCES  "Contact" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "CompositeSeqWIDBioSeqWID"
ADD CONSTRAINT "FK_CompositeSeqWIDBioSeqWID1"
FOREIGN KEY  ("CompositeSequenceWID")
REFERENCES  "DesignElement" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ComposSeqWIDRepoComposMapWID"
ADD CONSTRAINT "FK_ComposSeqWIDRepoComposMap1"
FOREIGN KEY  ("CompositeSequenceWID")
REFERENCES  "DesignElement" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ComposSeqWIDRepoComposMapWID"
ADD CONSTRAINT "FK_ComposSeqWIDRepoComposMap2"
FOREIGN KEY  ("ReporterCompositeMapWID")
REFERENCES  "BioEvent" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ComposSeqWIDComposComposMapWID"
ADD CONSTRAINT "FK_ComposSeqWIDComposComposM1"
FOREIGN KEY  ("CompositeSequenceWID")
REFERENCES  "DesignElement" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ComposSeqWIDComposComposMapWID"
ADD CONSTRAINT "FK_ComposSeqWIDComposComposM2"
FOREIGN KEY  ("CompositeCompositeMapWID")
REFERENCES  "BioEvent" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "FeatureWIDFeatureWID"
ADD CONSTRAINT "FK_FeatureWIDFeatureWID1"
FOREIGN KEY  ("FeatureWID1")
REFERENCES  "DesignElement" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "FeatureWIDFeatureWID"
ADD CONSTRAINT "FK_FeatureWIDFeatureWID2"
FOREIGN KEY  ("FeatureWID2")
REFERENCES  "DesignElement" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "FeatureWIDFeatureWID2"
ADD CONSTRAINT "FK_FeatureWIDFeatureWID21"
FOREIGN KEY  ("FeatureWID1")
REFERENCES  "DesignElement" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "FeatureWIDFeatureWID2"
ADD CONSTRAINT "FK_FeatureWIDFeatureWID22"
FOREIGN KEY  ("FeatureWID2")
REFERENCES  "DesignElement" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ReporterWIDBioSequenceWID"
ADD CONSTRAINT "FK_ReporterWIDBioSequenceWID1"
FOREIGN KEY  ("ReporterWID")
REFERENCES  "DesignElement" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ReporterWIDFeatureReporMapWID"
ADD CONSTRAINT "FK_ReporterWIDFeatureReporMa1"
FOREIGN KEY  ("ReporterWID")
REFERENCES  "DesignElement" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ReporterWIDFeatureReporMapWID"
ADD CONSTRAINT "FK_ReporterWIDFeatureReporMa2"
FOREIGN KEY  ("FeatureReporterMapWID")
REFERENCES  "BioEvent" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioAssayDimensioWIDBioAssayWID"
ADD CONSTRAINT "FK_BioAssayDimensioWIDBioAss1"
FOREIGN KEY  ("BioAssayDimensionWID")
REFERENCES  "BioAssayDimension" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioAssayDimensioWIDBioAssayWID"
ADD CONSTRAINT "FK_BioAssayDimensioWIDBioAss2"
FOREIGN KEY  ("BioAssayWID")
REFERENCES  "BioAssay" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioAssayMapWIDBioAssayWID"
ADD CONSTRAINT "FK_BioAssayMapWIDBioAssayWID1"
FOREIGN KEY  ("BioAssayMapWID")
REFERENCES  "BioEvent" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioAssayMapWIDBioAssayWID"
ADD CONSTRAINT "FK_BioAssayMapWIDBioAssayWID2"
FOREIGN KEY  ("BioAssayWID")
REFERENCES  "BioAssay" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BAssayMappingWIDBAssayMapWID"
ADD CONSTRAINT "FK_BAssayMappingWIDBAssayMap1"
FOREIGN KEY  ("BioAssayMappingWID")
REFERENCES  "BioAssayMapping" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BAssayMappingWIDBAssayMapWID"
ADD CONSTRAINT "FK_BAssayMappingWIDBAssayMap2"
FOREIGN KEY  ("BioAssayMapWID")
REFERENCES  "BioEvent" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ComposSeqDimensWIDComposSeqWID"
ADD CONSTRAINT "FK_ComposSeqDimensWIDComposS1"
FOREIGN KEY  ("CompositeSequenceDimensionWID")
REFERENCES  "DesignElementDimension" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ComposSeqDimensWIDComposSeqWID"
ADD CONSTRAINT "FK_ComposSeqDimensWIDComposS2"
FOREIGN KEY  ("CompositeSequenceWID")
REFERENCES  "DesignElement" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "DesnElMappingWIDDesnElMapWID"
ADD CONSTRAINT "FK_DesnElMappingWIDDesnElMap1"
FOREIGN KEY  ("DesignElementMappingWID")
REFERENCES  "DesignElementMapping" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "DesnElMappingWIDDesnElMapWID"
ADD CONSTRAINT "FK_DesnElMappingWIDDesnElMap2"
FOREIGN KEY  ("DesignElementMapWID")
REFERENCES  "BioEvent" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "FeatureDimensionWIDFeatureWID"
ADD CONSTRAINT "FK_FeatureDimensionWIDFeatur1"
FOREIGN KEY  ("FeatureDimensionWID")
REFERENCES  "DesignElementDimension" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "FeatureDimensionWIDFeatureWID"
ADD CONSTRAINT "FK_FeatureDimensionWIDFeatur2"
FOREIGN KEY  ("FeatureWID")
REFERENCES  "DesignElement" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "QuantTypeDimensWIDQuantTypeWID"
ADD CONSTRAINT "FK_QuantTypeDimensWIDQuantTy1"
FOREIGN KEY  ("QuantitationTypeDimensionWID")
REFERENCES  "QuantitationTypeDimension" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "QuantTypeDimensWIDQuantTypeWID"
ADD CONSTRAINT "FK_QuantTypeDimensWIDQuantTy2"
FOREIGN KEY  ("QuantitationTypeWID")
REFERENCES  "QuantitationType" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "QuantTypeMapWIDQuantTypeWID"
ADD CONSTRAINT "FK_QuantTypeMapWIDQuantTypeW1"
FOREIGN KEY  ("QuantitationTypeMapWID")
REFERENCES  "BioEvent" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "QuantTypeMapWIDQuantTypeWID"
ADD CONSTRAINT "FK_QuantTypeMapWIDQuantTypeW2"
FOREIGN KEY  ("QuantitationTypeWID")
REFERENCES  "QuantitationType" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "QuantTyMapWIDQuantTyMapWI"
ADD CONSTRAINT "FK_QuantTyMapWIDQuantTyMapWI1"
FOREIGN KEY  ("QuantitationTypeMappingWID")
REFERENCES  "QuantitationTypeMapping" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "QuantTyMapWIDQuantTyMapWI"
ADD CONSTRAINT "FK_QuantTyMapWIDQuantTyMapWI2"
FOREIGN KEY  ("QuantitationTypeMapWID")
REFERENCES  "BioEvent" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ReporterDimensWIDReporterWID"
ADD CONSTRAINT "FK_ReporterDimensWIDReporter1"
FOREIGN KEY  ("ReporterDimensionWID")
REFERENCES  "DesignElementDimension" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ReporterDimensWIDReporterWID"
ADD CONSTRAINT "FK_ReporterDimensWIDReporter2"
FOREIGN KEY  ("ReporterWID")
REFERENCES  "DesignElement" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "TransformWIDBioAssayDataWID"
ADD CONSTRAINT "FK_TransformWIDBioAssayDataW1"
FOREIGN KEY  ("TransformationWID")
REFERENCES  "BioEvent" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "TransformWIDBioAssayDataWID"
ADD CONSTRAINT "FK_TransformWIDBioAssayDataW2"
FOREIGN KEY  ("BioAssayDataWID")
REFERENCES  "BioAssayData" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioSourceWIDContactWID"
ADD CONSTRAINT "FK_BioSourceWIDContactWID1"
FOREIGN KEY  ("BioSourceWID")
REFERENCES  "BioSource" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioSourceWIDContactWID"
ADD CONSTRAINT "FK_BioSourceWIDContactWID2"
FOREIGN KEY  ("ContactWID")
REFERENCES  "Contact" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "LabeledExtractWIDCompoundWID"
ADD CONSTRAINT "FK_LabeledExtractWIDCompound1"
FOREIGN KEY  ("LabeledExtractWID")
REFERENCES  "BioSource" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "LabeledExtractWIDCompoundWID"
ADD CONSTRAINT "FK_LabeledExtractWIDCompound2"
FOREIGN KEY  ("CompoundWID")
REFERENCES  "Chemical" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioAssayWIDChannelWID"
ADD CONSTRAINT "FK_BioAssayWIDChannelWID1"
FOREIGN KEY  ("BioAssayWID")
REFERENCES  "BioAssay" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioAssayWIDChannelWID"
ADD CONSTRAINT "FK_BioAssayWIDChannelWID2"
FOREIGN KEY  ("ChannelWID")
REFERENCES  "Channel" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioAssayWIDFactorValueWID"
ADD CONSTRAINT "FK_BioAssayWIDFactorValueWID1"
FOREIGN KEY  ("BioAssayWID")
REFERENCES  "BioAssay" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "BioAssayWIDFactorValueWID"
ADD CONSTRAINT "FK_BioAssayWIDFactorValueWID2"
FOREIGN KEY  ("FactorValueWID")
REFERENCES  "FactorValue" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ChannelWIDCompoundWID"
ADD CONSTRAINT "FK_ChannelWIDCompoundWID1"
FOREIGN KEY  ("ChannelWID")
REFERENCES  "Channel" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ChannelWIDCompoundWID"
ADD CONSTRAINT "FK_ChannelWIDCompoundWID2"
FOREIGN KEY  ("CompoundWID")
REFERENCES  "Chemical" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "DerivBioAWIDDerivBioADataWID"
ADD CONSTRAINT "FK_DerivBioAWIDDerivBioAData1"
FOREIGN KEY  ("DerivedBioAssayWID")
REFERENCES  "BioAssay" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "DerivBioAWIDDerivBioADataWID"
ADD CONSTRAINT "FK_DerivBioAWIDDerivBioAData2"
FOREIGN KEY  ("DerivedBioAssayDataWID")
REFERENCES  "BioAssayData" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "DerivBioAssayWIDBioAssayMapWID"
ADD CONSTRAINT "FK_DerivBioAssayWIDBioAssayM1"
FOREIGN KEY  ("DerivedBioAssayWID")
REFERENCES  "BioAssay" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "DerivBioAssayWIDBioAssayMapWID"
ADD CONSTRAINT "FK_DerivBioAssayWIDBioAssayM2"
FOREIGN KEY  ("BioAssayMapWID")
REFERENCES  "BioEvent" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ImageWIDChannelWID"
ADD CONSTRAINT "FK_ImageWIDChannelWID1"
FOREIGN KEY  ("ImageWID")
REFERENCES  "Image" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ImageWIDChannelWID"
ADD CONSTRAINT "FK_ImageWIDChannelWID2"
FOREIGN KEY  ("ChannelWID")
REFERENCES  "Channel" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ImageAcquisitionWIDImageWID"
ADD CONSTRAINT "FK_ImageAcquisitionWIDImageW1"
FOREIGN KEY  ("ImageAcquisitionWID")
REFERENCES  "BioEvent" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ImageAcquisitionWIDImageWID"
ADD CONSTRAINT "FK_ImageAcquisitionWIDImageW2"
FOREIGN KEY  ("ImageWID")
REFERENCES  "Image" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "MeasBAssayWIDMeasBAssayDataWID"
ADD CONSTRAINT "FK_MeasBAssayWIDMeasBAssayDa1"
FOREIGN KEY  ("MeasuredBioAssayWID")
REFERENCES  "BioAssay" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "MeasBAssayWIDMeasBAssayDataWID"
ADD CONSTRAINT "FK_MeasBAssayWIDMeasBAssayDa2"
FOREIGN KEY  ("MeasuredBioAssayDataWID")
REFERENCES  "BioAssayData" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "HardwareWIDSoftwareWID"
ADD CONSTRAINT "FK_HardwareWIDSoftwareWID1"
FOREIGN KEY  ("HardwareWID")
REFERENCES  "Parameterizable" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "HardwareWIDSoftwareWID"
ADD CONSTRAINT "FK_HardwareWIDSoftwareWID2"
FOREIGN KEY  ("SoftwareWID")
REFERENCES  "Parameterizable" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "HardwareWIDContactWID"
ADD CONSTRAINT "FK_HardwareWIDContactWID1"
FOREIGN KEY  ("HardwareWID")
REFERENCES  "Parameterizable" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "HardwareWIDContactWID"
ADD CONSTRAINT "FK_HardwareWIDContactWID2"
FOREIGN KEY  ("ContactWID")
REFERENCES  "Contact" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ProtocolWIDHardwareWID"
ADD CONSTRAINT "FK_ProtocolWIDHardwareWID1"
FOREIGN KEY  ("ProtocolWID")
REFERENCES  "Parameterizable" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ProtocolWIDHardwareWID"
ADD CONSTRAINT "FK_ProtocolWIDHardwareWID2"
FOREIGN KEY  ("HardwareWID")
REFERENCES  "Parameterizable" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ProtocolWIDSoftwareWID"
ADD CONSTRAINT "FK_ProtocolWIDSoftwareWID1"
FOREIGN KEY  ("ProtocolWID")
REFERENCES  "Parameterizable" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ProtocolWIDSoftwareWID"
ADD CONSTRAINT "FK_ProtocolWIDSoftwareWID2"
FOREIGN KEY  ("SoftwareWID")
REFERENCES  "Parameterizable" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ProtocolApplWIDPersonWID"
ADD CONSTRAINT "FK_ProtocolApplWIDPersonWID1"
FOREIGN KEY  ("ProtocolApplicationWID")
REFERENCES  "ParameterizableApplication" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "ProtocolApplWIDPersonWID"
ADD CONSTRAINT "FK_ProtocolApplWIDPersonWID2"
FOREIGN KEY  ("PersonWID")
REFERENCES  "Contact" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "SoftwareWIDSoftwareWID"
ADD CONSTRAINT "FK_SoftwareWIDSoftwareWID1"
FOREIGN KEY  ("SoftwareWID1")
REFERENCES  "Parameterizable" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "SoftwareWIDSoftwareWID"
ADD CONSTRAINT "FK_SoftwareWIDSoftwareWID2"
FOREIGN KEY  ("SoftwareWID2")
REFERENCES  "Parameterizable" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "SoftwareWIDContactWID"
ADD CONSTRAINT "FK_SoftwareWIDContactWID1"
FOREIGN KEY  ("SoftwareWID")
REFERENCES  "Parameterizable" ("WID") ON DELETE CASCADE DEFERRABLE;

ALTER TABLE "SoftwareWIDContactWID"
ADD CONSTRAINT "FK_SoftwareWIDContactWID2"
FOREIGN KEY  ("ContactWID")
REFERENCES  "Contact" ("WID") ON DELETE CASCADE DEFERRABLE;



