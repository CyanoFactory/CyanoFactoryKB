Bio-SPICE BioWarehouse Version 4.6
June 17, 2009

For documentation on the Bio-SPICE BioWarehouse, open doc/index.html with the
browser of your choice.

CHANGES for release 4.6 (June 2009)
=======================

Schema:
  * Added index: BioSource.Strain.
  * Added column: Feature.Variant

All C loaders:
  * Accept YYYY/MM/DD format for release date.

KEGG loader:
  * Support for KEGG version 50.
  * Major change: the LIGAND:REACTION file is now translated. Most reactions
    are created by translating it rather than during translation of the
    LIGAND:ENZYME file.

CMR loader:
  * Changed logic to exclude non-original RNA locuses based on feat_type
    rather than naming conventions (which are not 100% reliable).
        
UniProt loader:
  * Updated to Uniprot 15.2 schema
  * Improved handling of feature references
  * Now store <variation>s in Feature.Variant


CHANGES for release 4.5 (December 2008)
=======================

Schema:
  * Modified column: DataSet.ReleaseDate to be Date type instead of String

All loaders:
  * Release date must be provided via the -r command line argument

All Java loaders:
  * Java loaders now require using Java 1.6+
  
All C loaders:
  * Fixed Makefile bug in which not all old files were cleaned out;
    this may have caused strange problems when compiling for different architectures
    
BioCyc Loader:
  * Extended allowed syntax for pathway links to gracefully ignore non-pathways

KEGG loader:
  * Support for KEGG version 46.0. NOTE: the loader gives several hundred spurious
    syntax errors when parsing the Ligand:Enzyme file upon encountering PRODUCT or
    SUBSTRATE entries that are not tagged with a KEGG Compound identifier; these are
    not loaded into the Chemical table.

NCBI Taxonomy loader:
  * Warns if a non-date or a date format other than YYYY-MM-DD is used to specify
    the dataset version

UniProt loader:
  * Updated to UniProt format version 14.2


CHANGES for release 4.4 (October 2008)
=======================

Uniprot Loader:
  * Updated loader to be compatible with latest Uniprot Schema (XSD), release 14.2


CHANGES for release 4.3 (September 2008)
=======================

Schema changes:
  * New table: DataSetHierarchy
    Datasets may now be hierarchically organized. To convert any query of a
    specific dataset into one that queries that dataset and/or any contained datasets,
    join with DataSetHierarchy and add
    'AND SomeTable.datasetwid=DataSetHierarchy.subwid AND DataSetHierarchy.superwid=X'
    to the SELECT clause where X is the dataset being queried.

New loader - ChIP-Chip Loader
  * Loads data for chIP-chip (chromatin immunoprecipitation with genome-tiling DNA chip expression) experiments
  
BioCyc loader:
  * New options: -l (link hierarchically to most recent "BioCyc" dataset)
                 -w <datasetwid> (add all data to an existing dataset)
                 -F (makes all input files optional)
  * Bug alert: Under Oracle, -w <datasetwid> hangs if datasetwid does not refer to
               a valid dataset.
		 
  
CHANGES for release 4.2 (April 2008)
=======================

Schema changes:
  * Added EvidenceType column to Support table and enumeration descriptions
  
BioCyc loader:
  * Loads regulation data from file regulation.dat:
    - Full support for enzymatic regulation
    - Minimal support for transcriptional regulation (see BioCyc loader manual)
  * Fixed bug in which no gene was linked to its transcription unit

UniProt Loader:
  * Domain and component names are now treated as synonyms of the protein.

CHANGES for release 4.1 (June 2007)
=======================

Schema changes:
  * New tables: TranscriptionUnit and TranscriptionUnitComponent.
  * New columns to Feature table: RegionOrPoint and PointType.

BioCyc loader:
  * Loads transcriptions units.
  * Loads promoters, terminators, and DNA binding sites.
  * Links MetaCyc terms to genes (if MetaCyc Ontology is loaded)
  * Links GO terms for proteins (if GO Ontology is loaded)
  * Fixed bug in loading publications - referent frames were sometimes
    clobbering the Citation columns of the publication referred to

MAGE loader:
  * Fixed bug where MAGE Description element was improperly mapped
  * MAGE OntologyEntry.category now mapped as parent Term of OntologyEntry.name
  

CHANGES for release 4.0 (November 2006)
=======================

Schema changes:
  * Extensions to the schema to accomodate MicroArray Gene Expression data
    Schema extensions were derived from the MAGE Object Model; see the MAGE
    loader documentation for more information (mage-loader/doc/index.html).
    109 new tables added to the schema.
  * Extensions to represent flow cytometry data:
    - new general-purpose table RelatedTerm for associatng a term with an object
    - new tables LightSource, FlowCytometryProbe, FlowCytometrySample
    - new BioSource columns Diseased, Disease, ATCCId
    - new ExperimentData column MageData
    - see the new tables and ExperimentData.html for representation details
  * Extensions to represent 2D protein gel data. Added tables ExperimentRelationship, GelLocation,
    ProteinWIDSpotWID, Spot, SpotIdMethod, and SpotWIDSpotIdMethodWID.
  * Added range of reserved WIDs to loaders with particular WID needs (used by eco2dbase loader)
    Special WID range - 1 to 999
    Reserved WID range - 1,000 to 999,999
    Regular WID range - 1,000,000 on
  * Replaced PointOfContact table with Contact table
  * New schema source file organization
  * New schema maintenance script (schema/build.xml)
  * See complete change list at head of schema/core-schema.xml

New loader - MicroArray Gene Expression (MAGE) Loader
  * Loads gene expression data in MAGE-ML format
  * See loader documentation at mage-loader/doc/index.html

New loader - BioPAX Loader
  * Loads data of BioPAX format
  * See loader documentation at biopax-loader/doc/index.html

New loader - Eco2dbase Loader
  * The Eco2dbase data set contains E. Coli 2D Protein Gel data
  * We provide MySQL and Oracle dump files of the eco2dbase data set that
    are directly loaded into the BioWarehouse
  * See loader documentation at eco2dbase-loader/doc/index.html

All Java loaders
  * Changed commit mechanism; should result in slight loader performance increases

New Warehouse Utility Programs
  * See doc/Utils.html for complete details
  * Master build and install script
  * Database summary utilities
  * Dataset deletion tool

New documentation:
  * How to Write a Loader (doc/HowToWriteALoader.html)
  * Schema design and modification guidelines (in doc/Schema.html)


CHANGES for release 3.8 (September 2006)
=======================

Schema changes:
  * New tables: Interaction and InteractionParticipant. These are not currently
    used, but will be used to represent molecular interactions.
  * New column: Chemical.Class ('T' if a chemical represents a class of related
    chemicals).

BioCyc loader:
  * Bug fix for missing pathway reactions when reaction has no predecessors.
  * Bug fix: proteins are now included in Reactant and Product tables when they
    occur as reaction substrates.
  * New feature: substrates occurring in reactions, but not defined as compounds,
    are assumed to be classes of compounds. Each class is added to the Chemical table.

    
CHANGES for release 3.7 (August 2006)
=======================

CMR Loader:
  * Support for CMR 19.0 - deprecate input files bcp_ident and bcp_nt_ident
    in favor of bcp_new_ident. Also, see the indexing change note below.

UniProt Loader:
  * Support for Uniprot 8.0 and later.

All C-based loaders:
  * Bug fix for system header files not included that cause errors
    on some platforms.

No schema table changes.

Indexing change:
  * For MySQL, indexing of the SequenceMatch table has been disabled by
    commenting out index statements in schema/mysql-index.sql. This should
    only affect users that load the CMR database, or who load their own data
    into this table. The reason is that the table size (180 million rows)
    after the load of CMR makes index computation prohibitively slow.
    We are working on a solution; if you need this table indexed, please contact us.
    

CHANGES for release 3.6 (April 2006)
=======================

New support email address:
  * For support and inquiries please contact support@biowarehouse.org

New loader - MetaCyc Ontology:
  * Loads the MetaCyc Chemical Compounds ontology, the Multifun Gene ontology,
    and the MetaCyc Pathway ontology into three separate datasets.
  * Each ontology term is stored in the Term table. 
    Superclass relationships are stored in the TermRelationship table.
    Synonyms and alternative identifiers for a term are stored in the SynonymTable.
  * See the loader documentation for full details about the schema mappings
    and how to build and run the loader (metacyc-ontology-loader/doc/index.html).
    
Major schema changes:
  * Added columns Term.Hierarchical and Term.Root to better represent hierarchcial
    ontologies.
  * Added DataSet.ChangeDate; it is NULL while a loader is running or if a fatal
    error has occurred, and can be used by applications to track change times.

KEGG Loader:
  * Support for versions 34.0 - 37.0. Loader now requires bison 1.875+.
  
UniProt Loader:
  * Support for verison 7.1 - 7.2 of UniProt.
  
CHANGES for release 3.5 (July 2005)
=======================

New loader - Gene Ontology:
  * Loads the Gene Ontology (termdb data only) into the BioWarehouse.  Associations
    (instance data) are not loaded at this time.
  * Each ontology term is stored in the Term table.  Relationships between terms 
    (is_a or part_of relationships) are stored in the TermRelationship table.
    Synonyms and alternative identifiers for a term are stored in the SynonymTable.
  * See the loader documentation for full details about the schema mappings
    and how to build and run the loader (go-loader/docs/index.html).
  * See also http://geneontology.org/.

MySQL changes:
  * Version 4.1.2+ of MySQL is now supported. Most recently tested version is 4.1.10.
    FOREIGN KEY constraint checking is now enforced in MySQL.

Oracle changes:
  * Upgraded the Oracle JDBC driver to be compatible with Oracle 10.1.0.4.0. 
    If you experience any difficulty running the Java loaders against your
    version of Oracle, please email us at support@biowarehouse.org.

Major schema changes:
  * Added columns LoadedBy, Application, and ApplicationVersion to
    DataSet table.  These are populated by our loader applications,
    with LoadedBy being the value of the $USER environment variable,
    Application being the name of the loader program, and
    ApplicationVersion being the version of the Biowarehouse release (e.g., 3.5).
  * Removed table TermWIDParentWID in favor of the new table TermRelationship.
    Added column Term.Definition.

Uniprot/SwissProt/TrEMBL Loader:
  * Updated the loader to work with UniProt schema v1.14 (UniProt release 5.0).
  * The contents of <name> are now stored in the SynonymTable instead of DBID.
  * Protein.MolecularWeightCalc is now stored as kDaltons, not Daltons.
  * BioSource.Strain is now populated by the contents of <strain> as the individual
    <name> elements under strain have been deprecated in the UniProt schema.
  
CMR Loader:
  * Support for another command line argument for the load variant - original
    (see cmr-loader/doc/index.html for usage, and cmr-loader/doc/cmr-manual.html
    for discussion): if a genome has been sequenced at a non-TIGR
    site, only the non-TIGR ORFs are loaded; all TIGR ORFs for the
    genome are excluded.  If a genome has been sequenced at TIGR,
    these TIGR ORFs are loaded.
  * Numerous bug fixes:
    - Strain name is now stripped out of Biosource.Name.
      (It still populates Biosource.Strain).
    - Gene.Direction, Gene.CodingRegionStart and Gene.CodingRegionEnd are now
      populated as documented in the schema.
    - NucleicAcid.MoleculeLength is now populated for both DNA and RNA.
  * Support for CMR version 16.0 (May 2005):
    - Read files named bcp_all_vs_all_00, bcp_all_vs_all_01, etc., as well as 
      those named using the old convention bcp_all_vs_all_1, bcp_all_vs_all_2, etc.

BioCyc Loader:
  * Fixed a bug (manifested in the MetaCyc DB) that caused many proteins and
    enzymatic reactions to be omitted from the load.
  * Support for version 9.1 of BioCyc:
    - If a subpathway occurs in the PREDECESSORS slot of a pathway, all its reactions
      are included in the reaction graph for that pathway that is stored in the
      PathwayReaction table.

NCBI Taxonomy Loader:
  * Removed extraneous white space from some names that also contain
    extraneous angle-bracketed information (e.g. the name of
    "Bacteria <sometag>" is now stored as "Bacteria").    

Enzyme Loader:
  * Application command-line parameters changed to accomodate setting of the
    data version and release date.


CHANGES for release 3.1
=======================

Bug fixes:
  * Some of the binary files in release 3.0 were corrupted.  They have been replaced with fresh copies.


CHANGES for release 3.0 (December 2004)
=======================

Major schema changes:
  * Organism table was removed in favor of BioSource and BioSubtype tables
  * NucleicAcid table was significantly generalized
  * Replicon table was removed in favor of generalized NucleicAcid table
  * Subsequence table was added
  * Renamed EnzReactionWIDChemicalWID table to EnzReactionAltCompound and added
    Cofactor column
  * Made Gene.Name nullable
  * Added GeneticCode.NCBIID
  * Added Description table, used for certain defining or distinguished comments
  * Generalized the Feature table
  * Clarified that GeneWIDNucleicAcidWID associates a gene with its
    nucleic acid product 
  * Added Version column to DBID and CrossReference tables
  
New tool - Enumeration definitions loader
  * Used when installing the Warehouse
  * Loads the Enumeration table with definitions for all terminology and 
    abbreviations used in various data columns of Warehouse tables

New loader - GenBank
  * Loads the XML version of the GenBank database.
  * Files need to be converted from the ASN1 format to XML, as described
    in the GenBank BioWarehouse documentation.
  * Currently only the BCT division has been tested.

Swiss-Prot loader:
  * The Swiss-Prot loader has been re-implemented as the "UniProt Loader",
    which reads files in the UniProt XML format.

Enzyme loader:
  * GUI version has been deprecated in favor of the command line version.
  * Simplified command line options.
    
CMR loader:
  * Updated for latest version of data (circa November 2004)
  * Now supports three variants - all, primary, and TIGR -
    for loading different subsets of Open Reading Frames
  * If NCBI Taxonomy database is loaded, populates NucleicAcid.GeneticCodeWID
  * If Enzyme database is loaded, populates EnzymaticReaction with
    EC number and associated reaction

KEGG loader:
  * Updated for latest version of data (32.0)
  * Improved parsing to remove spurious parse errors.
  * Converted manual to HTML format.

BioCyc loader:    
  * Updated for latest version of data (8.5)
  * Populates the EnzReactionCofactor, EnzReactionAltCompound, and
    EnzReactionInhibitorActivator tables
  * Stores all comments defined in the source database
  * Stores all DBLINKS defined in the source database in the CrossReference table

NCBI Taxonomy loader:
  * Updated for latest version of data (circa August 2004)

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
CHANGES for release 2.1 (2004)
=======================
Added mysql and oracle targets to schemadoc Makefile to construct the documentation.

Separated non-required table indexes from the schema proper:
  * Updated schemadoc tool to accept concatenated schema+index file.
  * Updated translate-schema.lisp to support separate file of indexes.

  
CHANGES for release 2.0.2
=========================
Support for KEGG version 27.0:
  * Support RELEASE property in genome file for version 27+.
  * Several fixes to compound and enzyme grammars.

Fixed bug in cmr-loader (Oracle): Gene.Direction is now always NULL,
since CMR apparently does not specify a transcription direction for its features.

Fixed bug in biocyc-loader: some syntax errors in genes, compounds, and
enxymatic reactions were not getting counted in error totals.

Fixed bug in miamexpress-loader/src/miamexpress-archive-loader.perl:
incorrect value for Experiment.ArchiveWID was being stored.


CHANGES for release 2.0.1
=========================
CMR loader - removed stale files from src/ that were causing a compile error.
NCBI loader (Oracle) - now allocates a special (low-valued) WID for its dataset.


CHANGES for release 2.0 (December 2003)
=======================

Extensive changes to support for the MySQL Database Management System.

Support for representation of experimental data within a Warehouse dataset.


New loaders
-----------

CMR (Comprehensive Microbial Resource)
NCBI Taxonomy
MIAMExpress archive

Versions of all loaders for MySQL.


Schema changes
--------------

The old schema file, warehouse-schema.sql, has been replaced with a
DBMS-independent template, schema-template.sql, as well as specific
schema for Oracle and MySQL, oracle-schema.sql and mysql-schema.sql.
The template is used to generate the specific schema files
programmatically.

New tables:

Warehouse
WIDTable
SpecialWIDTable
NucleicAcid
NucleicAcidWIDFeatureWID
Computation
Experiment
ExperimentData
Archive
PointOfContact

New columns:

Replicon.NucleicAcidWID

Gene.NucleicAcidWID


Loader changes
--------------
Numerous fixes to error detection and error count reporting have been made.


Release 1.0 (November 2002)
-----------
