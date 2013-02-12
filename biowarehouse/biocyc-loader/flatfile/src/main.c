/* *******************************************************************
 * The contents of this file are subject to the Mozilla Public License
 * Version 1.1 (the "License"); you may not use this file except in
 * compliance with the License. You may obtain a copy of the License at
 * http://www.mozilla.org/MPL/
 *
 * Software distributed under the License is distributed on an "AS IS"
 * basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See
 * the License for the specific language governing rights and
 * limitations under the License.
 *
 * The Original Code is the BioWarehouse.
 *
 * The Initial Developer of the Original Code is SRI International.
 * Portions created by SRI International are Copyright (C) 2004.
 * All Rights Reserved.
 ******************************************************************* */
/*
 * BioCyc loader for the BioSpice Data Warehouse
 * Thomas J Lee, SRI International, August 2002
 *
 * Liberally adapted from the KEGG loader of the Warehouse
 */

/* Control behavior of EXTERN macro */
#define MAIN_PGM 1

#include <stdio.h>
#include <getopt.h>
#include <dirent.h>
#include <limits.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>
#include "wh.h"
#include "main.h"
#include "db.h"

/* Globals related to input files */
extern int pubs_records, pubs_errors;
extern int compound_records, compound_errors;
extern int reaction_records, reaction_errors;
extern int protein_records, protein_errors;
extern int protseq_records, protseq_errors;
extern int gene_records, gene_errors;
extern int transunit_records, transunit_errors;
extern int promoter_records, promoter_errors;
extern int terminator_records, terminator_errors;
extern int dnabindsite_records, dnabindsite_errors;
extern int enzrxn_records, enzrxn_errors;
extern int regulation_records, regulation_errors;
extern int pathway_records, pathway_errors;
extern FILE *pubsin, *compoundin, *reactionin, *proteinin, *protseqin, *enzrxnin, *regulationin,
  *pathwayin, *transunitin, *promoterin, *terminatorin, *dnabindsitein, *genein;

extern char *organism_name;

char * min_version_supported = "11.5";
char * source_version;
char * source_db_name;
int link_hierarchically = 0;  // true => add row to DataSetHierarchy, parent dataset is latest "BioCyc"
int fail_on_missing_file = 1;  // -F option changes this
int loading_metacyc = 0;  // 0 unless loading organism named "MetaCyc"
char * release_date = "";

/* Global vars for the lexers */
char *lexedword;
char lexedchar;



/*
 * Take appropriate action for a missing file that is required
 * for all PGDBs except MetaCyc.
 */
void
missing_file(int file_is_required, char *filename) { 
  if (loading_metacyc)
    printf("Ignoring non-MetaCyc input file %s\n", filename);
  else if (file_is_required && fail_on_missing_file)
    fatal_error(filename);
  //--- else printf("*** Warning: missing file: %s\n", filename);
}


// Some false pathways may have been added to the Pathway table due to
// pathway links referencing non-pathways. This traverses the table of all
// defined pathway WIDs, querying the DBID table to see if pathway is present,
// and deletes that pathway if missing from DBID.
void
remove_nonpathways() {
  int slot;
  struct wid_object *tmp;
  for (slot=0; slot < WIDSLOTS; slot++) {
    for (tmp = pathway_wids[slot]; tmp; tmp = tmp->next) {
      // query this dataset's DBID table for pathway name
      if (!wh_select_dbid_otherwid(tmp->name, dataset_wid)) {
        printf("Removing nonpathway %s, WID %d\n", tmp->name, tmp->wid);
        wh_delete_object("Pathway", tmp->wid);
      }
    }
  }
}

void
init_loader() {
  /* Post-commandline processing Initializations */
  int biocyc_wid;
	       
  if (0 == strcasecmp(source_version, "unknown"))
    printf("*** Warning: no version specified, using 'unknown'\n");
  else if (atof(source_version) < atof(min_version_supported)) {
    fprintf(stderr, "Version %s specified in command line precedes earliest supported version %s; exiting\n",
            source_version, min_version_supported);
    exit(-1);
  }
     
  /*! Build full pathnames of input files */
  /*! Check that all files are present - abort if not */
	       
  /* Connect and initialize global WIDs */
  wh_connect();
  if (dataset_wid == -1) {
    dataset_wid = wh_select_dataset_wid("BioCyc");
    if (dataset_wid == 0) {
      printf("Warning: no dataset named BioCyc found in Warehouse, creating new dataset\n");
      db_make_biocyc_dataset("BioCyc", source_version, release_date);  // add entry to DataSet, remember global WID  
    }
    if (link_hierarchically) {
      biocyc_wid = dataset_wid;
      printf("Linking to BioCyc dataset - WID=%d\n", biocyc_wid);
      db_make_biocyc_dataset(source_db_name, source_version, release_date);  // add entry to DataSet, remember global WID
      wh_insert_named_linktable("DataSetHierarchy", "SuperWID", biocyc_wid, "SubWID", dataset_wid);
    }
    else
      printf("Merging into existing BioCyc dataset - WID=%d\n", dataset_wid);
  }
  else if (0 == dataset_wid)
    db_make_biocyc_dataset(source_db_name, source_version, release_date);  // add entry to Dataset, remember global WID
  else {
    if (!wh_dataset_exists(dataset_wid)) {
      fprintf(stderr, "Invalid dataset WID was specified with -w: %d\n", dataset_wid);
      exit(-1);
    }
    printf("Merging into existing dataset - WID=%d\n", dataset_wid);
  }
  
  // Add a reflexive entry to the dataset hierarchy, to indicate that this dataset 'contains' itself
  wh_insert_named_linktable("DataSetHierarchy", "SuperWID", dataset_wid, "SubWID", dataset_wid);
  
  /* Confirm presence of datasets used by this loader */
  multifun_dataset_wid = wh_select_dataset_wid("MultiFun Gene Ontology");  /* Find most recent Multifun load, if any */
  if (!multifun_dataset_wid)
    printf("Warning: MetaCyc Ontology MultiFun dataset is not loaded; no RelatedTerm rows for genes with TYPES attribibute referencing Multifun terms will be created\n");
  else 
    printf("Found Multifun DataSet WID=%d\n", multifun_dataset_wid);
  /* Confirm presence of datasets used by this loader */
  go_dataset_wid = wh_select_dataset_wid("Gene Ontology");  /* Find most recent Go load, if any */
  if (!go_dataset_wid)
    printf("Warning: Gene Ontology dataset is not loaded; no RelatedTerm rows for proteins with GO-TERMS attribibute referencing GO terms will be created\n");
  else 
    printf("Found Gene Ontology DataSet WID=%d\n", go_dataset_wid);
    
  /* Add entry to Organism, remember global WID */
  db_make_biocyc_organism(); 
  
  /* Init global name->WID mapping tables */
  widtable_init(citation_wids);
  widtable_init(protein_wids);
  widtable_init(pathway_wids);
  widtable_init(transunit_wids);
  widtable_init(promoter_positions);

  // Init predecessor subpathways
  predecessor_subpathways = NULL;
}


void
print_usage() {
#ifdef ORACLE
  printf("Required arguments:\n  -d datapath  -o \"OrganismName\" -u \"Userid/Password\" \n");
#elif DEF_MYSQL
  printf("Required arguments:\n  -d datapath -o \"OrganismName\" -u userId -h host -p password -b database_name \n");
#elif DEF_POSTGRES
  printf("Required arguments:\n  -d datapath -o \"OrganismName\" -u userId -h host -p password -b database_name \n");
#endif
  printf("Options:\n");
  printf("  -l datasetWID:\t Link dataset created by this loader to given dataset via the DataSetHierarchy table\n");
  printf("  -m:\t\t\t Merge loaded data into the most recently loaded dataset named BioCyc\n");
  printf("  -v version:\t\t [Default 'unknown'] Loaded into DataSet.Version, and checked against minimum version supported\n");
  printf("  -w datasetWID:\t Merge loaded data into the given dataset\n");
  printf("  -C commitFrequency:\t [Default 1000] Perform database commit every commitFrequency inputs\n");
  printf("  -F:\t\t\t Ignore any missing files\n");
}

int
main(int argc, char **argv) {
  /*
    Parses the flat file rep'n of a BsubCyc KB.
    Loads BioSPICE Warehouse with the resulting contents.

    Files input:
    pubs.dat
    compounds.dat
    reactions.dat
    proteins.dat
    protseq.dat
    enzrxns.dat
    regulation.dat
    genes.dat
    transunits.dat
    promoters.dat
    terminators.dat
    dnabindsites.dat
    pathways.dat
    classes.dat
  */  
  char *datadir = NULL;
  char pubsfilename[PATH_MAX + NAME_MAX];
  char compoundfilename[PATH_MAX + NAME_MAX];
  char reactionfilename[PATH_MAX + NAME_MAX];
  char proteinfilename[PATH_MAX + NAME_MAX];
  char protseqfilename[PATH_MAX + NAME_MAX];
  char enzrxnfilename[PATH_MAX + NAME_MAX];
  char regulationfilename[PATH_MAX + NAME_MAX];
  char genefilename[PATH_MAX + NAME_MAX];
  char transunitfilename[PATH_MAX + NAME_MAX];
  char promoterfilename[PATH_MAX + NAME_MAX];
  char terminatorfilename[PATH_MAX + NAME_MAX];
  char dnabindsitefilename[PATH_MAX + NAME_MAX];
  char pathwayfilename[PATH_MAX + NAME_MAX];
  
  FILE *pubsfile, *compoundfile, *reactionfile, *proteinfile, *protseqfile, *enzrxnfile, *regulationfile,
    *genefile, *transunitfile, *promoterfile, *terminatorfile, *dnabindsitefile, *pathwayfile;

  struct timeval time_started;
  long time_to_load; 
  int c;

  printf("BioSPICE Warehouse - BioCyc Dataset loader version %s\n",
	 WAREHOUSE_VERSION_STRING);
  if (argc <= 1) {
    print_usage();
    exit(0);
  }

  wh_init();
  commit_freq = 1000;
  source_db_name = "BioCyc";  /* default */
  source_version = "unknown"; /* default */
  dataset_wid = 0;            /* If not overridden, a new dataset is created */

  /* Follow getopt char with : if it has a required arg, :: if optional  */
  while ((c = getopt(argc, argv, "d:o:t:u:h:n:v:b:p:r:w:C:?Flm")) != -1)  
    switch (c) {
    case 'b':
      database_name = optarg;
      printf("Database name is %s\n", database_name);
      break;
    case 'C':
      commit_freq = atoi(optarg);
      break;
    case 'd':
      datadir = optarg;
      printf("Data directory is %s\n", datadir);
      break;
    case 'F':
      fail_on_missing_file = 0;
      printf("Missing files will cause warnings rather than failure \n");
      break;
    case 'h':
      host = optarg;
      printf("Host is %s\n", host);
      break;
    case 'l':
      dataset_wid = -1;  /* indicate merge into existing BioCyc dataset */
      link_hierarchically = 1;
      printf("Linking database to parent BioCyc dataset \n");
      break;
    case 'm':
      dataset_wid = -1;  /* indicate merge into existing BioCyc dataset */
      printf("Merging database into one BioCyc dataset \n");
      break;
    case 'n':
      source_db_name = optarg;
      if (0 == strcasecmp(source_db_name, "MetaCyc"))
        loading_metacyc = 1;
      printf("Source database name is %s\n", source_db_name);
      break;
    case 'o':
      organism_name = optarg;
      printf("Organism name is %s\n", organism_name);
      if (0 == strcasecmp(organism_name, "MetaCyc"))
        loading_metacyc = 1;
      break;
    case 'p':
      password = optarg;
      printf("Password is %s\n", password);
      break;
    case 'r':
      release_date = optarg;
      printf("Release date is %s\n", release_date);
      break;
    case 't':
      server_port = atoi(optarg);
      printf("Server port is %d\n", server_port);
      break;
    case 'u':
      userid = optarg;
      break;
    case 'v':
      source_version = optarg;
      printf("Source version is %s\n", source_version);
      break;
    case 'w':
      dataset_wid = atoi(optarg);
      // Merge notification is echoed to stdout later
      break;
    case '?':
      print_usage();
      exit(0);
    default:
      fprintf(stderr, "Run %s -? for help\n", argv[0]);
      exit(-1);
    }

  if (!organism_name) {
    fprintf(stderr, "%s: missing '-o OrganismName'.  Use -h for help.\n", argv[0]);
    exit(-1);
  }

  if (!datadir) {
    fprintf(stderr, "%s: missing '-d DATADIRECTORY'.  Use -h for help.\n", argv[0]);
    exit(-1);
  }

  if (strlen(datadir) > PATH_MAX) {
    fprintf(stderr, "Data directory name too long\n");
    exit(-1);
  }

  if (commit_freq)
    printf("Committing after every table load, and every %d WIDs\n", commit_freq);
    
  /* stores the time of the day when the loading of the database started */
  gettimeofday(&time_started, NULL); 

  init_loader();
  
  /*** Load DB a file at a time ***/
  
  widtable_init(undefined_words);
  sprintf(pubsfilename, "%s/%s\0", datadir, "pubs.dat");
  if (!(pubsfile = fopen(pubsfilename, "r")))
    missing_file(0, pubsfilename);
  else {
    printf("\nParsing PUBLICATIONS data from %s...\n", pubsfilename);
    pubsin = pubsfile;
    pubs_clear_entry(0);
    (void) pubsparse();
    printf("Found %d publications.  Encountered %d parse errors\n",
	   pubs_records, pubs_errors);
    if (commit_freq) wh_commit(); 
  }
  
  widtable_init(undefined_words);
  sprintf(compoundfilename, "%s/%s\0", datadir, "compounds.dat");
  if (!(compoundfile = fopen(compoundfilename, "r")))
    missing_file(1, compoundfilename);
  else {
    printf("\nParsing COMPOUND data from %s...\n", compoundfilename);
    compoundin = compoundfile;
    compound_clear_entry(0);
    (void) compoundparse();
    printf("Found %d compounds.  Encountered %d parse errors\n",
           compound_records, compound_errors);
    if (commit_freq) wh_commit(); 
  }

  widtable_init(undefined_words);
  sprintf(proteinfilename, "%s/%s\0", datadir, "proteins.dat");
  if (!(proteinfile = fopen(proteinfilename, "r")))
    missing_file(1, proteinfilename);
  else {
    printf("\nParsing PROTEIN data from %s...\n", proteinfilename);
    proteinin = proteinfile;
    protein_clear_entry(0);
    (void) proteinparse();
    printf("Found %d proteins.  Encountered %d parse errors\n",
           protein_records, protein_errors);
    if (commit_freq) wh_commit(); 
  }  
  widtable_init(undefined_words);
  sprintf(protseqfilename, "%s/%s\0", datadir, "protseq.fasta");
  if (!(protseqfile = fopen(protseqfilename, "r")))
    missing_file(0, protseqfilename);
  else {
    printf("\nParsing PROTEIN SEQUENCE data from %s...\n", protseqfilename);
    protseqin = protseqfile;
    protseq_clear_entry(0);
    (void) protseqparse();
    printf("Found %d protein sequences.  Encountered %d parse errors\n",
	   protseq_records, protseq_errors);
    if (commit_freq) wh_commit(); 
  }
  
  widtable_init(undefined_words);
  sprintf(transunitfilename, "%s/%s\0", datadir, "transunits.dat");
  if (!(transunitfile = fopen(transunitfilename, "r")))
    missing_file(0, transunitfilename);
  else {
    printf("\nParsing TRANSUNIT data from %s...\n", transunitfilename);
    transunitin = transunitfile;
    transunit_clear_entry(0);
    (void) transunitparse();
    printf("Found %d transcription units.  Encountered %d parse errors\n",
           transunit_records, transunit_errors);
    if (commit_freq) wh_commit();
  }

  widtable_init(undefined_words);
  sprintf(genefilename, "%s/%s\0", datadir, "genes.dat");
  if (!(genefile = fopen(genefilename, "r")))
    missing_file(1, genefilename);
  else {
    printf("\nParsing GENE data from %s...\n", genefilename);
    genein = genefile;
    gene_clear_entry(0);
    (void) geneparse();
    printf("Found %d genes.  Encountered %d parse errors\n",
           gene_records, gene_errors);
    if (commit_freq) wh_commit(); 
  }
  
  widtable_init(undefined_words);
  sprintf(promoterfilename, "%s/%s\0", datadir, "promoters.dat");
  if (!(promoterfile = fopen(promoterfilename, "r")))
    missing_file(1, promoterfilename);
  else {
    printf("\nParsing PROMOTER data from %s...\n", promoterfilename);
    promoterin = promoterfile;
    promoter_clear_entry(0);
    (void) promoterparse();
    printf("Found %d promoters.  Encountered %d parse errors\n",
           promoter_records, promoter_errors);
    if (commit_freq) wh_commit();
  }

  widtable_init(undefined_words);
  sprintf(terminatorfilename, "%s/%s\0", datadir, "terminators.dat");
  if (!(terminatorfile = fopen(terminatorfilename, "r")))
    missing_file(1, terminatorfilename);
  else {
    printf("\nParsing TERMINATOR data from %s...\n", terminatorfilename);
    terminatorin = terminatorfile;
    terminator_clear_entry(0);
    (void) terminatorparse();
    printf("Found %d terminators.  Encountered %d parse errors\n",
           terminator_records, terminator_errors);
    if (commit_freq) wh_commit(); 
  }

  widtable_init(undefined_words);
  sprintf(dnabindsitefilename, "%s/%s\0", datadir, "dnabindsites.dat");
  if (!(dnabindsitefile = fopen(dnabindsitefilename, "r")))
    missing_file(1, dnabindsitefilename);
  else {
    printf("\nParsing DNABINDSITE data from %s...\n", dnabindsitefilename);
    dnabindsitein = dnabindsitefile;
    dnabindsite_clear_entry(0);
    (void) dnabindsiteparse();
    printf("Found %d dnabindsites.  Encountered %d parse errors\n",
           dnabindsite_records, dnabindsite_errors);
    if (commit_freq) wh_commit(); 
  }

  widtable_init(undefined_words);
  sprintf(reactionfilename, "%s/%s\0", datadir, "reactions.dat");
  if (!(reactionfile = fopen(reactionfilename, "r")))
    missing_file(1, reactionfilename);
  else {
    printf("\nParsing REACTION data from %s...\n", reactionfilename);
    reactionin = reactionfile;
    reaction_clear_entry(0);
    (void) reactionparse();
    printf("Found %d reactions.  Encountered %d parse errors\n",
           reaction_records, reaction_errors);
    if (commit_freq) wh_commit(); 
  }
  
  widtable_init(undefined_words);
  sprintf(enzrxnfilename, "%s/%s\0", datadir, "enzrxns.dat");
  if (!(enzrxnfile = fopen(enzrxnfilename, "r")))
    missing_file(1, enzrxnfilename);
  else {
    printf("\nParsing ENZRXN data from %s...\n", enzrxnfilename);
    enzrxnin = enzrxnfile;
    enzrxn_clear_entry(0);
    (void) enzrxnparse();
    printf("Found %d enzrxns.  Encountered %d parse errors\n", enzrxn_records, enzrxn_errors);
    if (commit_freq) wh_commit(); 
  }

  widtable_init(undefined_words);
  sprintf(regulationfilename, "%s/%s\0", datadir, "regulation.dat");
  if (!(regulationfile = fopen(regulationfilename, "r")))
    missing_file(1, regulationfilename);
  else {
    printf("\nParsing REGULATION data from %s...\n", regulationfilename);
    regulationin = regulationfile;
    regulation_clear_entry(0);
    (void) regulationparse();
    printf("Found %d regulation entries.  Encountered %d parse errors\n", regulation_records, regulation_errors);
    if (commit_freq) wh_commit(); 
  }
  
  widtable_init(undefined_words);
  sprintf(pathwayfilename, "%s/%s\0", datadir, "pathways.dat");
  if (!(pathwayfile = fopen(pathwayfilename, "r")))
    missing_file(1, pathwayfilename);
  else {
    printf("\nParsing PATHWAY data from %s...\n", pathwayfilename);
    pathwayin = pathwayfile;
    pathway_clear_entry(0);
    (void) pathwayparse();
    printf("Found %d pathways.  Encountered %d parse errors\n",
           pathway_records, pathway_errors);
    if (commit_freq) wh_commit(); 
    db_insert_all_subpathway_predecessors();
    remove_nonpathways();
  }
  
  /* All done loading */
  time_to_load = GetTimeDiff(&time_started);
  printf("\nTotal load time: %ld minutes\n" , time_to_load/60000);

  /* Close up shop */
  wh_commit(); 
  wh_finalize_dataset();  // set time of completion of loader
  wh_disconnect();
  exit(0);
}

/*
 * fatal errors (e.g. out of memory)
 */
void
fatal_error(char *message) { 
  perror(message);
  exit(-1);
}

struct dblinks *
mk_dblink(char *database, char *xid, struct stringlist *description) {

  struct dblinks *d;

  d = (struct dblinks *) malloc (sizeof(struct dblinks));
  if (!d)
    fatal_error("mk_dblinks malloc");

  d->database    = database ;
  d->xid         = xid ;
  d->description = flatten_list_old(description);
  d->next        = NULL_DBLINKS;
  return(d);
}

struct dblinks *
mk_dblinks(char *database, struct stringlist *xids) {

  struct stringlist *xid = xids;
  struct dblinks *d = NULL_DBLINKS;

  for (;xid;xid=xid->next)
    d = prepend_dblink(mk_dblink(strdup(database), strdup(xid->string), NULL_LIST), d);

  free_stringlist(xids);
  return(d);
}

struct dblinks *
prepend_dblink(struct dblinks *dblinks1, struct dblinks *dblinks2) {
  
  struct dblinks *ptr = dblinks1;

  if (!dblinks1)
    return(dblinks2);
  for (;ptr->next;ptr=ptr->next);
  ptr->next = dblinks2;
  return(dblinks1);
}

void
free_dblinks(struct dblinks *dblink) {
  if (dblink) {
    free_dblinks(dblink->next);
    free(dblink->database);
    free(dblink->xid);
    free(dblink->description);
    free(dblink);
  }
}


void
insert_dblinks(int object_wid, struct stringlist * dbs, struct stringlist * ids) {
  /* Given corresponding lists of DB names and unique IDs, add a row to CrossReference for each */
  struct stringlist *dbptr = dbs;
  struct stringlist *idptr = ids;
  //- printf ("insert_dblinks\n"); 
  while (dbptr && idptr) {
    wh_insert_crossreference(object_wid, idptr->string, dbptr->string);
    //- printf("Crossref for %s\n", idptr->string); 
    dbptr = dbptr->next;
    idptr = idptr->next;
  }
}


int
find_promoter_position(char *name) {
  return(widtable_lookup(promoter_positions, name));
}

int
find_transunit(char *name) {
  return(widtable_lookup(transunit_wids, name));
}

int
find_pathway(char *name) {
  return(widtable_lookup(pathway_wids, name));
}

int
find_protein(char *name) {
  //- return(db_select_protein(name));
  return(widtable_lookup(protein_wids, name));
}

int
find_citation(char *citation) {
  // Returns the WID of the Citation matching the give id.
  // Matches if the UNIQUE-ID of a publication from pubs.dat:
  //  - matches id exactly
  //  - matches "PUB-" + id
  // where id is the citation given, with its optional enclosing []s removed.
			       
  char * id;
  int result;
  int len;
  
  //- printf("find_citation citation=%s\n", citation);
  if (!citation || 0 == strlen(citation)) return 0;
  if (citation[0] == ':' || strstr(citation, ":EV-")) {
    // Assume citation is a reference to evidence and igonre
    return 0;
  }
  
  len = strlen(citation);
  id = malloc(len + 5);  // leave room for PUB- prefix and null terminator 
  strcpy(id, "PUB-");
  if (citation[0] == '[' ) {
    strncpy(id+4, citation+1, len - 1); // Remove leading [ and last char (assumedly ]}
    id[4+len-2] = '\0';
  }
  else
    strcpy(id+4, citation);
  string_uppercase(id);
  //- printf("find_citation id=%s\n", id);
  
  result = widtable_lookup(citation_wids, id);
  if (!result)
    result = widtable_lookup(citation_wids, id+4); // Effectively remove PUBS- prefix
  
  free(id);
  return result;
}

short
insert_citation(int other_wid, char *citation) {
  // If publication for given citation has been parsed,
  // insert a row into CitationWIDOtherWID.
  // Else if it is an evidence citation, add a row to Support for it.
  // Else issue a warning and return an error indication.

  short error = 0;
  int citation_wid;

  assert(other_wid != 0);
  if (citation && strlen(citation) > 0) {
    if (citation[0] == ':' || strstr(citation, ":EV-")) {
      // Assume citation is a reference to evidence and insert it into the
      // Support table, along with the CitationWIDOtherWID table if the
      // evidence has a citation as well.
      error = insert_support(other_wid, citation);
      return error;
    }
    citation_wid = find_citation(citation);
    if (citation_wid)
      wh_insert_linktable("CitationWID", citation_wid, "OtherWID", other_wid);
    else {
      printf("Warning: no publication found for citation %s\n", citation);
      error = 1;
    }
  }
  return error;
}

/*
  Called by insert_citation when the citation is a reference to evidence
  Param:
  other_wid - the WID of the fact with the evidence
  citation - the citation from the BioCyc file
  Return:
  0 if no error, else 1
*/
short
insert_support(int other_wid, char *citation) {
  char *delim = ":"; // delimiter for the strtok function
  char *evidence; // evidence of an entity, extracted from citation
  char *citation_id; // id of a citation, extracted from citation
  int citation_wid; // WID of a citation
  int support_wid; // WID of the support table
  short evidence_ind; // Indicator is set depending on whether evidence is NULL
  char *buffer; // temporary buffer used in strtok
  short error = 0; // return value
  
  buffer = malloc(strlen(citation) + 1);
  strcpy(buffer, citation);
  if (citation[0] == ':') { // Evidence with no literature citation
    evidence = strtok(buffer, delim);
    evidence = string_column(evidence, 40, &evidence_ind);
    support_wid = wh_get_new_wid();
    wh_insert_support(support_wid, other_wid, evidence, evidence_ind);
    wh_insert_entry(support_wid, 0, 0, NULL);
  }
  else if (strstr(citation, ":EV-")) { // Evidence with literature citation
    citation_id = strtok(buffer, delim);
    evidence = strtok(NULL, delim);
    evidence = string_column(evidence, 40, &evidence_ind);
    citation_wid = find_citation(citation_id);
    if (citation_wid) {
      support_wid = wh_get_new_wid();
      wh_insert_support(support_wid, other_wid, evidence, evidence_ind);
      wh_insert_linktable("CitationWID", citation_wid, "OtherWID", support_wid);
      wh_insert_entry(support_wid, 0, 0, NULL);
    }
    else {
      printf("Warning: no publication found for citation %s\n", citation);
      error = 1;
    }      
  }
  else { // Not an evidence citation
    printf("Warning: expected evidence citation, but found citation without evidence: %s\n",
	   citation);
    error = 1;
  }
  
  free(buffer);
  return error;
}

int
obtain_enzrxn_compound(char *name, int allow_protein) {
  /* If name is defined in compounds.dat, return its WID by querying for it
     Else if allow_protein is true and name is defined in proteins.dat, return its WID.
     Else create a Chemical entry for the name, populating only the name & datasetwid,
     and return its WID.
     NOTE: In this case, no checking of name is done, even though
     it may be a class of chemicals, a loose description of some compounds, a protein, etc.
  */
  int compound_wid;
  compound_wid = db_select_chemical_wid(name);
  /*- printf("Found chemical %d\n", compound_wid);*/
  if (compound_wid) return compound_wid;
  
  if (allow_protein) {
    compound_wid = db_select_protein(name);
    /*- printf("Found protein %d\n", compound_wid);*/
    if (compound_wid) return compound_wid;
  }
  
  /* Compound is new; add it to Chemical table */
  compound_wid = wh_get_new_wid();
  wh_insert_chemical_simple(compound_wid, name);
  if (0) printf(">>> Chemical %s created\n", name);
  wh_insert_comment(compound_wid, ">>> Chemical added by BioCyc loader due to its presence in an enzymatic reaction or regulator");
  return compound_wid;
}

#ifdef IGNORE

#endif
