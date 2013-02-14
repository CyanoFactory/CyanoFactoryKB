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
#include "mysql.h"
#include "stdio.h"
#include <errno.h>
#include "wh.h"
#include "wh_mysql_util.h"

#define MAX_QUERY_SIZE 2000000

MYSQL mysql_handle;
char query[MAX_QUERY_SIZE];

char NullString[10] = "null";

void
wh_mysql_connect(void) {
  mysql_init(&mysql_handle);
  unsigned int timeout = 20;
  mysql_options(&mysql_handle, MYSQL_OPT_CONNECT_TIMEOUT, (const char *) &timeout);
  if (!mysql_real_connect(&mysql_handle, host, userid, password, database_name, server_port,  NULL, 0)) {
    fprintf(stderr, "Connection to MySQL server failed (host=%s userid=%s password=%s database_name=%s server_port=%d timeout=%d).\nError: %s\n",
	    host, userid, password, database_name, server_port, timeout, mysql_error(&mysql_handle));
    exit(-1);
  }
  if(-1 == snprintf(query, MAX_QUERY_SIZE, "Set AUTOCOMMIT=0"))
	wh_mysql_fatal_error("In wh_mysql_connect: Query string too long\n");
  wh_mysql_run_query(query, "Failed to turn auto commit off");
}

void
wh_mysql_finalize_dataset() {
  // Set ChangeDate to a nonNULL value to indicate normal termination of the loader.
  if (-1 == snprintf(query, MAX_QUERY_SIZE,
		     "UPDATE DataSet SET ChangeDate = NOW() WHERE WID = %d",
		     dataset_wid) >= MAX_QUERY_SIZE)
    wh_mysql_fatal_error("In wh_mysql_finalize_dataset : Query string too long\n");
    
  wh_mysql_run_query(query,"Failed to update DataSet table");
}

void
wh_mysql_commit(void) {
  if(-1 == snprintf(query, MAX_QUERY_SIZE, "COMMIT"))
	wh_mysql_fatal_error("In wh_mysql_commit: Query string too long\n");
  wh_mysql_run_query(query, "Failed to commit");
}

void
wh_mysql_disconnect(void) {
  wh_mysql_commit();   /* Added 11/22/2003 TJL */
  mysql_close(&mysql_handle);
}

char *
wh_mysql_malloc_and_create_escaped_string(char * from_string) {
 /* Maximum length of the escaped string assuming each character needs to be escaped.
    NEEDS TO BE FREED FROM THE FUNCTION THAT IT GETS USED FROM */
 char * escaped_string;
 char * temp_escaped_string;

 if(from_string == NULL)
 {
   return NULL;
 }
 escaped_string = (char *)(malloc(strlen(from_string) * sizeof(char) * 2 + 1));
 temp_escaped_string = (char *)(malloc(strlen(from_string) * sizeof(char) * 2 + 1));
 if((!escaped_string) || (!temp_escaped_string))
 {
    wh_mysql_fatal_error(" Error: Failed to allocate memory for escaped string\n");
 }
 mysql_real_escape_string(&mysql_handle, temp_escaped_string, from_string, strlen(from_string));
 sprintf(escaped_string, "'%s'", temp_escaped_string);
 free(temp_escaped_string);
 return escaped_string;
}

void
reset_wid_tables() {
  // Called when WIDTable and/or SpecialWIDTable have been corrupted, which
  // seems to happen when the mysql server goes down.

  MYSQL_RES *result;
  MYSQL_ROW row;
  int wid;

  printf("Resetting WID tables to legal values\n");  
  
  // Determine likely legal value for next special WID, and insert into SpecialWIDTable
  snprintf(query, MAX_QUERY_SIZE, "SELECT max(WID) FROM DataSet");
  wh_mysql_run_query(query, "Failed to get max WID from DataSet table\n");
  result = mysql_store_result(&mysql_handle);
  wid = 0;
  if (result && (row = mysql_fetch_row(result)) && row[0]) {
    errno = 0;
    wid = (int) strtol(row[0], NULL, 10);
    if (errno) wid = 0;
  }
  mysql_free_result(result);
  if (wid <= 1) {
    printf("Reset of wid table failed; WIDTable and SpecialWIDTable must be repaired manually");
    exit(-1);
  }
  snprintf(query, MAX_QUERY_SIZE, "ALTER TABLE SpecialWIDTable AUTO_INCREMENT=%d", wid+2);
  wh_mysql_run_query(query, "Failed to fix SpecialWIDTable\n");
  printf(" SpecialWIDTable counter set to %d\n", wid+2);

  // Determine likely legal value for next (regular) WID, and insert into WIDTable
  snprintf(query, MAX_QUERY_SIZE, "SELECT max(OtherWID) FROM Entry");
  wh_mysql_run_query(query, "Failed to get max OtherWID from Entry table\n");
  result = mysql_store_result(&mysql_handle);
  wid = 0;
  if (result && (row = mysql_fetch_row(result)) && row[0]) {
    errno = 0;
    wid = (int) strtol(row[0], NULL, 10);
    if (errno) wid = 0;
  }
  mysql_free_result(result);
  if (wid <= 1) {
    printf("Reset of wid table failed; WIDTable and SpecialWIDTable must be repaired manually");
    exit(-1);
  }
  snprintf(query, MAX_QUERY_SIZE, "ALTER TABLE WIDTable AUTO_INCREMENT=%d", wid*2);// Add factor of 2 to be conservative
  wh_mysql_run_query(query, "Failed to fix WIDTable\n");
  printf(" WIDTable counter set to %d\n", wid*2);
}

int
wh_mysql_get_new_wid(void) {
  int new_wid;
  if(-1 == snprintf(query, MAX_QUERY_SIZE, "DELETE FROM WIDTable"))
	wh_mysql_fatal_error("In wh_mysql_get_new_wid: Query string too long\n");
  wh_mysql_run_query(query, "Failed to delete from widtable");
  if(-1 == snprintf(query, MAX_QUERY_SIZE,"INSERT INTO WIDTable () Values ()"))
	wh_mysql_fatal_error("In wh_mysql_get_new_wid: Query string too long\n");
  wh_mysql_run_query(query, "Failed to insert into WIDTable to get new wid"); 
  new_wid = mysql_insert_id(&mysql_handle);
  
  /* Show progress */
  if (new_wid % 1000 == 0) {
    printf(".");
    fflush(stdout);
    fflush(stderr);
  }

  /* Optionally commit at intervals */
  if (commit_freq && (new_wid % commit_freq == 0)) {
    wh_mysql_commit();
  }

  return(new_wid);
}

int
wh_mysql_get_new_special_wid(void) {
  int new_wid;
  int max_special_wid = 0; 
  MYSQL_RES *result;
  MYSQL_ROW row;
  
  if(-1 == snprintf(query,MAX_QUERY_SIZE,"SELECT MaxSpecialWID FROM Warehouse"))
	wh_mysql_fatal_error("In wh_mysql_get_new_special_wid: Query string too long\n");
  wh_mysql_run_query(query, "Failed to get MaxSpecialWID");

  result = mysql_store_result(&mysql_handle);
  row = mysql_fetch_row(result);
  if(row)
  {
    max_special_wid = atoi((row[0]));
  } 

  if(-1 == snprintf(query, MAX_QUERY_SIZE, "DELETE FROM SpecialWIDTable"))
	wh_mysql_fatal_error("In wh_mysql_get_new_special_wid: Query string too long\n");
  wh_mysql_run_query(query, "Failed to delete from SpecialWIDtable");

  if(-1 == snprintf(query, MAX_QUERY_SIZE,"INSERT INTO SpecialWIDTable () Values ()"))
	wh_mysql_fatal_error("In wh_mysql_get_new_special_wid: Query string too long\n");
  wh_mysql_run_query(query, "Failed to insert into SpecialWIDTable to get new wid");

  new_wid = mysql_insert_id(&mysql_handle);
  if (0) printf(" The new wid inserted in the special table is %d\n", new_wid);

  if (new_wid <= 1) {
    // This is an error indicator, usually indicating the mysql server has gone down,
    //  and the WID counters in the WIDTable and SpecialWIDTable are no longer valid.
    // So reset these counters to a value very like to be higher than any previously
    // allocated WIDs.
    reset_wid_tables();
    return wh_mysql_get_new_special_wid();
  } else if (new_wid > max_special_wid)
    return( wh_mysql_get_new_wid() );

  return new_wid;
}

void
wh_mysql_insert_crossreference(int other_wid, char * xid, char * dataset_name) {
  char * escaped_xid = wh_mysql_malloc_and_create_escaped_string(xid);
  char * escaped_dataset_name = wh_mysql_malloc_and_create_escaped_string(dataset_name);
    
  if (-1 == snprintf(query, MAX_QUERY_SIZE,
		     "INSERT INTO CrossReference (OtherWID, XID, DataBaseName) Values(%d, %s, %s)",
		     other_wid, escaped_xid, escaped_dataset_name))
    wh_mysql_fatal_error("In wh_mysql_insert_crossreference: Query string too long\n");

  free(escaped_dataset_name);
  free(escaped_xid);
  wh_mysql_run_query(query, "Failed to insert into CrossReference table.");
}

void
wh_mysql_insert_entry(int object_wid, int load_error, int lineno) {
  /* Add object_wid to Entry table */
  if (load_error) {
    if(-1 == snprintf(query,MAX_QUERY_SIZE,
		      "INSERT INTO Entry (OtherWID, InsertDate, LoadError, LineNumber, DatasetWID) Values( %d , Now(), 'T', %d , %d )",
		      object_wid, lineno , dataset_wid))
	wh_mysql_fatal_error("In wh_mysql_insert_entry: Query string too long\n");
  }
  else {
    if(-1 == snprintf(query, MAX_QUERY_SIZE,
		      "INSERT INTO Entry (OtherWID, InsertDate, LoadError, LineNumber, DatasetWID) Values( %d , Now(), 'F' , %d , %d )",
		      object_wid, lineno , dataset_wid))
	wh_mysql_fatal_error("In wh_mysql_insert_entry: Query string too long\n");
  }
  wh_mysql_run_query(query, "Failed to insert into Entry Table.");
}

void
wh_mysql_insert_citation(int citation_wid, char * fulltext,
			  char * pubmed_id, short pubmed_id_ind) {
  char * escaped_fulltext = wh_mysql_malloc_and_create_escaped_string(fulltext);
  char * escaped_pubmed_id = wh_mysql_malloc_and_create_escaped_string(pubmed_id);
  short fulltext_ind = INDICATE_NONNULL;
  if (!fulltext || (0 == strlen(fulltext)))
    fulltext_ind = INDICATE_NULL;
    
  if (-1 == snprintf(query, MAX_QUERY_SIZE,
		     "INSERT INTO Citation (WID, Citation, PMID, DataSetWID) Values(%d, %s, %s, %d)",
		     citation_wid,
		     wh_mysql_string_column(escaped_fulltext, fulltext_ind),
		     wh_mysql_string_column(escaped_pubmed_id, pubmed_id_ind),
		     dataset_wid))
    wh_mysql_fatal_error("In wh_mysql_insert_citation: Query string too long\n");

  free(escaped_pubmed_id);
  free(escaped_fulltext);
  wh_mysql_run_query(query, "Failed to insert into Citation table.");
}

void
wh_mysql_update_citation(int citation_wid, char * fulltext,
			  char * pubmed_id, short pubmed_id_ind) {
  char * escaped_fulltext = wh_mysql_malloc_and_create_escaped_string(fulltext);
  char * escaped_pubmed_id = wh_mysql_malloc_and_create_escaped_string(pubmed_id);
  if (-1 == snprintf(query, MAX_QUERY_SIZE,
		     "UPDATE Citation SET Citation = %s, PMID=%s WHERE WID = %d",
		     escaped_fulltext,
		     wh_mysql_string_column(escaped_pubmed_id, pubmed_id_ind),
		     citation_wid) >= MAX_QUERY_SIZE)
    wh_mysql_fatal_error("In wh_mysql_update_citation: Query string too long\n");
    
  free(escaped_pubmed_id);
  free(escaped_fulltext);
  wh_mysql_run_query(query,"Failed to update Citation table");
}

void
wh_mysql_insert_description(int object_wid, char* table_name, char * description) {
  char * escaped_string;
  if (description && (0 != strlen(description)) ) {
    escaped_string = wh_mysql_malloc_and_create_escaped_string(description);
    if (-1 == snprintf(query, MAX_QUERY_SIZE,
                       "INSERT INTO Description (OtherWID, TableName, Comm) VALUES(%d, '%s', %s)",
                       object_wid, table_name, escaped_string ))
      wh_mysql_fatal_error("In wh_mysql_insert_description: Query string too long\n");
    free(escaped_string);
    wh_mysql_run_query(query, "Failed to insert into Description.");
  }
}

void
wh_mysql_insert_comment(int object_wid, char* comment) {
  char * escaped_string;
  if (comment && (0 != strlen(comment)) && (0 != strcmp(comment, "-")) ) {
    escaped_string = wh_mysql_malloc_and_create_escaped_string(comment);
    if (-1 == snprintf(query,MAX_QUERY_SIZE,
                       "INSERT INTO CommentTable (OtherWID, Comm) VALUES(%d, %s)",
                       object_wid, escaped_string ))
      wh_mysql_fatal_error("In wh_mysql_insert_comment: Query string too long\n");
    free(escaped_string);
    wh_mysql_run_query(query, "Failed to insert into CommentTable.");
  }
}

void
wh_mysql_insert_dbid(int other_wid,  char * xid) {
  char * escaped_xid = wh_mysql_malloc_and_create_escaped_string(xid);
  if(-1 == snprintf(query, MAX_QUERY_SIZE,"INSERT INTO DBID (OtherWID, XID) Values(%d, %s)",
		    other_wid, escaped_xid))
	wh_mysql_fatal_error("In wh_mysql_insert_into_dbid: Query string too long\n");
  wh_mysql_run_query(query, "Failed to insert into DBID table.");
  free(escaped_xid);
}

void
wh_mysql_insert_synonymtable(int wid, char * syn_name) {
  char * escaped_syn_name = wh_mysql_malloc_and_create_escaped_string(syn_name);

  if (-1 == snprintf(query, MAX_QUERY_SIZE,
		     "INSERT INTO SynonymTable (OtherWID, Syn) Values(%d, %s)",
		     wid, escaped_syn_name))
	wh_mysql_fatal_error("In wh_mysql_insert_into_synonymtable:Query string too long\n");

  free(escaped_syn_name);
  wh_mysql_run_query(query, "Failed to insert into synonym table.");
}


void
wh_mysql_insert_chemical_simple(int chemical_wid, char * name) {
  char * escaped_name = wh_mysql_malloc_and_create_escaped_string(name);
    
  if (-1 == snprintf(query, MAX_QUERY_SIZE,
		     "INSERT INTO Chemical (WID, Name, DataSetWID) Values(%d, %s, %d)",
		     chemical_wid, escaped_name, dataset_wid))
	wh_mysql_fatal_error("In wh_mysql_insert_chemical_simple: Query string too long\n");

  free(escaped_name);
  wh_mysql_run_query(query, "Failed to insert into Chemical table.");
}

void
wh_mysql_insert_enzymaticreaction(int enzrxn_wid, int reaction_wid, int protein_wid,
				  int complex_wid, short complex_ind,
				  char * reaction_direction, short reaction_direction_ind) {
  char complex_buf[MAX_FIXED_COLUMN_SIZE]; 
  char * escaped_reaction_direction = wh_mysql_malloc_and_create_escaped_string(reaction_direction);

  if (snprintf(query,MAX_QUERY_SIZE,
	       "INSERT INTO EnzymaticReaction (WID, ReactionWID, ProteinWID, ComplexWID, ReactionDirection, DataSetWID) VALUES (%d, %d, %d, %s, %s, %d)",
	       enzrxn_wid,
	       reaction_wid,
	       protein_wid,
	       wh_mysql_int_column(complex_wid, complex_ind, complex_buf),
	       wh_mysql_string_column(escaped_reaction_direction, reaction_direction_ind),
	       dataset_wid) >= MAX_QUERY_SIZE)
    wh_mysql_fatal_error("In wh_mysql_insert_enzymaticreaction : Query string too long\n");
  wh_mysql_run_query(query, "Failed to insert into EnzymaticReaction Table");
}

void
wh_mysql_insert_enzrxn_cofactor(int enzrxn_wid, int chemical_wid,
				char prosthetic, short prosthetic_ind) {
  char prosthetic_buf[MAX_FIXED_COLUMN_SIZE]; 

  if (snprintf(query,MAX_QUERY_SIZE,
	       "INSERT INTO EnzReactionCofactor (EnzymaticReactionWID, ChemicalWID, Prosthetic) VALUES (%d, %d, %s)",
	       enzrxn_wid,
	       chemical_wid,
	       wh_mysql_char_column(prosthetic, prosthetic_ind, prosthetic_buf)) >= MAX_QUERY_SIZE)
	wh_mysql_fatal_error("In wh_mysql_insert_enzrxn_cofactor : Query string too long\n");
  wh_mysql_run_query(query, "Failed to insert into EnzReactionCofactor Table");
}
							   
void
wh_mysql_insert_enzrxn_chemical(int enzrxn_wid, int chemical_wid,
				char mechanism, char inhibit_or_activate, char physiologically_relevant) {
  if (snprintf(query,MAX_QUERY_SIZE,
	       "INSERT INTO EnzReactionInhibitorActivator (EnzymaticReactionWID, CompoundWID, Mechanism, InhibitOrActivate, PhysioRelevant) VALUES (%d, %d, '%c', '%c', '%c')",
	       enzrxn_wid,
	       chemical_wid,
	       mechanism,
	       inhibit_or_activate,
	       physiologically_relevant) >= MAX_QUERY_SIZE)
	wh_mysql_fatal_error("In wh_mysql_insert_enzrxn_chemical: Query string too long\n");
  wh_mysql_run_query(query, "Failed to insert into EnzReactionInhibitorActivator Table");
}

void
wh_mysql_insert_enzrxn_alternate(int enzrxn_wid, int primary_wid, int alternate_wid, char cofactor) {
  if(snprintf(query, MAX_QUERY_SIZE,
	      "INSERT INTO EnzReactionAltCompound (EnzymaticReactionWID, PrimaryWID, AlternativeWID, Cofactor) VALUES (%d, %d, %d, '%c')",
	      enzrxn_wid,
	      primary_wid,
	      alternate_wid,
	      cofactor) >= MAX_QUERY_SIZE)
	wh_mysql_fatal_error("In wh_mysql_insert_enzrxn_alternate : Query string too long\n");
  wh_mysql_run_query(query, "Failed to insert into EnzReactionAltCompound table");
}

void
wh_mysql_insert_support(int support_wid, int other_wid, char *type, short type_ind) {
  char *escaped_type;

  escaped_type =  wh_mysql_malloc_and_create_escaped_string(type);
  if(snprintf(query, MAX_QUERY_SIZE,
	      "INSERT INTO Support (WID, OtherWID, Type, DataSetWID) VALUES (%d, %d, %s, %d)",
	      support_wid,
	      other_wid,
	      wh_mysql_string_column(escaped_type, type_ind),
	      dataset_wid) >= MAX_QUERY_SIZE)
    wh_mysql_fatal_error("In wh_mysql_insert_support : Query string too long\n");
  free(escaped_type);
  wh_mysql_run_query(query, "Failed to insert into Support table");
}

void
wh_mysql_insert_transcription_unit_component(int tu_wid, int other_wid, char * type) {
  char *escaped_type;

  escaped_type =  wh_mysql_malloc_and_create_escaped_string(type);
  
  if(snprintf(query, MAX_QUERY_SIZE,
	      "INSERT INTO TranscriptionUnitComponent (TranscriptionUnitWID, OtherWID, Type) VALUES (%d, %d, %s)",
	      tu_wid,
	      other_wid,
	      wh_mysql_string_column(escaped_type, INDICATE_NONNULL)) >= MAX_QUERY_SIZE)
    wh_mysql_fatal_error("In wh_mysql_insert_support : Query string too long\n");
  free(escaped_type);
  wh_mysql_run_query(query, "Failed to insert into TranscriptionUnitComponent table");
}


void
wh_mysql_insert_relatedterm(int term_wid, int other_wid, char *relationship, short relationship_ind) {
  char *escaped_relationship;

  escaped_relationship =  wh_mysql_malloc_and_create_escaped_string(relationship);
  if(snprintf(query, MAX_QUERY_SIZE,
	      "INSERT INTO RelatedTerm (TermWID, OtherWID, Relationship) VALUES (%d, %d, %s)",
	      term_wid,
	      other_wid,
	      wh_mysql_string_column(escaped_relationship, relationship_ind)) >= MAX_QUERY_SIZE)
    wh_mysql_fatal_error("In wh_mysql_insert_relatedterm : Query string too long\n");
  free(escaped_relationship);
  wh_mysql_run_query(query, "Failed to insert into RelatedTerm table");
}

void
wh_mysql_insert_transcription_unit(int transcription_unit_wid, char * name) {
  char * escaped_name = wh_mysql_malloc_and_create_escaped_string(name);
    
  if (-1 == snprintf(query, MAX_QUERY_SIZE,
		     "INSERT INTO TranscriptionUnit (WID, Name, DataSetWID) Values(%d, %s, %d)",
		     transcription_unit_wid,
		     escaped_name,
		     dataset_wid))
    wh_mysql_fatal_error("In wh_mysql_insert_transcription_unit: Query string too long\n");

  free(escaped_name);
  wh_mysql_run_query(query, "Failed to insert into TranscriptionUnit table.");
}

void
wh_mysql_insert_location(int protein_wid, char *location) {
  char *escaped_location = wh_mysql_malloc_and_create_escaped_string(location);
  
  if(snprintf(query, MAX_QUERY_SIZE,
	      "INSERT INTO Location (ProteinWID, Location) VALUES (%d, %s)",
	      protein_wid,
	      escaped_location) >= MAX_QUERY_SIZE)
    wh_mysql_fatal_error("In wh_mysql_insert_location : Query string too long\n");
  free(escaped_location);
  wh_mysql_run_query(query, "Failed to insert into Location table");
}



/*** Selects ***/

int
wh_mysql_dataset_exists(int dswid) {
  MYSQL_RES *result;
  MYSQL_ROW row;
  int wid = 0;
  
  if(snprintf(query, MAX_QUERY_SIZE,
	      "SELECT max(WID) FROM DataSet WHERE WID = %d",
	      dswid) >= MAX_QUERY_SIZE)
    wh_mysql_fatal_error("In wh_mysql_dataset_exists: Query string too long\n");
  /* printf("Query=%s\n", query);*/
  wh_mysql_run_query(query, "Failed to get WID from DataSet table\n");
  result = mysql_store_result(&mysql_handle);
  if (result && (row = mysql_fetch_row(result)) && row[0]) {
    errno = 0;
    wid = (int) strtol(row[0], NULL, 10);
    if (errno) wid = 0;
  }
  mysql_free_result(result);
  return wid;
}

int
wh_mysql_select_dataset_wid(char * db_name) {
  /* Returns the maximum (most recently loaded) DataSet.WID for the given dataset.
     Return 0 if not found.
  */
  MYSQL_RES *result;
  MYSQL_ROW row;
  int wid = 0;
  char * escaped_name = wh_mysql_malloc_and_create_escaped_string(db_name);
  
  if(snprintf(query, MAX_QUERY_SIZE,
	      "SELECT max(WID) FROM DataSet WHERE Name = %s",
	      escaped_name) >= MAX_QUERY_SIZE)
    wh_mysql_fatal_error("In wh_mysql_select_dataset_wid: Query string too long\n");
  free(escaped_name);
  /* printf("Query=%s\n", query);*/
  wh_mysql_run_query(query, "Failed to get WID from DataSet table\n");
  result = mysql_store_result(&mysql_handle);
  if (result && (row = mysql_fetch_row(result)) && row[0]) {
    errno = 0;
    wid = (int) strtol(row[0], NULL, 10);
    if (errno) wid = 0;
  }
  mysql_free_result(result);
  return wid;
}

int
wh_mysql_select_ecnumber_reaction(char * ecnumber, int enzyme_dataset_wid) {
  /* Returns the maximum (most recently loaded) Reaction.WID for the given Reaction.ECNumber,
     from the Enzyme database. Return 0 if not found.
  */
  MYSQL_RES *result;
  MYSQL_ROW row;
  int wid = 0;
  char * escaped_ec = wh_mysql_malloc_and_create_escaped_string(ecnumber);

  if(snprintf(query, MAX_QUERY_SIZE,
	      "SELECT max(WID) FROM Reaction WHERE ECNumber = %s AND DataSetWID = %d",
	      escaped_ec, enzyme_dataset_wid) >= MAX_QUERY_SIZE)
    wh_mysql_fatal_error("In wh_mysql_select_ecnumber_reaction: Query string too long\n");
  /* printf("wh_mysql_select_ecnumber_reaction Query=%s\n", query); */
  free(escaped_ec);
  wh_mysql_run_query(query, "Failed to get WID from Reaction table\n");
  result = mysql_store_result(&mysql_handle);
  if (result && (row = mysql_fetch_row(result)) && row[0]) {
    errno = 0;
    wid = (int) strtol(row[0], NULL, 10);
    if (errno) wid = 0;
  }
  mysql_free_result(result);
  return wid;
}

int
wh_mysql_select_geneticcode_wid (int gencode_id) {
  //!! Recent versions of KEGG don't provide a gencode to look up, so stub this out for now.
  return 0;
}

int
wh_mysql_select_term(char * term_name, int query_dataset_wid) {
  /* Returns the maximum (most recently loaded) Term.WID for the given Term.Name,
     from the given. Return 0 if not found.
  */
  MYSQL_RES *result;
  MYSQL_ROW row;
  int wid = 0;
  char * escaped_term = wh_mysql_malloc_and_create_escaped_string(term_name);

  if(snprintf(query, MAX_QUERY_SIZE,
	      "SELECT max(WID) FROM Term WHERE Name = %s AND DataSetWID = %d",
	      escaped_term, query_dataset_wid) >= MAX_QUERY_SIZE)
    wh_mysql_fatal_error("In wh_mysql_select_term: Query string too long\n");
  /* printf("wh_mysql_select_term Query=%s\n", query); */
  free(escaped_term);
  wh_mysql_run_query(query, "Failed to get WID from Term table\n");
  result = mysql_store_result(&mysql_handle);
  if (result && (row = mysql_fetch_row(result)) && row[0]) {
    errno = 0;
    wid = (int) strtol(row[0], NULL, 10);
    if (errno) wid = 0;
  }
  mysql_free_result(result);
  return wid;
}
                                                                

int
wh_mysql_select_dbid_otherwid(char * xid, int query_dataset_wid) {
  MYSQL_RES *result;
  MYSQL_ROW row;
  int wid = 0;
  char * escaped_xid = wh_mysql_malloc_and_create_escaped_string(xid);

  if(snprintf(query, MAX_QUERY_SIZE,
	      "SELECT max(DBID.OtherWID) FROM DBID,Entry WHERE DBID.XID = %s AND DBID.OtherWID = Entry.OtherWID AND Entry.DataSetWID = %d",
	      escaped_xid, query_dataset_wid) >= MAX_QUERY_SIZE)
    wh_mysql_fatal_error("In wh_mysql_select_term: Query string too long\n");
  free(escaped_xid);
  wh_mysql_run_query(query, "Failed to get OtherWID from DBID table\n");
  result = mysql_store_result(&mysql_handle);
  if (result && (row = mysql_fetch_row(result)) && row[0]) {
    errno = 0;
    wid = (int) strtol(row[0], NULL, 10);
    if (errno) wid = 0;
  }
  mysql_free_result(result);
  return wid;
}


void
wh_mysql_run_query(char * sql_query, char * error_string) {
  /*printf("Query string formed is %s \n", sql_query);*/
  if (mysql_real_query(&mysql_handle, sql_query, (unsigned int) strlen(sql_query)))
  {
    fprintf(stderr,"%s Error:%s\n",error_string,mysql_error(&mysql_handle));
  }
}

char *
wh_mysql_int_column(int num, short indicator, char * buffer)
{
   if(indicator >= 0)
   {
   	sprintf(buffer,"%d", num);
   }
   else
   {
        sprintf(buffer,"%s", NullString);
   }
   return buffer;
}

char *
wh_mysql_double_column(double num, short indicator, char * buffer)
{
   if(indicator >= 0)
   {
        sprintf(buffer,"%f", num);
   }
   else
   {
        sprintf(buffer,"%s", NullString);
   }
   return buffer;
}

char *
wh_mysql_char_column(char c, short indicator, char * buffer)
{
   if(indicator >= 0)
   {
        sprintf(buffer,"'%c'", c);
   }
   else
   {
        sprintf(buffer,"%s", NullString);
   }
   return buffer;
}

char *
wh_mysql_string_column(char * string, short indicator)
{
   if(indicator >= 0)
   {
     return string; 
   }
   else
   {
     return NullString;
   }
}

/*
 * fatal errors (e.g. out of memory)
 */
void
wh_mysql_fatal_error(char *message) {
  perror(message);
  exit(-1);
}

