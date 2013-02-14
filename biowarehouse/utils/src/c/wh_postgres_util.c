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
//#include "mysql.h"
#include <libpq-fe.h>
#include "stdio.h"
#include <errno.h>
#include "wh.h"
#include "wh_postgres_util.h"

#define MAX_QUERY_SIZE 2000000

PGconn* pg_handle;
char query[MAX_QUERY_SIZE];

char NullString[10] = "null";

void
wh_postgres_connect(void) {
  //mysql_init(&pg_handle);
  unsigned int timeout = 20;
  
  char _server_port[255];
  sprintf(_server_port, "%d", server_port);
  
  char* keys[] = { "host", "port", "dbname", "user", "password", NULL };
  char* values[] = {host, _server_port, database_name, userid, password, NULL};
  printf("%s, %s, %s, %s, %s\n", host, _server_port, database_name, userid, password);
  pg_handle = PQconnectdbParams(keys, values, 0);
  printf("connection\n");
  //pg_handle = PQsetdblogin(host, server_port, NULL, NULL, database_name, userid, password);
  //mysql_options(&pg_handle, MYSQL_OPT_CONNECT_TIMEOUT, (const char *) &timeout);
  if (PQstatus(pg_handle) != CONNECTION_OK) {
    fprintf(stderr, "Connection to Postgres server failed (host=%s userid=%s password=%s database_name=%s server_port=%d timeout=%d).\nError: %s\n",
	    host, userid, password, database_name, server_port, timeout, PQerrorMessage(pg_handle));
    exit(-1);
  }
  //if(-1 == snprintf(query, MAX_QUERY_SIZE, "Set AUTOCOMMIT=0"))
	//wh_postgres_fatal_error("In wh_postgres_connect: Query string too long\n");
  //wh_postgres_run_query(query, "Failed to turn auto commit off");
}

void
wh_postgres_finalize_dataset() {
  // Set ChangeDate to a nonNULL value to indicate normal termination of the loader.
  if (-1 == snprintf(query, MAX_QUERY_SIZE,
		     "UPDATE \"DataSet\" SET \"ChangeDate\" = NOW() WHERE \"WID\" = %d",
		     dataset_wid) >= MAX_QUERY_SIZE)
    wh_postgres_fatal_error("In wh_postgres_finalize_dataset : Query string too long\n");
    
  wh_postgres_run_query(query,"Failed to update DataSet table");
}

void
wh_postgres_commit(void) {
  /*if(-1 == snprintf(query, MAX_QUERY_SIZE, "COMMIT"))
	wh_postgres_fatal_error("In wh_postgres_commit: Query string too long\n");
  wh_postgres_run_query(query, "Failed to commit");*/
}

void
wh_postgres_disconnect(void) {
  PGresult   *res;
  res = PQexec(pg_handle, "END");
  PQfinish(pg_handle);
}

char *
wh_postgres_malloc_and_create_escaped_string(char * from_string) {
 /* Maximum length of the escaped string assuming each character needs to be escaped.
    NEEDS TO BE FREED FROM THE FUNCTION THAT IT GETS USED FROM */
 char * escaped_string;
 char * temp_escaped_string;

 if(from_string == NULL)
 {
   return NULL;
 }
 escaped_string = (char *)(malloc(strlen(from_string) * sizeof(char) * 2 + 1));
 if((!escaped_string))
 {
    wh_postgres_fatal_error(" Error: Failed to allocate memory for escaped string\n");
 }
 temp_escaped_string = PQescapeLiteral(pg_handle, from_string, strlen(from_string));
 //mysql_real_escape_string(&pg_handle, temp_escaped_string, from_string, strlen(from_string));
 sprintf(escaped_string, "%s", temp_escaped_string);
 PQfreemem(temp_escaped_string);
 return escaped_string;
}

void
reset_wid_tables() {
  // Called when WIDTable and/or SpecialWIDTable have been corrupted, which
  // seems to happen when the mysql server goes down.
/*
  PGresult *result;
  MYSQL_ROW row;
  int wid;

  printf("Resetting WID tables to legal values\n");  
  
  // Determine likely legal value for next special WID, and insert into SpecialWIDTable
  snprintf(query, MAX_QUERY_SIZE, "SELECT max(WID) FROM DataSet");
  wh_postgres_run_query(query, "Failed to get max WID from DataSet table\n");
  result = mysql_store_result(&pg_handle);
  wid = 0;
  if (result && (row = mysql_fetch_row(result)) && row[0]) {
    errno = 0;
    wid = (int) strtol(row[0], NULL, 10);
    if (errno) wid = 0;
  }
  PQclear(result);
  if (wid <= 1) {
    printf("Reset of wid table failed; WIDTable and SpecialWIDTable must be repaired manually");
    exit(-1);
  }
  snprintf(query, MAX_QUERY_SIZE, "ALTER TABLE SpecialWIDTable AUTO_INCREMENT=%d", wid+2);
  wh_postgres_run_query(query, "Failed to fix SpecialWIDTable\n");
  printf(" SpecialWIDTable counter set to %d\n", wid+2);

  // Determine likely legal value for next (regular) WID, and insert into WIDTable
  snprintf(query, MAX_QUERY_SIZE, "SELECT max(OtherWID) FROM Entry");
  wh_postgres_run_query(query, "Failed to get max OtherWID from Entry table\n");
  result = mysql_store_result(&pg_handle);
  wid = 0;
  if (result && (row = mysql_fetch_row(result)) && row[0]) {
    errno = 0;
    wid = (int) strtol(row[0], NULL, 10);
    if (errno) wid = 0;
  }
  PQclear(result);
  if (wid <= 1) {
    printf("Reset of wid table failed; WIDTable and SpecialWIDTable must be repaired manually");
    exit(-1);
  }
  snprintf(query, MAX_QUERY_SIZE, "ALTER TABLE WIDTable AUTO_INCREMENT=%d", wid*2);// Add factor of 2 to be conservative
  wh_postgres_run_query(query, "Failed to fix WIDTable\n");
  printf(" WIDTable counter set to %d\n", wid*2);*/
}

int
wh_postgres_get_new_wid(void) {
  int new_wid;
  if(-1 == snprintf(query, MAX_QUERY_SIZE, "DELETE FROM \"WIDTable\""))
	wh_postgres_fatal_error("In wh_postgres_get_new_wid: Query string too long\n");
  wh_postgres_run_query(query, "Failed to delete from widtable");
 
  if (-1 == snprintf(query, MAX_QUERY_SIZE, "SELECT nextval('\"WIDTable_PreviousWID_seq\"'::regclass)"))
	wh_postgres_fatal_error("In wh_postgres_get_new_wid: Query string too long\n");
  PGresult* res = wh_postgres_run_query(query, "Failed to insert into WIDTable to get new wid"); 
  
  if (wh_postgres_num_rows(res) == 0) _asm int 3;
  char* val = PQgetvalue(res, 0, 0);
  new_wid = atoi(val);
  //printf("NEW WID %d ", new_wid);
  if(-1 == snprintf(query, MAX_QUERY_SIZE,"INSERT INTO \"WIDTable\" Values (%d)", new_wid))
     wh_postgres_fatal_error("In wh_postgres_get_new_wid: Query string too long\n");
  wh_postgres_run_query(query, "Failed to insert into WIDTable to get new wid"); 
  //new_wid = mysql_insert_id(&pg_handle);
  
  //printf("The new wid inserted in the special table is %d\n", new_wid);
  
  /* Show progress */
  if (new_wid % 1000 == 0) {
    printf(".");
    fflush(stdout);
    fflush(stderr);
  }

  /* Optionally commit at intervals */
  if (commit_freq && (new_wid % commit_freq == 0)) {
    wh_postgres_commit();
  }

  return(new_wid);
}

int
wh_postgres_get_new_special_wid(void) {
  int new_wid;
  int max_special_wid = 0; 
  PGresult *result;
  
  if(-1 == snprintf(query,MAX_QUERY_SIZE,"SELECT \"MaxSpecialWID\" FROM \"Warehouse\""))
	wh_postgres_fatal_error("In wh_postgres_get_new_special_wid: Query string too long\n");
  result = wh_postgres_run_query(query, "Failed to get MaxSpecialWID");

  if (wh_postgres_num_rows(result) == 0) _asm int 3;
  char* val = PQgetvalue(result, 0, 0);
  max_special_wid = atoi(val);
  /*result = mysql_store_result(&pg_handle);
  row = mysql_fetch_row(result);
  if(row)
  {
    max_special_wid = atoi((row[0]));
  } */

  if(-1 == snprintf(query, MAX_QUERY_SIZE, "DELETE FROM \"SpecialWIDTable\""))
	wh_postgres_fatal_error("In wh_postgres_get_new_special_wid: Query string too long\n");
  wh_postgres_run_query(query, "Failed to delete from SpecialWIDtable");

  //if(-1 == snprintf(query, MAX_QUERY_SIZE,"INSERT INTO \"SpecialWIDTable\" () \"Values\" ()"))
  if (-1 == snprintf(query, MAX_QUERY_SIZE, "SELECT nextval('\"SpecialWIDTable_PreviousWID_seq\"'::regclass)"))
	wh_postgres_fatal_error("In wh_postgres_get_new_special_wid: Query string too long\n");
  result = wh_postgres_run_query(query, "Failed to insert into SpecialWIDTable to get new wid");

  
  //new_wid = mysql_insert_id(&pg_handle);
  //result = wh_postgres_run_query("SELECT lastval()", "SELECT lastval failed2");
  if (wh_postgres_num_rows(result) == 0) _asm int 3;
  val = PQgetvalue(result, 0, 0);
  new_wid = atoi(val);

  if(-1 == snprintf(query, MAX_QUERY_SIZE,"INSERT INTO \"SpecialWIDTable\" VALUES (%d)", new_wid)) {
	wh_postgres_fatal_error("In wh_postgres_get_new_special_wid: Query string too long\n");
  }

  /*if(-1 == snprintf(query, MAX_QUERY_SIZE,"INSERT INTO \"WIDTable\" Values (%d)", new_wid))
     wh_postgres_fatal_error("In wh_postgres_get_new_wid: Query string too long\n");*/
  wh_postgres_run_query(query, "Failed to insert into WIDTable to get new wid");
  
  //printf(" The new wid inserted in the special table is %d\n", new_wid);

  if (new_wid <= 1) { 
    // This is an error indicator, usually indicating the mysql server has gone down,
    //  and the WID counters in the WIDTable and SpecialWIDTable are no longer valid.
    // So reset these counters to a value very like to be higher than any previously
    // allocated WIDs.
    reset_wid_tables();
    return wh_postgres_get_new_special_wid();
  } else if (new_wid > max_special_wid)
    return( wh_postgres_get_new_wid() );

  return new_wid;
}

void
wh_postgres_insert_crossreference(int other_wid, char * xid, char * dataset_name) {
  char * escaped_xid = wh_postgres_malloc_and_create_escaped_string(xid);
  char * escaped_dataset_name = wh_postgres_malloc_and_create_escaped_string(dataset_name);
    
  if (-1 == snprintf(query, MAX_QUERY_SIZE,
		     "INSERT INTO \"CrossReference\" (\"OtherWID\", \"XID\", \"DataBaseName\") Values(%d, %s, %s)",
		     other_wid, escaped_xid, escaped_dataset_name))
    wh_postgres_fatal_error("In wh_postgres_insert_crossreference: Query string too long\n");

  free(escaped_dataset_name);
  free(escaped_xid);
  wh_postgres_run_query(query, "Failed to insert into CrossReference table.");
}

void
wh_postgres_insert_entry(int object_wid, int load_error, int lineno) {
  /* Add object_wid to Entry table */
  if (load_error) {
    if(-1 == snprintf(query,MAX_QUERY_SIZE,
		      "INSERT INTO \"Entry\" (\"OtherWID\", \"InsertDate\", \"LoadError\", \"LineNumber\", \"DataSetWID\") Values( %d , Now(), 'T', %d , %d )",
		      object_wid, lineno , dataset_wid))
	wh_postgres_fatal_error("In wh_postgres_insert_entry: Query string too long\n");
  }
  else {
    if(-1 == snprintf(query, MAX_QUERY_SIZE,
		      "INSERT INTO \"Entry\" (\"OtherWID\", \"InsertDate\", \"LoadError\", \"LineNumber\", \"DataSetWID\") Values( %d , Now(), 'F' , %d , %d )",
		      object_wid, lineno , dataset_wid))
	wh_postgres_fatal_error("In wh_postgres_insert_entry: Query string too long\n");
  }
  wh_postgres_run_query(query, "Failed to insert into Entry Table.");
}

void
wh_postgres_insert_citation(int citation_wid, char * fulltext,
			  char * pubmed_id, short pubmed_id_ind) {
  char * escaped_fulltext = wh_postgres_malloc_and_create_escaped_string(fulltext);
  char * escaped_pubmed_id = wh_postgres_malloc_and_create_escaped_string(pubmed_id);
  short fulltext_ind = INDICATE_NONNULL;
  if (!fulltext || (0 == strlen(fulltext)))
    fulltext_ind = INDICATE_NULL;
    
  if (-1 == snprintf(query, MAX_QUERY_SIZE,
		     "INSERT INTO \"Citation\" (\"WID\", \"Citation\", \"PMID\", \"DataSetWID\") Values(%d, %s, %s, %d)",
		     citation_wid,
		     wh_postgres_string_column(escaped_fulltext, fulltext_ind),
		     wh_postgres_string_column(escaped_pubmed_id, pubmed_id_ind),
		     dataset_wid))
    wh_postgres_fatal_error("In wh_postgres_insert_citation: Query string too long\n");

  free(escaped_pubmed_id);
  free(escaped_fulltext);
  wh_postgres_run_query(query, "Failed to insert into Citation table.");
}

void
wh_postgres_update_citation(int citation_wid, char * fulltext,
			  char * pubmed_id, short pubmed_id_ind) {
  char * escaped_fulltext = wh_postgres_malloc_and_create_escaped_string(fulltext);
  char * escaped_pubmed_id = wh_postgres_malloc_and_create_escaped_string(pubmed_id);
  if (-1 == snprintf(query, MAX_QUERY_SIZE,
		     "UPDATE \"Citation\" SET \"Citation\" = %s, \"PMID\"=%s WHERE \"WID\" = %d",
		     escaped_fulltext,
		     wh_postgres_string_column(escaped_pubmed_id, pubmed_id_ind),
		     citation_wid) >= MAX_QUERY_SIZE)
    wh_postgres_fatal_error("In wh_postgres_update_citation: Query string too long\n");
    
  free(escaped_pubmed_id);
  free(escaped_fulltext);
  wh_postgres_run_query(query,"Failed to update Citation table");
}

void
wh_postgres_insert_description(int object_wid, char* table_name, char * description) {
  char * escaped_string;
  if (description && (0 != strlen(description)) ) {
    escaped_string = wh_postgres_malloc_and_create_escaped_string(description);
    if (-1 == snprintf(query, MAX_QUERY_SIZE,
                       "INSERT INTO \"Description\" (\"OtherWID\", \"TableName\", \"Comm\") VALUES(%d, '%s', %s)",
                       object_wid, table_name, escaped_string ))
      wh_postgres_fatal_error("In wh_postgres_insert_description: Query string too long\n");
    free(escaped_string);
    wh_postgres_run_query(query, "Failed to insert into Description.");
  }
}

void
wh_postgres_insert_comment(int object_wid, char* comment) {
  char * escaped_string;
  if (comment && (0 != strlen(comment)) && (0 != strcmp(comment, "-")) ) {
    escaped_string = wh_postgres_malloc_and_create_escaped_string(comment);
    if (-1 == snprintf(query,MAX_QUERY_SIZE,
                       "INSERT INTO \"CommentTable\" (\"OtherWID\", \"Comm\") VALUES(%d, %s)",
                       object_wid, escaped_string ))
      wh_postgres_fatal_error("In wh_postgres_insert_comment: Query string too long\n");
    free(escaped_string);
    wh_postgres_run_query(query, "Failed to insert into CommentTable.");
  }
}

void
wh_postgres_insert_dbid(int other_wid,  char * xid) {
  char * escaped_xid = wh_postgres_malloc_and_create_escaped_string(xid);
  if(-1 == snprintf(query, MAX_QUERY_SIZE,"INSERT INTO \"DBID\" (\"OtherWID\", \"XID\") Values(%d, %s)",
		    other_wid, escaped_xid))
	wh_postgres_fatal_error("In wh_postgres_insert_into_dbid: Query string too long\n");
  wh_postgres_run_query(query, "Failed to insert into DBID table.");
  free(escaped_xid);
}

void
wh_postgres_insert_synonymtable(int wid, char * syn_name) {
  char * escaped_syn_name = wh_postgres_malloc_and_create_escaped_string(syn_name);

  if (-1 == snprintf(query, MAX_QUERY_SIZE,
		     "INSERT INTO \"SynonymTable\" (\"OtherWID\", \"Syn\") Values(%d, %s)",
		     wid, escaped_syn_name))
	wh_postgres_fatal_error("In wh_postgres_insert_into_synonymtable:Query string too long\n");

  free(escaped_syn_name);
  wh_postgres_run_query(query, "Failed to insert into synonym table.");
}


void
wh_postgres_insert_chemical_simple(int chemical_wid, char * name) {
  char * escaped_name = wh_postgres_malloc_and_create_escaped_string(name);
    
  if (-1 == snprintf(query, MAX_QUERY_SIZE,
		     "INSERT INTO \"Chemical\" (\"WID\", \"Name\", \"DataSetWID\") Values(%d, %s, %d)",
		     chemical_wid, escaped_name, dataset_wid))
	wh_postgres_fatal_error("In wh_postgres_insert_chemical_simple: Query string too long\n");

  free(escaped_name);
  wh_postgres_run_query(query, "Failed to insert into Chemical table.");
}

void
wh_postgres_insert_enzymaticreaction(int enzrxn_wid, int reaction_wid, int protein_wid,
				  int complex_wid, short complex_ind,
				  char * reaction_direction, short reaction_direction_ind) {
  char complex_buf[MAX_FIXED_COLUMN_SIZE]; 
  char * escaped_reaction_direction = wh_postgres_malloc_and_create_escaped_string(reaction_direction);

  if (snprintf(query,MAX_QUERY_SIZE,
	       "INSERT INTO \"EnzymaticReaction\" (\"WID\", \"ReactionWID\", \"ProteinWID\", \"ComplexWID\", \"ReactionDirection\", \"DataSetWID\") VALUES (%d, %d, %d, %s, %s, %d)",
	       enzrxn_wid,
	       reaction_wid,
	       protein_wid,
	       wh_postgres_int_column(complex_wid, complex_ind, complex_buf),
	       wh_postgres_string_column(escaped_reaction_direction, reaction_direction_ind),
	       dataset_wid) >= MAX_QUERY_SIZE)
    wh_postgres_fatal_error("In wh_postgres_insert_enzymaticreaction : Query string too long\n");
  wh_postgres_run_query(query, "Failed to insert into EnzymaticReaction Table");
}

void
wh_postgres_insert_enzrxn_cofactor(int enzrxn_wid, int chemical_wid,
				char prosthetic, short prosthetic_ind) {
  char prosthetic_buf[MAX_FIXED_COLUMN_SIZE]; 

  if (snprintf(query,MAX_QUERY_SIZE,
	       "INSERT INTO \"EnzReactionCofactor\" (\"EnzymaticReactionWID\", \"ChemicalWID\", \"Prosthetic\") VALUES (%d, %d, %s)",
	       enzrxn_wid,
	       chemical_wid,
	       wh_postgres_char_column(prosthetic, prosthetic_ind, prosthetic_buf)) >= MAX_QUERY_SIZE)
	wh_postgres_fatal_error("In wh_postgres_insert_enzrxn_cofactor : Query string too long\n");
  wh_postgres_run_query(query, "Failed to insert into EnzReactionCofactor Table");
}
							   
void
wh_postgres_insert_enzrxn_chemical(int enzrxn_wid, int chemical_wid,
				char mechanism, char inhibit_or_activate, char physiologically_relevant) {
  if (snprintf(query,MAX_QUERY_SIZE,
	       "INSERT INTO \"EnzReactionInhibitorActivator\" (\"EnzymaticReactionWID\", \"CompoundWID\", \"Mechanism\", \"InhibitOrActivate\", \"PhysioRelevant\") VALUES (%d, %d, '%c', '%c', '%c')",
	       enzrxn_wid,
	       chemical_wid,
	       mechanism,
	       inhibit_or_activate,
	       physiologically_relevant) >= MAX_QUERY_SIZE)
	wh_postgres_fatal_error("In wh_postgres_insert_enzrxn_chemical: Query string too long\n");
  wh_postgres_run_query(query, "Failed to insert into EnzReactionInhibitorActivator Table");
}

void
wh_postgres_insert_enzrxn_alternate(int enzrxn_wid, int primary_wid, int alternate_wid, char cofactor) {
  if(snprintf(query, MAX_QUERY_SIZE,
	      "INSERT INTO \"EnzReactionAltCompound\" (\"EnzymaticReactionWID\", \"PrimaryWID\", \"AlternativeWID\", \"Cofactor\") VALUES (%d, %d, %d, '%c')",
	      enzrxn_wid,
	      primary_wid,
	      alternate_wid,
	      cofactor) >= MAX_QUERY_SIZE)
	wh_postgres_fatal_error("In wh_postgres_insert_enzrxn_alternate : Query string too long\n");
  wh_postgres_run_query(query, "Failed to insert into EnzReactionAltCompound table");
}

void
wh_postgres_insert_support(int support_wid, int other_wid, char *type, short type_ind) {
  char *escaped_type;

  escaped_type =  wh_postgres_malloc_and_create_escaped_string(type);
  if(snprintf(query, MAX_QUERY_SIZE,
	      "INSERT INTO \"Support\" (\"WID\", \"OtherWID\", \"Type\", \"DataSetWID\") VALUES (%d, %d, %s, %d)",
	      support_wid,
	      other_wid,
	      wh_postgres_string_column(escaped_type, type_ind),
	      dataset_wid) >= MAX_QUERY_SIZE)
    wh_postgres_fatal_error("In wh_postgres_insert_support : Query string too long\n");
  free(escaped_type);
  wh_postgres_run_query(query, "Failed to insert into Support table");
}

void
wh_postgres_insert_transcription_unit_component(int tu_wid, int other_wid, char * type) {
  char *escaped_type;

  escaped_type =  wh_postgres_malloc_and_create_escaped_string(type);
  
  if(snprintf(query, MAX_QUERY_SIZE,
	      "INSERT INTO \"TranscriptionUnitComponent\" (\"TranscriptionUnitWID\", \"OtherWID\", \"Type\") VALUES (%d, %d, %s)",
	      tu_wid,
	      other_wid,
	      wh_postgres_string_column(escaped_type, INDICATE_NONNULL)) >= MAX_QUERY_SIZE)
    wh_postgres_fatal_error("In wh_postgres_insert_support : Query string too long\n");
  free(escaped_type);
  wh_postgres_run_query(query, "Failed to insert into TranscriptionUnitComponent table");
}


void
wh_postgres_insert_relatedterm(int term_wid, int other_wid, char *relationship, short relationship_ind) {
  char *escaped_relationship;

  escaped_relationship =  wh_postgres_malloc_and_create_escaped_string(relationship);
  if(snprintf(query, MAX_QUERY_SIZE,
	      "INSERT INTO \"RelatedTerm\" (\"TermWID\", \"OtherWID\", \"Relationship\") VALUES (%d, %d, %s)",
	      term_wid,
	      other_wid,
	      wh_postgres_string_column(escaped_relationship, relationship_ind)) >= MAX_QUERY_SIZE)
    wh_postgres_fatal_error("In wh_postgres_insert_relatedterm : Query string too long\n");
  free(escaped_relationship);
  wh_postgres_run_query(query, "Failed to insert into RelatedTerm table");
}

void
wh_postgres_insert_transcription_unit(int transcription_unit_wid, char * name) {
  char * escaped_name = wh_postgres_malloc_and_create_escaped_string(name);
    
  if (-1 == snprintf(query, MAX_QUERY_SIZE,
		     "INSERT INTO \"TranscriptionUnit\" (\"WID\", \"Name\", \"DataSetWID\") Values(%d, %s, %d)",
		     transcription_unit_wid,
		     escaped_name,
		     dataset_wid))
    wh_postgres_fatal_error("In wh_postgres_insert_transcription_unit: Query string too long\n");

  free(escaped_name);
  wh_postgres_run_query(query, "Failed to insert into TranscriptionUnit table.");
}

void
wh_postgres_insert_location(int protein_wid, char *location) {
  char *escaped_location = wh_postgres_malloc_and_create_escaped_string(location);
  
  if(snprintf(query, MAX_QUERY_SIZE,
	      "INSERT INTO \"Location\" (\"ProteinWID\", \"Location\") VALUES (%d, %s)",
	      protein_wid,
	      escaped_location) >= MAX_QUERY_SIZE)
    wh_postgres_fatal_error("In wh_postgres_insert_location : Query string too long\n");
  free(escaped_location);
  wh_postgres_run_query(query, "Failed to insert into Location table");
}



/*** Selects ***/

int
wh_postgres_dataset_exists(int dswid) {
  PGresult *result;
  //MYSQL_ROW row;
  int wid = 0;
  
  if(snprintf(query, MAX_QUERY_SIZE,
	      "SELECT max(\"WID\") FROM \"DataSet\" WHERE \"WID\" = %d",
	      dswid) >= MAX_QUERY_SIZE)
    wh_postgres_fatal_error("In wh_postgres_dataset_exists: Query string too long\n");
  /* printf("Query=%s\n", query);*/
  result = wh_postgres_run_query(query, "Failed to get WID from DataSet table\n");

  if (wh_postgres_num_rows(result) == 0)
	wid = 0;
  else
	wid  = (int) strtol(PQgetvalue(result, 0, 0), NULL, 10);
  //result = mysql_store_result(&pg_handle);
  //if (result && (row = mysql_fetch_row(result)) && row[0]) {
  //  errno = 0;
  //  wid = (int) strtol(row[0], NULL, 10);
  //  if (errno) wid = 0;
  //}
  PQclear(result);
  return wid;
}

int
wh_postgres_select_dataset_wid(char * db_name) {
  /* Returns the maximum (most recently loaded) DataSet.WID for the given dataset.
     Return 0 if not found.
  */
  PGresult *result;
  //MYSQL_ROW row;
  int wid = 0;
  char * escaped_name = wh_postgres_malloc_and_create_escaped_string(db_name);
  
  if(snprintf(query, MAX_QUERY_SIZE,
	      "SELECT max(\"WID\") FROM \"DataSet\" WHERE \"Name\" = %s",
	      escaped_name) >= MAX_QUERY_SIZE)
    wh_postgres_fatal_error("In wh_postgres_select_dataset_wid: Query string too long\n");
  free(escaped_name);
  /* printf("Query=%s\n", query);*/
  //wh_postgres_run_query(query, "Failed to get WID from DataSet table\n");
  
  
  result = wh_postgres_run_query(query, "Failed to get WID from DataSet table\n");

  if (wh_postgres_num_rows(result) == 0)
	wid = 0;
  else
	wid  = (int) strtol(PQgetvalue(result, 0, 0), NULL, 10);
  
  
  
  /*result = mysql_store_result(&pg_handle);
  if (result && (row = mysql_fetch_row(result)) && row[0]) {
    errno = 0;
    wid = (int) strtol(row[0], NULL, 10);
    if (errno) wid = 0;
  }*/
  PQclear(result);
  return wid;
}

int
wh_postgres_select_ecnumber_reaction(char * ecnumber, int enzyme_dataset_wid) {
  /* Returns the maximum (most recently loaded) Reaction.WID for the given Reaction.ECNumber,
     from the Enzyme database. Return 0 if not found.
  */
  PGresult *result;
  //MYSQL_ROW row;
  int wid = 0;
  char * escaped_ec = wh_postgres_malloc_and_create_escaped_string(ecnumber);

  if(snprintf(query, MAX_QUERY_SIZE,
	      "SELECT max(\"WID\") FROM \"Reaction\" WHERE \"ECNumber\" = %s AND \"DataSetWID\" = %d",
	      escaped_ec, enzyme_dataset_wid) >= MAX_QUERY_SIZE)
    wh_postgres_fatal_error("In wh_postgres_select_ecnumber_reaction: Query string too long\n");
  /* printf("wh_postgres_select_ecnumber_reaction Query=%s\n", query); */
  free(escaped_ec);
  
  
  
  result = wh_postgres_run_query(query, "Failed to get WID from Reaction table\n");
    
  if (wh_postgres_num_rows(result) == 0)
	wid = 0;
  else
	wid  = (int) strtol(PQgetvalue(result, 0, 0), NULL, 10);
  
  
  
  /*wh_postgres_run_query(query, "Failed to get WID from Reaction table\n");
  result = mysql_store_result(&pg_handle);
  if (result && (row = mysql_fetch_row(result)) && row[0]) {
    errno = 0;
    wid = (int) strtol(row[0], NULL, 10);
    if (errno) wid = 0;
  }*/
  PQclear(result);
  return wid;
}

int
wh_postgres_select_geneticcode_wid (int gencode_id) {
  //!! Recent versions of KEGG don't provide a gencode to look up, so stub this out for now.
  return 0;
}

int
wh_postgres_select_term(char * term_name, int query_dataset_wid) {
  /* Returns the maximum (most recently loaded) Term.WID for the given Term.Name,
     from the given. Return 0 if not found.
  */
  PGresult *result;
  //MYSQL_ROW row;
  int wid = 0;
  char * escaped_term = wh_postgres_malloc_and_create_escaped_string(term_name);

  if(snprintf(query, MAX_QUERY_SIZE,
	      "SELECT max(\"WID\") FROM \"Term\" WHERE \"Name\" = %s AND \"DataSetWID\" = %d",
	      escaped_term, query_dataset_wid) >= MAX_QUERY_SIZE)
    wh_postgres_fatal_error("In wh_postgres_select_term: Query string too long\n");
  /* printf("wh_postgres_select_term Query=%s\n", query); */
  free(escaped_term);
  
  result = wh_postgres_run_query(query, "Failed to get WID from Term table\n");
  
  if (wh_postgres_num_rows(result) == 0)
	wid = 0;
  else
	wid  = (int) strtol(PQgetvalue(result, 0, 0), NULL, 10);
  
  /*wh_postgres_run_query(query, "Failed to get WID from Term table\n");
  result = mysql_store_result(&pg_handle);
  if (result && (row = mysql_fetch_row(result)) && row[0]) {
    errno = 0;
    wid = (int) strtol(row[0], NULL, 10);
    if (errno) wid = 0;
  }*/
  PQclear(result);
  return wid;
}
                                                                

int
wh_postgres_select_dbid_otherwid(char * xid, int query_dataset_wid) {
  PGresult *result;
  //MYSQL_ROW row;
  int wid = 0;
  char * escaped_xid = wh_postgres_malloc_and_create_escaped_string(xid);

  if(snprintf(query, MAX_QUERY_SIZE,
	      "SELECT max(\"DBID\".\"OtherWID\") FROM \"DBID\",\"Entry\" WHERE \"DBID\".\"XID\" = %s AND \"DBID\".\"OtherWID\" = \"Entry\".\"OtherWID\" AND \"Entry\".\"DataSetWID\" = %d",
	      escaped_xid, query_dataset_wid) >= MAX_QUERY_SIZE)
    wh_postgres_fatal_error("In wh_postgres_select_term: Query string too long\n");
  free(escaped_xid);
  
  
  result = wh_postgres_run_query(query, "Failed to get OtherWID from DBID table\n");
  
  if (wh_postgres_num_rows(result) == 0)
	wid = 0;
  else
	wid  = (int) strtol(PQgetvalue(result, 0, 0), NULL, 10);
  
  
  /*result = mysql_store_result(&pg_handle);
  if (result && (row = mysql_fetch_row(result)) && row[0]) {
    errno = 0;
    wid = (int) strtol(row[0], NULL, 10);
    if (errno) wid = 0;
  }*/
  PQclear(result);
  return wid;
}


PGresult*
wh_postgres_run_query(char * sql_query, char * error_string) {
  /*printf("Query string formed is %s \n", sql_query);*/
  PGresult* res = PQexec(pg_handle, sql_query);
  if (PQresultStatus(res) != PGRES_COMMAND_OK && 
    PQresultStatus(res) != PGRES_TUPLES_OK)
  {
    fprintf(stderr,"%s Error:%s\n",error_string,PQerrorMessage(pg_handle));
	PQfinish(pg_handle);
	_asm{ int 3 };
	exit(-1);
  }
  return res;
}

char *
wh_postgres_int_column(int num, short indicator, char * buffer)
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
wh_postgres_double_column(double num, short indicator, char * buffer)
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
wh_postgres_char_column(char c, short indicator, char * buffer)
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
wh_postgres_string_column(char * string, short indicator)
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
wh_postgres_fatal_error(char *message) {
  perror(message);
  exit(-1);
}

int
wh_postgres_num_rows(const PGresult *res) {
	int rows = PQntuples(res);
	//if (rows == 0) {
	//	_asm{ int 3 };
	//}
	return rows;
}
