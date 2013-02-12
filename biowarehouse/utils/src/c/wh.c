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
#include "wh.h"
#include <limits.h>
#include <math.h>
#include <signal.h>
#include <stdio.h>

struct tm time_buffer;  // used by parse_date
WIDTABLE undefined_words;

//// Common code for all C Biowarehouse loaders (Oracle and MySQL)

/*
 * Initialiaztion and error-handling
 */

void
handle_termination_signal();

void
handle_segv();

void
handle_illegal_instruction();

// Handler for signals, to flush debug output and commit to DB.
void
exit_gracefully() {
  fflush(stdout);
  fflush(stderr);
  wh_commit();
  exit(-2);
}

void
handle_segv() {
  fprintf(stderr, "\nA segmentation violation signal was caught, terminating\n");
  exit_gracefully();
}

void
handle_illegal_instruction() {
  fprintf(stderr, "\nAn illegal instruction signal was caught, terminating\n");
  exit_gracefully();
}

void
handle_termination_signal() {
  fprintf(stderr, "\nA termination signal was caught, terminating\n");
  exit_gracefully();
}

void
wh_init() {
#ifdef ORACLE
  server_port = 1523; //!! Currently not used
#else
  server_port = 3306;
#endif

  // Trap signals to flush buffers and commit DB
/*   signal(SIGSEGV, handle_segv); */
/*   signal(SIGILL, handle_illegal_instruction); */
/*   signal(SIGHUP, handle_termination_signal); */
/*   signal(SIGINT, handle_termination_signal); */
/*   signal(SIGQUIT, handle_termination_signal); */
}

/*
 * Utility Numeric Functions
 */

int
ifloor(double x, short* conversion_error) {
  double xint = floor(x);
  if ((xint < INT_MIN) || (xint > INT_MAX)) {
    *conversion_error = 1;
    return(-1);
  }
  return( (int) xint );
}

int
iceil(double x, short* conversion_error){
  double xint = ceil(x);
  if ((xint < INT_MIN) || (xint > INT_MAX)) {
    *conversion_error = 1;
    return(-1);
  }
  return( (int) xint );
}

                                           
/*
 * Utility Functions for handling column values
 */

double
number_column(char *columnval, short *indicator, short *error) {
  /* Set indicator depending on whether columnval is NULL.
     Convert columnval to a number, cause parse error if invalid.
  */
  short silent = *indicator;  // kludge to make error msg report optional
  double result;
  if (columnval && strlen(columnval) > 0) { // Bug 533, added length check
    char *junk;
    *indicator = INDICATE_NONNULL;
    errno = 0;
    result = strtod(columnval, &junk);
    if (errno != 0 || *junk != '\0') {
      *error = 1;
      *indicator = INDICATE_NULL;
      if (!silent)
        printf("*** Conversion error on %s at '%s'\n", columnval, junk);
      return(0.0);
      /* printf("error columnval=%s number_column=%f\n", columnval, result);*/
    }
    /*printf("number_column columnval=%s result=%f\n", columnval, result);*/
    return(result);
  }
  else {
    *indicator = INDICATE_NULL;
    return(0.0);
  }
}

char *
string_column(char *columnval, int maxlen, short *indicator) {
  /* Set indicator depending on whether columnval is NULL.
     If needed, duplicate and truncate columnval. Return columnval.
  */
  char * result = string_truncate(columnval, maxlen);
  if (columnval)
    *indicator = INDICATE_NONNULL;
  else
    *indicator = INDICATE_NULL;
  return(result);
}

char
boolean_column(char *columnval, short *indicator) {
  /* Set indicator depending on whether columnval is NULL.
     Return 'F' if columnval is missing, "NIL" or "NO". Else return 'T'.
  */
  if (!columnval) {
    *indicator = INDICATE_NULL;
    return 'F';  
  }
  *indicator = INDICATE_NONNULL;  /* 2_0_2 fix, was INDICATE_NULL */
  if (0 == strcasecmp(columnval, "NO")) return 'F';
  if (0 == strcasecmp(columnval, "NIL")) return 'F';
  return 'T';
}


/*
 * Utility functions common to all database loaders.
 * Generally these functions simply dispatch on the DBMS variants.
 * They may also perform some DBMS-independent semantics like
 * truncating strings to their max lengths specified in the schema.
 */

/* TJL: reduced this from 2000000 since a large value caused seg fault sometimes */
#define MAX_QUERY_SIZE 20000
#define MAX_FIXED_COLUMN_SIZE 100
#define MAX_SQL_NAME_SIZE 256

void
wh_insert_bare_object(char * object_table, int wid) {
  /* Inserts into an object table with the given name.
     Insert given WID and global dataset WID. */

  char query[MAX_QUERY_SIZE];
  
  if (snprintf(query, MAX_QUERY_SIZE,
	       "INSERT INTO \"%s\" (\"WID\", \"DataSetWID\") VALUES (%d, %d)",
	       object_table, wid, dataset_wid)  >= MAX_QUERY_SIZE )
    fatal_error("wh_insert_bare_object: Query too long\n");
  
#ifdef ORACLE
  wh_oracle_run_query(query, "Failed to insert into link table");
#elif DEF_MYSQL
  wh_mysql_run_query(query, "Failed to insert into link table");
#elif DEF_POSTGRES
  wh_postgres_run_query(query, "Failed to insert into link table");
#endif
}

void
wh_insert_linktable(char * column1, int wid1, char* column2, int wid2) {
  /* Inserts into a table named after its (only) two columns.
     Inserts two WID values into those columns */

  char query[MAX_QUERY_SIZE];
  char table_name[MAX_SQL_NAME_SIZE];

  /* Construct table name from column names */
  if (snprintf(table_name, MAX_SQL_NAME_SIZE,
	       "%s%s", column1, column2) >= MAX_SQL_NAME_SIZE )
    fatal_error("Table name too long\n");
  
  if (snprintf(query, MAX_QUERY_SIZE,
	       "INSERT INTO \"%s\" (\"%s\", \"%s\") VALUES (%d, %d)",
	       table_name, column1, column2, wid1, wid2)  >= MAX_QUERY_SIZE )
    fatal_error("Query too long\n");
  //_asm int 3;
#ifdef ORACLE
  wh_oracle_run_query(query, "Failed to insert into link table");
#elif DEF_MYSQL
  wh_mysql_run_query(query, "Failed to insert into link table");
#elif DEF_POSTGRES
  wh_postgres_run_query(query, "Failed to insert into link table");
#endif
}

void
wh_insert_named_linktable(char * table_name, char * column1, int wid1, char* column2, int wid2) {
  /* Inserts two WID values into named columns of table */

  char query[MAX_QUERY_SIZE];

  if (snprintf(query, MAX_QUERY_SIZE,
	       "INSERT INTO \"%s\" (\"%s\", \"%s\") VALUES (%d, %d)",
	       table_name, column1, column2, wid1, wid2)  >= MAX_QUERY_SIZE )
    fatal_error("Query too long\n");
  
#ifdef ORACLE
  wh_oracle_run_query(query, "Failed to insert into link table");
#elif DEF_MYSQL
  wh_mysql_run_query(query, "Failed to insert into link table");
#elif DEF_POSTGRES
  wh_postgres_run_query(query, "Failed to insert into link table");
#endif
}

void
wh_insert_subunit(int complex_wid, int subunit_wid, int coefficient) {
  /* Inserts a row into the Subunit table */

  char query[MAX_QUERY_SIZE];
  
  if (snprintf(query, MAX_QUERY_SIZE,
	       "INSERT INTO \"Subunit\" (\"ComplexWID\", \"SubunitWID\", \"Coefficient\") VALUES (%d, %d, %d)",
	       complex_wid, subunit_wid, coefficient)  >= MAX_QUERY_SIZE )
    fatal_error("wh_insert_subunit: Query too long\n");
  
#ifdef ORACLE
  wh_oracle_run_query(query, "Failed to insert into Subunit");
#elif DEF_MYSQL
  wh_mysql_run_query(query, "Failed to insert into Subunit");
#elif DEF_POSTGRES
  wh_postgres_run_query(query, "Failed to insert into Subunit");
#endif
}

void
wh_commit(void) {
#ifdef ORACLE
  wh_oracle_commit();
#elif DEF_MYSQL
  wh_mysql_commit();
#elif DEF_POSTGRES
  wh_postgres_commit();
#endif
}

void
wh_connect(void) {
#ifdef ORACLE
  wh_oracle_connect();
#elif DEF_MYSQL
  wh_mysql_connect();
#elif DEF_POSTGRES
  wh_postgres_connect();
#endif
}

void
wh_disconnect(void) {
#ifdef ORACLE
  wh_oracle_disconnect();
#elif DEF_MYSQL
  wh_mysql_disconnect();
#elif DEF_POSTGRES
  wh_postgres_disconnect();
#endif
  if (dataset_wid > 1)
    printf("Load is committed, disconnected from database, DataSet.WID is %d\n", dataset_wid);
  else 
    printf("Load is committed, disconnected from database\n");
}

void
wh_finalize_dataset(void) {
#ifdef ORACLE
  wh_oracle_finalize_dataset();
#elif DEF_MYSQL
  wh_mysql_finalize_dataset();
#elif DEF_POSTGRES
  wh_postgres_finalize_dataset();
#endif
}

int
wh_get_new_wid(void) {
#ifdef ORACLE
  return wh_oracle_get_new_wid(); 
#elif DEF_MYSQL
  return wh_mysql_get_new_wid();
#elif DEF_POSTGRES
  return wh_postgres_get_new_wid();
#endif
}

int
wh_get_new_special_wid(void) {
#ifdef ORACLE
  return wh_oracle_get_new_special_wid(); 
#elif DEF_MYSQL
  return wh_mysql_get_new_special_wid();
#elif DEF_POSTGRES
  return wh_postgres_get_new_special_wid();
#endif
}

void
wh_sql_error (char *msg) {
#ifdef ORACLE
  wh_oracle_sql_error(msg);
/* #elif DEF_MYSQL */
/*   wh_mysql_sql_error(msg); */
#endif
}

void
wh_on_error(char * error_string) {
  printf("wh_on_error called : %s\n", error_string);
}

int
wh_dataset_exists(int dswid) {
#ifdef ORACLE
  return wh_oracle_dataset_exists (dswid);
#elif DEF_MYSQL
  return wh_mysql_dataset_exists (dswid);
#elif DEF_POSTGRES
  return wh_postgres_dataset_exists (dswid);
#endif  
}


/*** Schema-specific utilities ***/

void
wh_insert_citation(int citation_wid, char * fulltext,
		   char * pubmed_id, short pubmed_id_ind) {
#ifdef ORACLE
  wh_oracle_insert_citation (citation_wid, fulltext,
			     pubmed_id, pubmed_id_ind);
#elif DEF_MYSQL
  wh_mysql_insert_citation (citation_wid, fulltext,
#elif DEF_POSTGRES
  wh_postgres_insert_citation (citation_wid, fulltext,
			    pubmed_id, pubmed_id_ind);
#endif
}

void
wh_update_citation(int citation_wid, char * fulltext,
		   char * pubmed_id, short pubmed_id_ind) {
#ifdef ORACLE
  wh_oracle_update_citation (citation_wid, fulltext,
			     pubmed_id, pubmed_id_ind);
#elif DEF_MYSQL
  wh_mysql_update_citation (citation_wid, fulltext,
#elif DEF_POSTGRES
  wh_postgres_update_citation (citation_wid, fulltext,
			    pubmed_id, pubmed_id_ind);
#endif
}

void
wh_insert_comment(int object_wid, char* comment) {
  if (comment &&
      (0 != strlen(comment)) &&
      (0 != strcmp(comment, "-")) &&
      (0 != strcmp(comment, " ")) ) {
#ifdef ORACLE
    wh_oracle_insert_comment(object_wid, comment);
#elif DEF_MYSQL
    wh_mysql_insert_comment(object_wid, comment);
#elif DEF_POSTGRES
    wh_postgres_insert_comment(object_wid, comment);
#endif
  }
}

void
wh_insert_crossreference(int other_wid, char * xid, char * dataset_name) {
  if (xid &&
      (0 != strlen(xid)) &&
      (0 != strcmp(xid, "NIL"))) {
#ifdef ORACLE
    wh_oracle_insert_crossreference(other_wid, xid, dataset_name);
#elif DEF_MYSQL
    wh_mysql_insert_crossreference(other_wid, xid, dataset_name);
#elif DEF_POSTGRES
    wh_postgres_insert_crossreference(other_wid, xid, dataset_name);
#endif
  }
}

void
wh_insert_description(int object_wid, char* table_name, char* description) {
  if (description &&
      (0 != strlen(description)) &&
      (0 != strcmp(description, "-")) &&
      (0 != strcmp(description, " ")) ) {
#ifdef ORACLE
    wh_oracle_insert_description(object_wid, table_name, description);
#elif DEF_MYSQL
    wh_mysql_insert_description(object_wid, table_name, description);
#elif DEF_POSTGRES
    wh_postgres_insert_description(object_wid, table_name, description);
#endif
  }
}

void
wh_insert_dbid(int object_wid, char* name) {
#ifdef ORACLE
  wh_oracle_insert_dbid(object_wid, name);
#elif DEF_MYSQL
  wh_mysql_insert_dbid(object_wid, name);
#elif DEF_POSTGRES
  wh_postgres_insert_dbid(object_wid, name);
#endif
}

void
wh_insert_entry(int object_wid, int load_error, int lineno, char* change_date) {
#ifdef ORACLE
  wh_oracle_insert_entry(object_wid, load_error, lineno, change_date);
#elif DEF_MYSQL
     /* !Currently ignores change_date */
  wh_mysql_insert_entry(object_wid, load_error, lineno);
#elif DEF_POSTGRES
  wh_postgres_insert_entry(object_wid, load_error, lineno);
#endif
}

void
wh_insert_synonymtable(int object_wid, char* name) {
  if (name && strlen(name) > 0 &&
      (0 != strcmp(name, "NULL")) ) {
#ifdef ORACLE
    wh_oracle_insert_synonymtable(object_wid, name);
#elif DEF_MYSQL
    wh_mysql_insert_synonymtable(object_wid, name);
#elif DEF_POSTGRES
    wh_postgres_insert_synonymtable(object_wid, name);
#endif
  }
}

void
wh_insert_enzymaticreaction(int enzrxn_wid, int reaction_wid, int protein_wid,
			    int complex_wid, short complex_ind,
			    char * reaction_direction, short reaction_direction_ind) {
#ifdef ORACLE
  wh_oracle_insert_enzymaticreaction (enzrxn_wid, reaction_wid, protein_wid,
				      complex_wid, complex_ind,
				      reaction_direction, reaction_direction_ind);
#elif DEF_MYSQL
  wh_mysql_insert_enzymaticreaction (enzrxn_wid, reaction_wid, protein_wid,
#elif DEF_POSTGRES
  wh_postgres_insert_enzymaticreaction (enzrxn_wid, reaction_wid, protein_wid,
				     complex_wid, complex_ind,
				     reaction_direction, reaction_direction_ind);
#endif
}

void
wh_insert_enzrxn_cofactor(int enzrxn_wid, int chemical_wid,
			  char prosthetic, short prosthetic_ind) {
#ifdef ORACLE
  wh_oracle_insert_enzrxn_cofactor (enzrxn_wid, chemical_wid, prosthetic, prosthetic_ind);
#elif DEF_MYSQL
  wh_mysql_insert_enzrxn_cofactor (enzrxn_wid, chemical_wid, prosthetic, prosthetic_ind);
#elif DEF_POSTGRES
  wh_postgres_insert_enzrxn_cofactor (enzrxn_wid, chemical_wid, prosthetic, prosthetic_ind);
#endif
}

void
wh_insert_enzrxn_chemical(int enzrxn_wid, int chemical_wid,
			  char mechanism, char inhibit_or_activate, char physiologically_relevant) {
#ifdef ORACLE
  wh_oracle_insert_enzrxn_chemical (enzrxn_wid, chemical_wid, mechanism, inhibit_or_activate, physiologically_relevant);
#elif DEF_MYSQL
  wh_mysql_insert_enzrxn_chemical (enzrxn_wid, chemical_wid, mechanism, inhibit_or_activate, physiologically_relevant);
#elif DEF_POSTGRES
  wh_postgres_insert_enzrxn_chemical (enzrxn_wid, chemical_wid, mechanism, inhibit_or_activate, physiologically_relevant);
#endif
}

void
wh_insert_enzrxn_alternate(int enzrxn_wid, int primary_wid, int alternate_wid, char cofactor) {
#ifdef ORACLE
  wh_oracle_insert_enzrxn_alternate (enzrxn_wid, primary_wid, alternate_wid, cofactor);
#elif DEF_MYSQL
  wh_mysql_insert_enzrxn_alternate (enzrxn_wid, primary_wid, alternate_wid, cofactor);
#elif DEF_POSTGRES
  wh_postgres_insert_enzrxn_alternate (enzrxn_wid, primary_wid, alternate_wid, cofactor);
#endif	      
}


void
wh_insert_chemical_simple(int chemical_wid, char * name) {
  short name_ind;
  char * tname = string_column(name, 255, &name_ind);
  /*-printf("Inserting %s into Chemical WID=%d\n", name, chemical_wid);  */
  if (name_ind == INDICATE_NULL) return;
  
#ifdef ORACLE
  wh_oracle_insert_chemical_simple (chemical_wid, tname);
#elif DEF_MYSQL
  wh_mysql_insert_chemical_simple (chemical_wid, tname);
#elif DEF_POSTGRES
  wh_postgres_insert_chemical_simple (chemical_wid, tname);
#endif  
}

void
wh_insert_support(int support_wid, int other_wid, char *type, short type_ind) {
#ifdef ORACLE
  wh_oracle_insert_support(support_wid, other_wid, type, type_ind);
#elif DEF_MYSQL
  wh_mysql_insert_support(support_wid, other_wid, type, type_ind);
#elif DEF_POSTGRES
  wh_postgres_insert_support(support_wid, other_wid, type, type_ind);
#endif
}

void
wh_insert_transcription_unit(int transcription_unit_wid, char* name) {
#ifdef ORACLE
  wh_oracle_insert_transcription_unit(transcription_unit_wid, name);
#elif DEF_MYSQL
  wh_mysql_insert_transcription_unit(transcription_unit_wid, name);
#elif DEF_POSTGRES
  wh_postgres_insert_transcription_unit(transcription_unit_wid, name);
#endif
}

void
wh_insert_relatedterm(int term_wid, int other_wid, char *relationship, short relationship_ind) {
#ifdef ORACLE
  wh_oracle_insert_relatedterm (term_wid, other_wid, relationship, relationship_ind);
#elif DEF_MYSQL
  wh_mysql_insert_relatedterm (term_wid, other_wid, relationship, relationship_ind);
#elif DEF_POSTGRES
  wh_postgres_insert_relatedterm (term_wid, other_wid, relationship, relationship_ind);
#endif
}

void
wh_insert_location(int protein_wid, char *location) {
#ifdef ORACLE
  wh_oracle_insert_location(protein_wid, location);
#elif DEF_MYSQL
  wh_mysql_insert_location(protein_wid, location);
#elif DEF_POSTGRES
  wh_postgres_insert_location(protein_wid, location);
#endif
}

void
wh_insert_transcription_unit_component(int tu_wid, int other_wid, char * type) {
#ifdef ORACLE
  wh_oracle_insert_transcription_unit_component (tu_wid, other_wid, type);
#elif DEF_MYSQL
  wh_mysql_insert_transcription_unit_component  (tu_wid, other_wid, type);
#elif DEF_POSTGRES
  wh_postgres_insert_transcription_unit_component  (tu_wid, other_wid, type);
#endif                                                                                
}


/*** SQL SELECTs ***/

int
wh_select_dataset_wid(char * db_name) {
#ifdef ORACLE
  return wh_oracle_select_dataset_wid (db_name);
#elif DEF_MYSQL
  return wh_mysql_select_dataset_wid (db_name);
#elif DEF_POSTGRES
  return wh_postgres_select_dataset_wid (db_name);
#endif  
}

int
wh_select_ecnumber_reaction(char * ecnumber, int enzyme_dataset_wid) {
#ifdef ORACLE
  return wh_oracle_select_ecnumber_reaction (ecnumber, enzyme_dataset_wid);
#elif DEF_MYSQL
  return wh_mysql_select_ecnumber_reaction (ecnumber, enzyme_dataset_wid);
#elif DEF_POSTGRES
  return wh_postgres_select_ecnumber_reaction (ecnumber, enzyme_dataset_wid);
#endif  
}

int
wh_select_geneticcode_wid (char * gencode_id) {
#ifdef ORACLE
  return wh_oracle_select_geneticcode_wid (gencode_id);
#elif DEF_MYSQL
  return wh_mysql_select_geneticcode_wid (gencode_id);
#elif DEF_POSTGRES
  return wh_postgres_select_geneticcode_wid (gencode_id);
#endif  					      
}
						     
int
wh_select_term(char * term_name, int query_dataset_wid) {
#ifdef ORACLE
  return wh_oracle_select_term (term_name, query_dataset_wid);
#elif DEF_MYSQL
  return wh_mysql_select_term (term_name, query_dataset_wid);
#elif DEF_POSTGRES
  return wh_postgres_select_term (term_name, query_dataset_wid);
#endif  					      
}

int
wh_select_dbid_otherwid(char * xid, int query_dataset_wid) {
#ifdef ORACLE
  return wh_oracle_select_dbid_otherwid (xid, query_dataset_wid);
#elif DEF_MYSQL
  return wh_mysql_select_dbid_otherwid (xid, query_dataset_wid);
#elif DEF_POSTGRES
  return wh_postgres_select_dbid_otherwid (xid, query_dataset_wid);
#endif  					      
}


/* SQL DELETEs */

void
wh_delete_object(char * table_name, int wid) {
  char query[MAX_QUERY_SIZE];

  if (snprintf(query, MAX_QUERY_SIZE,
	       "DELETE FROM \"%s\" WHERE \"WID\" = %d", table_name, wid)  >= MAX_QUERY_SIZE )
    fatal_error("Query too long\n");
  
#ifdef ORACLE
  wh_oracle_run_query(query, "Failed to delete object");
#elif DEF_MYSQL
  wh_mysql_run_query(query, "Failed to delete object");
#elif DEF_POSTGRES
  wh_postgres_run_query(query, "Failed to delete object");
#endif
}



/*** Non-database Utilities ***/

void
wh_add_undefined(char * word, char * message_format) {
  /* If word is not yet in the undefined words table, add it.
     Increment the int associated with it to count the #of occurrences.
     Print a message in the given format.
     message_format should contain one %s for word.
  */
  int num_occurrences = wh_find_undefined(word);
  /*-- printf("wh_add_undefined %s\n", word); */
  if (num_occurrences == 0) {
    widtable_insert(undefined_words, 1, word);
    printf(message_format, word);
  }
  else
    widtable_update(undefined_words, num_occurrences+1, word);
}

int
wh_find_undefined(char * word) {
  /* Return 0 if word has not been passed to wh_add_undefined.
     Else return the #times word has been noted as undefined.
  */
  return(widtable_lookup(undefined_words, word));
}

/*---------------------------------------------------------------------
* GetTimeDiff
* This function takes the struct timeval, and returns the difference
* between the current time and this time in milliseconds.
*---------------------------------------------------------------------*/
long
GetTimeDiff(struct timeval *time_sent)
{
  struct timeval * curr_time;
  long time_diff;

  curr_time = (struct timeval *) malloc(sizeof(struct timeval));
  if (!curr_time)
  {
    fprintf(stderr, "\nERROR: Failed to allocate memory!\n");
    exit(1);
  }
  gettimeofday(curr_time, NULL);

  time_diff = ((curr_time->tv_sec - time_sent->tv_sec) * 1000) +
    ((curr_time->tv_usec - time_sent->tv_usec)/1000);

  free(curr_time);
  return time_diff;
}

// Return NULL if input is not a date of form YYYY-MM-DD or YYYY/MM/DD,
//  allowing for optional lead zeroes where valid
const char *
validate_date(const char *input, struct tm *time_buffer)
{
  const char *cp;
  /* First clear the result structure.  */
  memset(time_buffer, '\0', sizeof(*time_buffer));
  /* Try the ISO format first.  */
  cp = strptime(input, "%Y-%m-%d", time_buffer);
  if (cp != NULL)
    return cp;
  memset(time_buffer, '\0', sizeof(*time_buffer));
  cp = strptime(input, "%Y/%m/%d", time_buffer);
  return cp;
}
