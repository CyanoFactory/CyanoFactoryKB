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
 * Utility code common to all database loaders.
 * Define MAIN_PGM in main to control declaration of globals
 * Include this file after definition of MAIN_PGM
 */

#include <assert.h>
#include <errno.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include "widtable.h"
#include "string-util.h"

#ifndef _WH_H
#define _WH_H 1

#ifdef MAIN_PGM
#define EXTERN
#else
#define EXTERN extern
#endif

#define WAREHOUSE_VERSION_STRING "4.6"

/* Values used in _ind variables to indicate presence/absence of value for SQL.
   ProC uses these, but we use them for all DBMSes.
 */
#define INDICATE_NULL -1
#define INDICATE_NONNULL 0 


/* Command-line settings */
EXTERN int commit_freq;  /* commit every this many WIDs, or sometimes input records;
			    nonzero also causes commit after load of each table */
EXTERN char *userid; /* for Oracle: USERID/PASSWORD
			for MYSQL : USERID */
EXTERN char * password;
EXTERN char * host;
EXTERN char * database_name;
EXTERN unsigned int server_port;

/* Globals */
EXTERN int dataset_wid;    


/*** Non-database Utilities ***/

extern struct tm time_buffer;     // used by validate_date
extern WIDTABLE undefined_words;  // Must be initialized via widtable_init(undefined_words)

void
wh_init();

void
wh_add_undefined(char * word, char * message_format);

int
wh_find_undefined(char * word);

long
GetTimeDiff(struct timeval *);

const char *
validate_date(const char *input, struct tm *time_buffer);

int
ifloor(double x, short* conversion_error);

int
iceil(double x, short* conversion_error);
                                           
/* Column value handling */
char *
string_column(char *columnval, int maxlen, short *indicator);

double
number_column(char *columnval, short *indicator, short *error);

char
boolean_column(char *columnval, short *indicator);


/*** General Warehouse utilities ***/

void
wh_commit (void);

void
wh_connect (void);

void
wh_disconnect (void);

void
wh_finalize_dataset (void);

int
wh_get_new_wid (void);

int
wh_get_new_special_wid (void);

void
wh_sql_error (char *);

void
wh_on_error(char * error_string);

int
wh_dataset_exists(int dswid);


/*** SQL INSERTs and UPDATEs ***/

void
wh_insert_bare_object(char * object_table, int wid);

void
wh_insert_linktable(char * column1, int wid1, char* column2, int wid2);

void
wh_insert_named_linktable(char * table_name, char * column1, int wid1, char* column2, int wid2);

void
wh_insert_citation(int citation_wid, char * fulltext,
		   char * pubmed_id, short pubmed_id_ind);

void
wh_update_citation(int citation_wid, char * fulltext,
		   char * pubmed_id, short pubmed_id_ind);

void
wh_insert_comment(int object_wid, char * comment);

void
wh_insert_description(int object_wid, char* table_name, char * description);

void
wh_insert_dbid(int object_wid, char * name);

void
wh_insert_entry(int object_wid, int load_error, int lineno, char* change_date);

void
wh_insert_synonymtable(int object_wid, char * name);

void
wh_insert_crossreference(int other_wid, char * xid, char * dataset_name);

void
wh_insert_enzymaticreaction(int enzrxn_wid, int reaction_wid, int protein_wid,
			    int complex_wid, short complex_ind,
			    char * reaction_direction, short reaction_direction_ind);

void
wh_insert_enzrxn_cofactor(int enzrxn_wid, int chemical_wid,
			  char prosthetic, short prosthetic_ind);

void
wh_insert_enzrxn_chemical(int enzrxn_wid, int chemical_wid, char mechanism, char inhibit_or_activate, char physiologically_relevant);

void
wh_insert_enzrxn_alternate(int enzrxn_wid, int primary_wid, int alternate_wid, char cofactor);

void
wh_insert_support(int support_wid, int other_wid, char *type, short type_ind);

void
wh_insert_relatedterm(int term_wid, int other_wid, char *relationship, short relationship_ind);

void
wh_insert_location(int protein_wid, char *location);

void
wh_insert_transcription_unit_component(int tu_wid, int other_wid, char * type);


/*** SQL SELECTs ***/

int
wh_select_dataset_wid(char * db_name);

int
wh_select_ecnumber_reaction(char * ecnumber, int enzyme_dataset_wid);

int
wh_select_geneticcode_wid(char * gencode_id);
						     
int
wh_select_term(char * term_name, int query_dataset_wid);

int
wh_select_dbid_otherwid(char * xid, int query_dataset_wid);


/*** SQL DELETEs ***/

void
wh_delete_object(char * table_name, int wid);

#endif
