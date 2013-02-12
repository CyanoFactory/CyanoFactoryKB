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
 * MYSQL functions
 */

#ifndef _WH_MYSQL_UTIL_H
#define _WH_MYSQL_UTIL_H 1

/* Buffer sizes for converting args to strings or escaping them */
#define MAX_FIXED_COLUMN_SIZE 100
#define MAX_SQL_NAME_SIZE 250

void
wh_mysql_connect (void);

void
wh_mysql_commit (void);

void
wh_mysql_disconnect (void);

int
wh_mysql_get_new_wid (void);

int
wh_mysql_get_new_special_wid(void);

char *
wh_mysql_malloc_and_create_escaped_string(char * from_string); 

void
wh_mysql_insert_entry(int object_wid, int load_error, int lineno);

void
wh_mysql_insert_citation(int citation_wid, char * fulltext,
			 char * pubmed_id, short pubmed_id_ind);

void
wh_mysql_update_citation(int citation_wid, char * fulltext,
			 char * pubmed_id, short pubmed_id_ind);

void
wh_mysql_insert_dbid(int other_wid, char * xid);

void
wh_mysql_run_query(char * query, char * error_string);

int
wh_mysql_dataset_exists(int dswid);

void
wh_mysql_insert_synonymtable(int wid, char * syn_name);

void
wh_mysql_insert_description(int object_wid, char* table_name, char * description);

void
wh_mysql_insert_comment(int object_wid, char * comment);

void
wh_mysql_insert_chemical_simple(int chemical_wid, char * name);

void
wh_mysql_insert_enzymaticreaction(int enzrxn_wid, int reaction_wid, int protein_wid,
				  int complex_wid, short complex_ind,
				  char * reaction_direction, short reaction_direction_ind);

void
wh_mysql_insert_enzrxn_cofactor(int enzrxn_wid, int chemical_wid,
				char prosthetic, short prosthetic_ind);
							   
void
wh_mysql_insert_enzrxn_chemical(int enzrxn_wid, int chemical_wid,
				char mechanism, char inhibit_or_activate, char physiologically_relevant);

void
wh_mysql_insert_enzrxn_alternate(int enzrxn_wid, int primary_wid, int alternate_wid, char cofactor);

void
wh_mysql_insert_support(int support_wid, int other_wid, char *type, short type_ind);

void
wh_mysql_insert_relatedterm(int term_wid, int other_wid, char *relationship, short relationship_ind);

void
wh_mysql_insert_location(int protein_wid, char *location);


char *
wh_mysql_int_column(int num, short indicator, char * buffer);

char *
wh_mysql_double_column(double num, short indicator, char * buffer);

char *
wh_mysql_char_column(char c, short indicator, char * buffer);

char *
wh_mysql_string_column(char * string, short indicator);

void
wh_mysql_fatal_error(char *message);

int
wh_mysql_select_dataset_wid(char * db_name);

int
wh_mysql_select_ecnumber_reaction(char * ecnumber, int enzyme_dataset_wid);

int
wh_mysql_select_geneticcode_wid (int gencode_id);

int
wh_mysql_select_term(char * term_name, int query_dataset_wid);

int
wh_mysql_select_dbid_otherwid(char * xid, int query_dataset_wid);

#endif

