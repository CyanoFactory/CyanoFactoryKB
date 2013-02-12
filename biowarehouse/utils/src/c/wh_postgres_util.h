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
 * Postgres functions
 */

#ifndef _WH_POSTGRES_UTIL_H
#define _WH_POSTGRES_UTIL_H 1

/* Buffer sizes for converting args to strings or escaping them */
#define MAX_FIXED_COLUMN_SIZE 100
#define MAX_SQL_NAME_SIZE 250

void
wh_postgres_connect (void);

void
wh_postgres_commit (void);

void
wh_postgres_disconnect (void);

int
wh_postgres_get_new_wid (void);

int
wh_postgres_get_new_special_wid(void);

char *
wh_postgres_malloc_and_create_escaped_string(char * from_string); 

void
wh_postgres_insert_entry(int object_wid, int load_error, int lineno);

void
wh_postgres_insert_citation(int citation_wid, char * fulltext,
			 char * pubmed_id, short pubmed_id_ind);

void
wh_postgres_update_citation(int citation_wid, char * fulltext,
			 char * pubmed_id, short pubmed_id_ind);

void
wh_postgres_insert_dbid(int other_wid, char * xid);

PGresult*
wh_postgres_run_query(char * query, char * error_string);

int
wh_postgres_dataset_exists(int dswid);

void
wh_postgres_insert_synonymtable(int wid, char * syn_name);

void
wh_postgres_insert_description(int object_wid, char* table_name, char * description);

void
wh_postgres_insert_comment(int object_wid, char * comment);

void
wh_postgres_insert_chemical_simple(int chemical_wid, char * name);

void
wh_postgres_insert_enzymaticreaction(int enzrxn_wid, int reaction_wid, int protein_wid,
				  int complex_wid, short complex_ind,
				  char * reaction_direction, short reaction_direction_ind);

void
wh_postgres_insert_enzrxn_cofactor(int enzrxn_wid, int chemical_wid,
				char prosthetic, short prosthetic_ind);
							   
void
wh_postgres_insert_enzrxn_chemical(int enzrxn_wid, int chemical_wid,
				char mechanism, char inhibit_or_activate, char physiologically_relevant);

void
wh_postgres_insert_enzrxn_alternate(int enzrxn_wid, int primary_wid, int alternate_wid, char cofactor);

void
wh_postgres_insert_support(int support_wid, int other_wid, char *type, short type_ind);

void
wh_postgres_insert_relatedterm(int term_wid, int other_wid, char *relationship, short relationship_ind);

void
wh_postgres_insert_location(int protein_wid, char *location);


char *
wh_postgres_int_column(int num, short indicator, char * buffer);

char *
wh_postgres_double_column(double num, short indicator, char * buffer);

char *
wh_postgres_char_column(char c, short indicator, char * buffer);

char *
wh_postgres_string_column(char * string, short indicator);

void
wh_postgres_fatal_error(char *message);

int
wh_postgres_select_dataset_wid(char * db_name);

int
wh_postgres_select_ecnumber_reaction(char * ecnumber, int enzyme_dataset_wid);

int
wh_postgres_select_geneticcode_wid (int gencode_id);

int
wh_postgres_select_term(char * term_name, int query_dataset_wid);

int
wh_postgres_select_dbid_otherwid(char * xid, int query_dataset_wid);

int
wh_postgres_num_rows(const PGresult *res);

#endif

