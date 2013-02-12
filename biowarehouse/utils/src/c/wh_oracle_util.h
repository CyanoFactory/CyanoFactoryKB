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
 * Oracle functions
 */

#ifndef _WH_ORACLE_UTIL_H
#define _WH_ORACLE_UTIL_H 1

void
wh_oracle_connect (void);

void
wh_oracle_commit (void);

void
wh_oracle_disconnect (void);

int
wh_oracle_get_new_wid (void);

int
wh_oracle_get_new_special_wid(void); 

void
wh_oracle_sql_error (char *);


void
wh_oracle_run_query (char * query, char * failure_msg);

int
wh_oracle_dataset_exists(int dswid);

void
wh_oracle_insert_entry(int object_wid, int load_error, int lineno, char* change_date);

void
wh_oracle_insert_citation(int citation_wid, char * fulltext,
			  char * pubmed_id, short pubmed_id_ind);

void
wh_oracle_update_citation(int citation_wid, char * fulltext,
			  char * pubmed_id, short pubmed_id_ind);


void
wh_oracle_insert_comment(int object_wid, char * comment);

void
wh_oracle_insert_description(int object_wid, char* table_name, char * description);

void
wh_oracle_insert_dbid(int other_wid, char * xid);

void
wh_oracle_insert_synonymtable(int wid, char * syn_name);

void
wh_oracle_insert_crossreference(int other_wid, char * xid, char * dataset_name);

void
wh_insert_linktable(char * column1, int wid1, char* column2, int wid2);

void
wh_oracle_insert_chemical_simple(int chemical_wid, char * name);

void
wh_oracle_insert_enzymaticreaction(int enzrxn_wid, int reaction_wid, int protein_wid,
				   int complex_wid, short complex_ind,
				   char * reaction_direction, short reaction_direction_ind);

void
wh_oracle_insert_enzrxn_cofactor(int enzrxn_wid, int chemical_wid,
				 char prosthetic, short prosthetic_ind);
							   
void
wh_oracle_insert_enzrxn_chemical(int enzrxn_wid, int chemical_wid,
				 char mechanism, char inhibit_or_activate, char physiologically_relevant);

void
wh_oracle_insert_enzrxn_alternate(int enzrxn_wid, int primary_wid, int alternate_wid, char cofactor);

void
wh_oracle_insert_support(int support_wid, int other_wid, char *type, short type_ind);

void
wh_oracle_insert_relatedterm(int term_wid, int other_wid, char *relationship, short relationship_ind);

void
wh_oracle_insert_location(int protein_wid, char *location);


/*** Selects ***/

int
wh_oracle_select_dataset_wid(char * db_name);

int
wh_oracle_select_ecnumber_reaction(char * ecnumber, int enzyme_dataset_wid);

int
wh_oracle_select_geneticcode_wid(char * gencode_id);
						     						     
int
wh_oracle_select_term(char * term_name, int query_dataset_wid);

int
wh_oracle_select_dbid_otherwid(char * xid, int query_dataset_wid);

#endif

