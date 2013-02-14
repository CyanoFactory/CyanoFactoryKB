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

#ifndef DB_H
#define DB_H 1

#include "compound-parse.h"
#include "protein-parse.h"

void
db_make_biocyc_dataset(char * name, char * version, char * release_date);
void
db_make_biocyc_organism (void);

void
db_insert_all_subpathway_predecessors();

void
db_insert_into_chemical(int compound_wid,
			struct compound_entry *entry,
			char * systematic_name,
			char * cas_registry_numbers,
			char * smiles, char * formula);
void
db_insert_class_into_chemical(int compound_wid, char * class_name);

int
db_select_chemical_wid(char * name);

void
db_insert_into_protein(int wid, struct protein_entry *entry,
		       double molecular_weight, double molecular_weight_exp);
int
db_select_protein(char * name);

void
db_update_protein_full(int wid, struct protein_entry *entry,
		       double molecular_weight, double molecular_weight_exp, double pi_calc);
void
db_update_protein_aasequence(int wid, char * aa_sequence);

int
db_select_enzrxn(char *name);

void
db_insert_into_gene(int gene_wid,
		    char * common_name,
		    char * unique_id,
		    char direction,
		    char * startpos, short startpos_ind,
		    char * endpos, short endpos_ind,
		    char interrupted, short interrupted_ind);
int
db_select_gene(char * name);

void
db_insert_into_feature(int feature_wid, char * class, char * geometry, char * point_type,
                       char * common_name, short common_name_ind,
                       int startpos, short startpos_ind,
                       int endpos, short endpos_ind);

void
db_insert_into_reaction(int wid,
			char * delta_g0, short delta_g0_ind,
			char * ec_number, short ec_number_ind, short ec_number_proposed_ind,
			char spontaneous, short spontaneous_ind); 
void
db_insert_into_reactant(int reaction_wid,int chemical_wid,double coefficient);
void
db_insert_into_product(int reaction_wid,int chemical_wid,double coefficient);
int
db_select_reaction(char* name);

void
db_insert_into_pathway(int wid, char * name);
void
db_update_pathway(int wid, char * name);
void
db_insert_into_pathwaylink(int pathway1_wid, int pathway2_wid, int chemical_wid);
void
db_insert_into_superpathway(int sub_pathway_wid, int super_pathway_wid);
void
db_insert_into_pathwayreaction(int pathway_wid,
			       int reaction_wid,
			       char hypothetical,
			       int prior_reaction_wid,
			       short prior_reaction_wid_ind);
int
db_select_pathway(char * name);


#endif

