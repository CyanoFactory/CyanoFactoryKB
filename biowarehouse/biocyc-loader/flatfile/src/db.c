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
#include "db.h"
#include "stdio.h"
#include "main.h"

#define MAX_QUERY_SIZE 100000

void
db_make_biocyc_dataset(char * name, char * version, char * release_date) {
  if (!validate_date(release_date, &time_buffer)) {
    fprintf(stderr, "Illegal release date: %s\n", release_date);
    exit(-1);
  }
  printf("Creating dataset %s, version %s\n", name, version);
  
#ifdef ORACLE
  wh_oracle_make_biocyc_dataset(name, version, release_date);
#elif DEF_MYSQL
  wh_mysql_make_biocyc_dataset(name, version, release_date);
#endif
}
 
void
db_make_biocyc_organism(void) {
#ifdef ORACLE
  wh_oracle_make_biocyc_organism();
#elif DEF_MYSQL
  wh_mysql_make_biocyc_organism();
#endif

  /* Add a row to Entry table for BioSource */
  wh_insert_entry(organism_wid, 0, 0, "");
}

void
db_insert_all_subpathway_predecessors() {
#ifdef ORACLE
  wh_oracle_insert_all_subpathway_predecessors();
#elif DEF_MYSQL
  wh_mysql_insert_all_subpathway_predecessors();
#endif
}

void
db_insert_into_chemical(int compound_wid,
			struct compound_entry *entry,
			char * systematic_name,
			char * cas_registry_numbers,
			char * smiles,
			char * formula) {
#ifdef ORACLE
  wh_oracle_insert_into_chemical(compound_wid,entry,systematic_name ,cas_registry_numbers, smiles,formula);
#elif DEF_MYSQL
  wh_mysql_insert_into_chemical(compound_wid,entry,systematic_name ,cas_registry_numbers, smiles,formula);
#endif
}

void
db_insert_class_into_chemical(int compound_wid, char * class_name) {
#ifdef ORACLE
  wh_oracle_insert_class_into_chemical(compound_wid, class_name);
#elif DEF_MYSQL
  wh_mysql_insert_class_into_chemical(compound_wid, class_name);
#endif
}


int
db_select_chemical_wid(char * name) {
#ifdef ORACLE
  return wh_oracle_select_chemical_wid(name);
#elif DEF_MYSQL
  return wh_mysql_select_chemical_wid(name);
#endif
}

void
db_insert_into_pathway(int wid, char * name){
#ifdef ORACLE
  wh_oracle_insert_into_pathway(wid,name);
#elif DEF_MYSQL
  wh_mysql_insert_into_pathway(wid,name);
#endif
}


void
db_update_pathway(int wid, char * name){
#ifdef ORACLE
  wh_oracle_update_pathway(wid,name);
#elif DEF_MYSQL
  wh_mysql_update_pathway(wid,name);
#endif
}

void
db_insert_into_pathwaylink(int pathway1_wid, int pathway2_wid, int chemical_wid){
#ifdef ORACLE
  wh_oracle_insert_into_pathwaylink(pathway1_wid,pathway2_wid,chemical_wid);
#elif DEF_MYSQL
  wh_mysql_insert_into_pathwaylink(pathway1_wid,pathway2_wid,chemical_wid);
#endif
}

void
db_insert_into_superpathway(int sub_pathway_wid, int super_pathway_wid){
#ifdef ORACLE
  wh_oracle_insert_into_superpathway(sub_pathway_wid, super_pathway_wid);
#elif DEF_MYSQL
  wh_mysql_insert_into_superpathway(sub_pathway_wid, super_pathway_wid);
#endif
}

void
db_insert_into_pathwayreaction(int pathway_wid, int reaction_wid, char hypothetical,
			       int prior_reaction_wid, short prior_reaction_wid_ind) {
#ifdef ORACLE
  wh_oracle_insert_into_pathwayreaction(pathway_wid, reaction_wid, hypothetical,
				      prior_reaction_wid, prior_reaction_wid_ind);
#elif DEF_MYSQL
  wh_mysql_insert_into_pathwayreaction(pathway_wid, reaction_wid, hypothetical,
				     prior_reaction_wid, prior_reaction_wid_ind);
#endif
}

int
db_select_pathway(char * name) {
#ifdef ORACLE
  return wh_oracle_select_pathway(name);
#elif DEF_MYSQL
  return wh_mysql_select_pathway(name);
#endif
}


void
db_insert_into_reaction(int wid,
			char * delta_g0, short delta_g0_ind,
			char * ec_number, short ec_number_ind, short ec_number_proposed_ind,
			char spontaneous, short spontaneous_ind)
{
#ifdef ORACLE
  wh_oracle_insert_into_reaction(wid,delta_g0,delta_g0_ind,ec_number,ec_number_ind,ec_number_proposed_ind,spontaneous,spontaneous_ind);
#elif DEF_MYSQL
  wh_mysql_insert_into_reaction(wid,delta_g0,delta_g0_ind,ec_number,ec_number_ind,ec_number_proposed_ind,spontaneous,spontaneous_ind);
#endif
}

void
db_insert_into_reactant(int reaction_wid, int chemical_wid, double coefficient)
{
#ifdef ORACLE
  wh_oracle_insert_into_reactant(reaction_wid,chemical_wid,coefficient);
#elif DEF_MYSQL
  wh_mysql_insert_into_reactant(reaction_wid,chemical_wid,coefficient);
#endif
}

void
db_insert_into_product(int reaction_wid, int chemical_wid, double coefficient)
{
#ifdef ORACLE
  wh_oracle_insert_into_product(reaction_wid,chemical_wid,coefficient);
#elif DEF_MYSQL
  wh_mysql_insert_into_product(reaction_wid,chemical_wid,coefficient);
#endif
}

int
db_select_reaction(char* name)
{
#ifdef ORACLE
  return wh_oracle_select_reaction(name);
#elif DEF_MYSQL
  return wh_mysql_select_reaction(name);
#endif
}

void
db_insert_into_protein(int wid, struct protein_entry *entry,
		       double molecular_weight, double molecular_weight_exp) {
#ifdef ORACLE
  wh_oracle_insert_into_protein (wid, entry, molecular_weight, molecular_weight_exp);
#elif DEF_MYSQL
  wh_mysql_insert_into_protein (wid, entry, molecular_weight, molecular_weight_exp);
#endif
}

int
db_select_protein(char* name)
{
#ifdef ORACLE
  return wh_oracle_select_protein(name);
#elif DEF_MYSQL
  return wh_mysql_select_protein(name);
#endif
}

void
db_update_protein_full(int wid, struct protein_entry *entry,
		       double molecular_weight, double molecular_weight_exp, double pi_calc) {
#ifdef ORACLE
  wh_oracle_update_protein_full (wid, entry,
			       molecular_weight, molecular_weight_exp, pi_calc);
#elif DEF_MYSQL
  wh_mysql_update_protein_full (wid, entry,
			      molecular_weight, molecular_weight_exp, pi_calc);
#endif
}

void
db_update_protein_aasequence(int wid, char * aa_sequence)
{
#ifdef ORACLE
  wh_oracle_update_protein_aasequence(wid, aa_sequence);
#elif DEF_MYSQL
  wh_mysql_update_protein_aasequence(wid, aa_sequence);
#endif
}

int
db_select_enzrxn(char* name) {
#ifdef ORACLE
  return wh_oracle_select_enzrxn(name);
#elif DEF_MYSQL
  return wh_mysql_select_enzrxn(name);
#endif
}

void
db_insert_into_feature(int feature_wid, char * class, char * geometry, char * point_type,
                       char * common_name, short common_name_ind,
                       int startpos, short startpos_ind,
                       int endpos, short endpos_ind) {
#ifdef ORACLE
  wh_oracle_insert_into_feature (feature_wid, class, geometry, point_type, common_name, common_name_ind, startpos, startpos_ind, endpos, endpos_ind);
#elif DEF_MYSQL
  wh_mysql_insert_into_feature  (feature_wid, class, geometry, point_type, common_name, common_name_ind, startpos, startpos_ind, endpos, endpos_ind);
#endif 
}


void
db_insert_into_gene(int gene_wid, char * common_name, char * unique_id, char direction,
		    char * startpos, short startpos_ind,
		    char * endpos, short endpos_ind,
		    char interrupted, short interrupted_ind) {
#ifdef ORACLE
  wh_oracle_insert_into_gene(gene_wid,common_name,unique_id,direction,startpos,startpos_ind,endpos,endpos_ind,interrupted,interrupted_ind);
#elif DEF_MYSQL
  wh_mysql_insert_into_gene(gene_wid,common_name,unique_id,direction,startpos,startpos_ind,endpos,endpos_ind,interrupted,interrupted_ind);
#endif
}

int
db_select_gene(char * name) {
#ifdef ORACLE
  return wh_oracle_select_gene(name);
#elif DEF_MYSQL
  return wh_mysql_select_gene(name);
#endif
}

