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
#include <libpq-fe.h>
#include "stdio.h"
#include "main.h"
#include "wh_postgres.h"
#include "wh_postgres_util.h"

#define MAX_QUERY_SIZE 2000000

extern int  dataset_wid;
extern PGconn* pg_handle;
char query[MAX_QUERY_SIZE];
int  organism_wid;
char *organism_name;  /* set by cmdline option */

void
wh_postgres_make_biocyc_dataset(char * name, char * version, char * release_date) {
  char * user_running_loader = getenv("USER");
  char * warehouse_version_string = WAREHOUSE_VERSION_STRING;  
  
  dataset_wid = wh_postgres_get_new_special_wid();

  if(snprintf(query, MAX_QUERY_SIZE,
	      "INSERT INTO \"DataSet\" (\"WID\", \"Name\", \"Version\", \"ReleaseDate\", \"LoadDate\", \"LoadedBy\", \"Application\", \"ApplicationVersion\", \"HomeURL\", \"QueryURL\") Values(%d, '%s', '%s', '%s', NOW(), '%s', 'BioCyc Loader', '%s', 'http://www.biocyc.org', 'http://www.biocyc.org:1555')" ,
	      dataset_wid, name, version, release_date, user_running_loader, warehouse_version_string) >= MAX_QUERY_SIZE)
	wh_postgres_fatal_error("In wh_postgres_make_biocyc_dataset: Query string too long\n");
  wh_postgres_run_query(query, "Failed to insert into DataSet");
  
  printf("Made dataset, wid is %d\n", dataset_wid);
}

void
wh_postgres_make_biocyc_organism(void) {
   organism_wid = wh_postgres_get_new_wid();
   if(snprintf(query, MAX_QUERY_SIZE,
	      "INSERT INTO \"BioSource\" (\"WID\", \"Name\", \"DataSetWID\") VALUES(%d,'%s',%d)",
	       organism_wid, organism_name, dataset_wid) >= MAX_QUERY_SIZE)
	wh_postgres_fatal_error("In wh_postgres_make_biocyc_organism: Query string too long\n");
   wh_postgres_run_query(query, "Failed to insert into BioSource");
}

void
wh_postgres_insert_all_subpathway_predecessors() {
  // For each subpathway-predecessor - superpathway pair in global list:
  //   duplicate all subpathway PathwayReaction rows, except update the
  //   PathwayReaction.pathwaywid to reference the superpathway rather than
  //   the subpathway.
  // In this way the super inherits the reaction graph of the sub.

  WIDNAMELIST psp;
  int subwid, superwid;
    
  // Create temp table to copy PathwayReaction rows into
  wh_postgres_run_query("CREATE TABLE \"TempPathwayReaction\" ( LIKE \"PathwayReaction\" )",
		     "Failed to create table TempPathwayReaction");

  // For each subpathway-predecessor - superpathway pair in global list:
  for ( psp = predecessor_subpathways; psp; psp = psp->next) {
    superwid = psp->wid;
    subwid = find_pathway(psp->name);
    if (subwid == 0) {
      printf("Undefined predecessor pathway %s\n", psp->name);
      break;
    }
    
    // Select all rows for subpathway into temp table
    if (snprintf(query, MAX_QUERY_SIZE,
		 "INSERT INTO \"TempPathwayReaction\" SELECT * FROM \"PathwayReaction\" WHERE \"PathwayWID\" = %d", subwid)
	 >= MAX_QUERY_SIZE)
      wh_postgres_fatal_error("In wh_postgres_insert_all_subpathway_predecessors : Query string too long\n");
    wh_postgres_run_query(query, "Failed to insert into TempPathwayReaction");

    // Change all PathwayWID values from subwid to superwid
    if (snprintf(query, MAX_QUERY_SIZE,
		 "UPDATE \"TempPathwayReaction\" SET \"PathwayWID\" = %d", superwid)
	 >= MAX_QUERY_SIZE)
      wh_postgres_fatal_error("In wh_postgres_insert_all_subpathway_predecessors : Query string too long\n");
    wh_postgres_run_query(query, "Failed to update TempPathwayReaction");

    // Copy all rows back into PathwayReaction
    wh_postgres_run_query("INSERT INTO \"PathwayReaction\" SELECT * FROM \"TempPathwayReaction\"",
		       "Failed to copy from TempPathwayReaction into PathwayReaction");

    // Clear the temp table
    wh_postgres_run_query("DELETE FROM \"TempPathwayReaction\"",
		       "Failed to clear TempPathwayReaction");
  }

  // Remove temp table
  wh_postgres_run_query("DROP TABLE \"TempPathwayReaction\"",
		     "Failed to drop TempPathwayReaction");
}


void
wh_postgres_insert_into_chemical(int compound_wid, struct compound_entry *entry, char * systematic_name, char * cas_registry_numbers, char * smiles, char * formula) {
   
   char  * escaped_common_name = wh_postgres_malloc_and_create_escaped_string(entry->common_name);
   char  * escaped_systematic_name = wh_postgres_malloc_and_create_escaped_string(systematic_name);
   char  * escaped_registry_numbers = wh_postgres_malloc_and_create_escaped_string(cas_registry_numbers);
   char  * escaped_smiles = wh_postgres_malloc_and_create_escaped_string(smiles);
   char  * escaped_formula = wh_postgres_malloc_and_create_escaped_string(formula);

   if (snprintf(query, MAX_QUERY_SIZE,
		"INSERT INTO \"Chemical\" (\"WID\",\"Name\",\"SystematicName\",\"CAS\",\"MolecularWeightCalc\",\"Charge\",\"Smiles\",\"EmpiricalFormula\",\"PKA1\",\"PKA2\",\"PKA3\",\"DataSetWID\")VALUES (%d, %s, %s, %s, %s, %s, %s, %s, %s , %s, %s, %d)",
	       compound_wid,
	       escaped_common_name,
	       wh_postgres_string_column(escaped_systematic_name, entry->systematic_name_ind),
	       wh_postgres_string_column(escaped_registry_numbers, entry->cas_registry_numbers_ind),
	       wh_postgres_string_column(entry->molecular_weight, entry->molecular_weight_ind),
	       wh_postgres_string_column(entry->charge, entry->charge_ind),
	       wh_postgres_string_column(escaped_smiles, entry->smiles_ind),
	       wh_postgres_string_column(escaped_formula, entry->formula_ind),
	       wh_postgres_string_column(entry->pka1, entry->pka1_ind),
	       wh_postgres_string_column(entry->pka2, entry->pka2_ind),
	       wh_postgres_string_column(entry->pka3, entry->pka3_ind),
	       dataset_wid
	       ) >= MAX_QUERY_SIZE)
	wh_postgres_fatal_error("In wh_postgres_insert_into_chemical:Query string too long\n");

   wh_postgres_run_query(query, "Failed to insert into Chemical");
   free(escaped_common_name);
   free(escaped_systematic_name);
   free(escaped_registry_numbers);
   free(escaped_smiles);
   free(escaped_formula);
}

void
wh_postgres_insert_class_into_chemical(int compound_wid, char * class_name) {
   char  * escaped_common_name = wh_postgres_malloc_and_create_escaped_string(class_name);
   if (snprintf(query, MAX_QUERY_SIZE,
		"INSERT INTO \"Chemical\" (\"WID\",\"Name\",\"Class\",\"DataSetWID\") VALUES (%d, %s, 'T', %d)",
	       compound_wid,
	       escaped_common_name,
	       dataset_wid
	       ) >= MAX_QUERY_SIZE)
	wh_postgres_fatal_error("In wh_postgres_insert_into_chemical:Query string too long\n");

   wh_postgres_run_query(query, "Failed to insert into Chemical");
   free(escaped_common_name);
}

int
wh_postgres_select_chemical_wid(char * name) {
   PGresult *result;
   //MYSQL_ROW row;
   int wid = 0;
   char * escaped_name = wh_postgres_malloc_and_create_escaped_string(name);
  
   if(snprintf(query,MAX_QUERY_SIZE,
	       "SELECT \"WID\" FROM \"Chemical\" WHERE \"Name\" = %s AND \"DataSetWID\" = %d",
	       escaped_name, dataset_wid) >= MAX_QUERY_SIZE)
	wh_postgres_fatal_error("In wh_postgres_select_chemical_wid:Query string too long\n");
   result = wh_postgres_run_query(query, "Failed to get Wid from Chemical Table\n");

     
   
  if (wh_postgres_num_rows(result) > 0) {
	free(escaped_name);
	wid = atoi(PQgetvalue(result, 0, 0));
	PQclear(result);
	return wid;
  }

   PQclear(result);
   
   /* no rows returned - try synonyms */
   if(snprintf(query,MAX_QUERY_SIZE,
	       "SELECT \"WID\", \"OtherWID\" FROM \"DBID\", \"Chemical\" WHERE \"XID\" = %s AND \"OtherWID\" = \"WID\" AND \"DataSetWID\" = %d",
	       escaped_name, dataset_wid) >= MAX_QUERY_SIZE)
	wh_postgres_fatal_error("In wh_postgres_select_chemical_wid:Query string too long\n");

	
   result = wh_postgres_run_query(query, "Failed to get Wid from DBID and Chemical Table\n");
     wh_postgres_num_rows(result);

  if (wh_postgres_num_rows(result) > 0) {
	free(escaped_name);
	wid = atoi(PQgetvalue(result, 0, 0));
	PQclear(result);
	return wid;
  }
  PQclear(result);
	


   /* no rows returned - try synonyms . We need to just look up in the first row.*/
   if(snprintf(query,MAX_QUERY_SIZE,
	       "SELECT \"WID\" FROM \"Chemical\", \"SynonymTable\" WHERE \"Syn\" = %s AND \"OtherWID\" = \"WID\" AND \"DataSetWID\" = %d",
	       escaped_name, dataset_wid) >= MAX_QUERY_SIZE)
	wh_postgres_fatal_error("In wh_postgres_select_chemical_wid:Query string too long\n");
	
	
   result = wh_postgres_run_query(query, "Failed to get Wid from SynonymTable and Chemical Table\n");
     wh_postgres_num_rows(result);

  if (wh_postgres_num_rows(result) > 0) {
	wid = atoi(PQgetvalue(result, 0, 0));
  }
    
   PQclear(result);
   free(escaped_name);
   return wid;
}

void
wh_postgres_insert_into_pathway(int wid, char * name)
{
   char  * escaped_name = wh_postgres_malloc_and_create_escaped_string(name);
   if(snprintf(query,MAX_QUERY_SIZE,
	       "INSERT INTO \"Pathway\" (\"WID\",\"Name\",\"Type\",\"BioSourceWID\",\"DataSetWID\") VALUES (%d,%s,'O',%d,%d)",
	       wid, escaped_name, organism_wid, dataset_wid) >= MAX_QUERY_SIZE)
	wh_postgres_fatal_error("In wh_postgres_insert_into_pathway:Query string too long\n");
   wh_postgres_run_query(query, "Failed to insert into Pathway Table");
   free(escaped_name);
}

void
wh_postgres_update_pathway(int wid, char * name)
{
   char  * escaped_name = wh_postgres_malloc_and_create_escaped_string(name);
   if(snprintf(query,MAX_QUERY_SIZE,
	       "UPDATE \"Pathway\" SET \"Name\" = %s WHERE \"WID\" = %d",
	       escaped_name, wid) >= MAX_QUERY_SIZE)
     wh_postgres_fatal_error("In wh_postgres_update_pathway: Query string too long\n");
   wh_postgres_run_query(query, "Failed to Update Pathway Table");
   free(escaped_name);
}

void
wh_postgres_insert_into_pathwaylink(int pathway1_wid, int pathway2_wid, int chemical_wid){
  if(snprintf(query,MAX_QUERY_SIZE,
	      "INSERT INTO \"PathwayLink\" (\"Pathway1WID\",\"Pathway2WID\",\"ChemicalWID\")VALUES (%d, %d, %d)",
	      pathway1_wid, pathway2_wid, chemical_wid) >= MAX_QUERY_SIZE)
	wh_postgres_fatal_error("In wh_postgres_insert_into_pathwaylink:Query string too long\n");
  wh_postgres_run_query(query, "Failed to insert into PathwayLink Table");
}

void
wh_postgres_insert_into_superpathway(int sub_pathway_wid, int super_pathway_wid){
  if(snprintf(query,MAX_QUERY_SIZE,
	      "INSERT INTO \"SuperPathway\" (\"SubPathwayWID\",\"SuperPathwayWID\") VALUES (%d, %d)",
	      sub_pathway_wid, super_pathway_wid) >= MAX_QUERY_SIZE)
	wh_postgres_fatal_error("In wh_postgres_insert_into_superpathway:Query string too long\n");
  wh_postgres_run_query(query, "Failed to insert into SuperPathway Table");
}

void
wh_postgres_insert_into_pathwayreaction(int pathway_wid, int reaction_wid, char hypothetical,
				     int prior_reaction_wid, short prior_reaction_wid_ind){
  char int_to_char[100];
  if (snprintf(query,MAX_QUERY_SIZE,
	       "INSERT INTO \"PathwayReaction\" (\"PathwayWID\",\"ReactionWID\",\"PriorReactionWID\",\"Hypothetical\") VALUES (%d, %d, %s, '%c')",
	       pathway_wid,
	       reaction_wid,
	       wh_postgres_int_column(prior_reaction_wid, prior_reaction_wid_ind, int_to_char),
	       hypothetical) >= MAX_QUERY_SIZE)
	wh_postgres_fatal_error("In wh_postgres_insert_into_pathwayreaction:Query string too long\n");
  wh_postgres_run_query(query, "Failed to insert into PathWayReaction table");
}

int
wh_postgres_select_pathway(char * name) {
  PGresult *result;
  //MYSQL_ROW row;
  int wid = 0;
  char  * escaped_name = wh_postgres_malloc_and_create_escaped_string(name);
  
  if(snprintf(query,MAX_QUERY_SIZE,
	      "SELECT \"WID\" FROM \"Pathway\", \"DBID\" WHERE \"XID\" = %s AND \"OtherWID\" = \"WID\" AND \"DataSetWID\" = %d",
	      escaped_name, dataset_wid) >= MAX_QUERY_SIZE)
	wh_postgres_fatal_error("In wh_postgres_select_pathway: Query string too long\n");
  
  
  result = wh_postgres_run_query(query, "Failed to get Wid from Pathway Table");

  if (wh_postgres_num_rows(result) > 0) {
	wid = atoi(PQgetvalue(result, 0, 0));
  }
  PQclear(result);
  
  free(escaped_name);
  return wid;
}
   
void
wh_postgres_insert_into_reaction(int wid, char * delta_g0, short delta_g0_ind, char * ec_number, short ec_number_ind, short ec_number_proposed_ind, char spontaneous, short spontaneous_ind){
  char buffer[50]; 
  char  * escaped_delta_g0 = wh_postgres_malloc_and_create_escaped_string(delta_g0);
  char  * escaped_ec_number = wh_postgres_malloc_and_create_escaped_string(ec_number);

  if(snprintf(query,MAX_QUERY_SIZE,
	      "INSERT INTO \"Reaction\" (\"WID\",\"DeltaG\",\"ECNumber\",\"ECNumberProposed\",\"Spontaneous\",\"DataSetWID\")VALUES (%d, %s, %s, %s, %s, %d)",
	wid,
	wh_postgres_string_column(delta_g0, delta_g0_ind),
	wh_postgres_string_column(escaped_ec_number, ec_number_ind),
	wh_postgres_string_column(escaped_ec_number, ec_number_proposed_ind),
	wh_postgres_char_column(spontaneous,spontaneous_ind, buffer),
	dataset_wid) >= MAX_QUERY_SIZE)
	wh_postgres_fatal_error("In wh_postgres_insert_into_reaction: Query string too long\n");
  wh_postgres_run_query(query, "Failed to insert into Reaction Table");
  free(escaped_delta_g0);
  free(escaped_ec_number);
}

void
wh_postgres_insert_into_reactant(int reaction_wid, int chemical_wid, double coefficient){
  if(snprintf(query,MAX_QUERY_SIZE,
	      "INSERT INTO \"Reactant\" (\"ReactionWID\",\"OtherWID\",\"Coefficient\")VALUES (%d, %d, %d)",
	      reaction_wid, chemical_wid, (int) coefficient) >= MAX_QUERY_SIZE)
	wh_postgres_fatal_error("In wh_postgres_insert_into_reactant: Query string too long\n");
  wh_postgres_run_query(query, "Failed to insert into Reactants");
}

void
wh_postgres_insert_into_product(int reaction_wid,int chemical_wid,double coefficient){
  if(snprintf(query,MAX_QUERY_SIZE,
	      "INSERT INTO \"Product\" (\"ReactionWID\",\"OtherWID\",\"Coefficient\")VALUES (%d, %d, %d)",
	      reaction_wid, chemical_wid, (int) coefficient) >= MAX_QUERY_SIZE)
	wh_postgres_fatal_error("In wh_postgres_insert_into_product: Query string too long\n");
  wh_postgres_run_query(query, "Failed to insert into Product");
}

int
wh_postgres_select_reaction(char * name) {
  PGresult *result;
  //MYSQL_ROW row;
  int wid = 0;
  char  * escaped_name = wh_postgres_malloc_and_create_escaped_string(name);
  
  if(snprintf(query,MAX_QUERY_SIZE,
	      "SELECT \"WID\" FROM \"Reaction\", \"DBID\" WHERE \"XID\" = %s AND \"OtherWID\" = \"WID\" AND \"DataSetWID\" = %d",
	      escaped_name, dataset_wid) >= MAX_QUERY_SIZE)
	wh_postgres_fatal_error("In wh_postgres_select_reaction: Query string too long\n");
	
	
	
  result = wh_postgres_run_query(query, "Failed to get Wid from Reaction Table");

  if (wh_postgres_num_rows(result) > 0) {
	wid = atoi(PQgetvalue(result, 0, 0));
  }
  PQclear(result);
	
	
  free(escaped_name);
  return wid;
}

void
wh_postgres_insert_into_protein(int wid, struct protein_entry *entry,
			     double molecular_weight, double molecular_weight_exp)
{
  char  * escaped_common_name = wh_postgres_malloc_and_create_escaped_string(entry->common_name);
  char molecular_weight_buf[MAX_FIXED_COLUMN_SIZE];
  char molecular_weight_exp_buf[MAX_FIXED_COLUMN_SIZE];
  
  if(snprintf(query,MAX_QUERY_SIZE,
	      "INSERT INTO \"Protein\" (\"WID\",\"Name\",\"MolecularWeightCalc\",\"MolecularWeightExp\",\"PICalc\",\"DataSetWID\")VALUES (%d, %s, %s, %s, %s,%d)",
	wid,
	escaped_common_name,
	wh_postgres_double_column(molecular_weight,     entry->molecular_weight_ind, molecular_weight_buf),
	wh_postgres_double_column(molecular_weight_exp, entry->molecular_weight_exp_ind, molecular_weight_exp_buf),
	wh_postgres_string_column(entry->pi, entry->pi_ind),
	dataset_wid) >= MAX_QUERY_SIZE)
	wh_postgres_fatal_error("In wh_postgres_insert_into_protein:Query string too long\n");
  wh_postgres_run_query(query, "Failed to insert into Protein");
  free(escaped_common_name);
}
 
int
wh_postgres_select_protein(char * name) {
  PGresult *result;
  //MYSQL_ROW row;
  int wid = 0;
  char  * escaped_name = wh_postgres_malloc_and_create_escaped_string(name);

 if(snprintf(query,MAX_QUERY_SIZE,
	     "SELECT \"WID\" FROM \"Protein\" , \"DBID\" WHERE \"XID\" = %s AND \"OtherWID\" = \"WID\" AND \"DataSetWID\" = %d",
	     escaped_name, dataset_wid) >= MAX_QUERY_SIZE)
	wh_postgres_fatal_error("In wh_postgres_select_protein: Query string too long\n");
	
	
	
  result = wh_postgres_run_query(query, "Failed to get Wid from Protein Table");
 
  if (wh_postgres_num_rows(result) > 0) {
	wid = atoi(PQgetvalue(result, 0, 0));
  }
  PQclear(result);
	
	
  free(escaped_name);
  return wid;
}

void
wh_postgres_update_protein_full(int wid, struct protein_entry *entry,
			     double molecular_weight, double molecular_weight_exp, double pi_calc) {
  char * escaped_name = wh_postgres_malloc_and_create_escaped_string(entry->common_name);
  char mw_calc_buf[MAX_FIXED_COLUMN_SIZE];
  char mw_exp_buf[MAX_FIXED_COLUMN_SIZE];
  char pi_calc_buf[MAX_FIXED_COLUMN_SIZE];

  if(snprintf(query,MAX_QUERY_SIZE,
	      "UPDATE \"Protein\" SET \"Name\" = %s, \"MolecularWeightCalc\" = %s, \"MolecularWeightExp\" = %s , \"PICalc\" = %s WHERE \"WID\" = %d",
	      escaped_name,
	      wh_postgres_double_column(molecular_weight, entry->molecular_weight_ind, mw_calc_buf),
	      wh_postgres_double_column(molecular_weight_exp, entry->molecular_weight_exp_ind, mw_exp_buf),
	      wh_postgres_double_column(pi_calc, entry->pi_ind, pi_calc_buf),
	      wid) >= MAX_QUERY_SIZE)
	wh_postgres_fatal_error("In wh_postgres_update_protein_aasequence: Query string too long\n");
  wh_postgres_run_query(query,"Failed to update Protein full");
  free(escaped_name);
}

void
wh_postgres_update_protein_aasequence(int wid, char * aa_sequence)
{
  char * escaped_aa_sequence = wh_postgres_malloc_and_create_escaped_string(aa_sequence);
  int len = strlen(aa_sequence);
  if(snprintf(query,MAX_QUERY_SIZE,
	      "UPDATE \"Protein\" SET \"AASequence\" = %s, \"Length\"=%d WHERE \"WID\" = %d",
	      escaped_aa_sequence, len, wid) >= MAX_QUERY_SIZE)
	wh_postgres_fatal_error("In wh_postgres_update_protein_aasequence: Query string too long\n");
  wh_postgres_run_query(query,"Failed to update into Protein.AASequence");
  free(escaped_aa_sequence);
}

void
wh_postgres_insert_into_enzymatic_reaction(int enzrxn_wid,int reaction_wid, int protein_wid,int complex_wid, short complex_ind, char * reaction_direction, short reaction_direction_ind){
  char int_buffer[100];
  char * escaped_reaction_direction = wh_postgres_malloc_and_create_escaped_string(reaction_direction);
  if(snprintf(query,MAX_QUERY_SIZE,
	      "INSERT INTO \"EnzymaticReaction\" (\"WID\",\"ReactionWID\",\"ProteinWID\",\"ComplexWID\",\"ReactionDirection\",\"DataSetWID\")VALUES (%d, %d, %d, %s, %s, %d)", 
	enzrxn_wid,
	reaction_wid,
	protein_wid,
	wh_postgres_int_column(complex_wid,complex_ind,int_buffer),
	wh_postgres_string_column(escaped_reaction_direction, reaction_direction_ind),
	dataset_wid) >= MAX_QUERY_SIZE)
	wh_postgres_fatal_error("In wh_postgres_insert_into_enzymatic_reaction:Query string too long\n");
  wh_postgres_run_query(query,"Failed to insert into EnzymaticReaction");
  free(escaped_reaction_direction);
}
   
int
wh_postgres_select_enzrxn(char* name) {
  PGresult *result;
  //MYSQL_ROW row;
  int wid = 0;
  char  * escaped_name = wh_postgres_malloc_and_create_escaped_string(name);

  if(snprintf(query,MAX_QUERY_SIZE,
	      "SELECT \"WID\" FROM \"EnzymaticReaction\", \"DBID\" WHERE \"XID\" = %s AND \"OtherWID\" = \"WID\" AND \"DataSetWID\" = %d",
	      escaped_name, dataset_wid) >= MAX_QUERY_SIZE)
    wh_postgres_fatal_error("In wh_postgres_select_enzrxn:Query string too long\n");
  
  
  result = wh_postgres_run_query(query, "Failed to get Wid from EnzymaticReaction Table");
 
  if (wh_postgres_num_rows(result) > 0) {
	wid = atoi(PQgetvalue(result, 0, 0));
  }
  PQclear(result);
  
  
  free(escaped_name);
  return wid;
}


void
wh_postgres_insert_into_feature(int feature_wid, char * class, char * geometry, char * point_type,
                             char * common_name, short common_name_ind,
                             int startpos, short startpos_ind,
                             int endpos, short endpos_ind) {
  char  * escaped_common_name = wh_postgres_malloc_and_create_escaped_string(common_name);
  char  * escaped_class = wh_postgres_malloc_and_create_escaped_string(class);
  char  * escaped_point_type = wh_postgres_malloc_and_create_escaped_string(point_type);
  short point_type_ind = (0 == strcasecmp("point", geometry)) ? INDICATE_NONNULL : INDICATE_NULL;
  char start_buffer[100];
  char end_buffer[100];
  
  if(snprintf(query, MAX_QUERY_SIZE,
	      "INSERT INTO \"Feature\" (\"WID\", \"SequenceWID\", \"Type\", \"SequenceType\", \"Description\", \"Class\", \"RegionOrPoint\", \"PointType\", \"StartPosition\", \"EndPosition\", \"DataSetWID\") VALUES (%d, NULL, NULL, 'N', %s, %s, '%s', %s, %s, %s, %d)", 
	      feature_wid,
	      wh_postgres_string_column(escaped_common_name, common_name_ind),
	      escaped_class,
	      geometry,
	      wh_postgres_string_column(escaped_point_type, point_type_ind),  // NULL unless geometry='point'
	      wh_postgres_int_column(startpos, startpos_ind, start_buffer),
	      wh_postgres_int_column(endpos, endpos_ind, end_buffer),
	      dataset_wid) >= MAX_QUERY_SIZE)
	 wh_postgres_fatal_error("In wh_postgres_insert_into_feature:Query string too long\n");
  wh_postgres_run_query(query, "Failed to insert into Feature Table");
  free(escaped_point_type);
  free(escaped_class);
  free(escaped_common_name);
}

void
wh_postgres_insert_into_gene(int gene_wid, char * common_name, char * unique_id, char direction,
			       char * startpos, short startpos_ind,
			       char * endpos, short endpos_ind,
			       char interrupted, short interrupted_ind) {
  char interrupted_buffer[100];
  char  * escaped_common_name= wh_postgres_malloc_and_create_escaped_string(common_name);
  char * escaped_unique_id = wh_postgres_malloc_and_create_escaped_string(unique_id);
  
  if(snprintf(query, MAX_QUERY_SIZE,
		    "INSERT INTO \"Gene\" (\"WID\",\"Name\",\"GenomeID\",\"CodingRegionStart\",\"CodingRegionEnd\",\"Direction\",\"Interrupted\",\"DataSetWID\") VALUES (%d, %s, %s, %s, %s, '%c', %s, %d)", 
	gene_wid,
	escaped_common_name,
	escaped_unique_id,
	wh_postgres_string_column(startpos, startpos_ind),
	wh_postgres_string_column(endpos, endpos_ind),
	direction,
	wh_postgres_char_column(interrupted, interrupted_ind, interrupted_buffer),
	dataset_wid) >= MAX_QUERY_SIZE)
	 wh_postgres_fatal_error("In wh_postgres_insert_into_gene:Query string too long\n");
  wh_postgres_run_query(query, "Failed to insert into Gene Table");
  free(escaped_common_name);
  free(escaped_unique_id);
}

int
wh_postgres_select_gene(char * name){
  PGresult *result;
  //MYSQL_ROW row;
  int wid = 0;
  char  * escaped_name = wh_postgres_malloc_and_create_escaped_string(name);

 if(snprintf(query,MAX_QUERY_SIZE,
     "SELECT \"WID\" FROM \"Gene\", \"DBID\" WHERE \"XID\" = %s AND \"OtherWID\" = \"WID\" AND \"DataSetWID\" = %d",
		   escaped_name, dataset_wid) >= MAX_QUERY_SIZE)
	wh_postgres_fatal_error("In wh_postgres_select_gene: Query string too long\n");
 
  result = wh_postgres_run_query(query, "Failed to get Wid from Gene Table");

  if (wh_postgres_num_rows(result) > 0) {
	wid = atoi(PQgetvalue(result, 0, 0));
  }
  PQclear(result);
 
 
  free(escaped_name);
  return wid;
}
