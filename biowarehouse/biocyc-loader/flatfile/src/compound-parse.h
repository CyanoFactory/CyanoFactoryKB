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
 * Structures and function prototypes for the BsubCyc Compound loader
 * Thomas J Lee, SRI International, August 2002
 */

#ifndef _COMPOUND_PARSE_H
#define _COMPOUND_PARSE_H 1

/*
 * The compound entry 
 */

struct compound_entry {
  int               load_error;       
  char              *unique_id;       
  char              *common_name;
  char              *systematic_name;
  char              *cas_registry_numbers;
  char              *smiles;             /* eg. "c1(c(CC=O)nc[nH]1)" */
  char              *formula;            /* eg. "C12H22O11" */   
  char              *charge;             /* short */
  char              *molecular_weight;   /* double */
  char              *pka1;               /* double */
  char              *pka2;               /* double */
  char              *pka3;               /* double */
  short 	    systematic_name_ind;
  short        	    cas_registry_numbers_ind;
  short 	    smiles_ind;
  short             formula_ind;
  short		    charge_ind;
  short		    molecular_weight_ind;
  short		    pka1_ind;
  short		    pka2_ind;
  short		    pka3_ind;
  
  /* DBLINKS attribures - lists correspond */
  struct stringlist *dblinks_dbs;
  struct stringlist *dblinks_ids;
    
  struct stringlist *citations;        
  struct stringlist *comments;        
  struct stringlist *synonyms;        
};  


/*
 * Function prototypes
 */

void compound_next_entry (void);
void compound_clear_entry (short freeflag);
void compound_load_entry (struct compound_entry *);
void compound_parse_error (void);
void add_element_to_formula(char *element, char *number);
int find_chemical (char *);

#endif
