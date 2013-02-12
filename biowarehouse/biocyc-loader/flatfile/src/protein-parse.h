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
 * Structures and function prototypes for the BsubCyc Protein loader
 * Thomas J Lee, SRI International, August 2002
 */

#ifndef _PROTEIN_PARSE_H
#define _PROTEIN_PARSE_H 1

/*
 * The protein entry 
 */

struct protein_entry {
  /* AAsequence is loaded from another file */
  int               load_error;       
  char              *unique_id;          /* referenced by pathways */    
  char              *common_name;
  char              *gene;               /*! not used */
  char              *molecular_weight;
  char              *molecular_weight_exp;
  char              *pi;                /* pI, isoelectric point? */
  short 	    common_name_ind;
  short 	    molecular_weight_ind;
  short		    molecular_weight_exp_ind;
  short		    pi_ind;

  struct stringlist *components;        /* list of subunits of a complex */
  struct stringlist *coeffs;            /* corresponding list of #occurrences in complex, default 1 */
  
  struct stringlist *go_terms;           
  struct stringlist *locations;

  /* DBLINKS attribures - lists correspond */
  struct stringlist *dblinks_dbs;
  struct stringlist *dblinks_ids;
    
  struct stringlist *citations;        
  struct stringlist *synonyms;           
  struct stringlist *comments;        
};  


/*
 * Function prototypes
 */

void protein_next_entry (void);
void protein_clear_entry (short freeflag);
void protein_load_entry (struct protein_entry *);
void protein_parse_error (void);

#endif
