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
 * Structures and function prototypes for the BsubCyc Reaction loader
 * Thomas J Lee, SRI International, August 2002
 */

#ifndef _REACTION_PARSE_H
#define _REACTION_PARSE_H 1

/*
 * The reaction entry 
 */

enum side {LEFT = 0, RIGHT = 1, SIDE_ERROR = 2};

struct reaction_entry {
  int               load_error;       
  char              *unique_id;          /* referenced by pathways */    
  char              *ec_number;          /* eg. 4.1.1.71 or 4.1.1.- (partial EC) */
  char              *official_ec_flag;   /* YES or NO (or NIL?) */
  char              *delta_g0;           /* always NIL in BsubCyc */
  char              *spontaneous;        /* always NIL in BsubCyc */
  enum side         side_last_parsed;    /* used to match coefficient to either left or right compound */
  struct stringlist *left;               /* reactants */
  struct stringlist *left_coeffs;        /* corresponding reactant coefficients */
  struct stringlist *right;              /* products */
  struct stringlist *right_coeffs;       /* corresponding product coefficients */
  
  /* DBLINKS attribures - lists correspond */
  struct stringlist *dblinks_dbs;
  struct stringlist *dblinks_ids;
    
  struct stringlist *citations;        
  struct stringlist *synonyms;           /* SYNONYMS and COMMON-NAME */
  struct stringlist *comments;        
};  


/*
 * Function prototypes
 */

void reaction_next_entry (void);
void reaction_clear_entry (short freeflag);
void reaction_load_entry (struct reaction_entry *);
void reaction_parse_error (void);

#endif
