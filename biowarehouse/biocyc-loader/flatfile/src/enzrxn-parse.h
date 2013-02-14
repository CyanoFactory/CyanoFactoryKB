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
 * Structures and function prototypes for the Enzrxn loader
 * Thomas J Lee, SRI International, August 2002
 */

#ifndef _ENZRXN_PARSE_H
#define _ENZRXN_PARSE_H 1

/*
 * The enzrxn entry 
 */

struct enzrxn_entry {
  int               load_error;       
  char              *unique_id;
  char              *common_name;        
  char              *complex;               
  char              *enzyme;               
  char              *reaction;
  char              *reaction_direction;

  /* Cofactor-related attributes */
  struct stringlist *cofactors;              /* Non-prosthetic-group cofactors */
  struct stringlist *cofactor_prosthetics;   /* Prosthetic-group cofactors */
  struct stringlist *cofactor_unknowns;      /* Cofactors which may or may not be Prosthetic groups */ 
  struct stringlist *cofactor_comments;      /* COFACTOR-BINDING-COMMENT: Comments that apply to all cofactors, and to the ER itself */
  
  /* Alternate cofactor/substrate attributes - first in list is primary, rest are alternates to it */
  struct stringlist *alt_cofactors;              /* names of alternate cofactors */
  struct stringlist *alt_substrates;             /* names of alternate substrates */

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

void enzrxn_next_entry (void);
void enzrxn_clear_entry (short freeflag);
void enzrxn_load_entry (struct enzrxn_entry *);
void enzrxn_parse_error (void);

#endif
