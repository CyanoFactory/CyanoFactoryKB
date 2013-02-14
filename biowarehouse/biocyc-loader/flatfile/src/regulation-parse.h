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
 * Regulation loader
 * Thomas J Lee, SRI International, August 2007
 */

#ifndef _REGULATION_PARSE_H
#define _REGULATION_PARSE_H 1

/*
 * The regulation entry 
 */

struct regulation_entry {
  int               load_error;       
  char              *unique_id; 
  char              *type;       // either 'Regulation-of-Enzyme-Activity' or 'Regulation-of-Transcription-Initiation'
  char              *mechanism;  // eg/ ':ALLOSTERIC'
  char              *mode;  // +, -, or NIL
  char              *physiologically_relevant;  // T or NIL (generally all specified values are T)
  char              *regulated_entity;   // ID of enzymatic rxn or promoter.
  char              *regulator;          // ID of compound or protein.


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

void regulation_next_entry (void);
void regulation_clear_entry (short freeflag);
void regulation_load_entry (struct regulation_entry *);
void regulation_parse_error (void);

#endif
