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
 * Structures and function prototypes for the Dnabindsite loader
 * Thomas J Lee, SRI International, June 2008
 */

#ifndef _DNABINDSITE_PARSE_H
#define _DNABINDSITE_PARSE_H 1

/*
 * The dnabindsite entry 
 */

struct dnabindsite_entry {
  int               load_error;       
  char              *unique_id; 
  char              *common_name;
  short 	    common_name_ind;
  char              *abs_center_pos;   // can be fractional
  char              *regulated_promoter;    // used to convert relative to absolute pos

  struct stringlist *component_of; 

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

void dnabindsite_next_entry (void);
void dnabindsite_clear_entry (short freeflag);
void dnabindsite_load_entry (struct dnabindsite_entry *);
void dnabindsite_parse_error (void);

#endif
