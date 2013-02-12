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
 * Structures and function prototypes for the Terminator loader
 * Thomas J Lee, SRI International, June 2008
 */

#ifndef _TERMINATOR_PARSE_H
#define _TERMINATOR_PARSE_H 1

/*
 * The terminator entry 
 */

struct terminator_entry {
  int               load_error;       
  char              *unique_id; 
  char              *common_name;
  short 	    common_name_ind;
  char              *right_end_position;
  char              *left_end_position;    

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

void terminator_next_entry (void);
void terminator_clear_entry (short freeflag);
void terminator_load_entry (struct terminator_entry *);
void terminator_parse_error (void);

#endif
