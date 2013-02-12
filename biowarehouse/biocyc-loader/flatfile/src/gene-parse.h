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
 * Structures and function prototypes for the BsubCyc Gene loader
 * Thomas J Lee, SRI International, August 2002
 */

#ifndef _GENE_PARSE_H
#define _GENE_PARSE_H 1

/*
 * The gene entry 
 */

struct gene_entry {
  int               load_error;       
  char              *unique_id;
  char              *common_name;
  char              *evidence;
  char              *left_end_position;
  char              *right_end_position;
  char              *interrupted;
  char              *transcription_direction;
  struct stringlist *component_of;  // indicates transcription unit, as well as other ingnored links
  struct stringlist *products;           
  struct stringlist *types;           
  
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

void gene_next_entry (void);
void gene_clear_entry (short freeflag);
void gene_load_entry (struct gene_entry *);
void gene_parse_error (void);

#endif
