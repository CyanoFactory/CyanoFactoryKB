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
 * Structures and function prototypes for the BsubCyc Pathway loader
 * Thomas J Lee, SRI International, August 2002
 */

#ifndef _PATHWAY_PARSE_H
#define _PATHWAY_PARSE_H 1

struct pathway_link {
  char                *chemical;          /* links the pathways */    
  struct stringlist   *pathways;
  struct pathway_link *next;            
};  

/*
 * The pathway entry 
 */

struct pathway_entry {
  int                 load_error;       
  char                *unique_id;          /* referenced by pathways */    
  char                *common_name;
  struct stringlist   *reactions;         /* REACTION-LIST, compnents are reaction or pathway UNIQUE-IDs */
  /* successors & predecessors are corresponding lists. A null predecessor is indicated by a "" */
  struct stringlist   *successors;     
  struct stringlist   *predecessors;     
  struct stringlist   *temppreds;      /* holds predecessor rxns for a single successor while one predecessors is being parsed */
  struct stringlist   *sub_pathways;       
  struct stringlist   *super_pathways;           
  struct stringlist   *hypothetical_rxns; /* UNIQUE-IDs of reactions */
  struct stringlist   *templist;      /* holds pathway names while one pathway_link is being parsed */
  struct pathway_link *links;         /* each is a linking chemical and the list of pathways it links */

  /* DBLINKS attribures - lists correspond */
  struct stringlist *dblinks_dbs;
  struct stringlist *dblinks_ids;
    
  struct stringlist   *citations;        
  struct stringlist   *synonyms;
  struct stringlist   *comments;        
};  

/*
 * Function prototypes
 */

void pathway_next_entry (void);
void pathway_clear_entry (short freeflag);
void pathway_load_entry (struct pathway_entry *);
void pathway_parse_error (void);

#endif
