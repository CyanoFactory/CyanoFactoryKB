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
 * Structures and function prototypes for the BsubCyc Pubs loader
 * Thomas J Lee, SRI International, August 2002
 *
 * (C) 2002 Defense Advanced Research Projects Agency
 */

#ifndef _PUBS_PARSE_H
#define _PUBS_PARSE_H 1

/*
 * The pubs entry 
 */

struct pubs_entry {
  int               load_error;       
  char              *unique_id;
  char              *referent_frame;  // alias for another frame
  char              *medline_uid;
  char              *pubmed_id;
  char              *source;
  char              *title;
  char              *url;
  char              *year;
  struct stringlist *authors;           
  struct stringlist *comments;        
};  


/*
 * Function prototypes
 */

void pubs_next_entry (void);
void pubs_clear_entry (short freeflag);
void pubs_load_entry (struct pubs_entry *);
void pubs_parse_error (void);

#endif
