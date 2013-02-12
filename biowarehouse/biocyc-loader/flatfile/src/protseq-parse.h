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
 * Structures and function prototypes for the BsubCyc Protseq loader
 * Thomas J Lee, SRI International, August 2002
 */

#ifndef _PROTSEQ_PARSE_H
#define _PROTSEQ_PARSE_H 1

/*
 * The protseq entry 
 */

struct protseq_entry {
  int               load_error;
  char              *protein;          /* unique_id of the protein */
  char              *aa_sequence;      /* Amino Acid sequence of a protein */ 
};  


/*
 * Function prototypes
 */

void protseq_next_entry (void);
void protseq_clear_entry (short freeflag);
void protseq_load_entry (struct protseq_entry *);
void protseq_parse_error (void);

#endif
