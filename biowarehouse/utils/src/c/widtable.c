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
#include <stdlib.h>
#include <string.h>
#include "widtable.h"


int
hash_string(char *string) {
  int sx;
  int result = 0;
  for (sx=0; sx<(strlen(string)); sx++) {
    result += (int) string[sx] << (sx % 16);
  }
  return(result);
}

void
widtable_init(WIDTABLE hashtbl) {
  int i;
  for (i=0; i<WIDSLOTS; i++) {
    hashtbl[i] = NULL;
  }
}

int
widtable_size(WIDTABLE hashtbl) {
  int slot;
  struct wid_object *tmp;
  int size = 0;

  for (slot = 0; slot < WIDSLOTS; slot++) {
    tmp = hashtbl[slot];
    while (tmp) {
      size++;
      tmp = tmp->next;
    }
  }
  return size;
}

int
widtable_lookup(WIDTABLE hashtbl, char *name) {
  int slot;
  struct wid_object *tmp;

  if (name == NULL) return 0;
  slot = hash_string(name) % WIDSLOTS;
  tmp = hashtbl[slot];
  while (tmp) {
    if (0 == strcmp(tmp->name, name)) {
      return tmp->wid;
    }
    tmp = tmp->next;
  }
  return 0;
}

int
widtable_update(WIDTABLE hashtbl, int new_wid, char *name) {
  int slot;
  struct wid_object *tmp;

  slot = hash_string(name) % WIDSLOTS;
  tmp = hashtbl[slot];
  while (tmp) {
    if (0 == strcmp(tmp->name, name)) {
      tmp->wid = new_wid;
      return new_wid;
    }
    tmp = tmp->next;
  }
  return 0;
}

void
widtable_insert(WIDTABLE hashtbl, int wid, char *name) {
  /* Inserts name:wid pair into hash tbl. Makes copy of name,
     possibly causing memory leak.
  */
  int slot;
  struct wid_object *tmp;

  slot = hash_string(name) % WIDSLOTS;
  if ((tmp = malloc(sizeof(struct wid_object))) == NULL) {
    fatal_error("malloc failed in widtable_insert");
  }
  tmp->wid = wid;
  tmp->next = hashtbl[slot];
  tmp->name = strdup(name);
  hashtbl[slot] = tmp;
  return;
}

void
widtable_make_undefined(WIDTABLE hashtbl, char *name) {
  /* Marks name as being undefined by associating a distinguished WID with it.
   */
  widtable_insert(hashtbl, WIDERROR, name);
}


/////////////////////////////////////////////////////////////////////////////
///  WIDNAMELIST operations
/////////////////////////////////////////////////////////////////////////////

WIDNAMELIST
widnamelist_insert(WIDNAMELIST list, int wid, char *name) {
  // Inserts the pair (wid, name) at the front of WIDNAMELIST.
  // Returns the new head of the list.
  
  struct wid_object *tmp;
  if ((tmp = malloc(sizeof(struct wid_object))) == NULL) {
    fatal_error("malloc failed in widnamelist_insert");
  }
  tmp->wid = wid;
  tmp->name = strdup(name);
  tmp->next = list;

  return tmp;
}

