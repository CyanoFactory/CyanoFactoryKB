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
#ifndef _WIDTABLE_H
#define _WIDTABLE_H 1

#include <stdlib.h>

#define WIDSLOTS 19999

/* Use this value on widtable_insert() to remember that a name is undefined */
#define WIDERROR -1


/*
 *  Hash table item structure
 */
struct wid_object
{
  int               wid;
  char              *name;
  struct wid_object *next;
};

typedef struct wid_object * WIDTABLE[WIDSLOTS];

typedef struct wid_object * WIDNAMELIST;

void
widtable_init(WIDTABLE hashtbl);

int
widtable_size(WIDTABLE hashtbl);

int
widtable_lookup(WIDTABLE hashtbl, char *key);

int
widtable_update(WIDTABLE hashtbl, int new_wid, char *name);

void
widtable_insert(WIDTABLE hashtbl, int wid, char *name);

void
widtable_make_undefined(WIDTABLE hashtbl, char *name);

WIDNAMELIST
widnamelist_insert(WIDNAMELIST list, int wid, char *name);


#endif
