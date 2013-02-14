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
package com.sri.biospice.warehouse.util;

import com.sri.biospice.warehouse.schema.object.ObjectTable;
import com.sri.biospice.warehouse.schema.*;
import java.sql.SQLException;
import java.util.Vector;

/**
 * @author Valerie Wagner
 *         Date: Apr 9, 2005
 */
public class ObjectInspector {

    public static void inspectObject(long datasetWID, long wid, int tableName) throws SQLException {

        Vector tables;

        ObjectTable theObjectTable = TableFactory.getObjectTable(datasetWID, wid, tableName);
        theObjectTable.load();
        System.out.println(theObjectTable);

        // Entry
        tables = Entry.loadEntries(theObjectTable.getWID());
        printTables(tables);

        // DBID
        tables = DBID.loadDBIDs(theObjectTable.getWID());
        printTables(tables);

        // CrossReference
        tables = CrossReference.loadCrossReferences(theObjectTable.getWID());
        printTables(tables);

        // CommentTable
        tables = CommentTable.loadCommentTables(theObjectTable.getWID());
        printTables(tables);

        // Description
        tables = Description.loadDescriptions(theObjectTable.getWID());
        printTables(tables);

        // SynonymTable
        tables = SynonymTable.loadSynonyms(theObjectTable.getWID());
        printTables(tables);

        // TermRelationship
        tables = TermRelationship.loadTermRelationships(theObjectTable.getWID());
        printTables(tables);



    }

    private static void printTables(Vector tables) {
        for (int i = 0; i < tables.size(); i++) {
            System.out.println(tables.elementAt(i));
        }
    }
}
