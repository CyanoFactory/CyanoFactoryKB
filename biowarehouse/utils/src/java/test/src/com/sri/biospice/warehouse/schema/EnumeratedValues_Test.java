/* $Id: EnumeratedValues_Test.java,v 1.2 2006/07/07 01:06:36 valerie Exp $ */
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
package com.sri.biospice.warehouse.schema;

import com.sri.biospice.warehouse.WarehouseFunctionalTest;
import com.sri.biospice.warehouse.database.table.Column;
import com.sri.biospice.warehouse.util.WarehouseTestUtils;
import org.apache.log4j.Logger;
import java.sql.SQLException;
import java.util.Hashtable;
import java.util.Vector;

/**
 * @author Valerie.Wagner@sri.com
 *         Date: Aug 20, 2004
 */
public class EnumeratedValues_Test extends WarehouseFunctionalTest {
    private static Logger log = Logger.getLogger(EnumeratedValues_Test.class);


    public EnumeratedValues_Test(String testName) throws SQLException {
        super(testName);
    }


    public void testAllTables() {
        int[] tableTypes = TableFactory.getTabletypes();
        for (int i = 0; i < tableTypes.length; i++) {
            int tableType = tableTypes[i];
            if (tableType != TableFactory.DATASET) {
                Table table = TableFactory.createTestTable(tableType);
                assertNotNull("TableFactory return null table for type " + tableType, table);
                testOneTable(table);
            }
        }
    }

    private void testOneTable(Table table) {
        log.debug("Testing " + table.getTableName());

        Vector allColumns = table.getColumns();
        for (int i = 0; i < allColumns.size(); i++) {
            Column column = (Column)allColumns.elementAt(i);
            testOneTableOneColumn(table, column.getColumnName());
        }
    }

//    private void oldtestOneTableOneColumn(Table table, String columnName) {
//        Vector allowedTypes = table.getAllowedValues(columnName);
//        Hashtable dbAllowedTypes = table.getTableMetaData().getAllowedValues(columnName);
//
//        if (allowedTypes == null) {
//            assertNull(table.getTableName() + " class does not define types for column " + columnName +
//                    ".  But there are defined types in the database Enumeration table.", dbAllowedTypes);
////            log.debug( "Both are null" );
//        } else if (dbAllowedTypes == null) {
//            assertNull(table.getTableName() + " class defines type for column " + columnName +
//                    ", but no enumerations were found in the Enumeration table for this column",
//                    allowedTypes);
//        } else {
//            assertEquals("The number of enumerated values in the " + table.getTableName() + "." + columnName +
//                    " definition doesn't match what is the Enumeration table",
//                    dbAllowedTypes.size(), allowedTypes.size());
//
//            for (int i = 0; i < allowedTypes.size(); i++) {
//                String classValue = (String)allowedTypes.elementAt(i);
//                assertEquals("This value in the " + table.getTableName() + "class definition doesn't match the Enumeration table",
//                        classValue, (String)dbAllowedTypes.get(classValue));
////                log.debug( "Matched value: " + classValue );
//            }
//        }
//    }

    private void testOneTableOneColumn(Table table, String columnName) {
        Enum[] allowedTypes = table.getAllowedValues(columnName);
        Hashtable dbAllowedTypes = table.getTableMetaData().getAllowedValues(columnName);

        if (allowedTypes == null) {
            assertNull(table.getTableName() + " class does not define types for column " + columnName +
                    ".  But there are defined types in the database Enumeration table.", dbAllowedTypes);
//            log.debug( "Both are null" );
        } else if (dbAllowedTypes == null) {
            assertNull(table.getTableName() + " class defines type for column " + columnName +
                    ", but no enumerations were found in the Enumeration table for this column",
                    allowedTypes);
        } else {
            assertEquals("The number of enumerated values in the " + table.getTableName() + "." + columnName +
                    " definition doesn't match what is the Enumeration table",
                    dbAllowedTypes.size(), allowedTypes.length);

            for (int i = 0; i < allowedTypes.length; i++) {
                String classValue = allowedTypes[i].toString();
                assertEquals("This value in the " + table.getTableName() + "class definition doesn't match the Enumeration table",
                        classValue, (String)dbAllowedTypes.get(classValue));
//                log.debug( "Matched value: " + classValue );
            }
        }
    }
}
