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

import com.sri.bw.dbschema.SchemaDocument;
import com.sri.bw.dbschema.TableDocument;
import com.sri.bw.dbschema.ColumnDocument;
import com.sri.bw.dbschema.ForeignKeyDocument;
import com.sri.bw.LoaderStatistics;
import org.apache.xmlbeans.XmlCursor;
import org.apache.log4j.Logger;

/**
 * @author Valerie Wagner
 *         Date: Aug 27, 2006
 */
public class SchemaUtils {

    private static final Logger log = Logger.getLogger(SchemaUtils.class);

    public static final int NOT_FOUND = -100;

    public static void removeTable(SchemaDocument.Schema schema, String tableName) {
        int index = findTableIndex(schema, tableName);
        if (index != NOT_FOUND) {
            schema.removeTable(index);
        }
    }

    private static int findTableIndex(SchemaDocument.Schema schema, String tableName) {
        for (int i = 0; i < schema.getTableArray().length; i++) {
            if (tableName.equals(schema.getTableArray(i).getName())) {
                return i;
            }
        }
        log.warn("Table not found: " + tableName);
        LoaderStatistics.warningOccurred();
        return NOT_FOUND;
    }

    public static boolean hasTable(SchemaDocument.Schema schema, String tableName) {
        return findTableIndex(schema, tableName) != NOT_FOUND;
    }

    public static TableDocument.Table findTable(SchemaDocument.Schema schema, String tableName) {
        int index = findTableIndex(schema, tableName);
        if (index == NOT_FOUND) {
            return null;
        }
        return schema.getTableArray(index);
    }

    private static int findColumnIndex(TableDocument.Table table, String columnName) {
        for (int i = 0; i < table.getColumnArray().length; i++) {
            if (columnName.equals(table.getColumnArray(i).getName())) {
                return i;
            }
        }
        return NOT_FOUND;
    }

    public static boolean hasColumn(TableDocument.Table table, String columnName) {
        return findColumnIndex(table, columnName) != NOT_FOUND;
    }

    public static boolean hasPrimaryKey(TableDocument.Table table, String columnName) {
        return (table.isSetPrimaryKey() && columnName.equals(table.getPrimaryKey().getColumn()));
    }

    public static void removeAllColumns(SchemaDocument.Schema schema, String columnName) {
        for (TableDocument.Table table : schema.getTableArray()) {
            removeColumn(table, columnName);
        }
    }

    public static void removeColumn(TableDocument.Table table, String columnName) {
        int index = findColumnIndex(table, columnName);
        if (index != NOT_FOUND) {
            table.removeColumn(index);
        }

    }

    public static ColumnDocument.Column findColumn(SchemaDocument.Schema schema, String tableName, String columnName) {
        TableDocument.Table table = findTable(schema, tableName);
        if (table != null) {
            return findColumn(table, columnName);
        }
        return null;
    }

    public static ColumnDocument.Column findColumn(TableDocument.Table table, String columnName) {
        int index = findColumnIndex(table, columnName);
        if (index != NOT_FOUND) {
            return table.getColumnArray(index);
        }
        return null;
    }

    public static int getNumForeignKeys(TableDocument.Table table) {
        XmlCursor cursor = table.newCursor();
        cursor.selectPath(".//foreignKey");
        return cursor.getSelectionCount();
    }

    /**
     *
     * @param schema
     * @param targetTableName
     * @param newTableName
     * @param newColumnName
     */
    public static void changeForeignKeyTarget(SchemaDocument.Schema schema, String targetTableName, String newTableName, String newColumnName) {
        XmlCursor cursor = schema.newCursor();
        cursor.push();
        cursor.selectPath(".//foreignKey");
        while (cursor.toNextSelection()) {
            ForeignKeyDocument.ForeignKey fk = (ForeignKeyDocument.ForeignKey)cursor.getObject();
            if (fk.getToTable().equals(targetTableName)) {
                fk.setToTable(newTableName);
                fk.setToColumn(newColumnName);
            }
        }
        cursor.pop();
    }

    public static void removeForeignKeyTarget(SchemaDocument.Schema schema, String targetTableName) {
        XmlCursor cursor = schema.newCursor();
        cursor.push();
        cursor.selectPath(".//foreignKey");
        while (cursor.toNextSelection()) {
            ForeignKeyDocument.ForeignKey fk = (ForeignKeyDocument.ForeignKey)cursor.getObject();
            if (fk.getToTable().equals(targetTableName)) {
                cursor.removeXml();
            }
        }
        cursor.pop();
    }
}
