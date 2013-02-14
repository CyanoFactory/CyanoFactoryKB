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
package com.sri.biospice.warehouse.database.table;

import java.sql.PreparedStatement;
import java.util.Enumeration;
import java.util.Hashtable;

/**
 * @author Valerie Wagner
 *         Date: Apr 2, 2005
 */
public interface TableMetaData {
    /**
     * Returns the maximum size allowed for this column.  If the column
     * is a text column, this is the maximum length of the text.  If the
     * column is a numeric column, this is the precision of the number.
     * @param columnName The name of the column in this table
     * @return The maximum size allowed for an entry in this column
     */
    int getMaxColumnSize( String columnName );

    /**
     * @return The number of columns in the table
     */
    int getNumberOfColumns();

    PreparedStatement getPreparedInsertSatement();

    /**
     * @param columnName The column name
     * @return The metadata for this column
     */
    ColumnMetaData getColumn( String columnName );

    /**
     * @return A list of all columns in the table
     */
    Enumeration getColumns();

    /**
     * The INSERT String used to create the PreparedStatement for this class.
     * Used only for convience/logging/debugging.
     * @return A string representation of the PreparedStatement syntax
     */
    String getInsertString();

    boolean hasColumn( String name );

    String getTableName();

    Hashtable getAllowedValues( String columnName );
}
