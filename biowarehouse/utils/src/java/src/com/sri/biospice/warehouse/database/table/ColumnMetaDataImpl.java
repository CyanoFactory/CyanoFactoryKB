/* $Id: ColumnMetaDataImpl.java,v 1.2 2006/07/27 18:19:42 valerie Exp $ */
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

import com.sri.biospice.warehouse.schema.Table;
import java.sql.ResultSet;
import java.sql.SQLException;

/**
 * Holds metadata about one column in one table in a database
 * @author Valerie Wagner
 *         Date: Jun 29, 2004
 */
public class ColumnMetaDataImpl implements ColumnMetaData
{
    private String columnName;
    private int columnSQLType;
    private int maxColumnSize;
    private int parameterIndex;
    private String columnTypeName;
    private boolean valueIsSet = false;
    private Column column;


    /**
     * @param columnSet The result set that contains information about a table,
     *                  which must be currently pointing to a single column.
     * @throws SQLException
     */
    ColumnMetaDataImpl( ResultSet columnSet ) throws SQLException
    {
        columnName = columnSet.getString( "COLUMN_NAME" );
        columnSQLType = columnSet.getInt( "DATA_TYPE" );
        maxColumnSize = columnSet.getInt( "COLUMN_SIZE" );
        columnTypeName = columnSet.getString( "TYPE_NAME" );
    }


    public String getColumnName()
    {
        return columnName;
    }


    public int getColumnSQLType()
    {
        return columnSQLType;
    }


    public int getMaxColumnSize()
    {
        return maxColumnSize;
    }


    public int getParameterIndex()
    {
        return parameterIndex;
    }


    public void setParameterIndex( int parameterIndex )
    {
        this.parameterIndex = parameterIndex;
    }


    public String getColumnTypeName()
    {
        return columnTypeName;
    }


    public void valueSatisfied()
    {
        valueIsSet = true;
    }


    public boolean isSatisifed()
    {
        return valueIsSet;
    }


    public void clearValue()
    {
        valueIsSet = false;
    }

    public Column makeColumn(Table table) {
        switch (columnSQLType)
        {
            case -5:
                column= new LongColumn(columnName, table);
                break;
            case 12:
                column= new StringColumn(columnName, table);
                break;
        }
        return column;
    }


    public String toString()
    {
        String description = "";

        description += "Column: " + columnName;
        description += "\tSQL Type: " + columnTypeName;
        description += "\tSQL Type Code: " + columnSQLType;
        description += "\tMax Column Entry Size: " + maxColumnSize;
        description += "\tPreparedStatement Parameter Index: " + parameterIndex;

        return description;
    }
}
