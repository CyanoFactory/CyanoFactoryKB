/* $Id: ColumnBase.java,v 1.1 2006/07/07 15:03:36 valerie Exp $ */
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
import java.sql.PreparedStatement;
import java.sql.SQLException;

/**
 * @author Valerie Wagner
 *         Date: Aug 11, 2004
 */
abstract class ColumnBase implements Column
{
    protected String columnName;
    protected String fullColumnName;
    protected ColumnMetaData columnMetaData;


    // todo: we need to make it impossible to instantiate a column that doesn't exist in the
    // database instance
    ColumnBase( String columnName, Table table )
    {
        this.columnMetaData = table.getTableMetaData().getColumn( columnName );
        this.columnName = columnName;
        this.fullColumnName = table.getTableName() + "." + columnName;
        table.getColumns().add( this );
    }


    public String getColumnName()
    {
        return columnName;
    }


    protected void storeNull( PreparedStatement insertStmt ) throws SQLException
    {
        insertStmt.setNull( columnMetaData.getParameterIndex(),
                            columnMetaData.getColumnSQLType(),
                            columnMetaData.getColumnTypeName() );
    }
}
