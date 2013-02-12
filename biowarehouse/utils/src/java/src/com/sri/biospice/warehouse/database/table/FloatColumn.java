/* $Id: FloatColumn.java,v 1.1 2006/07/07 15:03:36 valerie Exp $ */
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
import java.sql.ResultSet;
import java.sql.SQLException;

/**
 * @author Valerie Wagner
 *         Date: Aug 12, 2004
 */
public class FloatColumn extends ColumnBase
{
    private Float floatValue;


//    public FloatColumn( String columnName, Vector allColumns, TableMetaData tableMetaData )
    public FloatColumn( String columnName, Table table )
    {
        super( columnName, table );
    }


    public void setValue( Float value )
    {
        floatValue = value;
    }


    public void setValue( float value )
    {
        floatValue = new Float( value );
    }


    public void setValue( String value )
    {
        if( value != null )
        {
            floatValue = new Float( value );
        }
    }


    public void store( PreparedStatement insertStmt ) throws SQLException
    {
        if( floatValue == null )
        {
            storeNull( insertStmt );
        }
        else
        {
            insertStmt.setFloat( columnMetaData.getParameterIndex(), floatValue.floatValue() );
        }
    }


    public void update( ResultSet row ) throws SQLException
    {
        if( floatValue != null )
        {
            row.updateFloat( columnName, floatValue.floatValue() );
        }
    }


    public void load( ResultSet row ) throws SQLException
    {
        // The ResultSet.getFloat() method returns zero if
        // the value is NULL.  So first we check if the value is
        // NULL so that we don't accidently get a value of zero.
        Object tempObject = row.getObject( columnName );
        if( tempObject != null )
        {
            floatValue = new Float( row.getFloat( columnName ) );
        }
        else
        {
            floatValue = null;
        }
    }


    public Float getValue()
    {
        return floatValue;
    }


    public String toString()
    {
        return floatValue == null ? "" : floatValue.toString();
    }
}
