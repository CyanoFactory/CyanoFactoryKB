/* $Id: LongColumn.java,v 1.1 2006/07/07 15:03:36 valerie Exp $ */
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
 *         Date: Aug 11, 2004
 */
public class LongColumn extends ColumnBase
{
    private Long longValue;


    public LongColumn( String columnName, Table table )
    {
        super( columnName, table );
    }


    public void setValue( Long value )
    {
        longValue = value;
    }


    public void setValue( long value )
    {
        longValue = new Long( value );
    }


    public void setValue( String value )
    {
        if( value != null )
        {
            longValue = new Long( value );
        }
    }


    public void store( PreparedStatement insertStmt ) throws SQLException
    {
        if( longValue == null )
        {
            storeNull( insertStmt );
        }
        else
        {
            insertStmt.setLong( columnMetaData.getParameterIndex(), longValue.longValue() );
        }
    }


    public void update( ResultSet row ) throws SQLException
    {
        if( longValue != null )
        {
            row.updateLong( columnName, longValue.longValue() );
        }
    }


    public void load( ResultSet row ) throws SQLException
    {
        // The ResultSet.getLong() method returns zero if
        // the value is NULL.  So first we check if the value is
        // NULL so that we don't accidently get a value of zero.
        Object tempObject = row.getObject( columnName );
        if( tempObject != null )
        {
            longValue = new Long( row.getLong( columnName ) );
        }
        else
        {
            longValue = null;
        }
    }


    public long getValue()
    {
        return longValue.longValue();
    }


    public Long getLongValue()
    {
        return longValue;
    }


    public String toString()
    {
        return longValue == null ? "" : longValue.toString();
    }


    public boolean equals( Object obj )
    {
        if( obj == this )
        {
            return true;
        }
        else if( obj == null || !(obj instanceof LongColumn) )
        {
            return false;
        }
        LongColumn col = (LongColumn)obj;
        if( longValue == null )
        {
            return col.longValue == null;
        }
        return longValue.equals( col.longValue );
    }


    public int hashCode()
    {
        if( longValue == null )
        {
            return 13;
        }
        return 13 + longValue.hashCode();
    }
}
