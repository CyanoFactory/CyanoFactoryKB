/* $Id: StringColumn.java,v 1.3 2008/10/06 21:31:36 valerie Exp $ */
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
import org.apache.log4j.Logger;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;

/**
 * @author Valerie Wagner
 *         Date: Aug 11, 2004
 */
public class StringColumn extends ColumnBase
{
    private static final Logger log = Logger.getLogger( StringColumn.class );
    private static final String EOL = System.getProperty("line.separator");

    protected String stringValue;


    public StringColumn( String columnName, Table table )
    {
        super( columnName, table );
    }


    public String getValue()
    {
        return stringValue;
    }


    public void setValue( String value )
    {
        stringValue = value;
    }

    public void setValue( Enum value )
    {
        if( value != null)
        {
            stringValue = value.toString();
        }
    }


    public void store( PreparedStatement insertStmt ) throws SQLException
    {
        if( stringValue == null )
        {
            storeNull( insertStmt );
        }
        else
        {
            checkStringLength();
            insertStmt.setString( columnMetaData.getParameterIndex(), stringValue );
        }
    }


    public void update( ResultSet row ) throws SQLException
    {
        checkStringLength();
        if( stringValue != null )
        {
            row.updateString( columnName, stringValue );
        }
    }


    public void load( ResultSet row ) throws SQLException
    {
        stringValue = row.getString( columnName );
    }


    private void checkStringLength()
    {
        int maxColumnLength = columnMetaData.getMaxColumnSize();

        // Check that the value doesn't exceed the maximum allowed size
        // If it does, truncate to fit the allowed size.

        if( stringValue != null && maxColumnLength > 0 &&  stringValue.length() > maxColumnLength )
        {
            final String margin = "    ";
            StringBuilder sb = new StringBuilder();
            sb.append("Value to be inserted has been truncated:");
            sb.append(EOL).append(margin).append("Column name: ").append(fullColumnName);
            sb.append(EOL).append(margin).append("Max column length: ").append(maxColumnLength);
            sb.append(EOL).append(margin).append("Value length: ").append(stringValue.length());
            sb.append(EOL).append(margin).append("Value: <").append(stringValue).append('>');
            log.warn(sb);

            stringValue = stringValue.substring( 0, maxColumnLength );
        }
    }


    public String toString()
    {
        return stringValue == null ? "" : stringValue;
    }


    public boolean equals( Object obj )
    {
        if( obj == this )
        {
            return true;
        }
        else if( obj == null || !(obj instanceof StringColumn) )
        {
            return false;
        }

        StringColumn col = (StringColumn)obj;
        if( stringValue == null )
        {
            return col.stringValue == null;
        }
        return stringValue.equals( col.stringValue );
    }


    public int hashCode()
    {
        if( stringValue == null )
        {
            return 11;
        }
        return 11 + stringValue.hashCode();
    }
}
