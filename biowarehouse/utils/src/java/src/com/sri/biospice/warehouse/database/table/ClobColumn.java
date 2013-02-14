/* $Id: ClobColumn.java,v 1.3 2006/11/02 22:06:58 kejariwa Exp $ */
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
import java.io.Reader;
import java.io.StringReader;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;

/**
 * Represents a value in a column of type CLOB.
 * The ClobColumn behaves the same as the StringColumn, except that it
 * does not check whether the string value exceeds the maximum allowed
 * value for this column.
 * @author Valerie Wagner
 *         Date: Aug 16, 2004
 */
public class ClobColumn extends StringColumn
{
    public ClobColumn( String columnName, Table table )
    {
        super( columnName, table );
    }


    public void store( PreparedStatement insertStmt ) throws SQLException
    {
        if( stringValue == null )
        {
            storeNull( insertStmt );
        }
        else
        {
            Reader clobReader = new StringReader( stringValue );
            int clobLength = stringValue.length();
            insertStmt.setCharacterStream( columnMetaData.getParameterIndex(), clobReader, clobLength );
        }
    }

    public void load( ResultSet row ) throws SQLException
    {
        Object tempObject = row.getObject( columnName );
        if (tempObject != null)
            stringValue = row.getClob(columnName ).toString();
        else
            stringValue = null;
    }



    public void update( ResultSet row ) throws SQLException
    {
        if( stringValue != null )
        {
            Reader clobReader = new StringReader( stringValue );
            int clobLength = stringValue.length();
            row.updateCharacterStream( columnName, clobReader, clobLength );
        }
    }

}
