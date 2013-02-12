/* $Id: TableMetaDataImpl.java,v 1.1 2006/07/07 15:03:36 valerie Exp $ */
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

import com.sri.biospice.warehouse.database.Warehouse;
import com.sri.biospice.warehouse.database.WarehouseManager;
import org.apache.log4j.Logger;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.Enumeration;
import java.util.Hashtable;

/**
 * Contains metadata about a single table, such as the size and type of each column.
 * @author Valerie Wagner
 *         Date: Jun 29, 2004
 */
public class TableMetaDataImpl  implements TableMetaData

{
    private static final Logger log = Logger.getLogger( TableMetaDataImpl.class );

    private Warehouse warehouse;
    private Hashtable columnMetaDataTable;
    private PreparedStatement prepStmt;
    private String tableName;
    private String sqlInsertLine;


    /**
     * This constructor accesses the Warehouse and obtains the metadata
     * for all columns associated with this table.  This operation has a large
     * overhead, so all attempts should be made to only do it once per table
     * per program.
     * @param tableName The name of the table to obtain metadata about
     */
    public TableMetaDataImpl( String tableName )
    {
        this.tableName = tableName;
        columnMetaDataTable = new Hashtable();
        ColumnMetaData colMetaData;

        // Read the metadata about this table from the Warehouse
        warehouse = WarehouseManager.getWarehouse();
        ResultSet columnSet = warehouse.getTableMetaData( tableName );

        if( columnSet != null )
        {

            // Read the metadata for each column
            try
            {
                while( columnSet.next() )
                {
                    colMetaData = new ColumnMetaDataImpl( columnSet );
                    columnMetaDataTable.put( colMetaData.getColumnName().toUpperCase(), colMetaData );
                }
            }
            catch( SQLException e )
            {
                log.error( "Error creating metadata for table " + tableName, e );
            }
        }

        try
        {
            createPreparedInsert();
        }
        catch( SQLException e )
        {
            log.error( "Error creating prepared statement for table " + tableName, e );
        }

    }


    /**
     * Creates a prepared statement to do an insert into this table.
     * @throws SQLException
     */
    private void createPreparedInsert() throws SQLException
    {
        sqlInsertLine = "INSERT INTO " + tableName + " (";
        int index = 1;

        Enumeration columns = columnMetaDataTable.elements();
        while( columns.hasMoreElements() )
        {
            ColumnMetaData columnMetaData = (ColumnMetaData)columns.nextElement();
            sqlInsertLine += columnMetaData.getColumnName();
            if( columns.hasMoreElements() )
            {
                sqlInsertLine += ", ";
            }
            columnMetaData.setParameterIndex( index++ );
        }
        sqlInsertLine += ") VALUES (";

        for( int i = 1; i < index; i++ )
        {
            sqlInsertLine += "?";
            if( i < index - 1 )
            {
                sqlInsertLine += ", ";
            }
        }
        sqlInsertLine += ")";

//        log.debug( "Prepared insert statment is: " + sqlInsertLine );

        prepStmt = warehouse.createPreparedStatement( sqlInsertLine );
//        prepInsertStmt = new PreparedInsertStatement( prepStmt, this );
    }


    /**
     * Returns the maximum size allowed for this column.  If the column
     * is a text column, this is the maximum length of the text.  If the
     * column is a numeric column, this is the precision of the number.
     * @param columnName The name of the column in this table
     * @return The maximum size allowed for an entry in this column
     */
    public int getMaxColumnSize( String columnName )
    {
        int size = -1;
        ColumnMetaData column = (ColumnMetaData)columnMetaDataTable.get( columnName.toUpperCase() );
        if( column != null )
        {
            size = column.getMaxColumnSize();
        }

        return size;
    }


    /**
     * @return The number of columns in the table
     */
    public int getNumberOfColumns()
    {
        return columnMetaDataTable.size();
    }


    public PreparedStatement getPreparedInsertSatement()
    {
        return prepStmt;
    }


    /**
     * @param columnName The column name
     * @return The metadata for this column
     */
    public ColumnMetaData getColumn( String columnName )
    {
        return (ColumnMetaData)columnMetaDataTable.get( columnName.toUpperCase() );
    }


    /**
     * @return A list of all columns in the table
     */
    public Enumeration getColumns()
    {
        return columnMetaDataTable.elements();
    }


    public String toString()
    {
        String description = "Table meta data for table: ";

        description += tableName + "\n";
        description += "Number of columns: " + getNumberOfColumns() + "\n";
        description += "Prepared insert statement format: " + sqlInsertLine + "\n";

        Enumeration columns = getColumns();
        while( columns.hasMoreElements() )
        {
            ColumnMetaData column = (ColumnMetaData)columns.nextElement();
            description += column.toString() + "\n";
        }

        return description;
    }


    /**
     * The INSERT String used to create the PreparedStatement for this class.
     * Used only for convience/logging/debugging.
     * @return A string representation of the PreparedStatement syntax
     */
    public String getInsertString()
    {
        return sqlInsertLine;
    }


    public boolean hasColumn( String name )
    {
        return columnMetaDataTable.containsKey( name.toUpperCase() );
    }


    public String getTableName()
    {
        return tableName;
    }


    public Hashtable getAllowedValues( String columnName )
    {
        String sqlLine = "Select Enumeration.Value from Enumeration where " +
                         "TableName='" + tableName + "' and " +
                         "ColumnName='" + columnName + "'";

        Hashtable allowedValues = null;

        try
        {
            ResultSet results = WarehouseManager.getWarehouse().executeQuery( sqlLine );
            if( results != null && results.next() )
            {
                allowedValues = new Hashtable();
                do
                {
                    String value = results.getString( "Value" );
                    allowedValues.put( value, value );
                }
                while( results.next() );
            }
        }
        catch( SQLException e )
        {
        }

        return allowedValues;
    }
}
