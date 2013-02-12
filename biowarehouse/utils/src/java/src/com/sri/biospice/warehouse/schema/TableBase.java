/* $Id: TableBase.java,v 1.1 2006/07/07 15:03:36 valerie Exp $ */
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

import com.sri.biospice.warehouse.database.table.Column;
import com.sri.biospice.warehouse.database.table.TableMetaData;
import com.sri.biospice.warehouse.database.table.TableMetaDataImpl;
import org.apache.log4j.Logger;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.Vector;

/**
 * @author Valerie Wagner
 *         Date: Aug 12, 2004
 */
public abstract class TableBase implements Table {
    private static final Logger log = Logger.getLogger(TableBase.class);

    protected String tableName;
    private Vector allColumns;


    public TableBase(String tableName) {
        this.tableName = tableName;
        allColumns = new Vector();

        if (getTableMetaData() == null) {
            setTableMetaData( new TableMetaDataImpl( tableName ) );
        }
    }


    protected abstract void setTableMetaData(TableMetaData tableMetaData);


    public void update() throws SQLException {
        ResultSet rowResultSet = null;
        log.debug("Updating " + tableName);
        try {
            rowResultSet = getRow();

            if (rowResultSet != null) {

                // Update each column value (columns will all check for null values)
                for (int i = 0; i < allColumns.size(); i++) {
                    ((Column) allColumns.elementAt(i)).update(rowResultSet);
                }

                // Actually update the row in the database
                rowResultSet.updateRow();
                rowResultSet.close();
            } else {
                log.error("Tried to update table " + tableName + ", but getRow() is not implemented for this table.");
            }
        } catch (SQLException e) {
            // Try to close the result set
            if (rowResultSet != null) {
                rowResultSet.close();
            }
            throw e;
        }
    }


    public void load() throws SQLException {
        ResultSet row = getRow();
        if (row != null) {
            for (int i = 0; i < allColumns.size(); i++) {
                ((Column) allColumns.elementAt(i)).load(row);
            }
        } else {
            log.error("Tried to load table " + tableName + ", but getRow() is not implemented for this table");
        }
    }


    public void store() throws SQLException {
        try {
            doInsert();
        } catch (SQLException e) {
            log.error("Error storing table " + tableName, e);
            log.error(toString());
            throw e;
        }
    }


    protected void doInsert() throws SQLException {
        if (log.isDebugEnabled()) {
            log.debug("Storing " + toString());
        }
        PreparedStatement insert = null;
        try {
            insert = getTableMetaData().getPreparedInsertSatement();
            insert.clearWarnings();
            insert.clearParameters();
            for (int i = 0; i < allColumns.size(); i++) {
                ((Column) allColumns.elementAt(i)).store(insert);
            }
            insert.execute();
        } catch (SQLException error1) {
            log.error("Error while storing table " + tableName, error1);
            SQLException error2 = new SQLException("Error storing table " + tableName);
            throw error2;
        }
    }


    public Vector getColumns() {
        return allColumns;
    }


    public String getTableName() {
        return tableName;
    }


    protected abstract ResultSet getRow() throws SQLException;


    public String toString() {
        StringBuffer description = new StringBuffer();
        description.append(tableName).append(": ");

        for (int i = 0; i < allColumns.size(); i++) {
            Column column = (Column) allColumns.elementAt(i);
            description.append(column.getColumnName())
                    .append("=")
                    .append(column.toString())
                    .append(", ");
        }
        return description.toString();
    }


    public Enum [] getAllowedValues(String columnName) {
        return null;
    }

}
