/* $Id: TimestampColumn.java,v 1.1 2006/07/07 15:03:36 valerie Exp $ */
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
import java.sql.Timestamp;
import java.util.Date;

/**
 * @author Valerie Wagner
 *         Date: Aug 12, 2004
 */
public class TimestampColumn extends ColumnBase {

    private Timestamp timestamp;


    public TimestampColumn(String columnName, Table table) {
        super(columnName, table);
    }


    public void setValue(Date date) {
        timestamp = new Timestamp(date.getTime());
    }


    public Date getValue() {
        return timestamp == null ? null : new Date(timestamp.getTime());
    }


    public void store(PreparedStatement insertStmt) throws SQLException {
        if (timestamp == null) {
            storeNull(insertStmt);
        } else {
            insertStmt.setTimestamp(columnMetaData.getParameterIndex(), timestamp);
        }
    }


    public void update(ResultSet row) throws SQLException {
        if (timestamp != null) {
            row.updateTimestamp(columnName, timestamp);
        }
    }


    public void load(ResultSet row) throws SQLException {
        timestamp = row.getTimestamp(columnName);
    }


    public String toString() {
        return timestamp == null ? "" : timestamp.toString();
    }

}
