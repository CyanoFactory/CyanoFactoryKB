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

import com.sri.biospice.warehouse.database.table.LongColumn;
import com.sri.biospice.warehouse.database.table.StringColumn;
import com.sri.biospice.warehouse.database.table.TableMetaData;
import com.sri.biospice.warehouse.schema.object.ObjectTable;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.Vector;

/**
 * @author Valerie Wagner
 *         Date: Mar 31, 2005
 */
public class Description extends TableBase {

    // Columns
    private LongColumn otherWidCol;
    private StringColumn tableNameCol;
    private StringColumn commCol;

    // Table metadata
    private static TableMetaData tableMetaData;

    public Description() {
        super("Description");
        init();
    }

    public Description(ObjectTable parentTable, String comment) {
        super("Description");
        init();
        otherWidCol.setValue(parentTable.getWID());
        tableNameCol.setValue(parentTable.getTableName());
        commCol.setValue(comment);
    }

    private void init() {
        otherWidCol = new LongColumn("OtherWID", this);
        tableNameCol = new StringColumn("TableName", this);
        commCol = new StringColumn("Comm", this);
    }

    protected void setTableMetaData(TableMetaData tableMetaData) {
        this.tableMetaData = tableMetaData;
    }

    /**
     * Not implemented yet.  Returns null;
     *
     * @return
     * @throws SQLException
     */
    protected ResultSet getRow() throws SQLException {
        return null;
    }

    public TableMetaData getTableMetaData() {
        return tableMetaData;
    }


    public long getOtherWid() {
        return otherWidCol.getValue();
    }


    /**
     * This function is named <code>getTableNameCol()</code> to distinguish
     * it from the <code>getTableName()</code> in <code>ObjectTable</code> that refers to the name of the
     * table ("Description").  Instead it returns the column in the Description
     * table called "TableName".
     *
     * @return Description.TableName
     */
    public String getTableNameCol() {
        return tableNameCol.getValue();
    }

    public String getComm() {
        return commCol.getValue();
    }


    public static Vector loadDescriptions(long otherWID) {
        return TableFactory.loadTables("select * from Description where OtherWID='" + otherWID + "'", TableFactory.DESCRIPTION);
    }
}
