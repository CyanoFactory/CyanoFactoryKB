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
package com.sri.bw;

import com.sri.biospice.warehouse.database.WarehouseManager;
import org.apache.log4j.Logger;
import java.sql.SQLException;
import java.util.Hashtable;

/**
 * @author Valerie Wagner
 *         Date: Jul 26, 2006
 */
public class TableInsert {

    private static final Logger log = Logger.getLogger(TableInsert.class);
    private GenericTable table;
    private Hashtable<String, Object> col2value = new Hashtable<String, Object>();

    public TableInsert(GenericTable table) {
        this.table = table;
    }

    public void put(String columnName, Object object) {
        if (columnName == null) {
            log.warn("put(): Column name is null for object " + object + " in table " + this);
            LoaderStatistics.warningOccurred();
        } else if (object != null) {
            col2value.put(columnName, object);
        }
    }

    public Object get(String columnName) {
        return col2value.get(columnName);
    }

    public String toString() {
        return table.getTableName() + " " + col2value;
    }

    public Hashtable<String, Object> getCol2value() {
        return col2value;
    }

    public void store() throws SQLException {
        try {
            table.store(this);
        } catch (SQLException e) {
            log.error("Error storing " + this.toString());
            log.error("Message is " + e.getMessage());
            throw e;
        }
    }

    /**
     * @param targetInsert Already-stored object row
     * @param columnName
     */
    public String linkFrom(TableInsert targetInsert, String columnName) throws SQLException {
//        log.info("Connecting " + targetInsert.getTableName() + "." + columnName + " to " + this.getTableName());
        targetInsert.put(columnName, this.getPrimaryKeyValue());
        if (targetInsert.table instanceof ObjectTable) {
            String update = "Update " + targetInsert.table.getTableName() + " set " +
                    columnName + " = " + this.getPrimaryKeyValue() + " where "
                    + targetInsert.table.getPrimaryKeyColumnName() + "=" + targetInsert.getPrimaryKeyValue();
//            log.info("update = " + update);
            WarehouseManager.getWarehouse().executeUpdate(update);
            return update;
        } else {
            log.error("Tried to link to non-object table");
            LoaderStatistics.errorOccurred();
        }
        return null;
    }

    public void linkTo(TableInsert targetInsert, String columnName) {
//        log.info("Connecting " + this.getTableName() + "." + columnName + " to " + targetInsert.getTableName());
        this.put(columnName, targetInsert.getPrimaryKeyValue());
    }

    public String linkTo(Long linkWID, String columnName) throws SQLException {
//        log.info("Connecting " + this.getTableName() + "." + columnName + " to " + linkWID);
        this.put(columnName, linkWID);
        String update = "Update " + this.table.getTableName() + " set " +
                columnName + " = " + linkWID + " where "
                + this.table.getPrimaryKeyColumnName() + "=" + this.getPrimaryKeyValue();
//        log.info("update = " + update);
        WarehouseManager.getWarehouse().executeUpdate(update);
        return update;
    }

    public Object getPrimaryKeyValue() {
        return get(table.getPrimaryKeyColumnName());
    }

    public String getTableName() {
        return table.getTableName();
    }

    public String getPrimaryKeyName(){
        return table.getPrimaryKeyColumnName();
    }

    public GenericTable getTable() {
        return table;
    }
}

