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
import com.sri.biospice.warehouse.database.table.TableMetaData;
import com.sri.biospice.warehouse.database.table.ColumnMetaData;
import com.sri.biospice.warehouse.schema.DataSet;
import com.sri.biospice.warehouse.schema.TableBase;
import org.apache.log4j.Logger;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.Enumeration;

/**
 * @author Valerie Wagner
 *         Date: Jul 26, 2006
 */
public class GenericTable extends TableBase implements Comparable {

    private static final Logger log = Logger.getLogger(GenericTable.class);

    private TableMetaData tableMetaData;
    private final PreparedStatement insert;
    protected final long datasetWID;
    protected final WIDGenerator widGenerator;
    protected final DataSet dataset;

    private static final int commitFrequency = 100;
    private static int commitCount = 0;

    public GenericTable(DataSet dataset, String tableName) throws SQLException {
        super(tableName);
        this.dataset = dataset;
        datasetWID = dataset.getWID();

        // The DataSet table needs the special wid generator
        if (dataset.getTableName().equals(tableName)) {
            widGenerator = WarehouseManager.getWarehouse().getSpecialWIDGenerator();
        } else {
            widGenerator = WarehouseManager.getWarehouse().getWIDGenerator();
        }

        insert = tableMetaData.getPreparedInsertSatement();
        clearValues();
    }


    protected void afterStore() throws SQLException {
        clearValues();
    }

    private void clearValues() throws SQLException {
//        log.info("clearing values");
        for (int i = 1; i <= tableMetaData.getNumberOfColumns(); i++) {
            insert.setObject(i, null);
        }
    }

    // todo: mechanism to store statistics

    public void store(TableInsert insertRow) throws SQLException {
        insert(insertRow);
        afterStore();
    }

    protected void insert(TableInsert insertRow) throws SQLException {
        Enumeration keys = insertRow.getCol2value().keys();
//        log.info("Storing " + insertRow);

        while (keys.hasMoreElements()) {
            String colName = (String)keys.nextElement();
            ColumnMetaData colMeta = tableMetaData.getColumn(colName);
            if (colMeta != null) {
                try {
                    Object value = insertRow.get(colName);
                    if (value instanceof String) {
                        String valStr = (String)value;
                        if (valStr.length() > colMeta.getMaxColumnSize()) {
                            value = valStr.subSequence(0, colMeta.getMaxColumnSize());
                        }
                    }
                    insert.setObject(colMeta.getParameterIndex(), value);
                } catch (SQLException e) {
                    log.error("Error setting object on prepared statement for column " + colName +", value " + insertRow.get(colName));
                    throw e;
                }
            } else {
                log.warn("No column named '" + colName + "' in table " + tableName);
                LoaderStatistics.warningOccurred();
            }
        }
        try {
            insert.executeUpdate();
            if (commitCount++ % commitFrequency == 0) {
                WarehouseManager.getWarehouse().commit();
            }
            LoaderStatistics.insertPerformed(this);
        } catch (SQLException e) {
            log.error("Error storing " + insertRow);
            LoaderStatistics.errorOccurred();
            throw e;
        }

    }


    public TableMetaData getTableMetaData() {
        return tableMetaData;
    }


    public void setTableMetaData(TableMetaData tableMetaData) {
        this.tableMetaData = tableMetaData;
    }

    protected ResultSet getRow() throws SQLException {
        return null;
    }

    public String toString() {
        return tableName;

    }

    public TableInsert newInsert() {
        return new TableInsert(this);
    }

    public String getPrimaryKeyColumnName() {
        return null;
    }

    public int compareTo(Object o) {
        if (o instanceof GenericTable) {
            GenericTable compare = (GenericTable)o;
            return tableName.compareTo(compare.tableName);
        }
        return 0;
    }
}