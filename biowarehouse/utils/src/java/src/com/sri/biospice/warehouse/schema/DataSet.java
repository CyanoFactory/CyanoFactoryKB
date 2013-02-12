/* $Id: DataSet.java,v 1.4 2008/10/15 21:48:44 valerie Exp $ */
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

import com.sri.biospice.warehouse.database.WarehouseManager;
import com.sri.biospice.warehouse.database.table.LongColumn;
import com.sri.biospice.warehouse.database.table.StringColumn;
import com.sri.biospice.warehouse.database.table.TableMetaData;
import com.sri.biospice.warehouse.database.table.TimestampColumn;
import com.sri.biospice.warehouse.schema.linking.DataSetHierarchy;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.Date;


/**
 * An entry in the DataSet table
 *
 * @author Priyanka Gupta
 * @author Valerie Wagner
 */
public class DataSet extends TableBase {
    // Columns in the DataSet table
    private LongColumn widCol;
    private StringColumn nameCol;
    private StringColumn versionCol;
    private TimestampColumn loadDateCol;
    private TimestampColumn changeDateCol;
    private TimestampColumn releaseDateCol;
    private StringColumn homeURLCol;
    private StringColumn queryURLCol;
    private StringColumn loadedByCol;
    private StringColumn applicationCol;
    private StringColumn applicationVersionCol;

    private static TableMetaData tableMetaData;
    private DataSetHierarchy datasetHierarchy;


    /**
     * Defines a new entry in the DataSet table.
     * Sets the load date to the current date/time.
     * Sets UserName to the system user name.
     *
     * @param datasetName The name of the dataset being loaded.
     */
    public DataSet(String datasetName) {
        super("DataSet");
        init();
        widCol.setValue(WarehouseManager.getWarehouse().getNextSpecialWID());
        nameCol.setValue(datasetName);
        loadDateCol.setValue(new Date());

        setLoadedBy(System.getProperty("user.name"));
        datasetHierarchy = new DataSetHierarchy(this);
    }

    /**
     * Represents a row in the DataSet table that already exists in the database
     * @param datasetWID WID of the existing DataSet entry
     */
    public DataSet(long datasetWID) {
        super("DataSet");
        init();
        widCol.setValue(datasetWID);
    }


    protected void setTableMetaData(TableMetaData tableMetaData) {
        this.tableMetaData = tableMetaData;
    }


    private void init() {
        widCol = new LongColumn("WID", this);
        nameCol = new StringColumn("Name", this);
        versionCol = new StringColumn("Version", this);
        loadDateCol = new TimestampColumn("LoadDate", this);
        changeDateCol = new TimestampColumn("ChangeDate", this);
        releaseDateCol = new TimestampColumn("ReleaseDate", this);
        homeURLCol = new StringColumn("HomeURL", this);
        queryURLCol = new StringColumn("QueryURL", this);
        loadedByCol = new StringColumn("LoadedBy", this);
        applicationCol = new StringColumn("Application", this);
        applicationVersionCol = new StringColumn("ApplicationVersion", this);
    }


    /**
     * Gets the dataset Wid
     *
     * @return Returns the dataset wid.
     */
    public long getWID() {
        return widCol.getValue();
    }


    public TableMetaData getTableMetaData() {
        return tableMetaData;
    }


    protected ResultSet getRow() throws SQLException {
        ResultSet entryRS = null;
        String sqlLine = "select " + tableName + ".* from " + tableName + " where WID=" + widCol.getValue();

        entryRS = WarehouseManager.getWarehouse().executeQuery(sqlLine);

        if (entryRS == null || entryRS.next() == false) {
            SQLException e = new SQLException("Entry not found in the " + tableName + " table with WID = " + widCol.getValue());
            throw e;
        }

        return entryRS;
    }


    public void setVersion(String version) {
        versionCol.setValue(version);
    }


    public void setReleaseDate(Date date) {
        releaseDateCol.setValue(date);
    }


    public void setHomeURL(String url) {
        homeURLCol.setValue(url);
    }


    public void setQueryURL(String url) {
        queryURLCol.setValue(url);
    }


    public void setChangeDate(Date changeDate) {
        changeDateCol.setValue(changeDate);
    }


    public String getLoadedBy() {
        return loadedByCol.getValue();
    }

    public void setLoadedBy(String loadedBy) {
        loadedByCol.setValue(loadedBy);
    }

    public String getApplication() {
        return applicationCol.getValue();
    }

    public void setApplication(String application) {
        applicationCol.setValue(application);
    }

    public String getApplicationVersion() {
        return applicationVersionCol.getValue();
    }

    public void setApplicationVersion(String applicationversion) {
        applicationVersionCol.setValue(applicationversion);
    }

    public Date getLoadDate() {
        return loadDateCol.getValue();
    }

    public void setName(String name) {
        nameCol.setValue(name);
    }


    @Override
    public void store() throws SQLException {
        super.store();
        if (datasetHierarchy == null) {
            datasetHierarchy = new DataSetHierarchy(this);
        }
        datasetHierarchy.store();
    }
}
