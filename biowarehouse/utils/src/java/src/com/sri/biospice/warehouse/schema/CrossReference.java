/* $Id: CrossReference.java,v 1.1 2006/07/07 15:03:36 valerie Exp $ */
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
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.Vector;


/**
 * An entry in the CrossReference table
 *
 * @author Priyanka Gupta
 * @author Valerie Wagner
 */
public class CrossReference extends TableBase {
    // Columns of the table
    private LongColumn otherWIDCol;
    private LongColumn crossWID;
    private StringColumn xidCol;
    private StringColumn typeCol;
    private StringColumn relationshipCol;
    private StringColumn versionCol;
    private StringColumn databaseNameCol;

    // Table metadata
    private static TableMetaData tableMetaData;

    // Enumerated values used by the Type column
//    public static final String ENUM_TYPE_ACCESSION = "Accession";
//    public static final String ENUM_TYPE_GUID = "GUID";
//    public static final String ENUM_TYPE_NAME = "Name";

    /**
     * Enumerated values allowed for column <i>CrossReference.Type</i>.
     */
    public enum Type {
        Accession,
        GUID
    }


    /**
     * Constructor for the class. Used to store a unique
     * identifier for an object from another database that
     * may or may not be loaded yet.
     *
     * @param xid          An external identifier.
     * @param type         The type of the external identifier.
     *                     i.e what the external identifier represents.
     * @param version      The version of the xid.
     * @param databaseName The name of the dataset this id belongs to.
     * @param otherWid     The wid of the warehouse object that this
     *                     gets linked from.
     * @deprecated todo Use typed constructor instead
     */
    public CrossReference(String xid, String type, String version, String databaseName, long otherWid) {
        super("CrossReference");
        init();
        xidCol.setValue(xid);
        typeCol.setValue(type);
        versionCol.setValue(version);
        databaseNameCol.setValue(databaseName);
        otherWIDCol.setValue(otherWid);
    }

    /**
     * Constructor for the class. Used to store a unique
     * identifier for an object from another database that
     * may or may not be loaded yet.
     *
     * @param xid          An external identifier.
     * @param type         The type of the external identifier.
     *                     i.e what the external identifier represents.
     * @param version      The version of the xid.
     * @param databaseName The name of the dataset this id belongs to.
     * @param otherWid     The wid of the warehouse object that this
     *                     gets linked from.
     */
    public CrossReference(String xid, Type type, String version, String databaseName, long otherWid) {
        super("CrossReference");
        init();
        xidCol.setValue(xid);
        typeCol.setValue(type);
        versionCol.setValue(version);
        databaseNameCol.setValue(databaseName);
        otherWIDCol.setValue(otherWid);
    }


    /**
     * Creates an empty CrossReference
     */
    public CrossReference() {
        super("CrossReference");
        init();
    }


    private void init() {
        otherWIDCol = new LongColumn("OtherWID", this);
        crossWID = new LongColumn("CrossWID", this);
        xidCol = new StringColumn("XID", this);
        typeCol = new StringColumn("Type", this);
        relationshipCol = new StringColumn("Relationship", this);
        versionCol = new StringColumn("Version", this);
        databaseNameCol = new StringColumn("DatabaseName", this);
    }


    public TableMetaData getTableMetaData() {
        return tableMetaData;
    }


    protected void setTableMetaData(TableMetaData tableMetaData) {
        this.tableMetaData = tableMetaData;
    }


    /**
     * todo not implemented
     *
     * @return
     * @throws SQLException
     */
    protected ResultSet getRow() throws SQLException {
        return null;
    }


    public Enum[] getAllowedValues(String columnName) {
        if (columnName.equalsIgnoreCase(typeCol.getColumnName())) {
            return Type.values();
        }

        return null;
    }

    public long getOtherWID() {
        return otherWIDCol.getValue();
    }

    public void setOtherWID(long otherWID) {
        this.otherWIDCol.setValue(otherWID);
    }

    public long getCrossWID() {
        return crossWID.getValue();
    }

    public void setCrossWID(long crossWID) {
        this.crossWID.setValue(crossWID);
    }

    public String getXid() {
        return xidCol.getValue();
    }

    public void setXid(String xid) {
        this.xidCol.setValue(xid);
    }

    public String getType() {
        return typeCol.getValue();
    }

    public void setType(Type type) {
        this.typeCol.setValue(type);
    }

    public String getVersion() {
        return versionCol.getValue();
    }

    public void setVersion(String version) {
        this.versionCol.setValue(version);
    }

    public String getDatabaseName() {
        return databaseNameCol.getValue();
    }

    public void setDatabaseName(String databaseName) {
        this.databaseNameCol.setValue(databaseName);
    }

    public static Vector loadCrossReferences(long otherWID) {
        String query = "Select * from CrossReference where OtherWID='" + otherWID + "'";
        return TableFactory.loadTables(query, TableFactory.CROSS_REFERENCE);
    }

}
