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
package com.sri.biospice.warehouse.schema.object;

import com.sri.biospice.warehouse.database.table.*;

import java.util.Date;
import java.sql.Timestamp;

/**
 * Created by IntelliJ IDEA.
 * User: kejariwa
 * Date: Oct 31, 2006
 * Time: 10:16:42 AM
 * To change this template use File | Settings | File Templates.
 */
public class ExperimentData extends ObjectTableBase {
    private LongColumn experimentWIDCol;
    private ClobColumn dataCol;
    private LongColumn mageDataCol;
    private StringColumn roleCol;
    private StringColumn kindCol;
    private TimestampColumn dateProducedCol;
    private LongColumn otherWIDCol;


    private static TableMetaData tableMetaData = null;

    public ExperimentData(long datasetWid) {
        super(datasetWid, "ExperimentData");
        init();
    }


    public enum Kind {
        O, //observation
        P, //parameter
        C, //computed from observation
        M  //metadata describing other data
    }


    /**
     * For an entry loaded from the Gene table in the database
     *
     * @param datasetWID
     * @param thisWID    The WID of the entry in the Gene table
     */
    public ExperimentData(long datasetWID, long thisWID) {
        super(datasetWID, "ExperimentData", thisWID);
        init();
    }

    private void init() {
        experimentWIDCol = new LongColumn("ExperimentWID", this);
        dataCol = new ClobColumn("Data",this);
        mageDataCol = new LongColumn("MageData",this);
        roleCol = new StringColumn("Role",this);
        kindCol = new StringColumn("Kind",this);
        dateProducedCol = new TimestampColumn("DateProduced",this);
        otherWIDCol = new LongColumn("OtherWID",this);
    }

    public Long getExperimentWID() {
        return experimentWIDCol.getValue();
    }

    public void setExperimentWID(Long experimentWID) {
        this.experimentWIDCol.setValue(experimentWID);
    }

    public String getDataCol() {
        return dataCol.getValue();
    }

    public void setData(String data) {
        this.dataCol.setValue(data);
    }

    public String getRole() {
        return roleCol.getValue();
    }

    public void setRole(String role) {
        this.roleCol.setValue(role);
    }

    public Date getDateProduced() {
        return dateProducedCol.getValue();
    }

    public void setDateProduced(Date dateProduced) {
        this.dateProducedCol.setValue(dateProduced);
    }

    public Long getOtherWID() {
        return otherWIDCol.getValue();
    }

    public void setOtherWID(Long otherWID) {
        this.otherWIDCol.setValue(otherWID);
    }


    public String getKind() {
        return kindCol.getValue();
    }

    public void setKind(Kind kind) {
        this.kindCol.setValue(kind);
    }

    public TableMetaData getTableMetaData() {
        return tableMetaData;
    }


    protected void setTableMetaData(TableMetaData tableMetaData) {
        this.tableMetaData = tableMetaData;
    }

    public Enum[] getAllowedValues(String columnName) {

        if (columnName.equalsIgnoreCase(kindCol.getColumnName())) {
            return Kind.values();
        }
        return null;
    }


    public Long getMageData() {
        return mageDataCol.getValue();
    }

    public void setMageData(Long mageData) {
        this.mageDataCol.setValue(mageData);
    }
}
