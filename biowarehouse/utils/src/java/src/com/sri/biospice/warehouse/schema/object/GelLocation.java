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

/**
 * @author Anish Kejariwal
 *         Date: Nov 1, 2006
 */
public class GelLocation extends ObjectTableBase{
    // Columns of the table
    private LongColumn spotWIDCol;
    private FloatColumn xCoordCol;
    private FloatColumn yCoordCol;
    private StringColumn refGelCol;
    private LongColumn experimentWIDCol;

    private static TableMetaData tableMetaData = null;    

    public GelLocation(long datasetWid) {
        super(datasetWid, "GelLocation");
        init();
    }


    /**
     * Represents an already-existing entry in the Feature table
     *
     * @param datasetWid
     * @param gelLocationWID the WID of the entry in the Feature table
     */
    public GelLocation(long datasetWid, long gelLocationWID) {
        super(datasetWid, "GelLocation", gelLocationWID);
        init();
    }



    private void init() {
        spotWIDCol = new LongColumn("SpotWID", this);
        xCoordCol = new FloatColumn("Xcoord", this);
        yCoordCol = new FloatColumn("Ycoord", this);
        refGelCol = new StringColumn("refGel", this);
        experimentWIDCol = new LongColumn("ExperimentWID", this);
    }


    public Long getSpotWID() {
        return spotWIDCol.getValue();
    }

    public void setSpotWID(Long spotWID) {
        this.spotWIDCol.setValue(spotWID);
    }

    public Float getxCoord() {
        return xCoordCol.getValue();
    }

    public void setxCoord(Float xCoord) {
        this.xCoordCol.setValue(xCoord);
    }

    public Float getyCoord() {
        return yCoordCol.getValue();
    }

    public void setyCoord(Float yCoord) {
        this.yCoordCol.setValue(yCoord);
    }

    public String getRefGel() {
        return refGelCol.getValue();
    }

    public void setRefGel(String refGel) {
        this.refGelCol.setValue(refGel);
    }

    public Long getExperimentWID() {
        return experimentWIDCol.getValue();
    }

    public void setExperimentWID(Long experimentWID) {
        this.experimentWIDCol.setValue(experimentWID);
    }

    public TableMetaData getTableMetaData() {
        return tableMetaData;
    }


    protected void setTableMetaData(TableMetaData tableMetaData) {
        this.tableMetaData = tableMetaData;
    }
}
