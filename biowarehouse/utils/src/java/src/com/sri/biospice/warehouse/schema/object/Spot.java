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
public class Spot extends ObjectTableBase {
    private StringColumn spotIdCol;
    private FloatColumn molecularWeightEstCol;
    private StringColumn pIEstCol;

    private static TableMetaData tableMetaData = null;


    /**
     * Creates a new entry in the Spot table
     *
     * @param datasetWid
     */
    public Spot(long datasetWid) {
        super(datasetWid, "Spot");
        init();
    }

    /**
     * For an entry loaded from the Gene table in the database
     *
     * @param datasetWID
     * @param thisWID    The WID of the entry in the Gene table
     */
    public Spot(long datasetWID, long thisWID) {
        super(datasetWID, "Spot", thisWID);
        init();
    }


    private void init() {
        spotIdCol = new StringColumn("SpotId", this);
        molecularWeightEstCol = new FloatColumn("MolecularWeightEst", this);
        pIEstCol = new StringColumn("PIEst", this);
    }


    public String getSpotId() {
        return spotIdCol.getValue();
    }

    public void setSpotId(String spotId) {
        this.spotIdCol.setValue(spotId);
    }

    public Float getMolecularWeightEst() {
        return molecularWeightEstCol.getValue();
    }

    public void setMolecularWeightEst(Float molecularWeightEst) {
        this.molecularWeightEstCol.setValue(molecularWeightEst);
    }

    public String getpIEst() {
        return pIEstCol.getValue();
    }

    public void setpIEst(String pIEst) {
        this.pIEstCol.setValue(pIEst);
    }

    public TableMetaData getTableMetaData() {
        return tableMetaData;
    }


    protected void setTableMetaData(TableMetaData tableMetaData) {
        this.tableMetaData = tableMetaData;
    }
}
