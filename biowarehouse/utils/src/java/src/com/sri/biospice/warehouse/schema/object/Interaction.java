/*
    #=========================================================================
    # Copyright 2006 SRI International.  All rights reserved.
    #
    # The material contained in this file is confidential and proprietary to SRI
    # International and may not be reproduced, published, or disclosed to others
    # without authorization from SRI International.
    #
    # DISCLAIMER OF WARRANTIES
    #
    # SRI International MAKES NO REPRESENTATIONS OR WARRANTIES ABOUT THE
    # SUITABILITY OF THE SOFTWARE, EITHER EXPRESS OR IMPLIED, INCLUDING BUT NOT
    # LIMITED TO THE IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
    # PARTICULAR PURPOSE, OR NON-INFRINGEMENT. SRI International SHALL NOT BE
    # LIABLE FOR ANY DAMAGES SUFFERED BY LICENSEE AS A RESULT OF USING, MODIFYING
    # OR DISTRIBUTING THIS SOFTWARE OR ITS DERIVATIVES
    #=========================================================================
*/
package com.sri.biospice.warehouse.schema.object;

import com.sri.biospice.warehouse.database.table.StringColumn;
import com.sri.biospice.warehouse.database.table.TableMetaData;

/**
 * Interaction
 *
 * @author David Dunkley
 * @version 0.1 Nov 10, 2006
 */
public class Interaction extends ObjectTableBase {
    private static final String TABLE_NAME = "Interaction";
    private static final String NAME_COL_NAME = "Name";

    private static TableMetaData tableMetaData;
    private StringColumn nameCol;
    private StringColumn typeCol;

    public Interaction(long datasetWID) {
        super(datasetWID, TABLE_NAME);
        init();
    }

    public Interaction(long datasetWID, long interactionWID) {
        super(datasetWID, TABLE_NAME, interactionWID);
        init();
    }


    private void init() {
        nameCol = new StringColumn(NAME_COL_NAME, this);
        typeCol = new StringColumn("Type", this);
    }


    public void setName(String name) {
        nameCol.setValue(name);
    }

    public String getName() {
        return nameCol.getValue();
    }

    public String getType() {
        return typeCol.getValue();
    }

    public void setType(String type) {
        typeCol.setValue(type);
    }

    protected void setTableMetaData(TableMetaData tableMetaData) {
        this.tableMetaData = tableMetaData;
    }

    public TableMetaData getTableMetaData() {
        return tableMetaData;
    }
}
