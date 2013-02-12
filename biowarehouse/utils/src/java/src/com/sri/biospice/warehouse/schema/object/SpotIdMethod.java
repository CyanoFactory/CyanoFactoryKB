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

import com.sri.biospice.warehouse.database.table.StringColumn;
import com.sri.biospice.warehouse.database.table.TableMetaData;

/**
 * @author Anish Kejariwal
 *         Date: Nov 1, 2006
 */
public class SpotIdMethod extends ObjectTableBase{
    private StringColumn methodNameCol;
    private StringColumn methodDescCol;
    private StringColumn methodAbbrevCol;

    private static TableMetaData tableMetaData = null;

    /**
      * Creates a new entry in the Gene table
      *
      * @param datasetWid
      */
    public SpotIdMethod(long datasetWid) {
        super(datasetWid, "SpotIdMethod");
        init();
    }

    /**
     * For an entry loaded from the Gene table in the database
     *
     * @param datasetWID
     * @param thisWID    The WID of the entry in the Gene table
     */
    public SpotIdMethod(long datasetWID, long thisWID) {
        super(datasetWID, "SpotIdMethod", thisWID);
        init();
    }


    private void init() {
        methodNameCol = new StringColumn("MethodName", this);
        methodDescCol = new StringColumn("MethodDesc", this);
        methodAbbrevCol = new StringColumn("MethodAbbrev", this);
    }

    public String getMethodName() {
        return methodNameCol.getValue();
    }

    public void setMethodName(String methodNameCol) {
        this.methodNameCol.setValue(methodNameCol);
    }

    public String getMethodDesc() {
        return methodDescCol.getValue();
    }

    public void setMethodDesc(String methodDesc) {
        this.methodDescCol.setValue(methodDesc);
    }

    public String getMethodAbbrev() {
        return methodAbbrevCol.getValue();
    }

    public void setMethodAbbrev(String methodAbbrev) {
        this.methodAbbrevCol.setValue(methodAbbrev);
    }

    public TableMetaData getTableMetaData() {
        return tableMetaData;
    }


    protected void setTableMetaData(TableMetaData tableMetaData) {
        this.tableMetaData = tableMetaData;
    }


}
