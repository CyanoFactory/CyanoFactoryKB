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

import com.sri.biospice.warehouse.schema.DataSet;
import java.sql.SQLException;
import java.util.Date;

/**
 * @author Valerie Wagner
 *         Date: Jul 27, 2006
 */
public class ObjectDescriptorTable extends GenericTable {

    public ObjectDescriptorTable(DataSet dataset, String tableName) throws SQLException {
        super(dataset, tableName);
    }

    // todo: mechanism to test for required values


    public TableInsert newInsert() {
        TableInsert insertRow = super.newInsert();
        return insertRow;
    }

    public String getPrimaryKeyColumnName() {
        return "OtherWID";
    }
}
