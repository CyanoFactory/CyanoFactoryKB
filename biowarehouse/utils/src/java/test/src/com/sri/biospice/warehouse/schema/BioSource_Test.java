/* $Id: BioSource_Test.java,v 1.3 2008/10/06 21:31:37 valerie Exp $ */
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

import com.sri.biospice.warehouse.WarehouseFunctionalTest;
import com.sri.biospice.warehouse.database.WarehouseManager;
import com.sri.biospice.warehouse.schema.object.BioSource;
import java.sql.SQLException;

/**
 * @author Valerie.Wagner@sri.com
 *         Date: Aug 18, 2004
 */
public class BioSource_Test extends WarehouseFunctionalTest {

    public BioSource_Test(String testName) throws SQLException {
        super(testName);
    }


    public void testTooLargeValues() throws SQLException {
        DataSet dataset = new DataSet("test");
        dataset.store();

        BioSource bioSource = new BioSource(dataset.getWID());

        // Create a cellLine vallue that is bigger than the max allowed length
        StringBuffer badCellLine = new StringBuffer();
        int maxLength = bioSource.getTableMetaData().getMaxColumnSize("CellLine");
        badCellLine.setLength(maxLength + 10);
        for (int i = 0; i < badCellLine.length(); i++) {
            badCellLine.setCharAt(i, 'x');
        }

        // Store
        bioSource.setCellLine(badCellLine.toString());
        bioSource.store();

        // Load the biosource from the DB
        BioSource storedBioSource = new BioSource(dataset.getWID(), bioSource.getWID());
        storedBioSource.load();
        assertEquals("CellLine not truncated", maxLength, storedBioSource.getCellLine().length());


        // Try updating to a new valid value
        String newCellLine = "aaa";
        bioSource.setCellLine(newCellLine);
        bioSource.update();
        storedBioSource.load();
        assertEquals("Wrong value for CellLine.", newCellLine, storedBioSource.getCellLine());

        // Try updating to a new invalid length value
        bioSource.setCellLine(badCellLine.toString());
        bioSource.update();
        storedBioSource.load();
        WarehouseManager.getWarehouse().commit();
        assertEquals("CellLine not truncated", maxLength, storedBioSource.getCellLine().length());
    }
}
