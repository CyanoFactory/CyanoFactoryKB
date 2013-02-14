/* $Id: Gene_Test.java,v 1.2 2006/07/07 01:06:36 valerie Exp $ */
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
import com.sri.biospice.warehouse.schema.object.Gene;
import java.sql.ResultSet;
import java.sql.SQLException;

/**
 * @author Valerie.Wagner@sri.com
 *         Date: Aug 11, 2004
 */
public class Gene_Test extends WarehouseFunctionalTest {


    public Gene_Test(String testName) throws SQLException {
        super(testName);
    }


    public void testLoadAndReload() throws SQLException {
        String name = "xxx";
        String newName = "yyy";
        String sqlLine;

        // Create a dataset
        DataSet dataset = new DataSet("test");
        dataset.store();

        // Create and store a nucleic acid
        Gene originalGene = new Gene(dataset.getWID());
        originalGene.setName("xxx");
        originalGene.store();
        WarehouseManager.getWarehouse().commit();

        // Check the values in the DB
        sqlLine = "Select * from Gene where WID=" + originalGene.getWID();
        ResultSet rs = WarehouseManager.getWarehouse().executeQuery(sqlLine);
        assertEquals("No result in Result Set", true, rs.next());
        assertEquals("Wrong type", name, rs.getString("Name"));
        rs.close();

        //Change some values and update
        com.sri.biospice.warehouse.schema.object.Gene storedGene = new Gene(dataset.getWID(), originalGene.getWID());
        storedGene.setName(newName);
        storedGene.update();

        WarehouseManager.getWarehouse().commit();
        rs = WarehouseManager.getWarehouse().executeQuery(sqlLine);
        assertEquals("No result in Result Set", true, rs.next());
        assertEquals("Wrong type", newName, rs.getString("Name"));
        rs.close();

        //Change some values and update
        storedGene = new Gene(dataset.getWID(), originalGene.getWID());
        storedGene.setName(null);
        storedGene.update();

        WarehouseManager.getWarehouse().commit();
        rs = WarehouseManager.getWarehouse().executeQuery(sqlLine);
        assertEquals("No result in Result Set", true, rs.next());
        assertEquals("Wrong type", newName, rs.getString("Name"));
        rs.close();
    }

}
