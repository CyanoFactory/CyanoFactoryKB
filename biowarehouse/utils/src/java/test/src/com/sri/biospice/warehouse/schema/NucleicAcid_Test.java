/* $Id: NucleicAcid_Test.java,v 1.2 2006/07/07 01:06:36 valerie Exp $ */
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
import com.sri.biospice.warehouse.schema.object.NucleicAcid;
import java.sql.ResultSet;
import java.sql.SQLException;

/**
 * @author Valerie.Wagner@sri.com
 *         Date: Aug 10, 2004
 */
public class NucleicAcid_Test extends WarehouseFunctionalTest {
    public NucleicAcid_Test(String testName) throws SQLException {
        super(testName);
    }


    public void testLoadAndReload() throws SQLException {
        NucleicAcid.Type type = NucleicAcid.Type.DNA;
        NucleicAcid.Strandedness strandedness = NucleicAcid.Strandedness.ss;
        NucleicAcid.Type newType = NucleicAcid.Type.RNA;
        NucleicAcid.Strandedness newStrand = NucleicAcid.Strandedness.ds;

        String sqlLine;

        // Create a dataset
        DataSet dataset = new DataSet("test");
        dataset.store();

        // Create and store a nucleic acid
        NucleicAcid originalNA = new NucleicAcid(dataset.getWID());
        originalNA.setType(NucleicAcid.Type.DNA);
        originalNA.setStrandedness(strandedness);
        originalNA.setClass(NucleicAcid.Class.plasmid);
        originalNA.setSequenceDerivation(NucleicAcid.SequenceDerivation.map);
        originalNA.setTopology(NucleicAcid.Topology.linear);
        originalNA.setMoleculeLength("100");
        originalNA.setMoleculeLengthApproximate("xxx");
//        originalNA.setGeneticCodeWID( "111" );
        originalNA.store();
        WarehouseManager.getWarehouse().commit();

        // Check the values in the DB
        sqlLine = "Select * from NucleicAcid where WID=" + originalNA.getWID();
        ResultSet rs = WarehouseManager.getWarehouse().executeQuery(sqlLine);
        assertEquals("No result in Result Set 1", true, rs.next());
        assertEquals("Wrong type", type.toString(), rs.getString("Type"));
        assertEquals("Wrong strandedness", strandedness.toString(), rs.getString("strandedness"));
        rs.close();

        //Change some values and update
        NucleicAcid storedNA = new NucleicAcid(dataset.getWID(), originalNA.getWID());
        storedNA.load();
        assertEquals("Wrong type after reload", type.toString(), storedNA.getType());
        assertEquals("Wrong strandedness after reload", strandedness.toString(), storedNA.getStrandedness());
        storedNA.setType(newType);
        storedNA.setStrandedness(newStrand);
        storedNA.update();

        WarehouseManager.getWarehouse().commit();
        rs = WarehouseManager.getWarehouse().executeQuery(sqlLine);
        assertEquals("No result in Result Set 2", true, rs.next());
        assertEquals("Wrong type", newType.toString(), rs.getString("Type"));
        assertEquals("Wrong strandedness", newStrand.toString(), rs.getString("strandedness"));
        rs.close();

        //Change some values and update
        storedNA = new NucleicAcid(dataset.getWID(), originalNA.getWID());
        storedNA.setType(null);
        storedNA.setStrandedness(null);
        storedNA.update();

        WarehouseManager.getWarehouse().commit();
        rs = WarehouseManager.getWarehouse().executeQuery(sqlLine);
        assertEquals("No result in Result Set 3", true, rs.next());
        assertEquals("Wrong type", newType.toString(), rs.getString("Type"));
        assertEquals("Wrong strandedness", newStrand.toString(), rs.getString("strandedness"));
        rs.close();
    }


    /**
     * Test that the minimum required columns are:
     * Type, Class, WID, DatasetWID
     */
    public void testMinimumEntry() throws SQLException {
        // Create a dataset
        DataSet dataset = new DataSet("test");
        dataset.store();

        // Create and store a nucleic acid
        NucleicAcid nucleicAcid = new NucleicAcid(dataset.getWID());
        nucleicAcid.setType(NucleicAcid.Type.DNA);
        nucleicAcid.setClass(NucleicAcid.Class.plasmid);
        nucleicAcid.store();
    }
}
