/* $Id: LinkingTables_Test.java,v 1.2 2006/07/07 01:06:36 valerie Exp $ */
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
import com.sri.biospice.warehouse.schema.linking.LinkingTables;
import com.sri.biospice.warehouse.schema.object.*;
import java.sql.ResultSet;
import java.sql.SQLException;

/**
 * @author Valerie.Wagner@sri.com
 *         Date: Aug 13, 2004
 */
public class LinkingTables_Test extends WarehouseFunctionalTest {

    public LinkingTables_Test(String testName) throws SQLException {
        super(testName);
    }


    public void testBioSourceWIDBioSubtypeWID() throws SQLException {
        DataSet dataset = new DataSet("test");
        dataset.store();
        LinkingTables linkingTables = new LinkingTables();

        com.sri.biospice.warehouse.schema.object.BioSource bioSource = new BioSource(dataset.getWID());
        bioSource.store();
        BioSubtype bioSubtype = new BioSubtype(dataset.getWID());
        bioSubtype.setValue("xxx");
        bioSubtype.store();

        linkingTables.addBioSourceWIDBioSubtypeWID(bioSource, bioSubtype);
        linkingTables.store();

        checkResult(bioSource, bioSubtype);
    }


    public void testProteinWIDFunctionWID() throws SQLException {
        DataSet dataset = new DataSet("test");
        dataset.store();
        LinkingTables linkingTables = new LinkingTables();

        Protein protein = new Protein(dataset.getWID());
        protein.store();
        Function function = new Function(dataset.getWID());
        function.store();

        linkingTables.addProteinWIDFunctionWID(protein, function);
        linkingTables.store();

        checkResult(protein, function);
    }


    public void testGeneWIDProteinWID() throws SQLException {
        DataSet dataset = new DataSet("test");
        dataset.store();
        LinkingTables linkingTables = new com.sri.biospice.warehouse.schema.linking.LinkingTables();

        Gene gene = new Gene(dataset.getWID());
        gene.store();
        Protein protein = new Protein(dataset.getWID());
        protein.store();

        linkingTables.addGeneWIDProteinWID(gene, protein);
        linkingTables.store();

        checkResult(gene, protein);
    }


    private void checkResult(ObjectTable table1, ObjectTable table2) throws SQLException {
        // We can recreate the names of the table and the columns because
        // by convention, all linking tables in the warehouse follow the form
        // table1WIDtable2WID.
        String wid1Name = table1.getTableName() + "WID";
        String wid2Name = table2.getTableName() + "WID";
        String tableName = wid1Name + wid2Name;
        String sqlLine = "Select * from " + tableName + " where " +
                wid1Name + "=" + table1.getWID() + " and " +
                wid2Name + "=" + table2.getWID();

        ResultSet results = WarehouseManager.getWarehouse().executeQuery(sqlLine);
        assertEquals("No entry found for " + tableName, true, results.next());
        assertEquals("Wrong WID for " + wid1Name, table1.getWID(), results.getLong(wid1Name));
        assertEquals("Wrong WID for " + wid2Name, table2.getWID(), results.getLong(wid2Name));
    }

}
