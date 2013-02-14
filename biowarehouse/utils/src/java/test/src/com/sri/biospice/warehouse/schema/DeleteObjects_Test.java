/* $Id: DeleteObjects_Test.java,v 1.2 2006/07/07 01:06:36 valerie Exp $ */
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
import com.sri.biospice.warehouse.schema.object.BioSource;
import com.sri.biospice.warehouse.schema.object.NucleicAcid;
import com.sri.biospice.warehouse.schema.object.Protein;
import com.sri.biospice.warehouse.schema.object.SubSequence;
import java.sql.ResultSet;
import java.sql.SQLException;

/**
 * @author Valerie.Wagner@sri.com
 *         Date: Sep 13, 2004
 */
public class DeleteObjects_Test extends WarehouseFunctionalTest {
    public DeleteObjects_Test(String testName) throws SQLException {
        super(testName);
    }


    public void testDeleteProtein() throws SQLException {
        DataSet dataset = new DataSet("test");
        dataset.store();

        int proteinTableCountBefore = getCount("Protein");
        int crossRefTableCountBefore = getCount("CrossReference");
        int linkTableCountBefore = getCount("BioSourceWIDProteinWID");
        Protein protein1 = new Protein(dataset.getWID());
        protein1.addCrossReference("test", "test", "test", "test");
        protein1.store();

        Protein protein2 = new Protein(dataset.getWID());
        protein2.store();

        BioSource bioSource = new BioSource(dataset.getWID());
        bioSource.store();

        LinkingTables linkingTables = new LinkingTables();
        linkingTables.addBioSourceWIDProteinWID(bioSource, protein1);
        linkingTables.store();

        int proteinCountAfterInsert = getCount("Protein");
        int crossRefCountAfterInsert = getCount("CrossReference");
        int linkCountAfterInsert = getCount("BioSourceWIDProteinWID");
        assertEquals("Protein count after insert should be one more than before", proteinTableCountBefore + 2, proteinCountAfterInsert);
        assertEquals("CrossReference count after insert should be one more than before", crossRefTableCountBefore + 1, crossRefCountAfterInsert);
        assertEquals("BioSourceWIDProteinWID count after insert should be one more than before", linkTableCountBefore + 1, linkCountAfterInsert);

        Protein proteinToDelete = new Protein(dataset.getWID(), protein1.getWID());
        proteinToDelete.delete();

        protein2.delete();

        int proteinCountAfterDelete = getCount("Protein");
        int crossRefCountAfterDelete = getCount("CrossReference");
        int linkCountAfterDelete = getCount("BioSourceWIDProteinWID");
        assertEquals("Protein count after delete should be equal to original count", proteinTableCountBefore, proteinCountAfterDelete);
        assertEquals("CrossReference count after delete should be equal to original count", crossRefTableCountBefore, crossRefCountAfterDelete);
        assertEquals("BioSourceWIDProteinWID count after delete should be equal to original count", linkTableCountBefore, linkCountAfterDelete);
    }


    public void testDeleteSubsequence() throws SQLException {
        DataSet dataset = new DataSet("test");
        dataset.store();

        int subsequenceTableCountBefore = getCount("Subsequence");
        SubSequence subsequence = new SubSequence(dataset.getWID());
        NucleicAcid nucleicAcid = new NucleicAcid(dataset.getWID());
        nucleicAcid.setType(NucleicAcid.Type.DNA);
        nucleicAcid.setClass(NucleicAcid.Class.plasmid);
        subsequence.setNucleicAcidWid(nucleicAcid.getWID());
        subsequence.setSequenceData("xxx");
        nucleicAcid.store();
        subsequence.store();

        int subsequenceCountAfterInsert = getCount("Subsequence");
        assertEquals("Subsequence count after insert should be one more than before", subsequenceTableCountBefore + 1, subsequenceCountAfterInsert);

        SubSequence subsequenceToDelete = new SubSequence(dataset.getWID(), subsequence.getWID());
        subsequenceToDelete.delete();

        int subsequenceCountAfterDelete = getCount("Subsequence");
        assertEquals("Subsequence count after delete should be equal to original count", subsequenceTableCountBefore, subsequenceCountAfterDelete);
    }


    private int getCount(String tableName) throws SQLException {
        String selectSqlCommand = "Select count(*) from " + tableName;
        ResultSet results = WarehouseManager.getWarehouse().executeQuery(selectSqlCommand);
        assertEquals("Result set from count(*) command contained no results!", true, results.next());
        return results.getInt(1);
    }
}
