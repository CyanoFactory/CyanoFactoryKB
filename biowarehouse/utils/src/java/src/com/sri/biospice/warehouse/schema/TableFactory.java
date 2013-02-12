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

import com.sri.biospice.warehouse.database.WarehouseManager;
import com.sri.biospice.warehouse.database.table.Column;
import com.sri.biospice.warehouse.schema.linking.*;
import com.sri.biospice.warehouse.schema.object.*;
import com.sri.biospice.warehouse.util.LoaderMain;
import com.sri.bw.dbschema.SchemaDocument;
import org.apache.log4j.Logger;
import org.apache.xmlbeans.XmlException;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.Vector;
import java.io.File;
import java.io.IOException;

/**
 * @author Valerie Wagner
 *         Date: Apr 8, 2005
 */
public class TableFactory {

    // DataSet
    public static final int DATASET = 1000;

    // Object tables
    public static final int BIOSOURCE = 0;
    public static final int BIOSUBTYPE = 1;
    public static final int CITATION = 2;
    public static final int ENZYMATIC_REACTION = 3;
    public static final int FEATURE = 4;
    public static final int FUNCTION = 5;
    public static final int GENE = 6;
    public static final int NUCLEIC_ACID = 7;
    public static final int PROTEIN = 8;
    public static final int REACTION = 9;
    public static final int SUBSEQUENCE = 10;
    public static final int TERM = 11;

    // Object description tables
    public static final int COMMENT_TABLE = 100;
    public static final int CROSS_REFERENCE = 101;
    public static final int DBID = 102;
    public static final int DESCRIPTION = 103;
    public static final int ENTRY = 104;
    public static final int SYNONYM_TABLE = 105;
    public static final int TERM_RELATIONSHIP = 106;

    // Linking tables
    public static final int BIOSOURCE_WID_BIOSUBTYPE_WID = 200;
    public static final int BIOSOURCE_WID_GENE_WID = 201;
    public static final int BIOSOURCE_WID_PROTEIN_WID = 202;
    public static final int CITATION_WID_OTHER_WID = 203;
    public static final int GENE_WID_NUCLEICACID_WID = 204;
    public static final int GENE_WID_PROTEIN_WID = 205;
    public static final int PROTEIN_WID_FUNCTION_WID = 206;

    private static final int[] tableTypes = new int[]
    {
        // DataSet
        DATASET,

        // Object tables
        BIOSOURCE,
        BIOSUBTYPE,
        CITATION,
        ENZYMATIC_REACTION,
        FEATURE,
        FUNCTION,
        GENE,
        NUCLEIC_ACID,
        PROTEIN,
        REACTION,
        SUBSEQUENCE,
        TERM,

        // Object description tables
        COMMENT_TABLE,
        CROSS_REFERENCE,
        DBID,
        DESCRIPTION,
        ENTRY,
        SYNONYM_TABLE,
        TERM_RELATIONSHIP,

        // Linking tables
        BIOSOURCE_WID_BIOSUBTYPE_WID,
        BIOSOURCE_WID_GENE_WID,
        BIOSOURCE_WID_PROTEIN_WID,
        CITATION_WID_OTHER_WID,
        GENE_WID_NUCLEICACID_WID,
        GENE_WID_PROTEIN_WID,
        PROTEIN_WID_FUNCTION_WID,
    };


    public static ObjectTable getObjectTable(long datasetWID, long wid, int tableType) {
        ObjectTable table = null;

        switch (tableType) {
            case BIOSOURCE:
                table = new BioSource(datasetWID, wid);
                break;
            case BIOSUBTYPE:
                table = new BioSubtype(datasetWID, wid);
                break;
            case CITATION:
                table = new Citation(datasetWID, wid);
                break;
            case ENZYMATIC_REACTION:
                table = new EnzymaticReaction(datasetWID, wid);
                break;
            case FEATURE:
                table = new Feature(datasetWID, wid);
                break;
            case FUNCTION:
                table = new Function(datasetWID, wid);
                break;
            case GENE:
                table = new Gene(datasetWID, wid);
                break;
            case NUCLEIC_ACID:
                table = new NucleicAcid(datasetWID, wid);
                break;
            case PROTEIN:
                table = new Protein(datasetWID, wid);
                break;
            case REACTION:
                table = new Reaction(datasetWID, wid);
                break;
            case SUBSEQUENCE:
                table = new SubSequence(datasetWID, wid);
                break;
            case TERM:
                table = new Term(datasetWID, wid);
                break;
        }

        return table;
    }

    public static Vector loadTables(String query, int tableType) {
        Vector refs = new Vector();
        ResultSet results = null;
        try {
            results = WarehouseManager.getWarehouse().executeQuery(query);
            if (results != null) {
                while (results.next()) {
                    Table oneTable = createTable(tableType);
                    Vector columns = oneTable.getColumns();
                    for (int i = 0; i < columns.size(); i++) {
                        Column column = (Column)columns.elementAt(i);
                        column.load(results);
                    }
                    refs.add(oneTable);
                }
                results.close();
            }
        } catch (SQLException e) {
            Logger.getLogger(TableFactory.class).error("Unable to load tables using query: " + query, e);
            if (results != null) {
                try {
                    results.close();
                } catch (SQLException e1) {
                }
            }
        }

        return refs;
    }

    private static Table createTable(int tableType) {
        Table table = null;

        table = getObjectTable(0, 0, tableType);
        if (table != null) {
            return table;
        }

        switch (tableType) {
            case DATASET:
                table = new DataSet("test");
                break;
            case CROSS_REFERENCE:
                table = new CrossReference();
                break;
            case COMMENT_TABLE:
                table = new CommentTable();
                break;
            case DBID:
                table = new DBID();
                break;
            case SYNONYM_TABLE:
                table = new SynonymTable();
                break;
            case TERM_RELATIONSHIP:
                table = new TermRelationship();
                break;
            case DESCRIPTION:
                table = new Description();
                break;
            case ENTRY:
                table = new Entry();
                break;
            default:
                table = createTestLinkingTable(tableType);
                break;

        }

        return table;
    }

    private static Table createTestLinkingTable(int tableType) {
        Table table = null;
        switch (tableType) {
            case BIOSOURCE_WID_BIOSUBTYPE_WID:
                table = new BioSourceWIDBioSubtypeWID((BioSource)createTestTable(BIOSOURCE),
                        (BioSubtype)createTestTable(BIOSUBTYPE));
                break;
            case BIOSOURCE_WID_GENE_WID:
                table = new BioSourceWIDGeneWID((BioSource)createTestTable(BIOSOURCE),
                        (Gene)createTestTable(GENE));
                break;
            case BIOSOURCE_WID_PROTEIN_WID:
                table = new BioSourceWIDProteinWID((BioSource)createTestTable(BIOSOURCE),
                        (Protein)createTestTable(PROTEIN));
                break;
            case CITATION_WID_OTHER_WID:
                table = new CitationWIDOtherWID((Citation)createTestTable(CITATION),
                        (Gene)createTestTable(GENE));
                break;
            case GENE_WID_NUCLEICACID_WID:
                table = new GeneWIDNucleicAcidWID((Gene)createTestTable(GENE),
                        (NucleicAcid)createTestTable(NUCLEIC_ACID));
                break;
            case GENE_WID_PROTEIN_WID:
                table = new GeneWIDProteinWID((Gene)createTestTable(GENE),
                        (Protein)createTestTable(PROTEIN));
                break;
            case PROTEIN_WID_FUNCTION_WID:
                table = new ProteinWIDFunctionWID((Protein)createTestTable(PROTEIN),
                        (Function)createTestTable(FUNCTION));
        }
        return table;
    }

    public static Table createTestTable(int tableType) {
        return createTable(tableType);
    }

    public static int[] getTabletypes() {
        return tableTypes;
    }

  

}
