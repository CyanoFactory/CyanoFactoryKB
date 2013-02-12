/* $Id: Gene.java,v 1.1 2006/07/07 15:03:37 valerie Exp $ */
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

import com.sri.biospice.warehouse.database.table.IntegerColumn;
import com.sri.biospice.warehouse.database.table.LongColumn;
import com.sri.biospice.warehouse.database.table.StringColumn;
import com.sri.biospice.warehouse.database.table.TableMetaData;

/**
 * An entry in the Gene table
 *
 * @author Priyanka Gupta
 * @author Valerie Wagner
 */
public class Gene extends ObjectTableBase {
    // Columns of the table
    private StringColumn nameCol;
    private LongColumn nucleicAcidWIDCol;
    private LongColumn subsequenceWIDCol;
    private StringColumn typeCol;
    private StringColumn genomeIDCol;
    private IntegerColumn codingRegionStartCol;
    private IntegerColumn codingRegionEndCol;
    private StringColumn codingRegionStartApproximateCol;
    private StringColumn codingRegionEndApproximateCol;
    private StringColumn directionCol;
    private StringColumn interruptedCol;

    private static TableMetaData tableMetaData = null;


    /**
     * Enumerated values allowed for column <i>Gene.Type</i>.
     */
    public enum Type {
        preRNA {
            public String toString() {
                return "pre-RNA";
            }
        },
        mRNA,
        rRNA,
        tRNA,
        snRNA,
        scRNA,
        polypeptide,
        snoRNA,
        other,
        unknown;
    }

    /**
     * Enumerated values allowed for column <i>Gene.Direction</i>.
     */
    public enum Direction {
        unknown,
        forward,
        reverse,
        forward_and_reverse,
        undefined_value;
    }
//    public static final String ENUM_DIRECTION_UNKOWN = "unknown";
//    public static final String ENUM_DIRECTION_FORWARD = "forward";
//    public static final String ENUM_DIRECTION_REVERSE = "reverse";
//    public static final String ENUM_DIRECTION_FORWARD_AND_REVERSE = "forward_and_reverse";
//    public static final String ENUM_DIRECTION_UNDEFINED_VALUE = "undefined_value";


    /**
     * Creates a new entry in the Gene table
     *
     * @param datasetWid
     */
    public Gene(long datasetWid) {
        super(datasetWid, "Gene");
        init();
    }


    /**
     * For an entry loaded from the Gene table in the database
     *
     * @param datasetWID
     * @param thisWID    The WID of the entry in the Gene table
     */
    public Gene(long datasetWID, long thisWID) {
        super(datasetWID, "Gene", thisWID);
        init();
    }


    private void init() {
        nameCol = new StringColumn("Name", this);
        nucleicAcidWIDCol = new LongColumn("NucleicAcidWID", this);
        subsequenceWIDCol = new LongColumn("SubsequenceWID", this);
        typeCol = new StringColumn("Type", this);
        genomeIDCol = new StringColumn("GenomeID", this);
        codingRegionStartCol = new IntegerColumn("CodingRegionStart", this);
        codingRegionEndCol = new IntegerColumn("CodingRegionEnd", this);
        codingRegionStartApproximateCol = new StringColumn("CodingRegionStartApproximate", this);
        codingRegionEndApproximateCol = new StringColumn("CodingRegionEndApproximate", this);
        directionCol = new StringColumn("Direction", this);
        interruptedCol = new StringColumn("Interrupted", this);
    }


    public void setNucleicAcidWid(long nucleicAcidWid) {
        nucleicAcidWIDCol.setValue(nucleicAcidWid);
    }


    public void setSubsequenceWid(long subsequenceWid) {
        subsequenceWIDCol.setValue(subsequenceWid);
    }


    public String getName() {
        return nameCol.getValue();
    }


    public String getType() {
        return typeCol.getValue();
    }


    public Integer getCodingRegionStart() {
        return codingRegionStartCol.getValue();
    }


    public String getCodingRegionEndApprox() {
        return codingRegionEndApproximateCol.getValue();
    }


    public String getCodingRegionStartApprox() {
        return codingRegionStartApproximateCol.getValue();
    }


    public Integer getCodingRegionEnd() {
        return codingRegionEndCol.getValue();
    }


    public String getDirection() {
        return directionCol.getValue();
    }


    public String getGenomeId() {
        return genomeIDCol.getValue();
    }


    public void setName(String name) {
        nameCol.setValue(name);
    }


//    public void setType(String type) {
//        typeCol.setValue(type);
//    }

    public void setType(Type type) {
        typeCol.setValue(type);
    }


    public void setCodingRegionStart(String start) {
        codingRegionStartCol.setValue(start);
    }


    public void setCodingRegionEndApprox(String endApprox) {
        codingRegionEndApproximateCol.setValue(endApprox);
    }


    public void setCodingRegionStartApprox(String startApprox) {
        codingRegionStartApproximateCol.setValue(startApprox);
    }


    public void setCodingRegionEnd(String end) {
        codingRegionEndCol.setValue(end);
    }


    public void setDirection(Direction direction) {
        directionCol.setValue(direction);
    }


    public void setGenomeId(String genomeId) {
        genomeIDCol.setValue(genomeId);
    }


    public TableMetaData getTableMetaData() {
        return tableMetaData;
    }


    public void setInterrupted(String interrupted) {
        interruptedCol.setValue(interrupted);
    }


    protected void setTableMetaData(TableMetaData tableMetaData) {
        this.tableMetaData = tableMetaData;
    }


    public Enum[] getAllowedValues(String columnName) {

        if (columnName.equalsIgnoreCase(typeCol.getColumnName())) {
            return Type.values();
        } else if (columnName.equalsIgnoreCase(directionCol.getColumnName())) {
            return Direction.values();
        }

        return null;
    }
}
