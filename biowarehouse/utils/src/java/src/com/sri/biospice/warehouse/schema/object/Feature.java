/* $Id: Feature.java,v 1.8 2009/06/15 22:08:41 valerie Exp $ */
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

import com.sri.biospice.warehouse.database.table.*;

/**
 * An entry in the Feature table
 *
 * @author Priyanka Gupta
 * @author Valerie Wagner
 */
public class Feature extends ObjectTableBase {
    // Columns in the Feature table
    private StringColumn descriptionCol;
    private StringColumn typeCol;
    private StringColumn classCol;
    private StringColumn sequenceTypeCol;
    private LongColumn sequenceWIDCol;
    private ClobColumn variantCol;
    private StringColumn regionOrPointCol;
    private StringColumn pointTypeCol;
    private IntegerColumn startPositionCol;
    private IntegerColumn endPositionCol;
    private StringColumn startPositionApproximateCol;
    private StringColumn endPositionApproximateCol;
    private StringColumn experimentalSupportCol;
    private StringColumn computationalSupportCol;

    private static TableMetaData tableMetaData = null;


    /**
     * Enumerated values for column <i>Feature.RegionOrPoint</i>
     */
    public enum RegionOrPoint {
        // todo: update API to use this enumeration

        /*
           <enumeration>
                <restriction value="region" description="Feature is specified by a start point and an end point on the sequence"/>
                <restriction value="point" description="Feature is specified by a single point on the sequence"/>
            </enumeration>
         */
        region,
        point
    }

    /**
     * Enumerated values for column <i>Feature.PointType</i>
     */
    public enum PointType {
        // todo: update API to use this enumeration

        /*
               <enumeration>
                   <restriction value="center" description="Feature is centered at location."/>
                   <restriction value="left" description="Feature extends to the left (decreasing position) of location."/>
                   <restriction value="right" description="Feature extends to the right (increasing position) of location."/>
               </enumeration>
        */
        center,
        left,
        right
    }

    /**
     * Enumerated values allowed for column <i>Feature.Class</i>.
     */
    public enum Class {
        /*
         <enumeration>
                <restriction value="binding site" description="Identifies the presence of a DNA binding site"/>
                <restriction value="promoter" description="Identifies the presence of a promoter"/>
                <restriction value="terminator" description="Identifies the presence of a terminator"/>
                <restriction value="pseudogene" description="Identifies a pseudogene, whether non-transcribed or processed"/>
                <restriction value="ORF" description="Identifies a truly unknown open reading frame according to warehouse definition (no strong evidence that a product is produced)"/>
                <restriction value="partial" description="Qualifier: States that feature is not complete"/>
                <restriction value="unknown product" description="Identifies an unspecified product is produced from this genomic location, as stated in dataset"/>
                <restriction value="notable" description="Qualifier: Characterizes the feature value as notable (same as 'Exceptional' in GB dataset)"/>
                <restriction value="by similarity" description="Qualifier: states that feature was derived by similarity analysis"/>
                <restriction value="potential" description="Qualifier: states that feature may be incorrect"/>
                <restriction value="probably" description="Qualifier: states that feature is probably correct"/>
                <restriction value="variant" description="Qualifier: states that the feature is a sequence variation. Feature.Variant will contain the alternate sequence."/>                
            </enumeration>
         */
        binding_site("binding site"),
        promoter("promoter"),
        terminator("terminator"),
        pseudogene("pseudogene"),
        ORF("ORF"),
        partial("partial"),
        unknown_product("unknown product"),
        notable("notable"),
        by_similarity("by similarity"),
        potential("potential"),
        probably("probably"),
        variant("variant");

        private String value;

        Class(String val) {
            this.value = val;
        }

        public String toString() {
            return value;
        }
    }

    public enum SequenceType {
//          <restriction value="P" description="Feature resides on a protein. Implies SequenceWID (if nonNULL) references a Protein"/>
//                <restriction value="S" description="Feature resides on a nucleic acid. Implies SequenceWID is nonNULL and references a Subsequence"/>
//                <restriction value="N" description="Feature resides on a nucleic acid. Implies SequenceWID (if nonNULL) references a NucleicAcid"/>
        P,
        S,
        N
    }


    // Enumerated values for the StartPositionApproximate
//    public static final String ENUM_START_POSITION_APPROXIMATE_GT = "gt";
//    public static final String ENUM_START_POSITION_APPROXIMATE_LT = "lt";
//    public static final String ENUM_START_POSITION_APPROXIMATE_NE = "ne";

    /**
     * Enumerated values allowed for columns <i>Feature.StartPositionApproximate</i> and
     * <i>Feature.EndPositionApproximate</i>
     */
    public enum PositionApproximate {
        gt,
        lt,
        ne
    }

    // Enumerated values for the EndPositionApproximate
//    public static final String ENUM_END_POSITION_APPROXIMATE_GT = "gt";
//    public static final String ENUM_END_POSITION_APPROXIMATE_LT = "lt";
//    public static final String ENUM_END_POSITION_APPROXIMATE_NE = "ne";


    /**
     * Creates a new entry in the Feature table
     *
     * @param datasetWid
     */
    public Feature(long datasetWid) {
        super(datasetWid, "Feature");
        init();
    }


    /**
     * Represents an already-existing entry in the Feature table
     *
     * @param datasetWid
     * @param featureWID the WID of the entry in the Feature table
     */
    public Feature(long datasetWid, long featureWID) {
        super(datasetWid, "Feature", featureWID);
        init();
    }


    private void init() {
        descriptionCol = new StringColumn("Description", this);
        typeCol = new StringColumn("Type", this);
        classCol = new StringColumn("Class", this);
        sequenceTypeCol = new StringColumn("SequenceType", this);
        sequenceWIDCol = new LongColumn("SequenceWID", this);
        variantCol = new ClobColumn("Variant", this);
        regionOrPointCol = new StringColumn("RegionOrPoint", this);
        pointTypeCol = new StringColumn("PointType", this);
        startPositionCol = new IntegerColumn("StartPosition", this);
        endPositionCol = new IntegerColumn("EndPosition", this);
        startPositionApproximateCol = new StringColumn("StartPositionApproximate", this);
        endPositionApproximateCol = new StringColumn("EndPositionApproximate", this);
        experimentalSupportCol = new StringColumn("ExperimentalSupport", this);
        computationalSupportCol = new StringColumn("ComputationalSupport", this);
    }


    public void setVariant(String variant) {
        variantCol.setValue(variant);
    }

    public void setSequenceType(String sequenceType) {
        sequenceTypeCol.setValue(sequenceType);
    }


    public void setDescription(String description) {
        descriptionCol.setValue(description);
    }

    public void setRegionOrPoint(String regionOrPoint) {
        regionOrPointCol.setValue(regionOrPoint);
    }

    public void setPointType(String pointType) {
        pointTypeCol.setValue(pointType);
    }

    public void setStartPosition(String start) {
        startPositionCol.setValue(start);
    }


    public void setEndPosition(String end) {
        endPositionCol.setValue(end);
    }


    public String getPositionEndApprox() {
        return endPositionApproximateCol.getValue();
    }


    public String getPositionStartApprox() {
        return startPositionApproximateCol.getValue();
    }


    public void setEndPositionApprox(PositionApproximate endApprox) {
        endPositionApproximateCol.setValue(endApprox);
    }


    public void setStartPositionApprox(PositionApproximate startApprox) {
        startPositionApproximateCol.setValue(startApprox);
    }


    public long getSubsequenceWID() {
        return sequenceWIDCol.getValue();
    }


    public void setSubsequenceWID(long subsequenceWID) {
        sequenceWIDCol.setValue(subsequenceWID);
    }


    public String getType() {
        return typeCol.getValue();
    }


    public void setType(String type) {
        typeCol.setValue(type);
    }


    public void setClass(Class classValue) {
        classCol.setValue(classValue);
    }


    public TableMetaData getTableMetaData() {
        return tableMetaData;
    }


    protected void setTableMetaData(TableMetaData tableMetaData) {
        this.tableMetaData = tableMetaData;
    }


    /**
     * This method clones data from this feature into the other feature
     *
     * @param anotherFeature Feature data is copied into
     */
    public void copyInto(Feature anotherFeature) {

        anotherFeature.descriptionCol.setValue(this.descriptionCol.getValue());
        anotherFeature.typeCol.setValue(this.typeCol.getValue());
        anotherFeature.classCol.setValue(this.classCol.getValue());
        anotherFeature.sequenceTypeCol.setValue(this.sequenceTypeCol.getValue());
        anotherFeature.sequenceWIDCol.setValue(this.sequenceWIDCol.getValue());
        anotherFeature.variantCol.setValue(this.variantCol.getValue());
        anotherFeature.regionOrPointCol.setValue(this.regionOrPointCol.getValue());
        anotherFeature.pointTypeCol.setValue(this.pointTypeCol.getValue());
        anotherFeature.startPositionCol.setValue(this.pointTypeCol.getValue());
        anotherFeature.endPositionCol.setValue(this.endPositionCol.getValue());
        anotherFeature.startPositionApproximateCol.setValue(this.startPositionApproximateCol.getValue());
        anotherFeature.endPositionApproximateCol.setValue(this.endPositionApproximateCol.getValue());
        anotherFeature.experimentalSupportCol.setValue(this.experimentalSupportCol.getValue());
        anotherFeature.computationalSupportCol.setValue(this.computationalSupportCol.getValue());
    }

    /**
     * For testing purposes
     *
     * @param columnName Name of the columne to get the allowed enumerated values for
     * @return A Vector of the allowed enumerated values
     */
    public Enum[] getAllowedValues(String columnName) {
        if (columnName.equalsIgnoreCase(classCol.getColumnName())) {
            return Class.values();
        } else if (columnName.equalsIgnoreCase(startPositionApproximateCol.getColumnName())) {
            return PositionApproximate.values();
        } else if (columnName.equalsIgnoreCase(endPositionApproximateCol.getColumnName())) {
            return PositionApproximate.values();
        } else if (columnName.equalsIgnoreCase(sequenceTypeCol.getColumnName())) {
            return SequenceType.values();
        } else if (columnName.equalsIgnoreCase(regionOrPointCol.getColumnName())) {
            return RegionOrPoint.values();
        } else if (columnName.equalsIgnoreCase(pointTypeCol.getColumnName())) {
            return PointType.values();
        }

        return null;
    }
}


