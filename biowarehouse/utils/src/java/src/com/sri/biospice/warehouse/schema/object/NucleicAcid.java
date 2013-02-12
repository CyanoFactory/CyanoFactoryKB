/* $Id: NucleicAcid.java,v 1.3 2008/10/06 21:18:57 valerie Exp $ */
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
 * An entry in the NucleicAcid table
 *
 * @author Priyanka Gupta
 * @author Valerie Wagner
 */
public class NucleicAcid extends ObjectTableBase {
    // Columns of the NucleicAcid table
    private StringColumn nameCol;
    private StringColumn typeCol;
    private StringColumn classCol;
    private StringColumn topologyCol;
    private StringColumn strandednessCol;
    private StringColumn sequenceDerivationCol;
    private StringColumn fragmentCol;
    private StringColumn fullySequencedCol;
    private IntegerColumn moleculeLengthCol;
    private StringColumn moleculeLengthApproximateCol;
    private IntegerColumn cumulativeLengthCol;
    private StringColumn cumulativeLengthApproximateCol;
    private LongColumn geneticCodeWIDCol;
    private LongColumn bioSourceWIDCol;


    // Metadata for the NucleicAcid table
    private static TableMetaData tableMetaData = null;

    // Enumerated values for the Class column
//    public static final String ENUM_Class_preRNA = "pre-RNA";
//    public static final String ENUM_Class_mRNA = "mRNA";
//    public static final String ENUM_Class_rRNA = "rRNA";
//    public static final String ENUM_Class_tRNA = "tRNA";
//    public static final String ENUM_Class_snRNA = "snRNA";
//    public static final String ENUM_Class_scRNA = "scRNA";
//    public static final String ENUM_Class_snoRNA = "snoRNA";
//    public static final String ENUM_Class_other = "other";
//    public static final String ENUM_Class_chromosome = "chromosome";
//    public static final String ENUM_Class_plasmid = "plasmid";
//    public static final String ENUM_Class_organelleChromosome = "organelle-chromosome";
//    public static final String ENUM_Class_transposon = "transposon";
//    public static final String ENUM_Class_virus = "virus";

    /**
     * Enumerated values allowed for column <i>NucleicAcid.Class</i>.
     */
    public enum Class {
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
        snoRNA,
        other,
        chromosome,
        plasmid,
        organelle_chromosome {
            public String toString() {
                return "organelle-chromosome";
            }
        },
        transposon,
        virus,
        unknown
    }

//    public static final String ENUM_Type_DNA = "DNA";
//    public static final String ENUM_Type_RNA = "RNA";
//    public static final String ENUM_Type_NA = "NA";

    /**
     * Enumerated values allowed for column <i>NucleicAcid.Type</i>.
     */
    public enum Type {
        DNA,
        RNA,
        NA
    }

    /**
     * Enumerated values allowed for column <i>NucleicAcid.MoleculeLengthApproximate</i>
     */
    public enum LengthApproximate {
        // todo: update API to use this enumeration
        /*
         <enumeration>
                <restriction value="gt" description="The length of the molecule's sequence is greater than the actual length specified"/>
                <restriction value="lt" description="The length of the molecule's sequence is less than the actual length specified"/>
                <restriction value="ne" description="The length of the molecule's sequence is less than or greater than the actual length. All we know is that its not the exact length"/>
            </enumeration>
         */
        gt,
        lt,
        ne
    }

//    public static final String ENUM_Topology_linear = "linear";
//    public static final String ENUM_Topology_circular = "circular";
//    public static final String ENUM_Topology_other = "other";

    /**
     * Enumerated values allowed for column <i>NucleicAcid.Topology</i>.
     */
    public enum Topology {
        linear,
        circular,
        other
    }

//    public static final String ENUM_LengthApproximate_gt = "gt";
//    public static final String ENUM_LengthApproximate_lt = "lt";
//    public static final String ENUM_LengthApproximate_ne = "ne";
//
//    public static final String ENUM_TotalLengthApproximate_gt = "gt";
//    public static final String ENUM_TotalLengthApproximate_lt = "lt";
//    public static final String ENUM_TotalLengthApproximate_ne = "ne";

//    public static final String ENUM_SequenceDerivation_virtual = "virtual";
//    public static final String ENUM_SequenceDerivation_raw = "raw";
//    public static final String ENUM_SequenceDerivation_seg = "seg";
//    public static final String ENUM_SequenceDerivation_reference = "reference";
//    public static final String ENUM_SequenceDerivation_constructed = "constructed";
//    public static final String ENUM_SequenceDerivation_consensus = "consensus";
//    public static final String ENUM_SequenceDerivation_map = "map";

    /**
     * Enumerated values allowed for column <i>NucleicAcid.SequenceDerivation</i>.
     */
    public enum SequenceDerivation {
        virtual,
        raw,
        seg,
        reference,
        constructed,
        consensus,
        map
    }

//    public static final String ENUM_Strandedness_ss = "ss";
//    public static final String ENUM_Strandedness_ds = "ds";
//    public static final String ENUM_Strandedness_mixed = "mixed";

    /**
     * Enumerated values allowed for column <i>NucleicAcid.Strandedness</i>.
     */
    public enum Strandedness {
        ss,
        ds,
        mixed
    }


    /**
     * Creates a new entry in the NucleicAcid table
     *
     * @param datasetWid
     */
    public NucleicAcid(long datasetWid) {
        super(datasetWid, "NucleicAcid");
        init();
    }


    /**
     * Used for loading an existing entry in the NucleicAcid table
     *
     * @param datasetWID
     * @param thisWID
     */
    public NucleicAcid(long datasetWID, long thisWID) {
        super(datasetWID, "NucleicAcid", thisWID);
        init();
    }


    private void init() {
        nameCol = new StringColumn("Name", this);
        classCol = new StringColumn("Class", this);
        geneticCodeWIDCol = new LongColumn("GeneticCodeWID", this);
        moleculeLengthCol = new IntegerColumn("MoleculeLength", this);
        moleculeLengthApproximateCol = new StringColumn("MoleculeLengthApproximate", this);
        sequenceDerivationCol = new StringColumn("SequenceDerivation", this);
        strandednessCol = new StringColumn("Strandedness", this);
        topologyCol = new StringColumn("Topology", this);
        typeCol = new StringColumn("Type", this);
        fragmentCol = new StringColumn("Fragment", this);
        fullySequencedCol = new StringColumn("FullySequenced", this);
        cumulativeLengthCol = new IntegerColumn("CumulativeLength", this);
        cumulativeLengthApproximateCol = new StringColumn("CumulativeLengthApproximate", this);
        bioSourceWIDCol = new LongColumn("BioSourceWID", this);
    }

    public void setName(String name) {
        nameCol.setValue(name);
    }
    
    public String getName() {
        return nameCol.getValue();
    }

    public void setClass(Class classValue) {
        classCol.setValue(classValue);
    }


    public void setGeneticCodeWID(String geneticCodeWID) {
        geneticCodeWIDCol.setValue(geneticCodeWID);
    }


    public void setSequenceDerivation(SequenceDerivation sequenceDerivation) {
        sequenceDerivationCol.setValue(sequenceDerivation);
    }


    public void setType(Type type) {
        typeCol.setValue(type);
    }


    public void setMoleculeLength(String moleculeLength) {
        moleculeLengthCol.setValue(moleculeLength);
    }


    public void setMoleculeLengthApproximate(String totalLengthApprox) {
        moleculeLengthApproximateCol.setValue(totalLengthApprox);
    }


    public String getMoleculeLengthApproximate(String lengthApprox) {
        return moleculeLengthApproximateCol.getValue();
    }


    public void setTopology(Topology topology) {
        topologyCol.setValue(topology);
    }


    public void setStrandedness(Strandedness strandedness) {
        strandednessCol.setValue(strandedness);
    }


    protected void setTableMetaData(TableMetaData tableMetaData) {
        this.tableMetaData = tableMetaData;
    }


    public TableMetaData getTableMetaData() {
        return tableMetaData;
    }


    public String getType() {
        return typeCol.getValue();
    }


    public String getStrandedness() {
        return strandednessCol.getValue();
    }


    public Enum[] getAllowedValues(String columnName) {

        if (columnName.equalsIgnoreCase(classCol.getColumnName())) {
            return Class.values();
        } else if (columnName.equalsIgnoreCase(typeCol.getColumnName())) {
            return Type.values();
        } else if (columnName.equalsIgnoreCase(topologyCol.getColumnName())) {
            return Topology.values();
        } else if (columnName.equalsIgnoreCase(sequenceDerivationCol.getColumnName())) {
            return SequenceDerivation.values();
        } else if (columnName.equalsIgnoreCase(strandednessCol.getColumnName())) {
            return Strandedness.values();
        } else if (columnName.equalsIgnoreCase(moleculeLengthApproximateCol.getColumnName())) {
            return LengthApproximate.values();
        } else if (columnName.equalsIgnoreCase(cumulativeLengthApproximateCol.getColumnName())) {
            return LengthApproximate.values();
        }

        return null;
    }


    public Long getGeneticCodeWID() {
        return geneticCodeWIDCol.getLongValue();
    }


    public void setBioSourceWID(long wid) {
        bioSourceWIDCol.setValue(wid);
    }

    public Long getBioSourceWID(long wid) {
        return bioSourceWIDCol.getLongValue();
    }
    
    public void setCumulativeLength(String length) {
        cumulativeLengthCol.setValue(length);
    }


    public void setCumulativeLengthApproximate(String qualifier) {
        cumulativeLengthApproximateCol.setValue(qualifier);
    }

}

