/* $Id: BioSource.java,v 1.3 2006/11/30 07:08:36 kejariwa Exp $ */
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

import com.sri.biospice.warehouse.database.table.LongColumn;
import com.sri.biospice.warehouse.database.table.StringColumn;
import com.sri.biospice.warehouse.database.table.TableMetaData;
import java.util.Vector;

/**
 * An entry in the BioSource table
 *
 * @author Valerie Wagner
 *         Date: May 12, 2004
 */
public class BioSource extends ObjectTableBase {

    // Columns in the BioSource table
    private LongColumn taxonWIDCol;
    private StringColumn nameCol;
    private StringColumn strainCol;
    private StringColumn organCol;
    private StringColumn organelleCol;
    private StringColumn tissueCol;
    private StringColumn cellTypeCol;
    private StringColumn cellLineCol;
    private StringColumn atccidCol;
    private StringColumn diseasedColumn;
    private StringColumn diseaseColumn;
    private StringColumn developmentStageCol;
    private StringColumn sexCol;
    private StringColumn mageClassCol;

    // Table metadata
    private static TableMetaData tableMetaData = null;

    // Enumerated values for the Sex column
//    public static final String ENUM_SEX_MALE = "Male";
//    public static final String ENUM_SEX_FEMALE = "Female";
//    public static final String ENUM_SEX_HERMAPHRODITE = "Hermaphrodite";
//    public static final String ENUM_SEX_DOES_NOT_APPLY = "DoesNotApply";

    /**
     * Enumerated values allowed for column <i>BioSource.Sex</i>.
     */
    public enum Sex {
        Male,
        Female,
        Hermaphrodite,
        DoesNotApply
    }


    /**
     * Defines a new entry in the BioSource table
     *
     * @param datasetWid
     */
    public BioSource(long datasetWid) {
        super(datasetWid, "BioSource");
        init();
    }


    /**
     * Represents an already-existing entry in the BioSource table
     *
     * @param datasetWid
     * @param bioSourceWID
     */
    public BioSource(long datasetWid, long bioSourceWID) {
        super(datasetWid, "BioSource", bioSourceWID);
        init();
    }


    private void init() {
        taxonWIDCol = new LongColumn("TaxonWID", this);
        nameCol = new StringColumn("Name", this);
        strainCol = new StringColumn("Strain", this);
        organCol = new StringColumn("Organ", this);
        organelleCol = new StringColumn("Organelle", this);
        tissueCol = new StringColumn("Tissue", this);
        cellTypeCol = new StringColumn("CellType", this);
        cellLineCol = new StringColumn("CellLine", this);
        atccidCol = new StringColumn("ATCCId",this);
        diseasedColumn = new StringColumn("Diseased",this);
        diseaseColumn = new StringColumn("Disease",this);
        developmentStageCol = new StringColumn("DevelopmentStage", this);
        sexCol = new StringColumn("Sex", this);
        mageClassCol = new StringColumn("MAGEClass",this);
    }


    public void setCellLine(String cellLine) {
        cellLineCol.setValue(cellLine);
    }


    public String getCellLine() {
        return cellLineCol.getValue();
    }


    public void setCellType(String cellType) {
        cellTypeCol.setValue(cellType);
    }


    public void setTissue(String tissue) {
        tissueCol.setValue(tissue);
    }


    public void setDevelopmentStage(String developmentStage) {
        developmentStageCol.setValue(developmentStage);
    }


    public void setStrain(String strain) {
        strainCol.setValue(strain);
    }


    public void setTaxonWID(String taxonWID) {
        taxonWIDCol.setValue(taxonWID);
    }


    public void setTaxonWID(Long taxonWID) {
        taxonWIDCol.setValue(taxonWID);
    }


    public TableMetaData getTableMetaData() {
        return tableMetaData;
    }


    public void setName(String name) {
        nameCol.setValue(name);
    }


    public void setOrgan(String organ) {
        organCol.setValue(organ);
    }


    public void setOrganelle(String organelle) {
        organelleCol.setValue(organelle);
    }


    public void setSex(Sex sex) {
        sexCol.setValue(sex);
    }


    protected void setTableMetaData(TableMetaData tableMetaData) {
        this.tableMetaData = tableMetaData;
    }


    public Enum[] getAllowedValues(String columnName) {
        if (columnName.equalsIgnoreCase(sexCol.getColumnName())) {
            return Sex.values();
        }
        return null;
    }


    public Long getTaxonWID() {
        return taxonWIDCol.getLongValue();
    }


    public String getName() {
        return nameCol.getValue();
    }


    public String getStrain() {
        return strainCol.getValue();
    }


    public String getOrgan() {
        return organCol.getValue();
    }


    public String getOrganelle() {
        return organelleCol.getValue();
    }


    public String getTissue() {
        return tissueCol.getValue();
    }


    public String getCellType() {
        return cellTypeCol.getValue();
    }


    public String getDevelopmentStage() {
        return developmentStageCol.getValue();
    }


    public String getSex() {
        return sexCol.getValue();
    }


    public String getMageClass() {
        return mageClassCol.getValue();
    }

    public void setMageClass(String mageClass) {
        this.mageClassCol.setValue(mageClass);
    }

   public String getAtccid() {
        return atccidCol.getValue();
    }

    public void setAtccid(String atccid) {
        this.atccidCol.setValue(atccid);
    }

    public String getDiseased() {
        return diseasedColumn.getValue();
    }

    public void setDiseased(String diseased) {
        this.diseasedColumn.setValue(diseased);
    }

    public String getDisease() {
        return diseaseColumn.getValue();
    }

    public void setDisease(String disease) {
        this.diseaseColumn.setValue(disease);
    }


}
