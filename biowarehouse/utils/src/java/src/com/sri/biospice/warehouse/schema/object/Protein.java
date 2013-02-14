/* $Id: Protein.java,v 1.2 2006/10/31 22:38:38 kejariwa Exp $ */
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
 * An entry in the Protein table
 * @author Priyanka Gupta
 */
public class Protein extends ObjectTableBase
{
    // Columns in the Protein table
    private StringColumn nameCol;
    private ClobColumn aaSequenceCol;
    private IntegerColumn lengthCol;
    private StringColumn lengthApproximateCol;
    private IntegerColumn chargeCol;
    private StringColumn fragmentCol;
    private FloatColumn molecularWeightCalcCol;
    private FloatColumn molecularWeightExpCol;
    private StringColumn piCalcCol;
    private StringColumn piExpCol;

    private static TableMetaData tableMetaData = null;


    private void init()
    {
        nameCol = new StringColumn( "Name", this );
        aaSequenceCol = new ClobColumn( "AASequence", this );
        lengthCol = new IntegerColumn( "Length", this );
        lengthApproximateCol = new StringColumn( "LengthApproximate", this );
        chargeCol = new IntegerColumn( "Charge", this );
        fragmentCol = new StringColumn( "Fragment", this );
        molecularWeightCalcCol = new FloatColumn( "MolecularWeightCalc", this );
        molecularWeightExpCol = new FloatColumn( "MolecularWeightExp", this );
        piCalcCol = new StringColumn( "PICalc", this );
        piExpCol = new StringColumn( "PIExp", this );
    }


    /**
     * Creates a new entry in the Protein table
     * @param datasetWid
     */
    public Protein( long datasetWid )
    {
        super( datasetWid, "Protein" );
        init();
    }


    /**
     * Represents an already-existing entry in the Protein table
     * @param datasetWid
     * @param proteinWID
     */
    public Protein( long datasetWid, long proteinWID )
    {
        super( datasetWid, "Protein", proteinWID );
        init();
    }


    public void setName( String name )
    {
        nameCol.setValue( name );
    }


    public Float getMolecularWeightCalc()
    {
        return molecularWeightCalcCol.getValue();
    }


    public void setMolecularWeightCalc( String molecularWeightCalc )
    {
        molecularWeightCalcCol.setValue( molecularWeightCalc );
    }


    public void setMolecularWeightCalc( float molecularWeightCalc )
    {
        molecularWeightCalcCol.setValue( molecularWeightCalc );
    }


    public TableMetaData getTableMetaData()
    {
        return tableMetaData;
    }


    public void setFragment( String fragment )
    {
        fragmentCol.setValue( fragment );
    }


    public void setLength( String length )
    {
        lengthCol.setValue( length );
    }

    public Integer getLength()
    {
        return lengthCol.getValue();
    }


    public void setAASequence( String aaSequence )
    {
        // todo: consider stripping off all whitespace
        aaSequenceCol.setValue( aaSequence );
    }


    public String getName()
    {
        return nameCol.getValue();
    }


    protected void setTableMetaData( TableMetaData tableMetaData )
    {
        this.tableMetaData = tableMetaData;
    }


    public void setLengthApproximate( String lengthApproximate )
    {
        lengthApproximateCol.setValue( lengthApproximate );
    }

    
    public String getAASequence() {
        return aaSequenceCol.getValue();
    }

    public String getLengthApproximate() {
        return lengthApproximateCol.getValue();
    }

    public Integer getCharge() {
        return chargeCol.getValue();
    }

    public String getFragment() {
        return fragmentCol.getValue();
    }

    public Float getMolecularWeightExp() {
        return molecularWeightExpCol.getValue();
    }

    public String getPiCalc() {
        return piCalcCol.getValue();
    }

    public String getPiExp() {
        return piExpCol.getValue();
    }


    public void setPiCalc(String piCalc) {
        piCalcCol.setValue(piCalc);
    }

    public void setPiExp(String piExp) {
        piExpCol.setValue(piExp);
    }
}
