/* $Id: SubSequence.java,v 1.1 2006/07/07 15:03:38 valerie Exp $ */
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
 * An entry in the SubSequence table
 * @author Priyanka Gupta
 * @author Valerie Wagner
 */
public class SubSequence extends ObjectTableBase
{
    // Columns in the SubSequence table
    private LongColumn nucleicAcidWIDCol;
    private StringColumn fullSequenceCol;
    private ClobColumn sequenceCol;
    private IntegerColumn lengthCol;
    private StringColumn lengthApproximateCol;
    private FloatColumn percentGCCol;
    private StringColumn versionCol;

    // Metadata for the Subsequence table
    private static TableMetaData tableMetaData = null;


    /**
     * Creates a new entry in the SubSequence table
     * @param datasetWid
     */
    public SubSequence( long datasetWid )
    {
        super( datasetWid, "Subsequence" );
        init();
    }


    /**
     * Represents a Subsequence already in the warehouse
     * @param datasetWid
     * @param subsequenceWID
     */
    public SubSequence( long datasetWid, long subsequenceWID )
    {
        super( datasetWid, "Subsequence", subsequenceWID );
        init();
    }


    private void init()
    {
        nucleicAcidWIDCol = new LongColumn( "NucleicAcidWID", this );
        fullSequenceCol = new StringColumn( "FullSequence", this );
        sequenceCol = new ClobColumn( "Sequence", this );
        lengthCol = new IntegerColumn( "Length", this );
        lengthApproximateCol = new StringColumn( "LengthApproximate", this );
        percentGCCol = new FloatColumn( "PercentGC", this );
        versionCol = new StringColumn( "Version", this );
    }


    /**
     * Sets the sequence, length, and percentGC for this sequence
     * This method is named setSequenceData because it does more than
     * just set the value for the Subsequence.Sequence column.
     * @param sequence
     */
    public void setSequenceData( String sequence )
    {
        sequenceCol.setValue( sequence );
        if( sequence != null )
        {
            lengthCol.setValue( sequence.length() );
            percentGCCol.setValue( calculatePercentGC( sequence ) );
        }
    }


    public void setNucleicAcidWid( long wid )
    {
        nucleicAcidWIDCol.setValue( wid );
    }


    public void setNucleicAcidWid( String wid )
    {
        nucleicAcidWIDCol.setValue( wid );
    }


    public String getSequence()
    {
        return sequenceCol.getValue();
    }


    public TableMetaData getTableMetaData()
    {
        return tableMetaData;
    }


    public void setFullSequence( boolean value )
    {
        if( value == true )
        {
            fullSequenceCol.setValue( "T" );
        }
        else
        {
            fullSequenceCol.setValue( "F" );
        }
    }


    /**
     * Calculates the percentage of Guanine and Cytosine
     * in this Sequence.
     * @return Returns the percent of GC in the subsequence.
     */
    private int calculatePercentGC( String sequence )
    {
        int lengthGC = 0;
        int percentGC = -1;
        int seqLength = sequence.length();
        for( int i = 0; i < seqLength; i++ )
        {
            char currChar = sequence.charAt( i );
            if( (currChar == 'G') ||
                (currChar == 'C') )
            {
                lengthGC++;
            }
        }

        if( seqLength != 0 )
        {
            percentGC = (lengthGC * 100) / seqLength;
        }

        return percentGC;
    }


    protected void setTableMetaData( TableMetaData tableMetaData )
    {
        this.tableMetaData = tableMetaData;
    }
}
