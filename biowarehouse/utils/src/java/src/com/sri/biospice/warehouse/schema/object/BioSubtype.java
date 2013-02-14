/* $Id: BioSubtype.java,v 1.1 2006/07/07 15:03:37 valerie Exp $ */
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

import com.sri.biospice.warehouse.database.table.StringColumn;
import com.sri.biospice.warehouse.database.table.TableMetaData;

/**
 * An entry in the BioSubtype table
 * @author Valerie Wagner
 *         Date: May 12, 2004
 */
public class BioSubtype extends ObjectTableBase
{
    // Columns in the BioSubtype table
    private StringColumn typeCol;
    private StringColumn valueCol;

    private static TableMetaData tableMetaData = null;


    /**
     * Creates a new entry in BioSubtype
     * @param datasetWid
     */
    public BioSubtype( long datasetWid )
    {
        super( datasetWid, "BioSubtype" );
        init();
    }


    /**
     * Used to represent a BioSubtype already in the warehouse
     * @param datasetWid
     */
    public BioSubtype( long datasetWid, long bioSubtypeWID )
    {
        super( datasetWid, "BioSubtype", bioSubtypeWID );
        init();
    }


    private void init()
    {
        typeCol = new StringColumn( "Type", this );
        valueCol = new StringColumn( "Value", this );
    }


    public TableMetaData getTableMetaData()
    {
        return tableMetaData;
    }


    public void setType( String type )
    {
        typeCol.setValue( type );
    }


    public void setValue( String value )
    {
        valueCol.setValue( value );
    }


    protected void setTableMetaData( TableMetaData tableMetaData )
    {
        this.tableMetaData = tableMetaData;
    }
}
