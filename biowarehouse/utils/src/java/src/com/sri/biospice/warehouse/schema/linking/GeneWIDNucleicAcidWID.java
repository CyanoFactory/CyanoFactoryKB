/* $Id: GeneWIDNucleicAcidWID.java,v 1.1 2006/07/07 15:03:37 valerie Exp $ */
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
package com.sri.biospice.warehouse.schema.linking;

import com.sri.biospice.warehouse.database.table.TableMetaData;
import com.sri.biospice.warehouse.schema.object.Gene;
import com.sri.biospice.warehouse.schema.object.NucleicAcid;

/**
 * @author Valerie.Wagner@sri.com
 *         Date: Sep 9, 2004
 */
public class GeneWIDNucleicAcidWID extends LinkingTable
{
    private static TableMetaData tableMetaData;


    public GeneWIDNucleicAcidWID( Gene gene, NucleicAcid nucleicAcid )
    {
        super( "GeneWID", "NucleicAcidWID" );
        this.setWID1( gene.getWID() );
        this.setWID2( nucleicAcid.getWID() );
    }


    protected void setTableMetaData( TableMetaData tableMetaData )
    {
        this.tableMetaData = tableMetaData;
    }


    public TableMetaData getTableMetaData()
    {
        return tableMetaData;
    }
}
