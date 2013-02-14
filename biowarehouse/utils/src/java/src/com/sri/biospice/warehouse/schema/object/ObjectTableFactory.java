/* $Id: ObjectTableFactory.java,v 1.2 2007/08/29 23:47:10 nguo Exp $ */
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


/**
 * Used to create instances of an ObjectTable of a specific type. E.g. to delete the object.
 * @author Valerie.Wagner@sri.com
 *         Date: Sep 28, 2004
 */
public class ObjectTableFactory
{

    public static ObjectTable createObjectTable( String objectTableName, long datasetWID, long objectTableWID )
    {
        if( objectTableName.equalsIgnoreCase( "protein" ) )
        {
            return new Protein( datasetWID, objectTableWID );
        }
        else if( objectTableName.equalsIgnoreCase( "biosource" ) )
        {
            return new BioSource( datasetWID, objectTableWID );
        }
        else if( objectTableName.equalsIgnoreCase( "biosubtype" ) )
        {
            return new BioSubtype( datasetWID, objectTableWID );
        }
        else if( objectTableName.equalsIgnoreCase( "citation" ) )
        {
            return new Citation( datasetWID, objectTableWID );
        }
        else if( objectTableName.equalsIgnoreCase( "feature" ) )
        {
            return new Feature( datasetWID, objectTableWID );
        }
        else if( objectTableName.equalsIgnoreCase( "gene" ) )
        {
            return new Gene( datasetWID, objectTableWID );
        }
        else if( objectTableName.equalsIgnoreCase( "function" ) )
        {
            return new Function( datasetWID, objectTableWID );
        }
        else if( objectTableName.equalsIgnoreCase( "nucleicacid" ) )
        {
            return new NucleicAcid( datasetWID, objectTableWID );
        }
        else if( objectTableName.equalsIgnoreCase( "subsequence" ) )
        {
            return new SubSequence( datasetWID, objectTableWID );
        }
        else if( objectTableName.equalsIgnoreCase( "support" ) )
        {
          return new Support( datasetWID, objectTableWID );
        }
        else if( objectTableName.equalsIgnoreCase( "term" ) )
        {
            return new Term( datasetWID, objectTableWID );
        }

        return null;
    }

}
