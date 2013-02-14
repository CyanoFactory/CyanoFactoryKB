/* $Id: TableUtils.java,v 1.1 2006/07/07 15:03:37 valerie Exp $ */
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

import java.sql.SQLException;
import java.util.Vector;

/**
 * @author Valerie Wagner
 *         Date: Aug 13, 2004
 */
public class TableUtils
{
    public static void storeTableVector( Vector tables )
    {
        if( tables != null )
        {
            for( int i = 0; i < tables.size(); i++ )
            {
                try
                {
                    ((Table)tables.elementAt( i )).store();
                }
                catch( SQLException e )
                {
                    // Already logged the error in the store method
                }
            }
        }
    }


    public static void storeTable( Table table )
    {
        if( table != null )
        {
            try
            {
                table.store();
            }
            catch( SQLException e )
            {
                // Already logged the error in the store method
            }
        }
    }
}
