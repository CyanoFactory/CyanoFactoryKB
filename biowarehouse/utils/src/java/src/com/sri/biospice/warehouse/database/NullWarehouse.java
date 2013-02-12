/* $Id: NullWarehouse.java,v 1.4 2006/11/07 17:51:17 valerie Exp $ */
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
package com.sri.biospice.warehouse.database;

import com.sri.bw.WIDGenerator;

/**
 * @author Valerie Wagner
 *         Date: Jul 19, 2004
 */
public class NullWarehouse extends NullDatabase implements Warehouse
{
    public long getNextWID()
    {
        return 0;
    }


    public long getNextSpecialWID()
    {
        return 0;
    }

    public WIDGenerator getWIDGenerator() {
        return null;
    }

    public WIDGenerator getSpecialWIDGenerator() {
        return null;
    }

    public WIDGenerator getReservedWIDGenerator()
    {
        return null;
    }


    public long getNextReservedWID() {
        return 0;
    }


    public void setWIDGenerator(WIDGenerator generator) {
        
    }

    public void setSpecialWIDGenerator(WIDGenerator generator) {
    }


}
