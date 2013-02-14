/* $Id: Table.java,v 1.1 2006/07/07 15:03:36 valerie Exp $ */
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

import com.sri.biospice.warehouse.database.table.TableMetaData;
import java.sql.SQLException;
import java.util.Vector;


/**
 * @author Valerie Wagner
 *         Date: Jun 14, 2004
 */
public interface Table
{
    public void store() throws SQLException;


    public TableMetaData getTableMetaData();


    public Vector getColumns();


    public String getTableName();


    public void load() throws SQLException;


    Enum [] getAllowedValues( String columnName );
}
