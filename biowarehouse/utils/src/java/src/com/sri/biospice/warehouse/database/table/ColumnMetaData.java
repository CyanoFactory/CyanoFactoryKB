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
package com.sri.biospice.warehouse.database.table;

import com.sri.biospice.warehouse.schema.Table;

/**
 * @author Valerie Wagner
 *         Date: Apr 2, 2005
 */
public interface ColumnMetaData {
    String getColumnName();

    int getColumnSQLType();

    int getMaxColumnSize();

    int getParameterIndex();

    void setParameterIndex( int parameterIndex );

    String getColumnTypeName();

    void valueSatisfied();

    boolean isSatisifed();

    void clearValue();

    Column makeColumn(Table table);
}
