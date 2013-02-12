/* $Id: Term.java,v 1.1 2006/07/07 15:03:38 valerie Exp $ */
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
 * @author Valerie Wagner
 *         Date: Aug 16, 2004
 */
public class Term extends ObjectTableBase {
    private static TableMetaData tableMetaData;
    private StringColumn nameCol;
    private StringColumn definitionCol;
    private StringColumn hierarchicalCol;
    private StringColumn rootCol;
    private StringColumn obsoleteCol;

    public static final String ENUM_OBSOLETE_TRUE = "T";
    public static final String ENUM_OBSOLETE_FALSE = "F";


    public Term(long datasetWID) {
        super(datasetWID, "Term");
        init();
    }


    public Term(long datasetWID, long termWID) {
        super(datasetWID, "Term", termWID);
        init();
    }


    private void init() {
        nameCol = new StringColumn("Name", this);
        definitionCol = new StringColumn("Definition", this);
        hierarchicalCol = new StringColumn("Hierarchical", this);
        rootCol = new StringColumn("Root", this);
        obsoleteCol = new StringColumn("Obsolete", this);
    }


    public void setName(String name) {
        nameCol.setValue(name);
    }

    public String getName() {
        return nameCol.getValue();
    }

    public String getDefinition() {
        return definitionCol.getValue();
    }

    public void setDefinition(String definition) {
        this.definitionCol.setValue(definition);
    }

    public String getObsolete() {
        return obsoleteCol.getValue();
    }

    public void setObsolete(String obsolete) {
        this.obsoleteCol.setValue(obsolete);
    }


    protected void setTableMetaData(TableMetaData tableMetaData) {
        this.tableMetaData = tableMetaData;
    }


    public TableMetaData getTableMetaData() {
        return tableMetaData;
    }

    public String getHierarchical() {
        return hierarchicalCol.getValue();
    }

    public void setHierarchical(String hierarchical) {
        this.hierarchicalCol.setValue(hierarchical);
    }

    public String getRoot() {
        return rootCol.getValue();
    }

    public void setRoot(String root) {
        this.rootCol.setValue(root);
    }
}
