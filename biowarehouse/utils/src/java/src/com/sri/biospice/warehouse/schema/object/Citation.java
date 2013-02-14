/* $Id: Citation.java,v 1.4 2008/02/05 23:24:05 valerie Exp $ */
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

import com.sri.biospice.warehouse.database.table.IntegerColumn;
import com.sri.biospice.warehouse.database.table.StringColumn;
import com.sri.biospice.warehouse.database.table.TableMetaData;

/**
 * @author Valerie Wagner
 *         Date: Aug 17, 2004
 */
public class Citation extends ObjectTableBase {
    private static TableMetaData tableMetaData;
    private StringColumn citationCol;
    private IntegerColumn pmidCol;
    private StringColumn titleCol;
    private StringColumn authorsCol;
    private StringColumn publicationCol;
    private StringColumn publisherCol;
    private StringColumn editorCol;
    private StringColumn yearCol;
    private StringColumn volumeCol;
    private StringColumn issueCol;
    private StringColumn pagesCol;
    private StringColumn uriCol;


    public Citation(long datasetWID) {
        super(datasetWID, "Citation");
        init();
    }


    /**
     * Represents a citation already in the warehouse
     *
     * @param datasetWID
     * @param citationWID
     */
    public Citation(long datasetWID, long citationWID) {
        super(datasetWID, "Citation", citationWID);
        init();
    }


    private void init() {
        citationCol = new StringColumn("Citation", this);
        pmidCol = new IntegerColumn("PMID", this);

        titleCol = new StringColumn("Title", this);
        authorsCol = new StringColumn("Authors", this);
        publicationCol = new StringColumn("Publication", this);
        publisherCol = new StringColumn("Publisher", this);
        editorCol = new StringColumn("Editor", this);
        yearCol = new StringColumn("Year", this);
        volumeCol = new StringColumn("Volume", this);
        issueCol = new StringColumn("Issue", this);
        pagesCol = new StringColumn("Pages", this);
        uriCol = new StringColumn("URI", this);
    }


    public void setCitation(String citation) {
        citationCol.setValue(citation);
    }


    public void setPMID(String pmid) {
        pmidCol.setValue(pmid);
    }


    protected void setTableMetaData(TableMetaData tableMetaData) {
        this.tableMetaData = tableMetaData;
    }


    public TableMetaData getTableMetaData() {
        return tableMetaData;
    }
    
    public String toString() {
        String nl = System.getProperty("line.separator");
        StringBuffer sb = new StringBuffer();

        String dec = "========";

        sb.append(dec).append(" [S] com.sri.biospice.warehouse.schema.object.Citation ").append(dec).append(nl);
        sb.append("#Super-Class of com.sri.biospice.warehouse.schema.object.Citation").append(nl);
        sb.append(super.toString()).append(nl);
        sb.append("tableMetaData=").append(tableMetaData).append(nl);
        sb.append("citationCol=").append(citationCol).append(nl);
        sb.append("pmidCol=").append(pmidCol).append(nl);
        sb.append("titleCol=").append(titleCol).append(nl);
        sb.append("authorsCol=").append(authorsCol).append(nl);
        sb.append("publicationCol=").append(publicationCol).append(nl);
        sb.append("publisherCol=").append(publisherCol).append(nl);
        sb.append("editorCol=").append(editorCol).append(nl);
        sb.append("yearCol=").append(yearCol).append(nl);
        sb.append("volumeCol=").append(volumeCol).append(nl);
        sb.append("issueCol=").append(issueCol).append(nl);
        sb.append("pagesCol=").append(pagesCol).append(nl);
        sb.append("uriCol=").append(uriCol).append(nl);

        sb.append(dec).append(" [E] com.sri.biospice.warehouse.schema.object.Citation ").append(dec);

        return sb.toString();
    }

    public Integer getPmid() {
        return pmidCol.getValue();
    }

    public String getCitation() {
        return citationCol.getValue();
    }
}
