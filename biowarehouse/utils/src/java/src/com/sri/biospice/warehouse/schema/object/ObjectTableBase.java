/* $Id */
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

import com.sri.biospice.warehouse.database.DeletionFactory;
import com.sri.biospice.warehouse.database.WarehouseManager;
import com.sri.biospice.warehouse.database.table.LongColumn;
import com.sri.biospice.warehouse.schema.*;
import org.apache.log4j.Logger;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.Date;
import java.util.Vector;

/**
 * A base class for all objects (tables) with WIDs. This class
 * stores DBIDs, CrossReferences, Synonyms, and Comments.
 *
 * @author Priyanka Gupta
 * @author Valerie Wagner
 */
public abstract class ObjectTableBase extends TableBase implements ObjectTable {
    private static final Logger log = Logger.getLogger(ObjectTableBase.class);

    protected long wid;
    protected long datasetWid;
    protected boolean alreadyStored = false;

    // References to other tables associated with this object
    // Every object table will have one entry in the Entry table
    private Entry entry;

    // And may have multiple entries in these tables:
    private Vector crossReferences = null;
    private Vector dbids = null;
    private Vector comments = null;
    private Vector synonyms = null;
    private Vector citations = null;
    private Vector descriptions = null;


    /**
     * Constructor for ObjectTableBase. Constructs a new object table with minimal required
     * initialization.
     *
     * @param datasetWID required for all objects in the database
     */
    public ObjectTableBase(long datasetWID, String tableName) {
        super(tableName);
        wid = WarehouseManager.getWarehouse().getNextWID();
        this.datasetWid = datasetWID;
        this.tableName = tableName;
        this.entry = new Entry(wid, datasetWID);

        // Set the values for the WID and DatasetWID
        // These columns are added to the this vector,
        // but won't ever be changed, so they are local variables here.
        LongColumn widCol = new LongColumn("WID", this);
        widCol.setValue(wid);

        LongColumn datasetCol = new LongColumn("DatasetWID", this);
        datasetCol.setValue(datasetWID);
    }


    /**
     * Constructor for an object table that has already been stored in the Warehouse
     *
     * @param datasetWID
     * @param tableName
     * @param thisWID
     */
    public ObjectTableBase(long datasetWID, String tableName, long thisWID) {
        super(tableName);
        alreadyStored = true;
        this.wid = thisWID;
        this.datasetWid = datasetWID;
        this.tableName = tableName;
        this.entry = new Entry(wid);
        // Don't add WID column here since WID is already
        // defined and we don't want to ever change it
        // Same goes for DatasetWID

        // todo !!! determine if this hurts anything
        LongColumn widCol = new LongColumn("WID", this);
        widCol.setValue(wid);

        LongColumn datasetCol = new LongColumn("DatasetWID", this);
        datasetCol.setValue(datasetWID);
    }


    /**
     * Gets the wid for this object.
     *
     * @return Returns the wid
     */
    public long getWID() {
        return wid;
    }


    public void addCrossReference(String xid, String type, String version, String databaseName) {
        if (crossReferences == null) {
            crossReferences = new Vector();
        }
        crossReferences.add(new CrossReference(xid, type, version, databaseName, wid));
    }

    public void addCrossReference(String xid, CrossReference.Type type, String version, String databaseName) {
        if (crossReferences == null) {
            crossReferences = new Vector();
        }
        crossReferences.add(new CrossReference(xid, type, version, databaseName, wid));
    }


    public void addDBID(String xid, String type, String version) {
        if (dbids == null) {
            dbids = new Vector();
        }
        DBID currDBId = new DBID(xid, type, version, wid);

        // Do not allow duplicate identical DBIDs
        if (!dbids.contains(currDBId)) {
            dbids.add(currDBId);
        }
    }

    public void addDBID(String xid, DBID.Type type, String version) {
        if (dbids == null) {
            dbids = new Vector();
        }
        DBID currDBId = new DBID(xid, type, version, wid);

        // Do not allow duplicate identical DBIDs
        if (!dbids.contains(currDBId)) {
            dbids.add(currDBId);
        }
    }

    /**
     * Adds a description to be stored in Description.Comm
     * If <code>description</code> is null, no Description
     * table entry is created.
     *
     * @param description Text of description to be added
     */
    public void addDescription(String description) {
        if (description != null) {
            if (descriptions == null) {
                descriptions = new Vector();
            }
            Description descr = new Description(this, description);
            descriptions.add(descr);
        }
    }


    public void addComment(String comment) {
        if (comment != null) {
            if (comments == null) {
                comments = new Vector();
            }
            comments.add(new CommentTable(this, comment));
        }
    }


    public void addSynonym(String synonym) {
        if (synonym != null && synonym.trim().length() > 0) {
            if (synonyms == null) {
                synonyms = new Vector();
            }
            synonyms.add(new SynonymTable(this, synonym));
        }
    }


    /**
     * @param syn A Vector of Strings containing synonyms
     */
    public void addSynonym(Vector syn) {
        if (synonyms == null) {
            synonyms = new Vector();
        }
        for (int i = 0; i < syn.size(); i++) {
            String synonym = (String) syn.elementAt(i);
            if (synonym.trim().length() > 0)
              synonyms.add(new SynonymTable(this, synonym));
        }
    }


    /**
     * Stores this object into the Database.
     */
    public void store() throws SQLException {
        try {
            if (alreadyStored) {
                update();
            } else {
                // Tell deriving classes to do their own store first
                doInsert();

                // Store the Entry
                entry.store();

                TableUtils.storeTableVector(dbids);
                TableUtils.storeTableVector(crossReferences);
                TableUtils.storeTableVector(comments);
                TableUtils.storeTableVector(synonyms);
                TableUtils.storeTableVector(citations);
                TableUtils.storeTableVector(descriptions);

                // Set already stored to true so the next time we just update
                // instead of insert.
                alreadyStored = true;
            }
        } catch (SQLException error) {
            log.error("Error storing table " + tableName, error);
            log.error(toString());
            throw error;
        }

    }


    public void setCreationDate() {
        entry.setCreationDate();
    }


    /**
     * @param CreationDate
     */
    public void setCreationDate(Date CreationDate) {
        entry.setCreationDate(CreationDate);
    }


    public void setErrorMessage(String ErrorMessage) {
        entry.setErrorMessage(ErrorMessage);
    }


    public void setLineNumber(String LineNumber) {
        entry.setLineNumber(LineNumber);
    }


    public void setLoadError() {
        entry.setLoadError();
    }


    public void setModifiedDate() {
        entry.setModifiedDate();
    }


    public void setModifiedDate(Date modifiedDate) {
        entry.setModifiedDate(modifiedDate);
    }


    protected ResultSet getRow() throws SQLException {
        ResultSet entryRS = null;
        String sqlLine = "select " + tableName + ".* from " + tableName + " where WID=" + wid;

        try {
            entryRS = WarehouseManager.getWarehouse().executeQuery(sqlLine);
        } catch (SQLException e) {
            log.error("Unable to execute query: " + sqlLine, e);
            throw e;
        }

        if (entryRS == null || entryRS.next() == false) {
            SQLException e = new SQLException("Entry not found in the " + tableName + " table with WID = " + wid);
            log.error(e.getMessage(), e);
            throw e;
        }

        return entryRS;
    }


    public void addCitation(Citation citation) {
        if (citations == null) {
            citations = new Vector();
        }
        citations.add(citation);
    }


    public void delete() {

        DeletionFactory.deletObjectTable(tableName, wid);
        DeletionFactory.deletOtherWIDTable("Entry", wid);
        DeletionFactory.deletOtherWIDTable("CrossReference", wid);
        DeletionFactory.deletOtherWIDTable("DBID", wid);
        DeletionFactory.deletOtherWIDTable("CommentTable", wid);
        DeletionFactory.deletOtherWIDTable("SynonymTable", wid);
        // todo: delete citations by looking up in CitationWIDOtherWID
        DeletionFactory.deleteLinkingTableEntries(wid);
    }


    public Vector getSynonyms() {
        return synonyms;
    }


    public Vector getDBIDs() {
        return dbids;
    }

    public Vector getComments() {
        return comments;
    }

    public Vector getCrossreferences() {
        return crossReferences;
    }

    public Vector getCitations() {
        return citations;
    }

    public Vector getDescriptions() {
        return descriptions;
    }

//    private void deleteRelated( String tableName, String columnName )
//    {
//        String sqlDeleteCommand = "delete from " + tableName + " where " + columnName + "=" + wid;
//        try
//        {
//            WarehouseManager.getWarehouse().executeUpdate( sqlDeleteCommand );
//        }
//        catch( SQLException e )
//        {
//            log.debug( "Could not perform delete: " + sqlDeleteCommand, e );
//        }
//    }

}
