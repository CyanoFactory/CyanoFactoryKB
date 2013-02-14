/* $Id: ObjectTable.java,v 1.1 2006/07/07 15:03:37 valerie Exp $ */
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

import com.sri.biospice.warehouse.schema.CrossReference;
import com.sri.biospice.warehouse.schema.Table;
import com.sri.biospice.warehouse.schema.DBID;
import java.util.Date;
import java.util.Vector;

/**
 * @author Valerie Wagner
 *         Date: Jul 29, 2004
 */
public interface ObjectTable extends Table {
    public void setCreationDate(Date creationDate);


    public void setModifiedDate(Date modifiedDate);


    /**
     * todo: deprecate this method
     * Adds a CrossReference.
     *
     * @param xid         An external identifier.
     * @param type        The type of the external identifier.
     *                    i.e what the external identifier represents.
     * @param version     The version of the xid.
     * @param datasetName The name of the dataset this id belongs to.
     * @deprecated Use typed version of this method instead
     */
    void addCrossReference(String xid, String type, String version, String datasetName);

    /**
     * Adds a CrossReference.
     *
     * @param xid         An external identifier.
     * @param type        The type of the external identifier.
     *                    i.e what the external identifier represents.
     * @param version     The version of the xid.
     * @param datasetName The name of the dataset this id belongs to.
     */
    void addCrossReference(String xid, CrossReference.Type type, String version, String datasetName);


    /**
     * Add a DBID.
     *
     * @param xid     An identifier for an object in this DataSet.
     * @param type    The type of the identifier.
     *                i.e what the identifier represents.
     * @param version The version of the xid.
     * @deprecated Use typed version instead
     */
    void addDBID(String xid, String type, String version);

    void addDBID(String xid, DBID.Type type, String version);


    /**
     * Adds a comment.
     *
     * @param line Comment to be added.
     */
    void addComment(String line);


    /**
     * Adds a synonym.
     *
     * @param syn Synonym to be added.
     */
    void addSynonym(String syn);

    /**
     * Adds a description to be stored in Description.Comm.
     * If <code>description</code> is null, no Description
     * table entry is created.
     *
     * @param description Text of description to be added
     */
    void addDescription(String description);


    /**
     * Adds a Vector of synonyms.
     *
     * @param syn Vector of Synonyms to be added.
     */
    void addSynonym(Vector syn);


    /**
     *
     */
    void setCreationDate();


    /**
     * @param ErrorMessage
     */
    void setErrorMessage(String ErrorMessage);


    /**
     * @param LineNumber
     */
    void setLineNumber(String LineNumber);


    /**
     *
     */
    void setLoadError();


    /**
     *
     */
    void setModifiedDate();


    long getWID();


    // todo: this should go away as a citation is an object table itself and may be linked to by many other objects
    public void addCitation(Citation citation);


    public void delete();


    public Vector getSynonyms();

    public Vector getDBIDs();

    public Vector getComments();

    public Vector getCrossreferences();

    public Vector getCitations();

    public Vector getDescriptions();

}

