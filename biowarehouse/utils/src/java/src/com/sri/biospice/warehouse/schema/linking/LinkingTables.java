/* $Id: LinkingTables.java,v 1.2 2008/07/16 22:05:34 valerie Exp $ */
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

import com.sri.biospice.warehouse.schema.TableUtils;
import com.sri.biospice.warehouse.schema.object.*;
import org.apache.log4j.Logger;
import java.util.Vector;

/**
 * A convince class for adding entires to the linking tables
 * @author Priyanka Gupta
 * @author Valerie Wagner
 */
public class LinkingTables
{
    private static final Logger log = Logger.getLogger( LinkingTables.class );

    private Vector proteinWIDFunctionWID;
    private Vector geneWIDProteinWID;
    private Vector bioSourceWIDBioSubtypeWID;
    private Vector citationWIDOtherWID;
    private Vector bioSourceWIDGeneWID;
    private Vector geneWIDNucleicAcidWID;
    private Vector bioSourceWIDProteinWID;


    public LinkingTables()
    {
        geneWIDProteinWID = new Vector();
        proteinWIDFunctionWID = new Vector();
        bioSourceWIDBioSubtypeWID = new Vector();
        citationWIDOtherWID = new Vector();
        bioSourceWIDGeneWID = new Vector();
        geneWIDNucleicAcidWID = new Vector();
        bioSourceWIDProteinWID = new Vector();
    }


    /**
     * Stores all linking tables into the Database.
     */
    public void store()
    {
        store( geneWIDProteinWID, "geneWIDProteinWID" );
        store( proteinWIDFunctionWID, "proteinWIDFunctionWID" );
        store( bioSourceWIDBioSubtypeWID, "bioSourceWIDBioSubtypeWID" );
        store( bioSourceWIDGeneWID, "bioSourceWIDGeneWID" );
        store( geneWIDNucleicAcidWID, "geneWIDNucleicAcidWID" );
        store( bioSourceWIDProteinWID, "bioSourceWIDProteinWID" );
        store( citationWIDOtherWID, "citationWIDOtherWID" );
    }


    private void store( Vector tableVector, String tableName )
    {
        if( tableVector.size() > 0 )
        {
            log.debug( "Storing table vector " + tableName );
            TableUtils.storeTableVector( tableVector );
        }
        tableVector.clear();
    }


    /**
     * Add an entry to the BioSourceWIDBioSubtypeWID table
     * @param bioSource  The BioSource to be linked
     * @param bioSubtype The BioSubtype to be linked
     */
    public void addBioSourceWIDBioSubtypeWID( BioSource bioSource, BioSubtype bioSubtype )
    {
        bioSourceWIDBioSubtypeWID.add( new BioSourceWIDBioSubtypeWID( bioSource, bioSubtype ) );
//        log.debug( "Added bioSourceWIDBioSubtypeWID: " + bioSource.getWID() + ", " + bioSubtype.getWID() );
    }


    /**
     * Add a entry to the ProteinWIDFunctionWID
     * @param protein  The Protein to be linked
     * @param function The Function to be linked
     */
    public void addProteinWIDFunctionWID( Protein protein, Function function )
    {
        proteinWIDFunctionWID.add( new ProteinWIDFunctionWID( protein, function ) );
//        log.debug( "Added ProteinWidFunctionWid: " + protein.getWID() + ", " + function.getWID() );
    }


    /**
     * Add an entry to the GeneWIDProteinWID table
     * @param gene    The Gene to be linked
     * @param protein The Protein to be linked
     */
    public void addGeneWIDProteinWID( Gene gene, Protein protein )
    {
        geneWIDProteinWID.add( new GeneWIDProteinWID( gene, protein ) );
    }


    public void addCitationWIDOtherWID( Citation citation, ObjectTable otherTable )
    {
        citationWIDOtherWID.add( new CitationWIDOtherWID( citation, otherTable ) );

    }


    public void addBioSourceWIDGeneWID( BioSource bioSource, Gene gene )
    {
        bioSourceWIDGeneWID.add( new BioSourceWIDGeneWID( bioSource, gene ) );
    }


    public void addGeneWIDNucleicAcidWID( Gene gene, NucleicAcid nucleicAcid )
    {
        geneWIDNucleicAcidWID.add( new GeneWIDNucleicAcidWID( gene, nucleicAcid ) );
    }


    public void addBioSourceWIDProteinWID( BioSource bioSource, Protein protein )
    {
        bioSourceWIDProteinWID.add( new BioSourceWIDProteinWID( bioSource, protein ) );
    }
}
