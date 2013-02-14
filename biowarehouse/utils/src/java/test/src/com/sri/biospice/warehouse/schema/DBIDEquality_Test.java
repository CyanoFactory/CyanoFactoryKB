/* $Id: DBIDEquality_Test.java,v 1.2 2006/07/07 01:06:36 valerie Exp $ */
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

import com.sri.biospice.warehouse.WarehouseFunctionalTest;
import java.sql.SQLException;
import java.util.Hashtable;
import java.util.Vector;

/**
 * @author Valerie.Wagner@sri.com
 *         Date: Oct 15, 2004
 */
public class DBIDEquality_Test extends WarehouseFunctionalTest
{
    DBID dbid1 = new DBID( "xid", "type", "version", 5 );
    DBID dbid2 = new DBID( "xid", "type", "version", 5 );
    DBID dbid3 = new DBID( "xid", "type", "version", 5 );
    DBID dbid4 = new DBID( "xid", "type", "Version", 5 );
    DBID dbid5 = new DBID( "xid", "type", null, 5 );
    DBID dbid6 = new DBID( "xid", "type", "version", 50 );


    public DBIDEquality_Test( String testName ) throws SQLException
    {
        super( testName );
    }


    public void testDBIDEquals()
    {
        // Test reflexivity
        assertEquals( "These should be equal 0", dbid1, dbid1 );

        // Test symmetry
        assertEquals( "These should be equal 1", dbid1, dbid2 );
        assertEquals( "These should be equal 2", dbid2, dbid1 );

        // Test transitivity
        assertEquals( "These should be equal 3", dbid2, dbid3 );
        assertEquals( "These should be equal 4", dbid1, dbid3 );

        // Test not equals
        assertEquals( "These should not be equal 5", true, !dbid1.equals( dbid4 ) );
        assertEquals( "These should not be equal 6", true, !dbid1.equals( dbid5 ) );
        assertEquals( "These should not be equal 7", true, !dbid1.equals( dbid6 ) );
    }


    public void testDBIDHashCodes()
    {
        assertEquals( "Hash code should be repeatable 1", dbid1.hashCode(), dbid1.hashCode() );
        assertEquals( "Same objects should have same hash codes 2", dbid1.hashCode(), dbid2.hashCode() );
        assertEquals( "Same objects should have same hash codes 3", dbid1.hashCode(), dbid3.hashCode() );

        assertEquals( "Different objects should have different hash codes 4", true, !(dbid1.hashCode() == dbid4.hashCode()) );
        assertEquals( "Different objects should have different hash codes 5", true, !(dbid1.hashCode() == dbid5.hashCode()) );
        assertEquals( "Different objects should have different hash codes 6", true, !(dbid1.hashCode() == dbid6.hashCode()) );
    }


    public void testDBIDEqualsInVector()
    {
        Vector dbids = new Vector();
        dbids.add( dbid1 );
        dbids.add( dbid4 );

        assertEquals( "Vector should contain this dbid 1", true, dbids.contains( dbid1 ) );
        assertEquals( "Vector should contain this dbid 2", true, dbids.contains( dbid4 ) );
        assertEquals( "Vector should not contain this dbid 3", true, !dbids.contains( dbid5 ) );
    }


    public void testDBIDEqualsInHashtable()
    {
        Hashtable dbids = new Hashtable();
        dbids.put( dbid1, dbid1 );
        dbids.put( dbid4, dbid4 );

        assertEquals( "Vector should contain this dbid 1", true, dbids.contains( dbid1 ) );
        assertEquals( "Vector should contain this dbid 2", true, dbids.contains( dbid4 ) );
        assertEquals( "Vector should contain this dbid 2a", true, dbids.contains( dbid2 ) );
        assertEquals( "Vector should not contain this dbid 3", true, !dbids.contains( dbid5 ) );

        assertEquals( "Vector should contain this dbid 4", true, dbids.containsKey( dbid1 ) );
        assertEquals( "Vector should contain this dbid 5", true, dbids.containsKey( dbid4 ) );
        assertEquals( "Vector should contain this dbid 5a", true, dbids.containsKey( dbid2 ) );
        assertEquals( "Vector should not contain this dbid 6", true, !dbids.containsKey( dbid5 ) );
    }
}

