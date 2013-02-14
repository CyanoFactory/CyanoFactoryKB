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
 ******************************************************************* *////* $Id: Protein_Test.java,v 1.1 2006/07/06 20:11:51 valerie Exp $ */
///******************************************************************************
// * Copyright 2004 SRI International.  All rights reserved.                     *
// *                                                                             *
// * The material contained in this file is confidential and proprietary to SRI  *
// * International and may not be reproduced, published, or disclosed to others  *
// * without authorization from SRI International.                               *
// *                                                                             *
// * DISCLAIMER OF WARRANTIES                                                    *
// *                                                                             *
// * SRI International MAKES NO REPRESENTATIONS OR WARRANTIES ABOUT THE          *
// * SUITABILITY OF THE SOFTWARE, EITHER EXPRESS OR IMPLIED, INCLUDING BUT NOT   *
// * LIMITED TO THE IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A         *
// * PARTICULAR PURPOSE, OR NON-INFRINGEMENT. SRI International SHALL NOT BE     *
// * LIABLE FOR ANY DAMAGES SUFFERED BY LICENSEE AS A RESULT OF USING, MODIFYING *
// * OR DISTRIBUTING THIS SOFTWARE OR ITS DERIVATIVES                            *
// *******************************************************************************/
//package com.sri.biospice.warehouse.schema;
//
//import com.sri.biospice.warehouse.util.WarehouseTestUtils;
//import junit.framework.TestCase;
//import java.sql.SQLException;
//
///**
// * @author Valerie.Wagner@sri.com
// *         Date: Aug 16, 2004
// */
//public class Protein_Test extends TestCase
//{
//    public Protein_Test( String testName ) throws SQLException
//    {
//        super( testName );
//        WarehouseTestUtils.initWarehouseTest();
//    }
//
//
//    public void testInsertHugeClob() throws SQLException
//    {
//        int sequenceLength = 40000;
//        DataSet dataset = new DataSet( "test" );
//        dataset.store();
//
//        Protein protein = new Protein( dataset.getWID() );
//
//        // Create a huge clob and store it
//        StringBuffer aaSequence = new StringBuffer();
//        aaSequence.setLength( sequenceLength );
//        for( int i = 0; i < sequenceLength; i++ )
//        {
//            aaSequence.setCharAt( i, 'x' );
//        }
//
//        protein.setAASequence( aaSequence.toString() );
//        protein.store();
//
//
//        // Change the clob a bit and update it
//        aaSequence.setCharAt( 0, 'a' );
//        protein.setAASequence( aaSequence.toString() );
//        protein.update();
//
//        WarehouseTestUtils.endWarehouseTest();
//
//    }
//
//
//}
