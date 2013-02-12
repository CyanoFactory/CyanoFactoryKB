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
 ******************************************************************* *////* $Id: LoaderProperties_Test.java,v 1.1 2006/07/06 20:11:50 valerie Exp $ */
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
//package com.sri.biospice.warehouse.database;
//
//import com.sri.biospice.warehouse.util.LoaderProperties;
//import junit.framework.TestCase;
//import org.apache.commons.cli.Option;
//import org.apache.commons.cli.OptionBuilder;
//import org.apache.commons.cli.Options;
//import org.apache.log4j.BasicConfigurator;
//
///**
// * @author valerie.wagner@sri.com
// *         Date: Jul 28, 2004
// */
//public class LoaderProperties_Test extends TestCase
//{
//    private Options commandLineOptions;
//
//
//    public LoaderProperties_Test( String testName )
//    {
//        super( testName );
//
//        BasicConfigurator.configure();
//
//        commandLineOptions = new Options();
//        Option propertiesFile = OptionBuilder.withArgName( "file" )
//            .hasArg()
//            .isRequired()
//            .withLongOpt( "properties" )
//            .withDescription( "Name of properties file" )
//            .create( 'p' );
//
//        Option inputFile = OptionBuilder.withArgName( "file" )
//            .hasArg()
//            .withLongOpt( "file" )
//            .withDescription( "Name of input data file" )
//            .create( 'f' );
//
//        commandLineOptions.addOption( propertiesFile );
//        commandLineOptions.addOption( inputFile );
//    }
//
//
//    public void testReadingFromPropertiesFile()
//    {
//        LoaderProperties properties = new LoaderProperties( commandLineOptions );
//
//        String[] args = new String[2];
//        args[0] = "-p";
//        args[1] = "data/sample.properties";
//
//        properties.parseCommandLine( args );
//        assertEquals( "Wrong input file", "testfile", properties.get( "file" ) );
//        assertEquals( "Expected that properties were valid", true, properties.validate() );
//    }
//
//
//    public void testCommandLineOverwritesPropertiesFile()
//    {
//        LoaderProperties properties = new LoaderProperties( commandLineOptions );
//
//        String[] args = new String[4];
//        args[0] = "-p";
//        args[1] = "data/sample.properties";
//        args[2] = "-f";
//        args[3] = "myfile";
//
//        properties.parseCommandLine( args );
//        assertEquals( "Wrong input file", "myfile", properties.get( "file" ) );
//        assertEquals( "Expected that properties were valid", true, properties.validate() );
//    }
//
//}
