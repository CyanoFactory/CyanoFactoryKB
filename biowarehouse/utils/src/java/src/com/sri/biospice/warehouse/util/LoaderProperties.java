/* $Id: LoaderProperties.java,v 1.3 2006/11/03 20:11:25 valerie Exp $ */
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
package com.sri.biospice.warehouse.util;

import org.apache.commons.cli.*;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;


import java.util.Enumeration;
import java.util.Iterator;
import java.util.Properties;
import java.util.Vector;

/**
 * @author Valerie Wagner
 *         Date: Jul 28, 2004
 */
public class LoaderProperties extends Properties
{
    private static final Logger log = Logger.getLogger( LoaderProperties.class );
    private Options cmdLineOptions;
    Vector requiredPropertyNames;
    public static final char PROPERTIES_FILE_CHAR = 'p';
    private String scriptName = "runLoader.sh";
    private String[] inputFiles;


    public LoaderProperties( Options commandLineOptions, Vector requiredProperties )
    {
        this.requiredPropertyNames = requiredProperties;
        this.cmdLineOptions = commandLineOptions;
    }

    public LoaderProperties() {

    }

    public void parseCommandLine( String[] args )
    {
        log.debug( "Parsing command line arguments" );
        CommandLineParser cliParser = new GnuParser();
        CommandLine line = null;

        try
        {
            // parse the command line arguments
            line = cliParser.parse( cmdLineOptions, args );

            // Check for help
            if( line.hasOption( 'h' ) )
            {
                printUsageAndExit();
            }

            // Create the properties object
            String propertiesFilename = line.getOptionValue( PROPERTIES_FILE_CHAR );
            loadPropertiesFile( propertiesFilename );

            Iterator parsedOptions = line.iterator();
            while( parsedOptions.hasNext() )
            {
                Option oneOption = (Option)parsedOptions.next();
                // If the file option takes multiple files, store them as an array
                if( oneOption.getLongOpt().equals( "file" ) && oneOption.hasArgs() )
                {
                    inputFiles = oneOption.getValues();
                }
                // All other options (but help) take only one argument
                if( !oneOption.getLongOpt().equals( "help" ) )
                {
                    log.debug( "Read command line option " + oneOption.getLongOpt() + "=" + oneOption.getValue() );
                    this.setProperty( oneOption.getLongOpt(), oneOption.getValue() == null ? "true" : oneOption.getValue() );
                }
            }
        }
        catch( ParseException exp )
        {
            log.error( "Parsing of command line options failed.  Reason: " + exp.getMessage() );
            printUsageAndExit();
        }
    }


    public void printUsageAndExit()
    {
        String usage = scriptName;
        String footer =
                new StringBuffer().append( "\nProperties may be set on the command line or in the properties " ).append( "file.  Values on the command line take precedence over those in a" )
                .append( " properties file. Properties in a property file are specified in name-value pairs." )
                .append( " For example: port=1234" )
                .toString();

        System.out.println();
        HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp( usage, null, cmdLineOptions, footer );
        System.exit( 1 );
    }


    private void loadPropertiesFile( String propertiesFilename )
    {


        if( propertiesFilename != null )
        {

            File propertiesFile = new File( propertiesFilename );
            if( !propertiesFile.exists() )
            {
                log.error( "The properties file specified (" + propertiesFilename + ") does not exist." );
                printUsageAndExit();
            }
            else
            {
                try
                {
                    log.debug( "Loading properties from property file: " + propertiesFilename );
                    this.load( new FileInputStream( propertiesFile ) );
                }
                catch( IOException e )
                {
                    log.error( "The properties file specified could not be read.  File was: " + propertiesFilename, e );
                    printUsageAndExit();
                }
            }
        }
    }


    public boolean validate()
    {
        // Report properties
        if( log.isDebugEnabled() )
        {
            Enumeration keys = keys();
            while( keys.hasMoreElements() )
            {
                String name = (String)keys.nextElement();
                if( !name.matches( ".*password.*" ) )
                {
                    log.debug( "Have property " + name + "=" + this.get( name ) );
                }
                else
                {
                    log.debug( "Have property " + name + "=[password value omitted]" );
                }
            }
        }

        boolean valid = true;
        for( int i = 0; i < requiredPropertyNames.size(); i++ )
        {
            String name = (String)requiredPropertyNames.elementAt( i );
            if( this.get( name ) == null )
            {
                log.error( "Missing required property: " + name );
                printUsageAndExit();
                valid = false;
            }
        }

        return valid;
    }


    public void setScriptName( String scriptName )
    {
        this.scriptName = scriptName;
    }


    public String getFile()
    {
        return getProperty( "file" );
    }


    public String getSourceDBMSType()
    {
        return getProperty( "sourcedbms" );
    }


    public String getSourceDatabaseHost()
    {
        return getProperty( "sourcehost" );
    }


    public String getSourceDatabasePort()
    {
        return getProperty( "sourceport" );
    }


    public String getSourceDatabaseUsername()
    {
        return getProperty( "sourceusername" );
    }


    public String getSourceDatabasePassword()
    {
        return getProperty( "sourcepassword" );
    }


    public String getSourceDatabaseName()
    {
        return getProperty( "sourcename" );
    }


        public String getDBMSType()
    {
        return getProperty( "dbms" );
    }


    public String getDatabaseHost()
    {
        return getProperty( "host" );
    }


    public String getDatabasePort()
    {
        return getProperty( "port" );
    }


    public String getDatabaseUsername()
    {
        return getProperty( "username" );
    }


    public String getDatabasePassword()
    {
        return getProperty( "password" );
    }


    public String getDatabaseName()
    {
        return getProperty( "name" );
    }


    public boolean isLoadAll()
    {
        return containsKey( "load-all" );
    }


    public String[] getMultipleFiles()
    {
        if( inputFiles != null )
        {
            return inputFiles;
        }
        else
        {
            // Split the "file" property based on one or more whitespace characters
            String allFiles = getFile();
            return allFiles.split( "\\s+" );
        }
    }


    public String getVersion()
    {
        return getProperty( "version" );
    }


    public String getReleaseDate()
    {
        return getProperty( "release" );
    }
}
