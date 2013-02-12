/* $Id: LoaderMain.java,v 1.11 2008/10/15 21:48:45 valerie Exp $ */
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

import com.sri.biospice.warehouse.database.Warehouse;
import com.sri.biospice.warehouse.database.WarehouseManager;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;
import org.apache.log4j.xml.DOMConfigurator;
import java.io.File;
import java.net.InetAddress;
import java.net.URL;
import java.net.UnknownHostException;
import java.sql.DatabaseMetaData;
import java.sql.SQLException;
import java.text.DateFormat;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.Date;
import java.util.GregorianCalendar;
import java.util.Vector;

/**
 * Common main starting point for all loaders
 *
 * @author Valerie Wagner
 *         Date: Jul 29, 2004
 */
public abstract class LoaderMain {
    private static final Logger log = Logger.getLogger(LoaderMain.class);
    protected LoaderProperties properties;
    protected String loaderName;

    public static final String STANDARD_LOG_CONFIG = "bw-log4j-config.xml";
    public static final String STANDARD_LOG_CONFIG_PATH = "/com/sri/biospice/warehouse/util/";

    /**
     * Time loader started
     */
    private Calendar startTime;
    private final int SECONDS_PER_MINUTE = 60;
    private final int MINUTES_PER_HOUR = 60;

    private Options commandLineOptions = new Options();
    private Vector requiredProperties = new Vector();

    private Date earliestSupportedDataRelease;
    private Date latestSupportedDataRelease;


    /**
     * Used for checking whether file user supplied is supported by our loader, or not
     */
    public enum DataSupport {
        SUPPORTED,
        NOT_SUPPORTED,
        UNKNOWN
    }

    public static final Option QUIT_OPTION = OptionBuilder.withLongOpt("quit-after")
            .hasArg()
            .withArgName("num entries")
            .withDescription("Quit after parsing X number of entries. (For testing purposes only.)")
            .create("q");

    public static final Option HELP_OPTION = OptionBuilder.withLongOpt("help")
            .withDescription("Print usage instructions")
            .create('h');

    public static final Option PROPERTIES_FILE_OPTION = OptionBuilder.withArgName("file")
            .hasArg()
            .withLongOpt("properties")
            .withDescription("Name of properties file")
            .create('p');

    public static final Option DB_HOST_OPTION = OptionBuilder.withArgName("host")
            .hasArg()
            .withLongOpt("host")
            .withDescription("Name or IP address of database server host")
            .create('s');

    public static final Option DB_NAME_OPTION = OptionBuilder.withArgName("name")
            .hasArg()
            .withLongOpt("name")
            .withDescription("Name or SID of database")
            .create('n');

    public static final Option DBMS_OPTION = OptionBuilder.withArgName("dbms")
            .hasArg()
            .withLongOpt("dbms")
            .withDescription("DBMS type (mysql or oracle)")
            .create('d');

    public static final Option DB_USERNAME_OPTION = OptionBuilder.withArgName("username")
            .hasArg()
            .withLongOpt("username")
            .withDescription("Username for connection to the database")
            .create('u');

    public static final Option DB_PASSWORD_OPTION = OptionBuilder.withArgName("password")
            .hasArg()
            .withLongOpt("password")
            .withDescription("Password for connection to the database")
            .create('w');

    public static final Option DB_PORT_OPTION = OptionBuilder.withArgName("port")
            .hasArg()
            .withLongOpt("port")
            .withDescription("Port database server is listening at")
            .create('t');

    public static final Option LOAD_INVALID_OPTION = OptionBuilder
            .withLongOpt("load-all")
            .withDescription("Load all data, including data found to be suspect")
            .create('l');

    public static final Option VERSION_OPTION = OptionBuilder
            .hasArg()
            .withArgName("version number")
            .withLongOpt("version")
            .withDescription("Version number of the input dataset")
            .create('v');

    public static final Option RELEASE_OPTION = OptionBuilder
            .hasArg()
            .withArgName("release date")
            .withLongOpt("release")
            .withDescription("Release date of the input dataset")
            .create('r');

    public static final Option MULTIPLE_INPUT_FILES_OPTION = OptionBuilder.withArgName("file")
            .hasArgs()
            .withLongOpt("file")
            .withDescription("Name of input data file(s)")
            .create('f');

    public static final Option SINGLE_INPUT_FILE_OPTION = OptionBuilder.withArgName("file")
            .hasArg()
            .withLongOpt("file")
            .withDescription("Name of input data file")
            .create('f');

    public static final Option SOURCE_DB_HOST_OPTION = OptionBuilder.withArgName("sourcehost")
            .hasArg()
            .withLongOpt("sourcehost")
            .withDescription("Name or IP address of souce database server host")
            .create('S');

    public static final Option SOURCE_DB_NAME_OPTION = OptionBuilder.withArgName("sourcename")
            .hasArg()
            .withLongOpt("sourcename")
            .withDescription("Name or SID of source database")
            .create('N');

    public static final Option SOURCE_DBMS_OPTION = OptionBuilder.withArgName("sourcedbms")
            .hasArg()
            .withLongOpt("sourcedbms")
            .withDescription("source DBMS type (mysql or oracle)")
            .create('D');

    public static final Option SOURCE_DB_USERNAME_OPTION = OptionBuilder.withArgName("sourceusername")
            .hasArg()
            .withLongOpt("sourceusername")
            .withDescription("Username for connection to the source database")
            .create('U');

    public static final Option SOURCE_DB_PASSWORD_OPTION = OptionBuilder.withArgName("sourcepassword")
            .hasArg()
            .withLongOpt("sourcepassword")
            .withDescription("Password for connection to the source database")
            .create('W');

    public static final Option SOURCE_DB_PORT_OPTION = OptionBuilder.withArgName("sourceport")
            .hasArg()
            .withLongOpt("sourceport")
            .withDescription("source Port database server is listening at")
            .create('T');

    /**
     * Constructs a LoaderMain, uses a single input file, builds standard options, and parses command-line arguments
     * Used by most loaders
     * @param args  Arguments passed in by user on command line
     * @param loaderName  Display-friendly name of the loader
     * @param log4jConfigFilename Name of the log4j config file
     * @param log4jConfigPath Path to the log4j config file in the JAR
     * @param scriptName Name of the script used to run the loader
     */
    public LoaderMain(String[] args, String loaderName, String scriptName) {
        this(args, loaderName, scriptName, false);
    }


     /**
     * Constructs a LoaderMain, choice of single or multiple input files, builds standard options, and parses command-line arguments
     * Used by loaders that may parse more than one input data file
     * @param args  Arguments passed in by user on command line
     * @param loaderName  Display-friendly name of the loader
     * @param log4jConfigFilename Name of the log4j config file
     * @param log4jConfigPath Path to the log4j config file in the JAR
     * @param scriptName Name of the script used to run the loader
     * @param multipleFiles True if parsing multiple input files; false if parsing a single input file
     */
    public LoaderMain(String[] args, String loaderName, String scriptName, boolean multipleFiles) {
        startTime = GregorianCalendar.getInstance();
        this.loaderName = loaderName;
        configureLogging();

        buildStandardOptions(multipleFiles);

        properties = new LoaderProperties(commandLineOptions, requiredProperties);
        properties.setScriptName(scriptName);
        properties.parseCommandLine(args);
    }


     /**
     * Constructs a LoaderMain, choice of single or multiple input files, builds standard options, but does not
      * parse command line options.  Clients should call parseCommandLine() separately.
     * Used by loaders that want the standard set of options, but also add more options via addOption().
     * @param loaderName  Display-friendly name of the loader
     * @param log4jConfigFilename Name of the log4j config file
     * @param log4jConfigPath Path to the log4j config file in the JAR
     * @param multipleFiles True if parsing multiple input files; false if parsing a single input file
     */
    public LoaderMain(String loaderName, boolean multipleFiles) {
        startTime = GregorianCalendar.getInstance();
        this.loaderName = loaderName;
        configureLogging();
        buildStandardOptions(multipleFiles);
    }


    /**
     * Constructor used by classes wishing to specify their own options (will not build
     * standard options or parse command line argments).
     * Users must call addOption() and parseCommandLine()
     * @param loaderName
     * @param log4jConfigFilename
     * @param log4jConfigPath
     */
    public LoaderMain(String loaderName, String log4jConfigFilename, String log4jConfigPath) {
        startTime = GregorianCalendar.getInstance();
        this.loaderName = loaderName;
        configureLogging();
    }


    /**
     * Parses command line options and verifies 
     * @param args
     * @param scriptName
     */
    public void parseCommandLine(String[] args, String scriptName) {
        properties = new LoaderProperties(commandLineOptions, requiredProperties);
        properties.setScriptName(scriptName);
        properties.parseCommandLine(args);
    }

    /**
     * Builds the standard set of command line options
     * @param multipleFiles Whether to use the option for multiple files (true) or a single input file (false)
     */
    private void buildStandardOptions(boolean multipleFiles) {
        addOption(HELP_OPTION, false);
        addOption(PROPERTIES_FILE_OPTION, false);
        addOption(DB_HOST_OPTION, true);
        addOption(DB_NAME_OPTION, true);
        addOption(DBMS_OPTION, true);
        addOption(DB_USERNAME_OPTION, true);
        addOption(DB_PASSWORD_OPTION, true);
        addOption(DB_PORT_OPTION, true);
        addOption(RELEASE_OPTION, true);
        addOption(VERSION_OPTION, true);
        addOption(LOAD_INVALID_OPTION, false);
        addOption(QUIT_OPTION, false);

        if (multipleFiles) {
            addOption(MULTIPLE_INPUT_FILES_OPTION, true);
        } else {
            addOption(SINGLE_INPUT_FILE_OPTION, true);
        }
    }

    /**
     * Add an option to the allowed command line options
     * @param option The option to add
     * @param required Whether or not the user is required to specify this option
     */
    public void addOption(Option option, boolean required) {
        commandLineOptions.addOption(option);
        if (required) {
            requiredProperties.add(option.getLongOpt());
        }
    }

    /**
     * Configure the log4j system with the config file.  If the user
     * has placed the file at the same location as the jar file, use it.
     * Otherwise, use the one packaged up with the jar.
     */
    private void configureLogging() {
        String log4jConfigFilename = STANDARD_LOG_CONFIG;
        String log4jConfigPath = STANDARD_LOG_CONFIG_PATH;
        if (log4jConfigFilename != null) {
            File localFile = new File(log4jConfigFilename);
            if (localFile.exists()) {
                DOMConfigurator.configure(log4jConfigFilename);
            } else {
                String packagedConfigFilename = log4jConfigPath + log4jConfigFilename;
                URL packagedURL = this.getClass().getResource(packagedConfigFilename);

                if (packagedURL != null) {
                    DOMConfigurator.configure(packagedURL);
                } else {
                    BasicConfigurator.configure();
                    log.error("Could not find log4j configuration file: " + log4jConfigFilename);
                    System.exit(1);
                }
            }
        } else {
            BasicConfigurator.configure();
        }
    }


    /**
     * Establish the connection to the Warehouse.
     */
    public void connectToWarehouse() {
        String databaseHost = properties.getDatabaseHost();
        String databasePort = properties.getDatabasePort();
        String databaseName = properties.getDatabaseName();
        String databaseDBMSType = properties.getDBMSType();
        String databaseUsername = properties.getDatabaseUsername();
        String databasePassword = properties.getDatabasePassword();


        Warehouse theWarehouse = WarehouseManager.initWarehouse(databaseDBMSType, databaseHost, databaseName, databasePort);

        try {
            theWarehouse.connectToDatabase(databaseUsername, databasePassword);
            log.info("***************");
            log.info("databaseDBMSType = " + databaseDBMSType);
            log.info("databaseHost = " + databaseHost);
            log.info("databasePort = " + databasePort);
            log.info("databaseName = " + databaseName);
            log.info("databaseUsername = " + databaseUsername);
            DatabaseMetaData meta = theWarehouse.getMetaData();
            if(meta != null ){
                log.info("URL = " + meta.getURL());
                log.info("RDBMS = " + meta.getDatabaseProductName() + ", version: " + meta.getDatabaseProductVersion());
                log.info("JDBC driver = " + meta.getDriverName() + ", version: " + meta.getDriverVersion());
            }
            log.info("Current system user = " + System.getProperty("user.name"));
            try {
                log.info("Running client on = " + InetAddress.getLocalHost().getHostName() +
                        " [" + InetAddress.getLocalHost().getHostAddress() + "]");
            } catch (UnknownHostException e) {
                log.debug(e);
            }
            log.info("***************");
        } catch (SQLException e) {
            log.error("Could not connect to the Warehouse.  Reason: " + e.getMessage());
            System.exit(1);
        }
    }


    /**
     * Close connection to warehouse
     */
    public void disconnectFromWarehouse() {
        try {
            WarehouseManager.getWarehouse().commit();
            WarehouseManager.getWarehouse().close();
        } catch (SQLException e) {
            log.error("Could not close connection to the Warehouse.  Reason: " + e.getMessage());
        }
    }


    /**
     * Checks if the input file exists and is readable
     *
     * @return true if the file exists, and is a file (not a dir), and is readable.
     */
    protected boolean fileIsReadable() {
        String filepath = properties.getFile();
        return fileIsReadable(filepath);
    }


    /**
     * Checks if the input file exists and is readable
     *
     * @param filepath Path to the input file
     * @return true if the file exists, and is a file (not a dir), and is readable
     */
    public static boolean fileIsReadable(String filepath) {

        boolean fileIsGood = true;

        File file = new File(filepath);
        if (!file.exists()) {
            log.error("File does not exist: [" + filepath + "]");
            fileIsGood = false;
        } else if (!file.isFile()) {
            log.error("This is not a file: [" + filepath + "]");
            fileIsGood = false;
        } else if (!file.canRead()) {
            log.error("This file is not readable: [" + filepath + "]");
            fileIsGood = false;
        }

        return fileIsGood;
    }


    /**
     * Checks that all files listed in multiple file inputs are readable
     *
     * @return True if all input files are readable
     */
    protected boolean filesAreReadable() {
        boolean readable = true;
        String[] filepaths = properties.getMultipleFiles();
        for (int i = 0; i < filepaths.length; i++) {
            String filepath = filepaths[i];
            readable = readable & fileIsReadable(filepath);
        }
        return readable;
    }


    /**
     * Print loader info to log
     */
    protected void reportVersion() {
        log.info("Warehouse loader: " + loaderName);
        log.info("Loader version: " + getLoaderVersionNumber());
        log.info("Loader build: " + getLoaderBuildNumber());
    }


    /**
     * Print time elapsed during loader run
     */
    public void reportTimeElapsed() {
        // Calculate how long the run took
        // todo: Shouldn't there be a simpler way to do this?
        Calendar endTime = GregorianCalendar.getInstance();
        long startTimeMS = startTime.getTimeInMillis();
        long endTimeMS = endTime.getTimeInMillis();
        long differenceMS = endTimeMS - startTimeMS;
        long differenceS = differenceMS / 1000;
        long seconds = differenceS % SECONDS_PER_MINUTE;
        long differenceM = differenceS / SECONDS_PER_MINUTE;
        long minutes = differenceM % MINUTES_PER_HOUR;
        long differenceH = differenceM / MINUTES_PER_HOUR;

//        SimpleDateFormat dateFormate = new SimpleDateFormat("HH:mm:ss");
//        Calendar cal = new GregorianCalendar();
//        cal.set(Calendar.HOUR, (int)differenceH);
//        cal.set(Calendar.MINUTE, (int)minutes);
//        cal.set(Calendar.SECOND, (int)seconds);
        log.info("Time elapsed: " + differenceH + ":" +
                minutes + ":" + seconds);
//        log.info("Time elapsed: " + dateFormate.format(cal.getTime()));
    }


    /**
     * Called in the event of an unrecoverable error.  Program will log the message and exit.     *
     *
     * @param message Message to use explaining why program failed
     */
    public static void fail(String message) {
        log.fatal(message);
        System.exit(1);
    }


    /**
     * Determine if we should try to parse input data
     *
     * @return true if data version is supported or unknown
     */
    public boolean dataVersionIsSupported() {
        DataSupport supported = checkDataVersionSupported();
        return supported == DataSupport.SUPPORTED || supported == DataSupport.UNKNOWN;
    }


    /**
     * Determine if we should try to parse input data
     *
     * @return true if data release is supported or unknown
     */
    public boolean dataReleaseIsSupported() {
        DataSupport supported = checkDataReleaseSupported();
        return supported == DataSupport.SUPPORTED || supported == DataSupport.UNKNOWN;
    }

    /**
     * Determine if data version is supported
     *
     * @return Whether the data version is supported (within known range), unsupported (earlier than known range), or unknown (later than known range)
     */
    public DataSupport checkDataVersionSupported() {
        float earliest = getEarliestSupportedDataVersion();
        float latest = getLatestSupportedDataVersion();
        float dataVersion = 0;
        log.info("Verifying that supplied data version number is supported.  Supported versions are " +
                earliest + " - " + latest + ".  Version supplied is " + properties.getVersion());
        try {
            dataVersion = Float.parseFloat(properties.getVersion());
        } catch (NumberFormatException e) {
            log.error("data version is not an number", e);
            log.fatal("The version of the data supplied on the command line is not a number.  Version supplied=" + properties.getVersion());
            properties.printUsageAndExit();
        }

        if (dataVersion < earliest) {
            log.error("The version of the data supplied is not supported.  This loader supports data versions " +
                    earliest + " - " + latest);
            return DataSupport.NOT_SUPPORTED;
        } else if (dataVersion > latest) {
            log.warn("The version of the data supplied is newer than the latest version of data tested with this loader.");
            log.warn("This loader supports data versions " + earliest + " - " + latest);
            log.warn("The loader will attempt to parse the data.");
            return DataSupport.UNKNOWN;
        }
        return DataSupport.SUPPORTED;
    }


    public Date getReleaseDate() {
        SimpleDateFormat dateFormat;
        String dateSupplied = properties.getReleaseDate();
        Date dataDate = null;

        try {
            dateFormat = new SimpleDateFormat("yyyy/MM/dd");
            dataDate = dateFormat.parse(dateSupplied);
        } catch (ParseException e) {
            try {
                dateFormat = new SimpleDateFormat("yyyy/MM");
                dataDate = dateFormat.parse(dateSupplied);
            } catch (ParseException e1) {
                log.error("Release date supplied (" + dateSupplied + ") is not a valid date. Please supply a release date in either yyyy/MM/dd or yyyy/MM format.");
                properties.printUsageAndExit();
            }
        }
        return dataDate;
    }


    /**
     * Determine if data release is supported
     *
     * @return Whether the data release is supported (within known range), unsupported (earlier than known range), or unknown (later than known range)
     */
    public DataSupport checkDataReleaseSupported() {
        Date earliest = getEarliestSupportedDataRelease();
        Date latest = getLatestSupportedDataRelease();
        if (earliest == null || latest == null) {
            return DataSupport.SUPPORTED;
        }
        String dateSupplied = properties.getReleaseDate();
        DateFormat dateFormat;
        Date dataDate = getReleaseDate();

        if(dataDate == null ){
                log.error("Release date supplied (" + dateSupplied + ") is not a valid date. Please supply a release date in either yyyy/MM/dd or yyyy/MM format.");
                properties.printUsageAndExit();
                return DataSupport.NOT_SUPPORTED;
        }

        if (earliest.compareTo(dataDate) > 0) {
            log.error("This data release is not supported. " +
                    "Supported dates are " + earliest + " - " + latest + ".  " +
                    "Supplied data release is " + dateSupplied);
            return DataSupport.NOT_SUPPORTED;
        } else if (latest.compareTo(dataDate) < 0) {
            log.warn("This data is newer than the supported range.  Loader will attempt to parse data.  " +
                    "Supported dates are " + earliest + " - " + latest + ".  " +
                    "Supplied data release is " + dateSupplied);
            return DataSupport.UNKNOWN;
        } else {
            return DataSupport.SUPPORTED;
        }
    }


    /**
     * @return Earliest supported data release date (day, month, and year)
     */
    public final Date getEarliestSupportedDataRelease() {
        return earliestSupportedDataRelease;
    }


    /**
     * @return Latest supported data release date (day, month, and year)
     */
    public final Date getLatestSupportedDataRelease() {
        return latestSupportedDataRelease;
    }

    /**
     * Set earliest data release the loader supports
     * Inheriting classes should call this method in their constructors
     * if they support data release checks
     *
     * @param date "yyyy/MM/dd"
     */
    protected void setEarliestSupportedDataRelease(String date) {
        DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd");
        try {
            this.earliestSupportedDataRelease = dateFormat.parse(date);
        } catch (ParseException e) {
            fail("Error parsing date " + date + " " + e.getMessage());
        }
    }


    /**
     * Set latest data release the loader supports
     * Inheriting classes should call this method in their constructors
     * if they support data release checks
     *
     * @param date "yyyy/MM/dd"
     */
    protected void setLatestSupportedDataRelease(String date) {
        DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd");
        try {
            this.latestSupportedDataRelease = dateFormat.parse(date);
        } catch (ParseException e) {
            fail("Error parsing date " + date + " " + e.getMessage());
        }
    }


    protected void setEarliestSupportedDataRelease(Date earliestSupportedDataRelease) {
        this.earliestSupportedDataRelease = earliestSupportedDataRelease;
    }

    protected void setLatestSupportedDataRelease(Date latestSupportedDataRelease) {
        this.latestSupportedDataRelease = latestSupportedDataRelease;
    }

    /**
     * Gives earliest suported data version
     * Inheriting classes should override this method if they support data version checks
     *
     * @return Earliest supported data version
     */
    public float getEarliestSupportedDataVersion() {
        return 0;
    }


    /**
     * Gives latest suported data version
     * Inheriting classes should override this method if they support data version checks
     *
     * @return Latest supported data version
     */
    public float getLatestSupportedDataVersion() {
        return 0;
    }


    /**
     * Inheriting classes should implement this method
     *
     * @return Loader version number
     */
    public String getLoaderVersionNumber() {
        return "";
    }


    /**
     * Inheriting classes should implement this method
     *
     * @return Loader build number
     */
    public String getLoaderBuildNumber() {
        return "";
    }

}

