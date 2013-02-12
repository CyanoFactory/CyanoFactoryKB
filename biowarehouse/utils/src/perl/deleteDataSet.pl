#!/usr/bin/perl

use DBI;

$oracle = "oracle";
$mysql = "mysql";

$argCount = @ARGV;
$dbWIDSQLLine = "select name from dataset where";
$dbWID = 0;
@blank = { " " };
$appendSQL = "";
print $argCount;
if( $argCount < 4)
{
        die "usage: For mysql: perl deleteDataSetWID.pl mysql <datasetWid> <userid> <password> <host> <database_name>\n For oracle: perl deleteDataSetWID.pl oracle <datasetWid> <userid> <password>\n"
}
else
{
        $db = shift(@ARGV);
	$dbWID = shift( @ARGV );
        $userid = shift( @ARGV );
        $password = shift( @ARGV );
}

print "db is $db\n";
    
if($db eq $oracle)
{
  $dbspec = $ENV{ORACLE_SID};
  $datasrc = "Oracle:host=localhost;sid=$dbspec";
}
else{
        if($db eq $mysql){
                $host = shift( @ARGV );
                $database = shift( @ARGV );
                $datasrc = "mysql:database=$database:$host";
        }
        else
        {
        die "usage: For mysql: perl deleteDataSetWID.pl mysql <datasetWid> <userid> <password> <host> <database_name>\n For oracle: perl deleteDataSetWID.pl oracle <datasetWid> <userid> <password>\n"
        }
}

 #$dbh = DBI->connect( "dbi:$datasrc", $userid, $password) or die DBI->errstr;

 $dbh = DBI->connect( qq{dbi:$datasrc}, $userid, $password, { PrintError => 1, RaiseError => 1, AutoCommit => 0} ) or die $DBI::errstr;

    $dbWIDSQLLine2 = join(" ", $dbWIDSQLLine, "wid =", $dbWID);
    print "line: $dbWIDSQLLine2 \n";
    $sel = $dbh->prepare( $dbWIDSQLLine2 );
    $sel->execute();
    if( @row = $sel->fetchrow_array ) {
	$dbany = @row[0];
	print "WID Found: $dbWID\n";
    }
    else {
	die "No Data Sets match parameters";
    }

    $sel->finish;

    #list of table names
    @lookUpTables = ("DBID", "SynonymTable", "CommentTable", "CrossReference","Support", "ToolAdvice" );

    $sqlLine = join(" ", "select OtherWID from Entry where datasetWID =", $dbWID);
    print "$sqlLine \n";
    $selA = $dbh->prepare( $sqlLine );
    $selA->execute(); 

    # get all the wid's for this object
    while( @row = $selA->fetchrow_array ) {
	$var = @row[0];
	
	### take care of lookup tables
	foreach $lookup ( @lookUpTables ) {
	    $sqlLine = join( " ", "Delete From ", $lookup, " where otherwid = ", $var );
	    $selB = $dbh->prepare( $sqlLine );
	    print "Execute: $sqlLine \n";
	    $selB->execute();
	}
        $dbh->commit;
        print "commit.\n";
    }

    # First delete  Computation, so that Sequence Match can be deleted with it
    # trying to delete everything by cascaded deletes causes Oracle to run out of rollback segment
    $sqlLine = join(" ", "Delete From Computation where datasetWID =", $dbWID);
    print "$sqlLine \n";
    $selA = $dbh->prepare( $sqlLine );
    $selA->execute();
    $dbh->commit;
    print "commit.\n"; 

    $sqlLine = join(" ", "Delete From Dataset where wid = ", $dbWID );
    $selB = $dbh->prepare( $sqlLine );
    print "Execute: $sqlLine \n";
    $selB->execute();

    $selB->finish;

    $dbh->commit;
    print "Commited.\n";

    $dbh->disconnect;
    print "Disconnected.\n";

    exit;
