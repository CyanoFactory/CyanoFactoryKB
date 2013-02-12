#!/usr/bin/perl

use DBI;

$argCount = @ARGV;
if( $argCount < 3 || $argCount > 5) {
	die "usage: perl $0 dbvendor userid password [db] [host]
       dbvendor = oracle | mysql
       db = MySQL database or Oracle SID; defaults to \$ORACLE_SID
       host = IP address or name of database server machine; defaults to localhost\n"
}

# Process command line args
$host=localhost;
$dbvendor = shift( @ARGV );
$userid = shift( @ARGV );
$password = shift( @ARGV );
if ( $argCount >= 4 ) {
    $dbspec = shift( @ARGV );
}
else {
    $dbspec = $ENV{ORACLE_SID};
}
if ( $argCount >= 5 ) {
    $host = shift( @ARGV );
}
if ( $dbvendor eq "oracle" ) {
	$datasrc = "Oracle:host=$host;sid=$dbspec";
}
elsif ( $dbvendor eq "mysql" ) {
	$datasrc = "mysql:database=$dbspec:$host";
}
else { die "Invalid dbvendor: $dbvendor\n"; }

# Connect to the database
print "Connecting $userid to $datasrc\n";
$dbh = DBI->connect("dbi:$datasrc",
		    #- qq{dbi:Oracle:$db},
		    $userid,
		    $password,
		    { PrintError => 1, RaiseError => 1, AutoCommit => 0 } )
	or die $DBI::errstr;

# Query for all datasets
print "Looking up information for all datasets in $userid...\n\n";
$sel = $dbh->prepare( "select * from DataSet" );
$sel->execute();
$selB = $dbh->prepare( "select count(OtherWID) from Entry where DataSetwid = ?" );

local $^W = 0;

@names = @{$sel->{NAME}};

while( @row = $sel->fetchrow_array ) {
    $tempName = @names[0];
    $tempValue = @row[0];
    print "$tempName: $tempValue \n";
    $tempName = @names[1];
    $tempValue = @row[1];
    print "$tempName: $tempValue \n";
    $tempName = @names[2];
    $tempValue = @row[2];
    print "$tempName: $tempValue \n";
    $tempName = @names[3];
    $tempValue = @row[3];
    print "$tempName: $tempValue \n";
    $tempName = @names[4];
    $tempValue = @row[4];
    print "$tempName: $tempValue \n";
    $tempName = @names[5];
    $tempValue = @row[5];
    print "$tempName: $tempValue \n";
    $tempName = @names[6];
    $tempValue = @row[6];
    print "$tempName: $tempValue \n";

    $WID = @row[0];
    $selB->execute($WID);
    if ( @row2 = $selB->fetchrow_array ) {
	$tempValue = @row2[0];
	print "Number of Entries: $tempValue \n";
    }
    print "\n\n";
}

$sel->finish;
$selB->finish;

$dbh->disconnect;
exit;
