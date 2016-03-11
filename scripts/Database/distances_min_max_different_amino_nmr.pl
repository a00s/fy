#!/usr/bin/perl
##########################################################################
# Create min and max data for results of NMR of different protein
##########################################################################
use DBI;
$my_db ="a00s_230";
$my_user = "a00s_230";
$my_pass = "testando";
$my_host="127.0.0.1";
$pdbid = @ARGV[0];
chomp $pdbid;
if($pdbid eq ""){
    print "Use \n./distances_min_max_different_amino_nmr.pl pdbid\n";
    exit;
}

$dbh = DBI->connect("DBI:mysql:$my_db:$my_host", $my_user, $my_pass, {AutoCommit => 1});
$query = "SELECT amino1, label1, amino2, label2, count(*), MIN(distance) FROM (SELECT atomo1.i_306337 amino1, atomo1.i_331770 label1, atomo2.i_306337 amino2, atomo2.i_331770 label2, ROUND((SQRT(POW((atomo1.i_306299 - atomo2.i_306299),2) + POW((atomo1.i_306307 - atomo2.i_306307),2) + POW((atomo1.i_306315 - atomo2.i_306315),2))),2) AS distance FROM a_306280 atomo1 LEFT JOIN a_306280 atomo2 ON atomo1.i_307676=atomo2.i_307676 WHERE atomo1.i_307676='$pdbid' AND atomo1.i_306393=1 AND atomo1.i_331770 IS NOT NULL AND atomo2.i_331770 IS NOT NULL AND atomo1.id<>atomo2.id AND atomo1.i_306323=atomo2.i_306323 AND atomo1.i_306401 = atomo2.i_306401 AND atomo1.i_306344 <> atomo2.i_306344) tabtemp GROUP BY amino1, amino2, label1, label2";
$sqlQuery  = $dbh->prepare($query);
$sqlQuery->execute;
while (@row= $sqlQuery->fetchrow_array()) {
    print "@row[0] @row[1] @row[2] @row[3] @row[4] @row[5]\n";
    my $sth = $dbh->prepare('INSERT INTO m_380480 SET id=NULL');
    $sth->execute();
    $last_id = $dbh->{'mysql_insertid'};
    $sth = $dbh->prepare("INSERT INTO a_380484 SET id=$last_id, i_380488='@row[0]'/*aminofrom*/, i_380494='@row[1]'/*atomfrom*/, i_380500='@row[2]'/*aminoto*/,i_380506='@row[3]'/*atomto*/,i_380512=NULL/*sameamino*/,i_380517='@row[5]'/*mindistance*/,i_380529='@row[4]', i_393242='$pdbid'");
    $sth->execute();	    
}
