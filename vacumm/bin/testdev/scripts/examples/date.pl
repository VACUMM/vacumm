#!/usr/bin/perl -w

use Date::Calc (Add_Delta_Days);

$arg= @ARGV;
if ( $arg == 0 ) {
$nb=0;
}
else {
$nb= $ARGV[0];
$ou= $ARGV[1];
}
($j,$m,$a)=(localtime)[3,4,5];
$m++;
$a +=1900;

$m=8;
$a=2008;
$j=14;

($aa,$mm,$jj)= Add_Delta_Days($a,$m,$j,$nb);

if ($mm<10) {$mm = "0" . $mm;};
if ($jj<10) {$jj = "0" . $jj;};  

if ( $ou == 1 )  {
print "$jj\n";
}
elsif ( $ou == 2 ) {
print "$mm\n";
}
elsif ( $ou == 3 ) {
print "$aa\n";
}
