<?php
$sIn="../repeatmodeler/Tba-families.fa";
$sOut = "Tba-te-families.fa";
$sSpecies = "Tridentiger_barbatus";

$h = fopen($sIn, 'r');
$hO = fopen($sOut, 'w');

while(false !== ($sLn = fgets($h)) ) {
	if ($sLn[0] != '>') {
		fwrite($hO, $sLn);
		continue;
	}

	list($s1, $s2)  =  explode(' ' , $sLn);

	fwrite($hO, "$s1 @$sSpecies\n");
}
?>
