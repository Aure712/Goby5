<?php
$sOrtho = "genespace.orthogroups.txt";
$sGFF = "/data/projects/rcui/macropodus_compare/improve_genemodels_UPhO/mop/mop.longest.addedUPhO.genesymbol.spgeneid.gff3";
$nCol = 3; //0-based column
$sOut = "groupid2coord.txt";

$h = fopen($sGFF, 'r');
$arrCoords = array();
while(false !== ($sLn = fgets($h))) {
	$sLn = trim($sLn);
	if ($sLn == '' || $sLn[0] =='#') continue;

	$arrF = explode("\t", $sLn);
	if ($arrF[2] != "mRNA") {
		continue;
	}

	preg_match('/ID=([^;]+)/', $arrF[8], $arrM);
	if (count($arrM)!=2) {
		continue;
	}

	$sID = trim($arrM[1]);
	$sID = str_replace('.', '_', $sID);
	$arrCoords[$sID] = array($arrF[0], $arrF[3]-1, $arrF[4]);
}

$hO = fopen($sOut, 'w');

$hOrth = fopen($sOrtho, 'r');

while(false !== ($sLn = fgets($hOrth))) {
	$sLn = trim($sLn);
	if ($sLn == '') continue;

	$arrF = explode("\t", $sLn);
	if ($arrF[0] == 'OrthoID') {
		continue;
	}
	list($s1, $s2, $sRNAID) = explode("|", $arrF[$nCol]);

	if (!array_key_exists($sRNAID, $arrCoords)) {
		die("$sRNAID not found in GFF3\n");
	}

	fwrite($hO, implode("\t", $arrCoords[$sRNAID])."\t".$arrF[0]."\n");
}
?>
