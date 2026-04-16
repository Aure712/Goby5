<?php

$sIn = "/data2/projects/yshao/Goby5Compare/synteny2/goby2sp_tbi/results/TridentigerBifasciatus_pangenomeDB.txt.gz";
$sOut = "synorthos2sp_tbi.txt";

$arrGenomes = fnGetGenomes($sIn); 
$hIn = popen("zcat $sIn" , 'r');
$hO = fopen($sOut, 'w');

$arrOrths = array();
while(false !== ($sLn = fgets($hIn)) ) {
	$sLn = trim($sLn);
	if ($sLn == '') continue;
	$arrF = explode("\t", $sLn);
	if ($arrF[0] == 'pgChr') {
		continue;
	}

	$nPGID = $arrF[2];
	if (!array_key_exists($nPGID , $arrOrths)) {
		$arrOrths[$nPGID] = array('pgChr' => $arrF[0], 'pgOrd' => $arrF[1], 'directSynOrth' => array() , 'NSOrths' => array() );
		foreach($arrGenomes as $sSample) {
			$arrOrths[$nPGID]['directSynOrth'][$sSample] = "NA";
			$arrOrths[$nPGID]['NSOrths'][$sSample] = array();
		}
	}

	$sGenome = $arrF[6];
	$bIsDirectSyn = ($arrF[7]=='TRUE');
	if ($bIsDirectSyn) {
		$arrOrths[$nPGID]['directSynOrth'][$sGenome] = $arrF[11];
	}

	$bIsNSOrtho = ($arrF[10] == 'TRUE');
	if ($bIsNSOrtho) {
		$arrOrths[$nPGID]['NSOrths'][$sGenome][$arrF[11]] = true ;
	}
}

fwrite($hO, "pgChr\tpgOrd\tpgID\t".implode("\t",$arrGenomes)."\t".implode("\t",$arrGenomes)."\n");
foreach($arrOrths as $nPGID => $arrOr) {
	fwrite($hO, $arrOr['pgChr']."\t".$arrOr['pgOrd']."\t$nPGID\t".implode("\t",$arrOr['directSynOrth'])."\t");
	$arrNSOrths = array();
	foreach($arrOr['NSOrths'] as $arrNSOr ) {
		if (count($arrNSOr)==0) {
			$arrNSOrths[] = 'NA';
		} else {
			$arrNSOrths[] = implode(';', array_keys($arrNSOr) );
		}
	}
	fwrite($hO, implode("\t",$arrNSOrths )."\n");
}


function fnGetGenomes($sIn) {

	$h = popen("zcat $sIn | cut -f7 | grep -v genome | sort | uniq", 'r');
	$arrRet = array();
	while(false !== ($sLn = fgets($h)) ) {
		$sLn = trim($sLn);
		if ($sLn == '') continue;
		$arrRet[] = $sLn;
	}

	pclose($h);
	return $arrRet;
}
?>
