<?php
$sIn = "merged.gff";
$sOut = "merged.fixed.gff";

$h = fopen($sIn, 'r');
$hO = fopen($sOut, 'w');

$sPrefix = "om";
$arrLocusMap = array();
$arrGeneMap = array();

while(false !== ($sLn = fgets($h)) ) {
  $sLn = trim($sLn);
  if ($sLn == '') continue;
  if ($sLn[0] == '#') {
    fwrite($hO, $sLn."\n");
    continue;
  }
  
  $arrF = explode("\t", $sLn);
  $arrAnn = fnParseFields($arrF[8]);
  if ($arrF[2] == 'locus') {
    continue;
  }
  
  if ($arrF[2] == 'gene') {
    $sLocus = $arrAnn['locus'];
    if (array_key_exists($sLocus, $arrLocusMap)) {
      $arrGeneMap[$arrAnn['ID']] = $sLocus;
      continue;
    }
    
    $arrLocusMap[$sLocus] = $arrAnn['ID'];
    $arrGeneMap[$arrAnn['ID']] = $sLocus;
    fwrite($hO, $sLn."\n");
    continue;
  }
  
  if ($arrF[2] == 'mRNA') {
    $sParent = $arrAnn['Parent'];
    $sLocus = $arrGeneMap[$sParent];
    $sNewParent = $arrLocusMap[$sLocus];
    $arrAnn['Parent'] = $sNewParent;
    $arrF[8] = fnArr2Annot($arrAnn);
    fwrite($hO, implode("\t",$arrF)."\n");
    continue;
  }
  
  fwrite($hO, $sLn."\n");
  
}

function fnParseFields($s) {
        $arrF1 = explode(";" , $s);
        $arrRet = array();
        foreach($arrF1 as $sF) {
                $arrF2 = explode("=" , $sF);
                $arrRet[trim($arrF2[0]) ] = trim($arrF2[1]);
        }

        return $arrRet;
}

function fnArr2Annot($arr) {
        $s = "";

        foreach($arr as $sKey => $sVal) {
                if ($sVal === false) $sVal=0;
                $s .= "$sKey=$sVal;";
        }

        return substr($s, 0, strlen($s)-1);
}


?>
