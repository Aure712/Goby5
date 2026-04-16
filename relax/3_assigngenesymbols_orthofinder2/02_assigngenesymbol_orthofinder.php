<?php
$sOrthoGroups = "/data2/projects/yshao/relax/2_genespace2/rundir/orthofinder/Results_Aug19/Phylogenetic_Hierarchical_Orthogroups/N0.tsv";
$sOut = "orthofinder_genesymbols.tsv";

$arrGeneSymbolMap = glob("*.id2genesymbol.map.txt");
$arrUnknownKeyWord = array( "/unknown/");

$arrGeneSymbolMaps = array();
$arrGeneDescriptionMaps = array();

$hOut = fopen($sOut , "w");

foreach($arrGeneSymbolMap as $sMap) {
	$arrM = explode('.', $sMap);
	$sSp = $arrM[0];
	list( $arrGeneSymbolMaps[$sSp] , $arrGeneDescriptionMaps[$sSp]) = fnParseMap( $sMap ); 
}

//print_r($arrGeneSymbolMaps);

$hOrthoGroups = fopen($sOrthoGroups , 'r');
$arrSpCol = array();
$arrRefSpp = array();
while( false !== ($sLn = fgets($hOrthoGroups) )) {
	$sLn = trim($sLn, "\n\r ");
	if ($sLn == '') continue;
	$arrF = explode("\t"  , $sLn);
	
	if ($arrF[0] == 'HOG') {
		//header line
		$arrSpCol = array_slice(array_flip($arrF) , 3 );
		$arrGeneSymbolMaps = array_intersect_key($arrGeneSymbolMaps, $arrSpCol );
		$arrGeneDescriptionMaps = array_intersect_key($arrGeneDescriptionMaps, $arrSpCol );
		$arrRefSpp = array_keys($arrGeneSymbolMaps);
		//print_r($arrSpCol);
		//die();
		fwrite($hOut, implode("\t", array_slice($arrF, 0, 3) )."\tGeneSymbol\tDescription\t" );
		fwrite($hOut , implode("\t", array_slice($arrF, 3))."\tAllGeneSymbols\tAllDescriptions\n" );

		continue;
	}
	
	$sGroupID = $arrF[0];
	$arrSymbols = array();
	$arrDescriptions = array();
	$arrDescriptions['unknown'] = "unknown";
	foreach($arrRefSpp as $sSp ) {
		$arrIDs = explode(",", $arrF[$arrSpCol[$sSp]] );
		//print_r($arrIDs);
		foreach($arrIDs as $sID) {
			$sID = trim($sID);
			list($sID) = explode('.', $sID);
			if (!array_key_exists($sID , $arrGeneSymbolMaps[$sSp]) ) { //if doesn't exist
				continue;
			}

			$sSymbol = $arrGeneSymbolMaps[$sSp][$sID];			
			$sSymbol = trim(preg_replace("/\(\d+ of .+/", "", $sSymbol));
			$sDes = $arrGeneDescriptionMaps[$sSp][$sID];

			foreach($arrUnknownKeyWord as $sUnknownKeyWord) {
				if (preg_match($sUnknownKeyWord , $sSymbol) ===1 ) {
					continue 2;
				}
			}
			$sSymbol = strtolower($sSymbol);
			if (!array_key_exists($sSymbol , $arrSymbols)) {
				$arrSymbols[$sSymbol] = 0;
			}
			$arrSymbols[$sSymbol] += 1;
				
			if (!array_key_exists($sSymbol,  $arrDescriptions)) {
				$arrDescriptions[$sSymbol] = $sDes;
			}

			$arrDescriptions[$sSymbol] = (strlen($arrDescriptions[$sSymbol]) < strlen($sDes) && strlen($arrDescriptions[$sSymbol]) > 10 )? $arrDescriptions[$sSymbol] : $sDes;

		
		}
	}
	

	arsort( $arrSymbols , SORT_NUMERIC);
	$arrSymbolsFlip = array_keys($arrSymbols);

	if (count($arrSymbolsFlip) > 1) {
		echo("warning: more than one possible gene symbols found for group $sGroupID : " . implode(', ' , $arrSymbolsFlip ) . "\n");
		
	}
	
	$sFinalSymbol = 'unknown';
	if (count($arrSymbolsFlip) > 0) {
		foreach($arrSymbolsFlip as $sSym) {
			if (strpos($sSym, ':')!==false) continue;
			if (substr($sSym, 0,3) == 'loc') continue;
			$sFinalSymbol = $sSym; 
		}

		if ($sFinalSymbol == 'unknown') {
			$sFinalSymbol = $arrSymbolsFlip[0];
		}
	}

	fwrite($hOut, implode("\t", array_slice($arrF, 0, 3) )."\t".$sFinalSymbol."\t".$arrDescriptions[$sFinalSymbol]."\t" );
	fwrite($hOut , implode("\t", array_slice($arrF, 3))."\t".implode('|', $arrSymbolsFlip)."\t".implode('|', array_unique(array_slice($arrDescriptions,1) ))."\n" );
}


function fnParseMap( $sMap ) {
	$hMap = fopen($sMap, 'r');
	$arrRet = array(array(), array());
	while( false !== ($sLn = fgets($hMap) )) {
		$sLn = trim($sLn);
		if ($sLn == '') continue;

		$arrF = explode("\t" , $sLn);
		$arrRet[0][$arrF[0]] = $arrF[1];
		$arrRet[1][$arrF[0]] = preg_replace("/%2C transcript variant.*/" , '', $arrF[2]);
	}

	return $arrRet;
}

?>
