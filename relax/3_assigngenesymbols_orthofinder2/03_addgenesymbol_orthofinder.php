<?php
//assign gene symbol
$sOrthoGroup = "/data2/projects/zma/Rhinogobius_spp/DNA/hyphy/assigngenesymbols_orthofinder/orthofinder_genesymbols.tsv";

$sGFF = "Rform.gff3";
$sOut = "Rform.genesymbol.gff3";
/*
$sGFF = "/data2/projects/zma/Rhinogobius_spp/DNA/hyphy/ref/Rform.softmask_longest_isoform.gff3";
$sOut = "Rform.softmasked.longest_isoform.genesymbol.gff3";
 */
$sThisSp = 'Rform';
$sPredictor = "maker"; //replace column 2 with this value

$hOut = fopen($sOut , 'w');

$arrGeneSymbolCounts = array();

$arrGene2RNAMap = array();
$hGFF = fopen($sGFF, 'r');

$arrAvailableRNAIDs = array();
$arrRNAIDsWithSymbols = array();
$arrRNAIDsWithSymbols2 = array();

	while( false !== ($sLn = fgets($hGFF ) ) ) {
		if ($sLn == '#') continue;
		$sLn = trim($sLn);
		if ($sLn == '') continue;

		$arrF = explode("\t" , $sLn);
		

		if (count($arrF) != 9) continue;
		
		if ($arrF[2] == 'mRNA') {
			$arrAnnot = fnParseFields($arrF[8]);
			if (!array_key_exists("ID" , $arrAnnot) ) continue;
			if (!array_key_exists("Parent" , $arrAnnot) ) continue;
			if (!array_key_exists($arrAnnot['Parent'], $arrGene2RNAMap)) {
				$arrGene2RNAMap[$arrAnnot['Parent'] ] = array();
			}
			$arrGene2RNAMap[$arrAnnot['Parent'] ][$arrAnnot['ID']] = true;
			$sRNAID =trim($arrAnnot['ID']);
			$arrAvailableRNAIDs[$sRNAID] = true;

			continue;
		}

	}
fclose($hGFF);

$arrID2Symbol = fnLoadOrthoGroups($sOrthoGroup);

//print_r($arrID2Symbol);


$hGFF = fopen($sGFF, 'r');

	while( false !== ($sLn = fgets($hGFF ) ) ) {
		if ($sLn == '#') {
			fwrite( $hOut , $sLn);
			continue;
		}

		$sLn = trim($sLn);
		if ($sLn == '') continue;

		$arrF = explode("\t" , $sLn);
		$arrF[1] = $sPredictor; //change predictor

		if (count($arrF) != 9) {
			fwrite( $hOut , $sLn."\n");
			continue;
		}
		
		if ($arrF[2] == 'gene' || $arrF[2] == 'mRNA') {
			$arrAnnot = fnParseFields($arrF[8]);
			if (!array_key_exists("ID" , $arrAnnot) ) continue;
			if ($arrF[2] == 'gene') {
				if (!array_key_exists($arrAnnot['ID'] , $arrGene2RNAMap) ) {
					fwrite( $hOut , implode("\t", $arrF)."\n");
					continue;
				}
			}
			
			$sRNAID = "";
			if ($arrF[2] == 'gene') {
				foreach($arrGene2RNAMap[$arrAnnot['ID']] as $sIsoformID => $b) {
					if (array_key_exists($sIsoformID, $arrID2Symbol)) {
						$sRNAID = $sIsoformID;
						break;
					}
				}
				
			} else {
				$sRNAID = $arrAnnot['ID'];
			}
			//$sRNAID = ($arrF[2] == 'gene')? $arrGene2RNAMap[$arrAnnot['ID']] : $arrAnnot['ID'];
			//$sRNAID =str_replace(".", "_", $sRNAID);



			if (!array_key_exists($sRNAID , $arrID2Symbol) ) {
				$arrAnnot['gene'] = 'unknown';
				$arrAnnot['description'] = 'unknown';
			} else {
				$arrAnnot['gene'] = $arrID2Symbol[$sRNAID]['genesymbol'];
				$arrAnnot['description'] = $arrID2Symbol[$sRNAID]['description'];
				//$arrAnnot['ensemblorthologs'] = $arrID2Symbol[$sRNAID]['ensembl_orthologs'];possible_gene_symbols
				$arrAnnot['possible_gene_symbols'] = $arrID2Symbol[$sRNAID]['possible_gene_symbols'];
				$arrAnnot['possible_gene_desc'] = $arrID2Symbol[$sRNAID]['possible_gene_desc'];

			}
			

			$arrF[8] = fnArr2Annot($arrAnnot);
			fwrite( $hOut , implode("\t", $arrF)."\n");
			continue;
		} else {
			fwrite( $hOut , implode("\t", $arrF)."\n");
			continue;
		}



	}
fclose($hGFF);



function fnLoadOrthoGroups($sOrthoGroup) {
	global $sThisSp,  $arrNum2ID , $arrGeneSymbolCounts, $arrAvailableRNAIDs, $arrRNAIDsWithSymbols, $arrRNAIDsWithSymbols2;

	$hF = fopen($sOrthoGroup , 'r');
	$arrRet = array(); //key is gene id
	$arrGeneSymbolCounts = array();
	$nSpCol = -1;
	$nAllGeneSymbolCol = -1;
	$nAllDescCol = -1;
	while( false !== ($sLn = fgets($hF) )) {
		$sLn = trim($sLn, "\n\r ");
		if ($sLn == '') continue;
		
		$arrF = explode("\t", $sLn);
		if ($arrF[0] == 'HOG') {
			$arrFFlip = array_flip($arrF);
			if (!array_key_exists($sThisSp, $arrFFlip)) {
				die("Error: species $sThisSp not defined in header of $sOrthoGroup\n");
			}
			
			$nSpCol = $arrFFlip[$sThisSp];
			$nAllGeneSymbolCol = $arrFFlip['AllGeneSymbols'];
			$nAllDescCol = $arrFFlip['AllDescriptions'];
			continue;
		}
		
		list($sGroupID, $sOGID, $sCladeID, $sGeneSymbol, $sDescription) = $arrF; 
		$arrGenes = explode(',' , trim($arrF[$nSpCol]) );

		$arrEnsemblOrthologs = array();
		$arrAddedGeneIDs = array();
		foreach( $arrGenes as $sGene) {

			$sGeneID = trim($sGene);

			//$sGeneID =str_replace(".", "_", $sGeneID);

			$arrAddedGeneIDs[] = $sGeneID;
			$arrRet[$sGeneID] = array();
			$arrRet[$sGeneID]['genesymbol'] = $sGeneSymbol;
			$arrRet[$sGeneID]['description'] = $sDescription;
			if (!array_key_exists($sGeneSymbol , $arrGeneSymbolCounts)) {
				$arrGeneSymbolCounts[$sGeneSymbol] = 0;
			}
			if (array_key_exists($sGeneID , $arrAvailableRNAIDs) && (!array_key_exists($sGeneID, $arrRNAIDsWithSymbols) )) {
				$arrGeneSymbolCounts[$sGeneSymbol] += 1;
				$arrRNAIDsWithSymbols[$sGeneID] = true;
			}
			
		}


		foreach($arrAddedGeneIDs as $sGeneID) {
			$arrRet[$sGeneID]['possible_gene_symbols'] = trim($arrF[$nAllGeneSymbolCol]);
			$arrRet[$sGeneID]['possible_gene_desc'] = trim($arrF[$nAllDescCol]);
		}
	}

	$arrGeneSymbolCounter = array();
	foreach($arrRet as $sGeneID => &$arrInfo) {
		if ($arrGeneSymbolCounts[$arrInfo['genesymbol']] > 1) { //
			if (!array_key_exists($arrInfo['genesymbol'] , $arrGeneSymbolCounter) ) {
				$arrGeneSymbolCounter[$arrInfo['genesymbol']] = 0;
			}
			if (array_key_exists($sGeneID , $arrAvailableRNAIDs) && (!array_key_exists($sGeneID, $arrRNAIDsWithSymbols2))) {
				$arrGeneSymbolCounter[$arrInfo['genesymbol']]  += 1;
				$arrRNAIDsWithSymbols2[$sGeneID] = true;
			}
			$arrInfo['genesymbol'] = $arrInfo['genesymbol']." (".$arrGeneSymbolCounter[$arrInfo['genesymbol']]." of ".$arrGeneSymbolCounts[$arrInfo['genesymbol']].")";

		}
	}

	return $arrRet;
}


function fnParseFields($s) {
	$arrF1 = explode(";" , $s);
	$arrRet = array();
	foreach($arrF1 as $sF) {
		$sF = trim($sF);
		if ($sF == '') continue;
		$arrF2 = explode("=" , $sF);
		if (count($arrF2)!=2) {
			echo("parse error: $sF\n");
		}
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
