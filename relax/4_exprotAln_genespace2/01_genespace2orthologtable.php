<?php
$sGeneSpace = "/data2/projects/yshao/relax/2_genespace2/rundir/results/gffWithOgs.txt.gz";
$sGeneSymbolMap = "/data2/projects/yshao/relax/3_assigngenesymbols_orthofinder2/orthofinder_genesymbols.tsv";
$sRNAID2SpGeneIDMap = "/data2/projects/yshao/relax/4_exprotAln_genespace2/preProt/allproteins.fa";
$arrTaxTranslate = array('Rhinogobiusformosanus' => 'Rhinogobius_formosanus',
'Gobiopsismacrostoma' => 'Gobiopsis_macrostoma',
'Glossogobiusgiuris' => 'Glossogobius_giuris',
'Tridentigerbifasciatus' => 'Tridentiger_bifasciatus',
'Tridentigerbarbatus' => 'Tridentiger_barbatus',
'Tridentigerradiatus' => 'Tridentiger_radiatus',
'Boleophthalmuspectinirostris' => 'Boleophthalmus_pectinirostris',
'Neogobiusmelanostomus' => 'Neogobius_melanostomus',
'Periophthalmusmodestus' => 'Periophthalmus_modestus',
'Rhinogobiussimilis' => 'Rhinogobius_similis',
'Siphamiatubifer' => 'Siphamia_tubifer',
'Taenioidessp' => 'Taenioides_sp',
'Odontamblyopusrebecca' => 'Odontamblyopus_rebecca',
'Oxyeleotrismarmorata' => 'Oxyeleotris_marmorata',
'Periophthalmusmagnuspinnatus' => 'Periophthalmus_magnuspinnatus'
);
$arrTaxTranslateBack = array_flip($arrTaxTranslate);
$sGeneSymbolIndexSp = "Tridentigerbarbatus";

$sOut = "genespace.orthogroups.txt";

$h = popen("zcat $sGeneSpace", 'r');
$hO = fopen($sOut, 'w');

$arrRNAID2SpGeneID = fnLoadSpGeneIdMap($sRNAID2SpGeneIDMap);

$arrFam = array();
$arrSpecies = array();

while(false !== ($sLn = fgets($h)) ) {
						$sLn = trim($sLn);
						if ($sLn == "") continue;
						
					 $arrF = explode("\t", $sLn);
					 if ($arrF[0] == 'genome') continue;
					 $sSpecies = $arrF[0];
					 $sOrthoFam = $arrF[13];
					 $bIsRepresent = ($arrF[11]=='TRUE');

					 $sGeneID = $arrF[1];

					 $sSpecies = (array_key_exists($sSpecies, $arrTaxTranslate) )?  $arrTaxTranslate[$sSpecies] : $sSpecies;
					 if (!$bIsRepresent) {
					 		continue;
					 }

					 $arrSpecies[$sSpecies] = true;
					 if (!array_key_exists($sOrthoFam, $arrFam) ) {
					 		$arrFam[$sOrthoFam] = array();
					 }

					 if (!array_key_exists($sSpecies, $arrFam[$sOrthoFam])) {
					 		$arrFam[$sOrthoFam][$sSpecies] = $sGeneID;
					 } else {
					 	  echo("Warning: duplicated gene in orthogroup $sOrthoFam in $sSpecies: $sGeneID ".$arrFam[$sOrthoFam][$sSpecies]."\n");
					 }
} 

$arrSpecies = array_keys($arrSpecies);

$arrID2Symbol = fnLoadGeneSymbolMap($sGeneSymbolMap , $sGeneSymbolIndexSp);
//print_r($arrID2Symbol );
//die();

fwrite($hO, "OrthoID\tGeneSymbol\tGeneName\t");
fwrite($hO, implode("\t", $arrSpecies)."\n" );

$sGeneSymbolIndexSp = (array_key_exists($sGeneSymbolIndexSp, $arrTaxTranslate) )?  $arrTaxTranslate[$sGeneSymbolIndexSp] : $sGeneSymbolIndexSp;

foreach($arrFam as $sOrthoFam => $arrSpp) {
					$sOrthoFam = "Group_".$sOrthoFam."_1";
					 if (!array_key_exists($sGeneSymbolIndexSp, $arrSpp)) {
					 		continue; //if the focal species is missing then ignore
					 }
					 if (!array_key_exists($arrSpp[$sGeneSymbolIndexSp],$arrID2Symbol )) {
					 		fwrite($hO, "$sOrthoFam\tunknown\tunknown\t");
					 } else {
					 	 fwrite($hO, "$sOrthoFam\t".implode("\t", $arrID2Symbol[$arrSpp[$sGeneSymbolIndexSp]] ) . "\t");
					 }
					 $arrSpIDs = array();
					 foreach($arrSpecies as $sSp) {
					 				 $sSpGeneId  = $sRNAID = array_key_exists($sSp, $arrSpp)? str_replace('.','_', $arrSpp[$sSp]):'NA';	
									 if (array_key_exists($sSp, $arrRNAID2SpGeneID ) &&  array_key_exists($sRNAID, $arrRNAID2SpGeneID[$sSp] ) ) {
									 		$sSpGeneId = $sSp. "|". $arrRNAID2SpGeneID[$sSp][$sRNAID] . "|".$sRNAID;
									 }
									  
					 				 $arrSpIDs[] = $sSpGeneId;					 
					 }

					 fwrite($hO, implode("\t", $arrSpIDs)."\n" );
}

function fnLoadGeneSymbolMap($sGeneSymbolMap , $sGeneSymbolIndexSp) {
				 $h = fopen($sGeneSymbolMap, 'r');
				 $nIdx = -1;
				 $arrRet = array();
				 while(false !== ($sLn = fgets($h)) ) {
						$sLn = trim($sLn);
						if ($sLn == "") continue;
						$arrF = explode("\t", $sLn);
						if ($arrF[0] =='HOG') {
							 $arrFlip = array_flip($arrF);
							 if (!array_key_exists($sGeneSymbolIndexSp , $arrFlip)) {
							 		die("Index species $sGeneSymbolIndexSp not found in $sGeneSymbolMap\n");
							 }
							 $nIdx = $arrFlip[$sGeneSymbolIndexSp];
							 continue;
						}

						$sSymbol = $arrF[3];
						$sGeneName = $arrF[4];
						$arrIDs = explode(',', $arrF[$nIdx]);
						foreach($arrIDs as $sID) {
									 $arrRet[trim($sID)] = array($sSymbol, $sGeneName);
						}
				}

				return $arrRet;
}

function fnLoadSpGeneIdMap($sRNAID2SpGeneIDMap) {
				 $h = popen("grep '>' $sRNAID2SpGeneIDMap", 'r');
				 $arrRet = array();
				 while(false !== ($sLn = fgets($h)) ) {
						$sLn = trim($sLn);
						if ($sLn == "") continue;

						$sLn = substr($sLn, 1);
						list($sSp, $sSpGeneId, $sRNAID) = explode("|", $sLn);
						if (!array_key_exists($sSp, $arrRet)) {
							 $arrRet[$sSp] = array();
						}

						$arrRet[$sSp][$sRNAID] = $sSpGeneId;
				 }

				 return $arrRet;
}
?>
