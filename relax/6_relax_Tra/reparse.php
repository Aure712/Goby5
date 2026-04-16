<?php
$sIn = "sum_BTP.txt";
$sOut = "sum_BTP_reparse.txt";

$h = fopen($sIn ,'r');
$hO = fopen($sOut, 'w');
while(false !== ($sLn = fgets($h) ) ) {
	$sLn = trim($sLn);
	if ($sLn == '') continue;
	$arrF = explode("\t", $sLn);
	if ($arrF[1] != 'Success') {
		fwrite($hO, $sLn."\n");
		continue;
	}
	$sTmpDir = $arrF[count($arrF)-1];
	$sF1 = fnReadText(  $sTmpDir."/in.nex.RELAX.json" ); //read output.
	$arrJSONret = json_decode($sF1 , true);
	$arrRet = array();

        $arrRet["MG94xREV"] = $arrJSONret["fits"]["MG94xREV with separate rates for branch sets"]["Log Likelihood"]; //loglik of MG94XREV
        $arrRet["NULL"] = $arrJSONret["fits"]["RELAX null"]["Log Likelihood"]; //Loglik of null model assuming K=1 (no relax or intensification of selection)
        $arrRet["Alternative"] = $arrJSONret["fits"]["RELAX alternative"]["Log Likelihood"]; //loglik of alt model, k estimated.
        $arrRet["K"] = $arrJSONret["test results"]["relaxation or intensification parameter"]; //estimated k 
        $arrRet["LR"] = $arrJSONret["test results"]["LRT"]; //likelihood ratio
        $arrRet["P"] = $arrJSONret["test results"]["p-value"]; //p value with 1 degree of freedom. chi sq.
        $arrRet["R_omega0"] = $arrJSONret["fits"]["RELAX alternative"]["Rate Distributions"]["Reference"][0]["omega"]; //omega 0 (<1) of reference branches.
        $arrRet["R_omega0_prop"] = $arrJSONret["fits"]["RELAX alternative"]["Rate Distributions"]["Reference"][0]["proportion"]; // proportion of omega 0 (<1) of reference branches.
        $arrRet["R_omega1"] = $arrJSONret["fits"]["RELAX alternative"]["Rate Distributions"]["Reference"][1]["omega"]; //omega 1 (near 1) of reference branches.
        $arrRet["R_omega1_prop"] = $arrJSONret["fits"]["RELAX alternative"]["Rate Distributions"]["Reference"][1]["proportion"]; // proportion of omega 1 (near 1) of reference branches.
        $arrRet["R_omega2"] = $arrJSONret["fits"]["RELAX alternative"]["Rate Distributions"]["Reference"][2]["omega"]; //omega 2 (>1) of reference branches.
        $arrRet["R_omega2_prop"] = $arrJSONret["fits"]["RELAX alternative"]["Rate Distributions"]["Reference"][2]["proportion"]; // proportion of omega 2 (>1) of reference branches.
        $arrRet["T_omega0"] = $arrJSONret["fits"]["RELAX alternative"]["Rate Distributions"]["Test"][0]["omega"]; //omega 0 (<1) of Test branches.
        $arrRet["T_omega0_prop"] = $arrJSONret["fits"]["RELAX alternative"]["Rate Distributions"]["Test"][0]["proportion"]; // proportion of omega 0 (<1) of Test branches.
        $arrRet["T_omega1"] = $arrJSONret["fits"]["RELAX alternative"]["Rate Distributions"]["Test"][1]["omega"]; //omega 1 (near 1) of Test branches.
        $arrRet["T_omega1_prop"] = $arrJSONret["fits"]["RELAX alternative"]["Rate Distributions"]["Test"][1]["proportion"]; // proportion of omega 1 (near 1) of Test branches.
        $arrRet["T_omega2"] = $arrJSONret["fits"]["RELAX alternative"]["Rate Distributions"]["Test"][2]["omega"]; //omega 2 (>1) of Test branches.
        $arrRet["T_omega2_prop"] = $arrJSONret["fits"]["RELAX alternative"]["Rate Distributions"]["Test"][2]["proportion"]; // proportion of omega 2 (>1) of Test branches.                
                     
	$arrHyphyRet = $arrRet;
	$sGeneName = $arrF[0];
	$sGeneSymbol= $arrF[2];
	$nSeqLen = $arrF[3]; 
	$sGeneFullName = $arrF[22];
	
	fwrite($hO, "$sGeneName\tSuccess\t$sGeneSymbol\t".($nSeqLen)
				."\t".$arrHyphyRet["P"] 
				."\t".$arrHyphyRet["MG94xREV"] 
				."\t".$arrHyphyRet["NULL"] 
				."\t".$arrHyphyRet["Alternative"] 
				."\t".$arrHyphyRet["K"] 
				."\t".$arrHyphyRet["LR"] 
				."\t".$arrHyphyRet["R_omega0"] 
				."\t".$arrHyphyRet["R_omega0_prop"] 
				."\t".$arrHyphyRet["R_omega1"] 
				."\t".$arrHyphyRet["R_omega1_prop"] 
				."\t".$arrHyphyRet["R_omega2"] 
				."\t".$arrHyphyRet["R_omega2_prop"] 
				."\t".$arrHyphyRet["T_omega0"] 
				."\t".$arrHyphyRet["T_omega0_prop"] 
				."\t".$arrHyphyRet["T_omega1"] 
				."\t".$arrHyphyRet["T_omega1_prop"] 
				."\t".$arrHyphyRet["T_omega2"] 
				."\t".$arrHyphyRet["T_omega2_prop"] 
				."\t$sGeneFullName"
				."\t$sTmpDir"
			."\n");



}
	function fnReadText($file) {
	/*
		$fh = fopen($file, 'r');
		if (!$fh) {
			return false;
		}
		
		$theData = fgets($fh);
		fclose($fh);
		return $theData;
		*/
		if (file_exists($file)) {
			return file_get_contents($file);
		}
		else {
			return false;
		}
	}

?>
