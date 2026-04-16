<?php
$sRunDir = "rundir2";

/*
$sScfPrefix = "GmaScf_"; //to be removed from scaffold names, leaving only numbers
$sSp = "Gobiopsis_macrostoma";
$sVersion = "1.0";
$sLongestIsoformProt = "/data2/projects/yshao/Goby5RNA/annotations/longest_iso/Gma_longest.M.fa";

$sScfPrefix = "GgiScf_"; 
$sSp = "Glossogobius_giuris";
$sVersion = "1.0";
$sLongestIsoformProt = "/data2/projects/yshao/Goby5RNA/annotations/longest_iso/Ggi_longest.M.fa";

$sScfPrefix = "TbiScf_"; 
$sSp = "Tridentiger_bifasciatus";
$sVersion = "1.0";
$sLongestIsoformProt = "/data2/projects/yshao/Goby5RNA/annotations/longest_iso/Tbi_longest.M.fa";

$sScfPrefix = "TbaScf_"; 
$sSp = "Tridentiger_barbatus";
$sVersion = "1.0";
$sLongestIsoformProt = "/data2/projects/yshao/Goby5RNA/annotations/longest_iso/Tba_longest.M.fa";

$sScfPrefix = "TraScf_"; 
$sSp = "Tridentiger_radiatus";
$sVersion = "1.0";
$sLongestIsoformProt = "/data2/projects/yshao/Goby5RNA/annotations/longest_iso/Tra_longest.M.fa";

$sSp = "Human";
$sVersion = "107";
$sLongestIsoformProt = "/data/projects/rcui/bahaha_assembly/annotations/other_refs/Homo_sapiens/longest_isoform.prot.titlecleaned.fa";

$sSp = "Zebrafish";
$sVersion = "107";
$sLongestIsoformProt = "/data/projects/rcui/bahaha_assembly/annotations/other_refs/Danio_rerio/longest_isoform.prot.titlecleaned.fa";
 */

$sScfPrefix = "GmaScf_"; //to be removed from scaffold names, leaving only numbers
$sSp = "Gobiopsis_macrostoma";
$sVersion = "1.0";
$sLongestIsoformProt = "/data2/projects/yshao/Goby5RNA/annotations/longest_iso/Gma_longest.M.fa";

$sWD = "$sRunDir/rawGenomes/$sSp/$sVersion/annotation";
exec("mkdir -p $sWD");

$hGFF = false;
if ($sScfPrefix !='') {
	$hGFF = popen("sort -k1,1n -k4,4n | gzip -c > $sWD/$sSp.gene.gff.gz", 'w');
} else {
	$hGFF = popen("sort -k1,1 -k4,4n | gzip -c > $sWD/$sSp.gene.gff.gz", 'w');
}

$hFas = popen("gzip -c > $sWD/$sSp.pep.fa.gz", 'w');

$hIn = popen("zcat -f $sLongestIsoformProt", 'r');

while(false!== ($sLn = fgets($hIn) )) {
	$sLn = trim($sLn);
	if ($sLn == '') {
		continue;
	}

	if ($sLn[0] == '>') {
		$arrF = explode(' ', substr($sLn, 1));
		$sSeqName = $arrF[0];
		$sSeqName = preg_replace('/[:|]/', '_', $sSeqName);
		$sPos = $arrF[count($arrF)-1];
		preg_match('/([^:]+):([^-]+)-([^(]+)\\((\\S)\\)/', $sPos, $arrM);
		if (count($arrM) != 5) {
			die("Failed to parse header: $sLn\n");
		}

		list($sChr, $nStart, $nEnd, $sOrient) = array_slice($arrM, 1);
		$sChr = str_replace($sScfPrefix, '', $sChr);
		fwrite($hGFF, "$sChr\tfunannotate\tgene\t$nStart\t$nEnd\t\t$sOrient\t\tlocus=$sSeqName\n");
		fwrite($hFas, ">$sSeqName\n");
	} else {
		$sLn = str_replace('*', '', $sLn);
		fwrite($hFas, $sLn."\n");
	}

	
}

?>
