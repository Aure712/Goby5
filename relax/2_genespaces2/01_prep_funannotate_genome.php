<?php

$sRunDir = "rundir";

/*
$sScfPrefix = "";
$sSp = "Rhinogobius_formosanus";
$sVersion = "1";
$sLongestIsoformProt = "/data2/projects/zma/Rhinogobius_spp/DNA/hyphy/ref/Rform_longest_isoform.prot.fa";

$sScfPrefix = "";
$sSp = "Gobiopsis_macrostoma";
$sVersion = "1";
$sLongestIsoformProt = "/data2/projects/yshao/Goby5RNA/annotations/longest_iso/Gma_longest.prot.fa";

$sScfPrefix = ""; 
$sSp = "Glossogobius_giuris";
$sVersion = "1";
$sLongestIsoformProt = "/data2/projects/yshao/Goby5RNA/annotations/longest_iso/Ggi_longest.prot.fa";

$sScfPrefix = ""; 
$sSp = "Tridentiger_bifasciatus";
$sVersion = "1";
$sLongestIsoformProt = "/data2/projects/yshao/Goby5RNA/annotations/longest_iso/Tbi_longest.prot.fa";

$sScfPrefix = ""; 
$sSp = "Tridentiger_barbatus";
$sVersion = "1";
$sLongestIsoformProt = "/data2/projects/yshao/Goby5RNA/annotations/longest_iso/Tba_longest.prot.fa";

$sScfPrefix = ""; 
$sSp = "Tridentiger_radiatus";
$sVersion = "1";
$sLongestIsoformProt = "/data2/projects/yshao/Goby5RNA/annotations/longest_iso/Tra_longest.prot.fa";

$sScfPrefix = "";
$sSp = "Boleophthalmus_pectinirostris";
$sVersion = "1";
$sLongestIsoformProt = "/data2/projects/yshao/relax/1_sp/Boleophthalmus_pectinirostris/longest_isoform.prot.fa";

$sScfPrefix = "";
$sSp = "Neogobius_melanostomus";
$sVersion = "1";
$sLongestIsoformProt = "/data2/projects/yshao/relax/1_sp/Neogobius_melanostomus/longest_isoform.prot.fa";

$sScfPrefix = "";
$sSp = "Periophthalmus_modestus";
$sVersion = "1";
$sLongestIsoformProt = "/data2/projects/yshao/relax/1_sp/Periophthalmus_modestus/longest_isoform.prot.fa";

$sScfPrefix = "";
$sSp = "Rhinogobius_similis";
$sVersion = "1";
$sLongestIsoformProt = "/data2/projects/yshao/relax/1_sp/Rhinogobius_similis/longest_isoform.prot.fa";

$sScfPrefix = "";
$sSp = "Siphamia_tubifer";
$sVersion = "1";
$sLongestIsoformProt = "/data2/projects/yshao/relax/1_sp/Siphamia_tubifer/longest_isoform.prot.fa";

$sScfPrefix = "";
$sSp = "Taenioides_sp";
$sVersion = "1";
$sLongestIsoformProt = "/data2/projects/yshao/relax/1_sp/Taenioides_sp/longest_isoform.prot.fa";

$sScfPrefix = "";
$sSp = "Oxyeleotris_marmorata";
$sVersion = "1";
$sLongestIsoformProt = "/data2/projects/yshao/relax/1_sp/Oxyeleotris_marmorata/longest_isoform.prot.fa";

$sScfPrefix = "";
$sSp = "Odontamblyopus_rebecca";
$sVersion = "1";
$sLongestIsoformProt = "/data2/projects/yshao/relax/1_sp/Odontamblyopus_rebecca/longest_isoform.prot.fa";

$sScfPrefix = "";
$sSp = "Periophthalmus_magnuspinnatus";
$sVersion = "1";
$sLongestIsoformProt = "/data2/projects/yshao/relax/1_sp/Periophthalmus_magnuspinnatus/longest_isoform.prot.fa";

*/

$sScfPrefix = "";
$sSp = "Periophthalmus_magnuspinnatus";
$sVersion = "1";
$sLongestIsoformProt = "/data2/projects/yshao/relax/1_sp/Periophthalmus_magnuspinnatus/longest_isoform.prot.fa";


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
