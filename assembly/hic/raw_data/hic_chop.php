<?php
$sRd1Files = "/public3/group_crf/home/g21shaoy23/HiC/raw_data/Tra/63.32G/HIC924_R1.fq.gz /public3/group_crf/home/g21shaoy23/HiC/raw_data/Tra/7.03G/HIC924_R1.fq.gz"; //separate each file by a space, multiple files are allowed. for example: "test_1_1.fq.gz test_2_1.fq.gz"
$sRd2Files = "/public3/group_crf/home/g21shaoy23/HiC/raw_data/Tra/63.32G/HIC924_R2.fq.gz /public3/group_crf/home/g21shaoy23/HiC/raw_data/Tra/7.03G/HIC924_R2.fq.gz"; //separate each file by a space, example "test_1_2.fq.gz test_2_2.fq.gz"

$sO1 = "chopped_1.fq.gz"; // output file for chopped read1, gzip format
$sO2 = "chopped_2.fq.gz"; // output file for chopped read2, gzip format

$nMinOutputLen = 25; //do not output shorter than 25bp, both mates will be discarded.
$sCutSite = "GATC"; //the enzyme site for digestion.

$nCutSiteLen = strlen($sCutSite); //length of the enzyme site.



$h1 = popen("zcat -f $sRd1Files", 'r'); // read from the decompression stream of the zcat command, this command will concatenate multiple gzip files into one stream.
$h2 = popen("zcat -f $sRd2Files", 'r');

$hO1 = popen("gzip -c > $sO1", 'w'); // open the write stream, piping the output to the gzip command for compression. 
$hO2 = popen("gzip -c > $sO2", 'w');

$nLn = -1;            // current line number of the input stream.
$arrLns1 = array();   // buffer for file 1
$arrLns2 = array();   // buffer for file 2
$nCutLen1 = 0;        // computed cut length of read 1
$nCutLen2 = 0;        // computed cut length of read 2
$nKeptPairs = 0;      // counter for kept pairs
$nDiscardPairs = 0;   // counter for discarded pairs.

while(true) {
	$nLn++;										// increase line number counter, line 1 = 0
	$sLn1 = fgets($h1);							// get the next line from the read 1 stream
	$sLn2 = fgets($h2);							// get the next line from the read 2 stream
	
	if ($sLn1 === false || $sLn2 === false) { // if end of file, sLn1 and sLn2 become false
		if ($sLn1 !== $sLn2) {
			echo("Warning: two read files do not reach the end at the same point\n"); //warn if the two input streams don't reach the end at the same time.
		}
		break;
	}
	
	$sLn1 = trim($sLn1); // remove trailing "\n"
	$sLn2 = trim($sLn2);
	
	if ($nLn % 4 == 0) { // title line
		//check title lines:
		if ($sLn1[0] != '@' || $sLn2[0] != '@') { //check if title lines start with @
			die("Error: unexpected title lines\n");			
		}
		list($sRd1Name, $tmp) = explode('/', $sLn1); // remove the trailing "/1, /2" from read names
		list($sRd2Name, $tmp) = explode('/', $sLn2);
		
		if ($sRd1Name != $sRd2Name) { // check if two read names are identical 
			die("Error: on line $nLn, two read names differ. $sRd1Name, $sRd2Name\n"); // if two read names differ, errors out.
		}
		
		//clear buffers:
		$arrLns1 = array();
		$arrLns2 = array();
		
		//push lines to buffer;
		$arrLns1[] = $sLn1;
		$arrLns2[] = $sLn2; 
		continue; //next line;

	}
	
	if ($nLn % 4 == 1) { // dna sequence line
		$nEnzymeSite1 = strpos($sLn1, $sCutSite); // search for the enzyme site; if not found, function returns false;
		$nEnzymeSite2 = strpos($sLn2, $sCutSite); // search for the enzyme site; if not found, function returns false;
		
		$nCutLen1 = ($nEnzymeSite1 === false)? strlen($sLn1) : ($nEnzymeSite1 + $nCutSiteLen ); // if not found, set cut length to the length of the full DNA, otherwise, set it to the cut point.
		$nCutLen2 = ($nEnzymeSite2 === false)? strlen($sLn2) : ($nEnzymeSite2 + $nCutSiteLen );
		
		$arrLns1[] = substr($sLn1, 0, $nCutLen1); // perform the cutting, and pushes the cut line into the buffer.
		$arrLns2[] = substr($sLn2, 0, $nCutLen2); // perform the cutting, and pushes the cut line into the buffer.
		
		continue;
		
	}
	
	if ($nLn % 4 == 2) { // + line
		if ($sLn1[0] != '+' || $sLn2[0] != '+') {
			die("+ is not found on expected line.\n"); //check if this line starts with "+", if not, errors out.
			
		}
		
		$arrLns1[] = '+';
		$arrLns2[] = '+'; //push line to buffer;
		
		continue;
	}
	
	if ($nLn % 4 == 3) { // quality line
		if ($nCutLen1 < $nMinOutputLen || $nCutLen2 < $nMinOutputLen) { // too short after cutting, discard by not outputting.
			$nDiscardPairs++; //counter +1
			continue;			
		}
		
		//if retained lengths ok, then write to file
		
		$arrLns1[] = substr($sLn1, 0, $nCutLen1); //cut quality line, and push to buffer;
		$arrLns2[] = substr($sLn2, 0, $nCutLen2); //cut quality line, and push to buffer;
		
		//write to output streams
		fwrite($hO1, implode("\n", $arrLns1) . "\n" );
		fwrite($hO2, implode("\n", $arrLns2) . "\n" );
		$nKeptPairs++; // counter +1
	}
	
}

echo("Discarded $nDiscardPairs pairs, trimmed and kept $nKeptPairs\n");

?>
