#!usr/bin/perl -w 

$ArrSize = @ARGV;

$IndexColumns = 0;
if ($ArrSize > 1) {
	if ($ARGV[1] > 0) { $IndexColumns = $ARGV[1]; }
}

@Files = <$ARGV[0]_*.dat>;

$FileCounter = 0;
$Rows = 0;

foreach $File (@Files) {
	
	open FH, $File or die $!;
	
	if ($FileCounter == 0) {
		
		while(1) {
			$l = <FH>;
			last unless(defined $l);
			$l =~ s/,//g;
			@lSplit = split(' ', $l);
			$lSplitLength = @lSplit;
			for ($i = 0; $i < $IndexColumns; $i++) { push(@Index, $lSplit[$i]); }
			for ($i = $IndexColumns; $i < $lSplitLength; $i++) { push(@Data, $lSplit[$i]) }
			$Rows++;
		}
		
	}
	else {
		
		$DataPos = 0;
		
		while(1) {
			$l = <FH>;
			last unless(defined $l);
			$l =~ s/,//g;
			@lSplit = split(' ', $l);
			for ($i = $IndexColumns; $i < $lSplitLength; $i++) { $Data[$DataPos++] += $lSplit[$i] }
		}
		
	}
	
	close(FH);
	
	$FileCounter++;
}

$IndexPos = 0;
$DataPos = 0;

open FH, ">$ARGV[0].dat" or die $!;

for ($j = 0; $j < $Rows; $j++) {
	for ($i = 0; $i < $IndexColumns; $i++) { printf FH "%f ", $Index[$IndexPos++]; }
	for ($i = $IndexColumns; $i < $lSplitLength; $i++) { printf FH "%f ", $Data[$DataPos++] / $FileCounter }
	printf FH "\n";
}

close(FH)