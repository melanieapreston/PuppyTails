use strict;

my @files = glob "*.txt";

my $key = "SampleKey.out";
open KEY, "$key", or die "Can't open $key\n";

my %names;

while (my $line = <KEY>) {
	chomp $line;
	my @split = split("\t", $line);
	$names{$split[0]}[0] = $split[0];
	$names{$split[0]}[1] = $split[1];
}

my @position;
my $countelement = 0;
my @firstline;
#my $sample;

my $out3 = "UnweightedTailCompSummaryForAllSamples.out";
chomp $out3;
open OUT3, ">$out3", or die "Can't open $out3\n";
print OUT3 "Sample ID\tSum Total unique tails\tUnweighted Average A incorporation\tUnweighted Average G incorporation\tUnweighted Average C incorporation\tUnweighted Average T incorporation\tAverage Tail Length\tPercent tailed\n";

my $out4 = "WeightedTailCompSummaryForAllSamples.out";
chomp $out4;
open OUT4, ">$out4", or die "Can't open $out4\n";
print OUT4 "Sample ID\tSum Total unique tails\tTotal number of nucleotides added\tWeighted Average A incorporation\tWeighted Average G incorporation\tWeighted Average C incorporation\tWeighted Average T incorporation\tAverage Tail length\tPercent tailed\n";



foreach my $file(@files) {
	chomp $file;
	open IN, "$file", or die "can't open $file\n";
	my $prefix = $file;
	$prefix =~ s/\.txt//;
	print "\n\nNow processing the file: $file\n\n";
	$countelement=0;

	
	while (my $line = <IN>) {
		$position[$countelement] = $line; 
		
		if ($countelement eq 8) {
			my $sample = $position[0];
			$sample =~ /MP(.+?)_/;
			$sample = "MP".$1; 
			#print "$position[0]\n";
			my $row2count = 0;
			my $row3count = 0;
			my $row4count = 0;
			my $row5count = 0;
			my $A =0;
			my $G =0;
			my $C =0;
			my $T =0;
			my $bases = 0;
			my $averagetail = 0;
			my $tailsum=0;
				
			print "\nFound the sample ID: $sample\n";
				if (exists $names{$sample}){
					print "Sample ID: $sample matched the key: $names{$sample}[0]\n";
					@firstline = split("\t", $position[0]);
					$firstline[0] = "$names{$sample}[0]-$names{$sample}[1]";
					$position[0] = join("\t", @firstline);
				}
			
			my $out = "$firstline[0].out";
			open OUT, ">$out", or die "can't open $out\n";
			print OUT $position[0].$position[1].$position[2].$position[3].$position[4].$position[6].$position[7];
		
			my @row1 = split ("\t", $position[0]);
			my @row2 = split ("\t", $position[1]);
			my @row3 = split ("\t", $position[2]);
			my @row4 = split ("\t", $position[3]);
			my @row5 = split ("\t", $position[4]);
			my @row6 = split ("\t", $position[5]);
			my @row7 = split ("\t", $position[6]);
			my @row8 = split ("\t", $position[7]);
			
			my $row1sum = 0;
			my $row2sum = 0;
			my $row3sum = 0;
			my $row4sum = 0;
			my $row5sum = 0;
			my $row6sum = 0;
			my $row7sum = 0;
			my $row8sum = 0;
		
		#$bases = 0;
			for (my $i=6; $i<36; $i++) {
				$row1sum = $row1sum + $row1[$i];
				
				
				if ($row1[$i] ne 0) {
					$row2sum = $row2sum + $row2[$i];
					$row2count++;
					$row3sum = $row3sum + $row3[$i];
					$row3count++;
					$row4sum = $row4sum + $row4[$i];
					$row4count++;
					$row5sum = $row5sum + $row5[$i];
					$row5count++;
					
					$A = $A + ($row1[$i]*$row2[$i]*$row7[$i]);
					$G = $G + ($row1[$i]*$row3[$i]*$row7[$i]);
					$C = $C + ($row1[$i]*$row4[$i]*$row7[$i]);
					$T = $T + ($row1[$i]*$row5[$i]*$row7[$i]);
					$bases = $bases + ($row1[$i]*$row7[$i]);
					
				}
			}
			print "\n#####$bases\n";
			my $row2ave = 0;
			my $row3ave = 0;
			my $row4ave = 0;
			my $row5ave = 0;
			
			my $percentA = 0;
			my $percentG = 0;
			my $percentC = 0;
			my $percentT = 0;
			my $tailed = 0;
		
			
			if ($row2count ne 0) {$row2ave = $row2sum / $row2count*100;} else {$row2ave = "NA";}
			if ($row3count ne 0) {$row3ave = $row3sum / $row3count*100;} else {$row3ave = "NA";}
			if ($row4count ne 0) {$row4ave = $row4sum / $row4count*100;} else {$row4ave = "NA";}
			if ($row5count ne 0) {$row5ave = $row5sum / $row5count*100;} else {$row5ave = "NA";}
			
			if ($bases ne 0) {
				$percentA = $A / $bases * 100;
				$percentG = $G / $bases * 100;
				$percentC = $C / $bases * 100;
				$percentT = $T / $bases * 100;
				$averagetail = $bases/$row1sum;
			}
			else {
				$percentA = "NA";
				$percentG = "NA";
				$percentC = "NA";
				$percentT = "NA";
				$averagetail = "NA";
			}
			
			if ($row1[1] ne 0) {$tailed = $row1sum/($row1[1] + $row1sum)*100;} else {$tailed = "NA";}
			
			my $out2 = "$firstline[0]-averages.out";
			open OUT2, ">$out2", or die "can't open $out2\n";
			print OUT2 "Sum Total unique tails\t$row1sum\n";
			print OUT2 "Unweighted Average A incorporation\t$row2ave\n";
			print OUT2 "Unweighted Average G incorporation\t$row3ave\n";
			print OUT2 "Unweighted Average C incorporation\t$row4ave\n";
			print OUT2 "Unweighted Average T incorporation\t$row5ave\n";
			print OUT2 "Weighted Average A incorporation\t$percentA\n";
			print OUT2 "Weighted Average G incorporation\t$percentG\n";
			print OUT2 "Weighted Average C incorporation\t$percentC\n";
			print OUT2 "Weighted Average T incorporation\t$percentT\n";
			print OUT2 "Total A nucleotides added\t$A\n";
			print OUT2 "Total G nucleotides added\t$G\n";
			print OUT2 "Total C nucleotides added\t$C\n";
			print OUT2 "Total T nucleotides added\t$T\n";
			print OUT2 "Total number of nucleotides added\t$bases\n";
			
			print OUT3 "$firstline[0]\t$row1sum\t$row2ave\t$row3ave\t$row4ave\t$row5ave\t$averagetail\t$tailed\n";
			print OUT4 "$firstline[0]\t$row1sum\t$bases\t$percentA\t$percentG\t$percentC\t$percentT\t$averagetail\t$tailed\n";
	
			undef @position;
			undef @firstline;
			$countelement = -1;
			#undef $sample;
		}
		
		$countelement++;
	}	
}



exit;

