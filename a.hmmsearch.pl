#!/usr/bin/perl
####handling hmmsearch results, trying to convert these files to pfam-like signature
 
$pep_file = $ARGV[0];
while(my $file = glob "hmmprofile/*.hmm"){
	$lable = $file;
	$lable =~ s/.hmm//g;
	$lable =~ s/hmmprofile\///g;
#	system("rm hmm.txt");
	$cmd = "hmmsearch -E 1e-3 $file $pep_file > hmm.txt";
	system($cmd);
	
	open(IN, "hmm.txt") or die"";
	while(<IN>){
		chomp;
		next if(/#/);
		next if(/Query:/);
		next if(/Scores/);
		next if(/full sequence/);
		next if(/E-value/);
		next if(/-------/);
		last if(/annotation/);
		@data = split(/\s+/,$_);
		$gene = $data[9];
		$infordb{$gene} .= $lable."	";
		}
	close IN;
	
	}
system("rm hmm.txt");

#open(OUT, "> lable.out") or die"";
foreach $gene(sort keys %infordb){
	next if($gene eq "");
	print "$gene	$infordb{$gene}\n";
	}
#close OUT;
