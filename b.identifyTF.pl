#!/usr/bin/perl 

system("rm -rf TFfamily");
mkdir TFfamily;

my %classdb; 

while(my $file = glob "TMP_rules/*.tmp"){
	my $family = $file;
	$family =~ s/.tmp//g;
	$family =~ s/TMP_rules\///g;
	my $outfile = $family.".txt";
	open(my $out, ">TFfamily/$outfile") or die"";
	open(my $fh, $file) or die"";
	$/='>';
	<$fh>;
	while(<$fh>){
		chomp;
		my ($gene,$in_domain,$out_domain) = split(/\n/,$_);
		my @data_in  = split(/\s+/,$in_domain);
		my @data_out = split(/\s+/,$out_domain);
		my $judgement_A = & judgeA(@data_in);
		my $judgement_B = & judgeB(@data_out);
		if(($judgement_A eq "OK") and ($judgement_B eq "OK")){
			#print $out ">$gene\n$in_domain\n$out_domain\n";  #####################For debug###################
			print $out "$gene\n";
			$classdb{$family}++;
			}
		}
	close $fh;
	close $out;
	}

open(OUT, "> member_num.txt") or die"";
foreach $family(sort keys %classdb){
	print OUT "$family	$classdb{$family}\n";
	}
close OUT;


sub judgeA{
	my @data = @_;
	my $num_of_domain = @data - 2;
	my $num_of_yes    = 0;
	$flag = $data[1];
	foreach $i(2..$#data){
		$num_of_yes++ if($data[$i] =~ /yes/);
		}
  return "OK" if(($flag eq "and") and ($num_of_yes == $num_of_domain));
  return "OK" if(($flag eq "or") and ($num_of_yes != 0));
  return "OK" if(($flag eq "only") and ($num_of_yes != 0));
	}

sub judgeB{
	my @data = @_;
	my $num_of_domain = @data - 2;
	my $num_of_yes    = 0;	
	return "OK" if($num_of_yes == 0);
	}