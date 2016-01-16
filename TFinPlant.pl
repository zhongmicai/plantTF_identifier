#!/usr/bin/perl 

use Getopt::Std;
getopts "i:p:";



if ((!defined $opt_i)|| (!defined $opt_p) ) {
    die "************************************************************************
    Usage: perl TFinPlant.pl -i pfam.out -p protein.fasta
      -h : help and usage.
      -i : pfam_scan output
      -p : protein fasta file 
************************************************************************\n";
}else{
  print "************************************************************************\n";
  print "Version 1.1\n";
  print "Copyright to Tanger\n";
  print "RUNNING...\n";
  print "************************************************************************\n";
        
        }

my $all_script;
while(my $file = glob "*.pl"){
	$all_script .= $file."-";
	}

if(($all_script =~ /a.hmmsearch.pl/) and ($all_script =~ /b.identifyTF.pl/)){
	print "Programming checking: looks good!\n";
}else{
	die "Programming checking: not complete!\nPlease make sure to copy all the scripts in your working path\n";
	}

my $protein     = $opt_p;
my $pfam_out    = $opt_i;    
my %familydb;   ###for storing rules 
my %infordb;    ###reading domain signatures for each gene

###0. reading rules for each TF family
while(<DATA>){
	next if(/#/);
	my @data = split(/,/,$_);
	$family     = $data[0];
	$domain_in  = $data[1];
	$domain_in  =~ s/\s+//g;
	$domain_out = $data[2];
	$domain_out =~ s/\s+//g;
	$familydb{$family}->{"IN"}  = $domain_in;
	$familydb{$family}->{"OUT"} = $domain_out;
	}
close DATA;
print "Reading rules finished\n";
###1. identify domain signatures for each gene, including pfam domain (from Pfam_scan) and hmm domain (from hmmsearch)

open(IN, $pfam_out) or die"";
while(<IN>){
	chomp;
	next if(/#/);
	@data    = split(/\s+/,$_);
	my $num = @data;
	next if($num == 0);
	my $gene    = $data[0];
	my $pfam_id = $data[5];
	$pfam_id =~ s/\..*//g;
	$infordb{$gene} .= $pfam_id."	"; 
	}
close IN;

open(IN, "perl a.hmmsearch.pl $protein|") or die"";
while(<IN>){
	chomp;
	my($gene, $lable) = split(/\s+/,$_,2);
	$infordb{$gene} .= $lable;
	}
close IN;

my $num_of_gene = keys %infordb;
open(LB, "> gene.domain.txt") or die"";
foreach $gene(sort keys %infordb){
	print LB "$gene:	$infordb{$gene}\n";
	}
close LB;

print "Finished scaning $num_of_gene genes in $protein file\n";

###2. generate rules for each family:
system("rm -rf TMP_rules");
mkdir TMP_rules;
foreach my $family(sort keys %familydb){
	my $tmp_out = $family.".tmp";
	open(my $fh, ">TMP_rules/$tmp_out") or die"$!\n";
	print $fh "##########################$family########################\n";
	my @in_domain;
	my @out_domain;
  my $flag; 	
	if($familydb{$family}->{'IN'} =~ /\//){
		$flag = "or";
		@in_domain = split(/\//,$familydb{$family}->{'IN'});
	}elsif($familydb{$family}->{'IN'} =~ /\+/){
		$flag = "and";
		@in_domain = split(/\+/,$familydb{$family}->{'IN'});
	}else{
		$flag = "only";
		$in_domain[0] = $familydb{$family}->{'IN'};
		}
	
	if($familydb{$family}->{'OUT'} =~ /\//){
		@out_domain = split(/\//,$familydb{$family}->{'OUT'});
	}elsif($familydb{$family}->{'OUT'} =~ /\+/){
		@out_domain = split(/\+/,$familydb{$family}->{'OUT'});
	}else{
		$out_domain[0] = $familydb{$family}->{'OUT'};
		}
###</For debug only, check rules
#  print ">$family\n";
#  print "in_domain	";
#  foreach my $in(@in_domain){
#  	print "$in	";
#  	}
#  print "\nout_domain	";
#  foreach my $out(@out_domain){
#  	print "$out	";
#  	}
#  print "\n";
###For debug only/>

	foreach my $gene(sort keys %infordb){
		my %tmp_pfam = ();
		@pfamdb = split(/\s+/,$infordb{$gene});
		foreach $pfam(@pfamdb){
			$tmp_pfam{$pfam} = 0;
			}
		
		print $fh ">$gene\nin_domain:	$flag	";
		###check in_domain
    foreach my $in(@in_domain){
	  	print $fh "$in,yes	" if(exists($tmp_pfam{$in}));
	  	print $fh "$in,no	" if(!exists($tmp_pfam{$in}));
	  	}
		###check out_domain
		print $fh "\nout_domain:	.	";
		foreach my $out(@out_domain){
	  	print $fh "$out,yes	" if(exists($tmp_pfam{$in}));
	  	print $fh "$out,no	" if(!exists($tmp_pfam{$in}));			
			}
		print $fh "\n";
		
		}
	close $fh;
	}

print "Generating rules for each gene\n";

###3. identify genes that meet the rules

print "Identifying transcription factor gene families\n";

`perl b.identifyTF.pl`;

system("rm -rf TMP_rules");

print "Pipeline finished\n";

__DATA__
#domain,should possess,should NOT possess
ABI,PF02362,PF00847/PF06507
Alfin-like,Alfin-like,PF00628/PF00046/PF02373/PF02375
AP2-EREBP,PF00847,
ARF,PF06507,
ARID,PF01388,
AUX_IAA,PF02309,PF06507/PF02362
BBR-BPC,PF06217,
BES,PF05687,
bHLH,PF00010,
BSD,PF03909,
bZIP,PF00170/PF07716/PF03131,PF00010
C2C2-CO-like,PF06203+PF00643,PF00320/PF04640
C2C2-Dof,PF02701,
C2C2-GATA,PF00320,
C2H2,PF00096,PF02373/PF02375/PF00628
C3H,PF00642,PF00628/PF00176/PF00096
CAMTA,PF03859/PF00612,
CCAAT,PF02045/PF00808/CCAAT-Dr1/HF-YB/NF-YC,
CPP,PF03638,
CSD,PF00313,
DBP,DNC+PF00481,
E2F-DP,PF02319,
EIL,PF04873,
FAR1,PF03101,
FHA,PF00498,
G2-like,G2-like,PF00072
GeBP,PF04504,
GRAS,PF03514,
GRF,PF08880+PF08879,
HB,PF00046/PF03790/PF03791,
HRT,HRT,
HSF,PF00447,
LFY,PF01698,
LIM,PF00412,
LOB,PF03195,
MADS,PF00319,
mTERF,PF02536,
MYB,PF00249,PF01388/PF00072/PF00176 
NAC,PF02365,
NOZZLE,NOZZLE_Angio,
OFP,PF04844,
PBF-2-like,PF08536,
PLATZ,PF04640,
RWP-RK,PF02042,
S1Fa-like,PF04689,
SAP,STER_AP ,
SBP,PF03110,
Sigma70-like,PF04542/PF04539/PF04545,
SRS,PF05142,
TAZ,PF02135,PF00628
TCP,PF03634,
Tify,PF06200,PF00320
TIG,PF01833,
Trihelix,trihelix,
TUB,PF01167,
ULT,ULT,
VARL,VARL,
VOZ,VOZ,
WRKY,PF03106,
ZF-HD,PF04770,
Zn-clus,PF00172,
Coactivator_p15,PF02229,
DDT,PF02791,Alfin-like/PF00046
GNAT,PF00583,PF00628
HMG,PF00505,PF01388/PF04690
IWS1,PF08711,
JUMONJI,PF02373/PF02375,
LUG,LUFS,
MBF,PF08523,
MED6,PF04934,
MED7,PF05983,
PHD,PF00628,PF02791/PF00046/PF02373/PF02375
Pseudo_ARR-B,PF06203+PF00072,
RB,PF01857,
Rcd1-like,PF04078,
SET,PF00856,PF03638/PF00628/PF00096
SNF2,PF00176,PF00847/PF00628
SOH1,PF05669,
SWI_SNF-BAF60b,PF02201,
SWI_SNF-SWI3,PF04433,PF00249
TRAF,PF00651,PF07707/PF00917/PF03000/PF02135