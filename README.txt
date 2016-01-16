Usage:

1. Prerequisites
   a. install hmmsearch
   b. run pfam_scan to identify pfam domain
2. Install plantTF_identifier
   git clone https://github.com/tangerzhang/plantTF_identifier.git
   mv plantTF_identifier /to/your/working/path/
3. run program
   perl TFinPlant.pl -i pfam.out -p protein.fasta
4. read results
   TFfamily/ directory contain members of each TF family
   member_num.txt is general statistics for each TF family
   gene.domain.txt contains domain signatures for each gene
