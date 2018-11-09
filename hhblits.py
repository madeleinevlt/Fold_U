import Bio
import subprocess
from conkit.applications import HHblitsCommandline

#run hhblits (input = query fasta, database, output = a3m)
hhblits_cline = HHblitsCommandline(
	cmd = "/bin/hh-suite/build/bin/hhblits", input ="Q92630.fasta", database ="uniprot20_29Feb2012",
	cpu = "3", matflix = "100000", oa3m = "Q92630.a3m", show_all
)
hhblits_cline()

#reformate a3m into fasta
reformat = subprocess.Popen(["./bin/hh-suite/scripts/reformat.pl","-r","a3m","fas","test.a3m","test.fasta"], stdout=subprocess.PIPE).communicate()[0] 
 
#format fasta into aln
convert_alignment = subprocess.Popen(["./bin/CCMpred/scripts/convert_alignment.py","test.fasta","fasta","test.aln"], stdout=subprocess.PIPE).communicate()[0]  

#extract top coupling
subprocess.call("./bin/CCMpred/scripts/top_couplings.py output.mat > tops_outputs.mat", shell=True) 



