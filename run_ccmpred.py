#! coding: utf-8
# Meet-U

#Script to compute fonction for getting from an query seq an top couplings co-ev scores

#Librairies imported
import Bio
import subprocess

from conkit.applications import CCMpredCommandline
from conkit.applications import HHblitsCommandline



def calculate_co_evolution_score(query,database_ref):
 	"""
        Calculates co-evolution score of the query based on the MSA alignment
        Co-evolution score measures co-occurence of a pair of amino acid in
	ortholog sequences. Two amino acid have co-evoluated if the occurence
	of one of this amino never occur whithout the other.

        Args:
             query:file in .fasta
             database_ref: database for hhblits
        Returns:
            matrix containing scores associated with top co-evoluated pairs
            of amino acids in the sequence
    """
    # conkit command hhblits produces Multiple Sequence Alignment from fasta
        #filter the sequences at a 90% identity threshold and E-value cutoff for
        # inclusion in result alignment= 0.001 are included in the tool
    

    #Getting multiple Alignment from a query sequence (FASTA format)
    hhblits_cline = HHblitsCommandline(cmd = "/bin/hh-suite/build/bin/hhblits", 
    	input =query,
	 	database = database_ref ,cpu = "3", matflix = "100000",
	 	oa3m = "query.a3m", show_all)
    hhblits_cline()

    #reformate a3m into fasta
    reformat = subprocess.Popen(
        ["./bin/hh-suite/scripts/reformat.pl","-r","a3m",
        "fas","query.a3m","query.fasta"], stdout=subprocess.PIPE).communicate()[0]  #why comunicate & others juste need to run the process and gives a file !

    
    #format fasta into aln
    convert_alignment = subprocess.Popen(
        ["./bin/CCMpred/scripts/convert_alignment.py","query.fasta","fasta",
        "query.aln"], stdout=subprocess.PIPE).communicate()[0]					#same !

    #CCMpred
    ccmpred_cline = CCMpredCommandline(
        cmd ='bin/CCMpred/bin/ccmpred', alnfile= "query.aln", matfile= "query.mat"
    )
    ccmpred_cline()
    #extract top coupling
    subprocess.call(
    "./bin/CCMpred/scripts/top_couplings.py query.mat > tops_outputs.mat",
     shell=True
     )
