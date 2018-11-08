#a mettre dans alignements (avec les autres socires)

import Bio
import sys
import subprocess
sys.path.append("bin/CCMpred/scripts/")
from top_couplings import get_top_pairs

#script allowing to convert fasta.hhr to aln

#script allowing to extract N residues which the higest probability to be in contact

from conkit.applications import CCMpredCommandline

def calculate_co_evolution_score(file_aln):
    """
        Calcule the coevolution score of the query based on the MSA alignment with a chosen database
        The co-evolution score represente the co-occurence of a pair of amino acid in a maximum of species. Two amino acid have co-evoluate if the occurence of one of this amino never occur whithout the other

        Args:
        file_aln obtained after conversion from file.fasta.hhr
        output_mat define name of the matrix output
        Returns:
            matrix containing all scores associated with each pair of amino acids in the sequence
    """
    ccmpred_cline = CCMpredCommandline(cmd ='bin/CCMpred/bin/ccmpred', alnfile= file_aln, matfile= "output.mat")
    ccmpred_cline()
	subprocess.call("./bin/CCMpred/scripts/top_couplings.py output.mat > tops_outputs.mat", shell=True)


        
calculate_co_evolution_score("bin/CCMpred/example/1atzA.aln")





