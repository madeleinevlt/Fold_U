#a mettre dans alignements (avec les autres socires)

import Bio

#script allowing to convert fasta.hhr to aln
import convert_alignment
#script allowing to extract N residues which the higest probability to be in contact
import top_couplings
from conkit.applications import CCMpredCommandline


def calculate_co_evolution_score(file_aln,output_mat):
    """
        Calcule the coevolution score of the query based on the MSA alignment with a chosen database
        The co-evolution score represente the co-occurence of a pair of amino acid in a maximum of species. Two amino acid have co-evoluate if the occurence of one of this amino never occur whithout the other

        Args:
		file_aln obtained after conversion from file.fasta.hhr
		output_mat define name of the matrix output
        Returns:
        	matrix containing all scores associated with each pair of amino acids in the sequence
    """
    ccmpred_cline = CCMpredCommandline(cmd ='/bin/ccmpred', alnfile= file_aln, matfile= output_mat)
    return ccmpred_cline()
	


calculate_co_evolution_score('1atza.aln','1atza.mat')






