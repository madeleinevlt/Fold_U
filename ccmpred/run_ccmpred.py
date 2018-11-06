#a mettre dans alignements (avec les autres socires)

import Bio
import conkit
#script allowing to convert fasta.hhr to aln
import convert_alignment
#script allowing to extract N residues which the higest probability to be in contact
import top_coupling
from conkit.applications import CCMpredCommandline


def calculate_co_evolution_score(file.aln):
    """
        Calcule the coevolution score of the query based on the MSA alignment with a chosen database
        The co-evolution score represente the co-occurence of a pair of amino acid in a maximum of species. Two amino acid have co-evoluate if the occurence of one of this amino never occur whithout the other

        Args:
            file.aln
            database
        Returns:
            matrix
    """
    ccmpred_cline = CCMpredCommandline(
        alnfile="file.aln", matfile="output.mat"
        )
    return output.mat


