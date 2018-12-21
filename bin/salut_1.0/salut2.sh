#!/bin/bash

#####Création des répertoires necessaires#####
# rm -Rf pssm_query pssm_homstrad tmp_query tmp_homstrad alignments
# mkdir Query pssm_query tmp_query alignments



# #####Déplacement du fasta dans Query####
# cp "$1" ./Query


# #####Traitements sur les queries#####
rm -f log/query.log
rm -f log/querypgp.log


fastaquery=`echo "$1"`
query=`echo ${fastaquery##*/} | cut -d. -f1`
echo "QUERY : "$query

# echo $query >> log/query.log

# echo "Psi-blast contre $3"
# #./blast-legacy-2.2.26-./bin/salut_1.0/src/blastpgp -i Query/$query.fasta -d db/"$3".fasta -j "$2" -o Query/$query.blast_out
# legacy_blast blastpgp -i "$1" -d "$3" -j "$2" -o data/queries/$query/$query.blast_out &>> log/querypgp.log

# ./bin/salut_1.0/src/blasttofaa.py data/queries/$query/$query "$2"

# # echo "Alignement multiple avec MUSCLE"
# ./bin/muscle3.8.31_i86linux64 -in data/queries/$query/$query.faa -out data/queries/$query/$query.mfasta &>> log/query.log

# echo "Création de la pssm"
# ./bin/salut_1.0/src/pssm_query.py data/queries/$query/$query $query &>> log/query.log

# echo "Prediction de la structure secondaire avec psipred"
# ./bin/psipred.4.02/psipred/runpsipred_single "$1" >> log/query.log
# mv -f $query.horiz tmp_query
# rm -f $query.ss2 $query.ss


# ######Création du Foldrec#####

# mquery=`echo pssm_query/$query.aamtx`

# for mtemplate in pssm_homstrad/* ; do
#     nametemplate=`echo ${mtemplate##*/} | cut -d. -f1`
#     echo "Comparaison "$query" <---> "$nametemplate
#     ./bin/salut_1.0/src/pssm_comparison.py $mquery $mtemplate >> log/comparaison.log
# done

# ./bin/salut_1.0/src/foldrec.py data/queries/$query/$query.aln data/queries/$query/$query.fasta $query
# cat data/foldrec/$query*part[0-9].foldrec > data/queries/$query/$query.foldrec
# rm -f data/foldrec/$query*part?.foldrec


# rm -Rf alignments tmp_query tmp_homstrad pssm_query pssm_homstrad
# rm -f log/tmp.log dsspAll_out.txt
