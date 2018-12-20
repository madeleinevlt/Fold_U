#!/bin/bash

#####Création des répertoires necessaires#####
rm -Rf Query pssm_query pssm_homstrad tmp_query tmp_homstrad suivi alignments foldrec
mkdir Query pssm_query pssm_homstrad tmp_query tmp_homstrad suivi alignments foldrec


#####Déplacement du fasta dans Query####

cp "$1" ./Query



#####Traitements sur les queries#####

rm -f suivi/query.txt


fastaquery=`echo "$1"`
query=`echo ${fastaquery##*/} | cut -d. -f1`
echo "QUERY : "$query

echo $query >> suivi/query.log

echo "Psi-blast contre $3"
#./blast-legacy-2.2.26-2/bin/blastpgp -i Query/$query.fasta -d db/"$3".fasta -j "$2" -o Query/$query.blast_out
legacy_blast blastpgp -i Query/$query.fasta -d db/"$3".fasta -j "$2" -o Query/$query.blast_out &>> suivi/querypgp.log

python3 bin/blasttofaa.py $query "$2" Query

echo "Alignement multiple avec MUSCLE"
./muscle3.8.31_i86linux64 -in Query/$query.faa -out Query/$query.mfasta &>> suivi/query.log

echo "Création de la pssm"
python3 bin/pssm_query.py $query.mfasta $query &>> suivi/query.log

echo "Prediction de la structure secondaire avec psipred"
psipred.4.02/psipred/runpsipred_single Query/$query.fasta >> suivi/query.log
mv -f $query.horiz tmp_query
rm -f $query.ss2 $query.ss


######Création du Foldrec#####

mquery=`echo pssm_query/$query.aamtx`

for mtemplate in pssm_homstrad/* ; do
    nametemplate=`echo ${mtemplate##*/} | cut -d. -f1`
    echo "Comparaison "$query" <---> "$nametemplate
    python3 bin/pssm_comparison.py $mquery $mtemplate >> suivi/comparaison.log
done

python3 bin/foldrec.py alignments/$query.aln Query/$query.fasta $query
cat foldrec/$query*part[0-9].foldrec > foldrec/$query.foldrec
rm -f foldrec/$query*part?.foldrec


rm -Rf alignments tmp_query tmp_homstrad pssm_query pssm_homstrad
rm -f suivi/tmp.log dsspAll_out.txt
