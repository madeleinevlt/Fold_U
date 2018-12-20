#!/bin/bash

#####Création des répertoires necessaires#####
rm -Rf Query pssm_query pssm_homstrad tmp_query tmp_homstrad suivi alignments foldrec
mkdir Query pssm_query pssm_homstrad tmp_query tmp_homstrad suivi alignments foldrec


#####Déplacement des fasta des queries#####

for query in "$1"* ; do
    cp $query ./Query
done


#####Traitements sur les queries#####

rm -f suivi/query.txt

for fastaquery in Query/*.fasta ; do

    query=`echo ${fastaquery##*/} | cut -d. -f1`
    echo "QUERY : "$query

    echo $query >> suivi/query.log

    echo "Psi-blast contre $3"
    legacy_blast blastpgp -i Query/$query.fasta -d db/"$3".fasta -j "$2" -o Query/$query.blast_out

    python3 bin/blasttofaa.py $query "$2" Query

    echo "Alignement multiple avec MUSCLE"
    ./muscle3.8.31_i86linux64 -in Query/$query.faa -out Query/$query.mfasta &>> suivi/query.log

    echo "Création de la pssm"
    python3 bin/pssm_query.py $query.mfasta $query &>> suivi/query.log

    echo "Prediction de la structure secondaire avec psipred"
    psipred.4.02/psipred/runpsipred_single Query/$query.fasta >> suivi/query.log
    mv -f $query.horiz tmp_query
    rm -f $query.ss2 $query.ss

done

#####PSSM sur HOMSTRAD#####

rm -f suivi/pssm_homstrad.log
echo "Création des pssm pour les templates"
for template in HOMSTRAD/* ; do
    nametemplate=`echo ${template##*/}`
    python3 bin/pssm_homstrad.py $template/*.map $nametemplate &>> suivi/pssm_homstrad.log

done


######Assignation structure secondaire sur HOMSTRAD#####

echo "Assignation de la structure secondaire pour les templates"
pdbMETAFOLD=`cut -d\  -f2 METAFOLD.txt`
for f in HOMSTRAD/*/*.atm ; do
    template=`echo $f | cut -d/ -f2`
    pdb=`echo ${f##*/} | cut -d. -f1`
    for item in ${pdbMETAFOLD[*]} ; do 
        if [ "$pdb" == "${item%%.*}" ] ; then       
            ./dssp-2.0.4-linux-amd64 -i ./$f -o dssp$template.out &>> suivi/tmp.log
            if [ -f "dssp$template.out" ] ; then
                python3 bin/parse_dssp.py dssp$template.out &>> suivi/tmp.log
            fi
        fi
    done
done

rm -f dssp*.out


######Création du Foldrec#####

for mquery in pssm_query/* ; do
    namequery=`echo ${mquery##*/} | cut -d. -f1`

    for mtemplate in pssm_homstrad/* ; do 
        nametemplate=`echo ${mtemplate##*/} | cut -d. -f1`
        echo "Comparaison "$namequery" <---> "$nametemplate
        python3 bin/pssm_comparison.py $mquery $mtemplate >> suivi/comparaison.log 
    done

    python3 bin/foldrec.py alignments/$namequery.aln Query/$namequery.fasta $namequery 
    cat foldrec/$namequery*part[0-9].foldrec > foldrec/$namequery.foldrec
    rm -f foldrec/$namequery*part?.foldrec

done

rm -Rf alignments tmp_query tmp_homstrad pssm_query pssm_homstrad
rm -f suivi/tmp.log dsspAll_out.txt
rm -f psitmp*.mtx
