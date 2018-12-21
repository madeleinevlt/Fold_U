#!/bin/bash

rm -Rf log
mkdir log

#####PSSM sur HOMSTRAD#####

rm -f log/pssm_homstrad.log
echo "Création des pssm pour les templates"
for template_path in "$1"/* ; do
    template_name=`echo ${template_path##*/}`
    python3 bin/salut_1.0/src/pssm_homstrad.py $template_path/$template_name.map $template_name &>> log/pssm_homstrad.log
done


######Assignation structure secondaire sur HOMSTRAD#####

echo "Assignation de la structure secondaire pour les templates"
pdbMETAFOLD=`cut -d\  -f2 "$2"`
for f in "$1"/*/*.atm ; do
    template=`echo $f | cut -d/ -f2`
    pdb=`echo ${f##*/} | cut -d. -f1`
    for item in ${pdbMETAFOLD[*]} ; do
        if [ "$pdb" == "${item%%.*}" ] ; then
            ./bin/dssp-2.0.4-linux-amd64 -i ./$f -o dssp$template.out &>> log/dssp.log
            if [ -f "dssp$template.out" ] ; then
                python3 bin/salut_1.0/src/parse_dssp.py dssp$template.out &>> log/dssp.log
            fi
        fi
    done
done

rm -f dssp*.out
