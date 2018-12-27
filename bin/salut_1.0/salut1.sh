#!/bin/bash

rm -Rf log
mkdir log

#####PSSM sur TEMPLATES#####

rm -f log/pssm_template.log
echo "PSSM files creation for each template"
for template_path in "$1"/* ; do
    template_name=`echo ${template_path##*/}`
    file_name=`echo "$template_path/$template_name.map"`
    ./bin/salut_1.0/src/pssm_template.py $template_path/$template_name $template_name &>> log/pssm_template.log
done


######Assignation structure secondaire sur TEMPLATES#####

mkdir tmp_dssp
echo "Secondary structures assignment for each template"
pdbMETAFOLD=`cut -d\  -f2 "$2"`
for template_path in "$1"/*; do
    template=`echo "$template_path" | rev | cut -d'/' -f1 | rev`
    for f in $template_path/*.atm ; do
        pdb=`echo ${f##*/} | cut -d. -f1`
        for item in ${pdbMETAFOLD[*]} ; do
            if [ "$pdb" == "${item%%.*}" ] ; then
                ./bin/dssp-2.0.4-linux-amd64 -i ./$f -o tmp_dssp/dssp$template.out &>> log/dssp.log
                if [ -f "tmp_dssp/dssp$template.out" ] ; then
                    ./bin/salut_1.0/src/parse_dssp.py tmp_dssp/dssp$template.out &>> log/dssp.log
                fi
            fi
        done
    done
done
rm -fr tmp_dssp
