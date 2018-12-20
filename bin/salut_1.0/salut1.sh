

#####PSSM sur HOMSTRAD#####

rm -f log/pssm_homstrad.log
echo "Création des pssm pour les templates"
for template in "$1"/* ; do
    nametemplate=`echo ${template##*/}`
    ./bin/salut_1.0/src/pssm_homstrad.py $template/*.map $nametemplate &>> log/pssm_homstrad.log

done


######Assignation structure secondaire sur HOMSTRAD#####

echo "Assignation de la structure secondaire pour les templates"
pdbMETAFOLD=`cut -d\  -f2 "$2"`
for f in "$1"/*/*.atm ; do
    template=`echo $f | cut -d/ -f2`
    pdb=`echo ${f##*/} | cut -d. -f1`
    for item in ${pdbMETAFOLD[*]} ; do
        if [ "$pdb" == "${item%%.*}" ] ; then
            ./dssp-2.0.4-linux-amd64 -i ./$f -o dssp$template.out &>> log/tmp.log
            if [ -f "dssp$template.out" ] ; then
                ./bin/salut_1.0/src/parse_dssp.py dssp$template.out &>> log/tmp.log
            fi
        fi
    done
done

rm -f dssp*.out
