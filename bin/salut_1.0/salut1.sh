

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
