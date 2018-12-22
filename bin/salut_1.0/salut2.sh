#!/bin/bash

#####Traitements sur les queries#####

fasta_query=`echo "$1"`
query=`echo ${fasta_query##*/} | cut -d. -f1`
path_query=`echo "data/queries/$query/$query"`

echo "Psi-blast against $3"
if [ -f $path_query.blast_out ] ; then
    echo "File $path_query.blast.out already exists"
else
    legacy_blast blastpgp -i "$1" -d "$3" -j "$2" -o $path_query.blast_out &>> log/querypgp.log
fi

echo "Convert blast to faa"
./bin/salut_1.0/src/blasttofaa.py $path_query "$2"

echo "Multiple alignment calculated with MUSCLE"
./bin/muscle3.8.31_i86linux64 -in $path_query.faa -out $path_query.mfasta &>> log/query.log
./bin/salut_1.0/src/mfasta.py $query

echo "PSSM file creation"
./bin/salut_1.0/src/pssm_query.py $path_query $query &>> log/query.log

echo "Secondary structure prediction using Psipred"
./bin/psipred.4.02/psipred/runpsipred_single "$1" >> log/query.log
mv -f $query.horiz data/queries/$query
rm -f $query.ss2 $query.ss


echo "PSSM comparison between the query and the templates"
mquery=`echo $path_query.aamtx`
if [ -f $path_query.aln ] ; then
    echo "File $path_query.aln already exists"
else
    for mtemplate in data/HOMSTRAD/*/*.aamtx ; do
        nametemplate=`echo ${mtemplate##*/} | cut -d. -f1`
        echo $query" <---> "$nametemplate" comparaison"
        ./bin/salut_1.0/src/pssm_comparison.py $mquery $mtemplate >> log/comparaison.log
    done
fi

# echo "Write foldrec file"
# ./bin/salut_1.0/src/foldrec.py $path_query.aln $path_query.fasta $query
# cat $path_query*part[0-9].foldrec > $path_query.foldrec
# rm -f $path_query*part?.foldrec
