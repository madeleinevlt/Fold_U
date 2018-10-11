
# Pseudo-code ~ Fold-U

## Ligne de commande
```
options :
	-i, --input : FILE foldrec
	-o, --output : FILE filname.out
	-n, --nbPDB : INT
______________________________________________
(scores implémentés)
	--coev : true/false
	--physics : true/false
	--acc : true/false
	...
```

## Classes
```
class Alignment:
	id_SCOP : string
	ppalign_score : int
	query : liste d'instances de la classe residue
	template : liste d'instances de la classe residue
	
	def threading_score(self):

class Residue:
	self.resNum = resNum
	self.resName = chain[resNum].get_resname()
	self.CA = chain[self.resNum]['CA'].get_vector()
	group : string
```

## Prototypes
```
def parsing_foldrec:
	entrée : Fichier foldrec, nombre de templates à générer
	sortie : liste d'instances de la classe alignment

	boucle des residus dans le fichier pdb: (voir récupérer données pdb/residus)
		si gap dans template:
			résidue = '-'
		sinon:
			récupérer les données du résidu i
	
	boucle 



```

## Récupérer données pdb/residus
```python
from Bio.PDB import *

p = PDBParser(QUIET=True) # QUIET=T : Warnings issued are suppressed
pdb = p.get_structure(opt.input,opt.input)

resList = [] # List of Residue instances

for chain in pdb.get_chains():
        # First residue of the current chain
        first = next(res.id[1] for res in chain.get_residues() if (res.id[1] < 99999))

        # Last residue of the current chain
        for res in chain.get_residues():
            if (res.get_id()[0] != ' '):
                break
			last = res.get_id()[1]

			for resNum in range(first,last+1):
				r = classes.Residue(chain, resNum)
				resList.append(r)
```