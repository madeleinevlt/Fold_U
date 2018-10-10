
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
class alignment :
	query : list of residue class
	template : list of residue class
	score : int

	def threading_score(self):

class residue :
	CAcoords = np.array[(x,y,z)]
	resName : string
	resNum : int
	group : string

```

## Prototypes
```
def parsing_foldrec:
	entrée : Fichier foldrec, nombre de templates à générer
	sortie : liste d'instances de la classe alignment
```
