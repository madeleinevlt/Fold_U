
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

def parsing_DOPE:
	entrée : DOPE file
	sortie : matrice DOPE = matrice de listes des valeurs DOPE pour chaque couple d'acides aminés (Carbones alpha (CA))

def calc_threading_score:
	1) matrice séquence cible / séquence cible
	2) calcul des distances entre chaque aa
	3) sélection des distances entre 5 et 10 A
	4) transformer les distances en energie (utiliser matrice DOPE)
	sortie : somme de la matrice

```
