
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

def parsing_metafold:
	entrée : metafold file
	sortie: dictionnaire { metafold_1 : pdb , metafold_2 : pdb, ... }


def calc_threading_score:
	1) matrice séquence cible / séquence cible
	2) calcul des distances entre chaque aa
	3) sélection des distances entre 5 et 10 A
	4) transformer les distances en energie (utiliser matrice DOPE)
	sortie : somme de la matrice

```

## Documentation

* Nous allons générer la documentation avec Sphinx avec le thème Read The Docs (exemple: https://docs.readthedocs.io/en/latest/)


### Docstrings


* Les docstrings doivent être écrites en reStructuredText et formatées de la façon suivante pour être *parsées* automatiquement par Sphinx:
```python
def fonction_utile(arg1, arg2, ...):
    """
        Petite description de la fonction

        Args:
    	    arg1: description de l'argument 1
    	    arg2: description de l'argument 2
    	    ...

        Returns:
	    type_de_la_variable_retournée: description de la variable retournée
    """
```
