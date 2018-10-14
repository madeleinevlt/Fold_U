# Meeting 1 : lundi 1/10  

## Répartition des rôles  

* **Manager** ~ Hélène
 - Organisation des meetings
 - Gestion de l'avancement du projet
 - Gestion du temps

* **Expert code/github** ~ Gabriel/Arnold
 - Répartition du code entre contributeurs
 - Gestion des branches dans GitHub
 - Vérification du code
 - Choix des classes / fichiers

* **Expert biblio** ~ Flora/Arnold
 - Recherche d'articles intéressants

* **Expert rapport/présentation** ~ Tom
 - Bon niveau en anglais / correction anglais
 - Gestion du rapport en latex (overleaf)
 - Gestion présentation diapos

## Commandes GitHub

#### Cloner le repository
```shell
git clone https://github.com/meetU-MasterStudents/2018---2019-Equipe-4.git
```

#### Ajout / Modification fichiers
```shell
cd repository
git add nom_du_fichier
git commit -m "commentaire"
git push
```

#### "Tirer" le repository
```shell
cd repository
git pull
```

### Create a new repository  
```shell
git init
git add README.md
git commit -m "first commit"
git remote add origin https://github.com/...
git push -u origin master
```
### Astuce : Sauvegarde des log de GitHub pendant 8h

Dans le ~/.gitconfig :
```
[user]
	email = mail@gmail.com
	name = pseudo
[core]
	pager = less
[url "https://pseudo@github.com"]
	insteadOf = https://github.com
[credential]
helper = cache --timeout=28800
```


### Pour générer la documentation

Si vous avez codé des fonctions et que vous avez fais des jolies *docstrings*,
ou même si vous avez modifié le README ou que vous voulez expliquer des choses
dans la doc il suffit de la mettre à jour en exécutant dans un terminal:
```
cd docs
make html
google-chrome-stable build/html/index.html
ou
firefox build/html/index.html
```
