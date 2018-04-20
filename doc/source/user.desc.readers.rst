.. _user.desc.readers:


Data readers
************



.. _user.desc.readers.info:

Grille des lecteurs et autres infos
====================================

Il serait utile d'ajouter des méthodes génériques
de récupération d'information sur les fichiers associés aux lecteurs.
Une des méthodes est la lecture de la grille : ::
    
    mars = MARS(config, dates)
    grid = mars.get_grid()
    
Elle peut-être sauvegardée en cache.
La date ne sert ici qu'a sélectionner les fichiers à lire.
La récupération de la grille peut alors se faire en cherchant
la première variable qui contient une grille : ::
    
    def get_grid(self):
        f = cdms2.open(self.ncfiles[0])
        for varname, varspec in f.variables.items():
            grid  = var.getGrid()
            if grid is not None: return grid




.. _user.desc.readers.sigma:


Sigma coordinates reader and converters
=======================================

Le but est d'adapter les formulations extantes de gestion des coordonnées sigma
afin de les intégrer dans les "lecteurs" de données (de type classe MARS).
Il est par exemple nécessaire de pouvoir s'adapter à la lecture de plusieurs fichiers.

Il existe le module :mod:`actimar.assimilation.sigma` (lui même adapté de
:mod:`pyinterpsig2z.sigma`) dont il faut s'inspirer.
Les classes de ce module sont explicitement attachée à un fichier,
dont elle tire le type de coordonnées sigma (via :meth:`~actimar.assimilation.sigma.Sigma.factory`), 
puis les différents paramètres
ou variables constantes dans le temps.
Le calcul des profondeurs se fait alors à un instant donné en lisant
eta, et les paramètres de profondeur pour une sélection de l'espace donné
(via par exemple :meth:`~actimar.assimilation.sigma.SigmaStandard.sigma_to_depths`).

Le but est alors de lire une seule fois le premier fichier pour
initialiser une classe de sigma en appelant :meth:`~actimar.assimilation.sigma.Sigma.factory`.
Ensuite, il doit être possible de fournir un autre nom (descripteur) de fichier
lors de la conversion de sigma vers les profondeurs.

Exemple : ::
    
    # Initialisation de l'instance sigma
    f1 = cdms2.open('mymarsfile1.nc')
    sigma = Sigma.factory(f1)
    f1.close()
    
    # Lecture des profondeurs
    from cdms2.selectors import Selector
    selector = Selector(lon=(-10, 0), lat=(43, 48))
    f2 = cdms2.open('mymarsfile2.nc')
    depths = sigma.sigma_to_depths(seletor, f=f2)
    
Le descripteur de fichier est sauvegardé dans l'attribut :attr:`Sigma.f`.
Il suffit simplement de le mettre à jour lors de l'appel de :meth:`~~actimar.assimilation.sigma.SigmaStandard.sigma_to_depths`,
et le tour devrait être joué.

L'intégration dans les classes de lecteur peut se faire par exemple
en ajoutant les méthodes suivantes :

- :meth:`load_sigma` : création d'une instance :class:`Sigma` avc la méthode :meth:`factory`.
- :meth:`load_sigma_depths` : afin de récupérer les profondeurs.
