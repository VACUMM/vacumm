.. _user.faq.graphics:

Graphiques
==========

Comment modifier les paramètres de tracé ?
------------------------------------------

Modifiez le fichier de configuration de Matplotlib.
Il est typiquement situé ici : :file:`~/.matplotlib/matplotlibrc`.

Pour en savoir plus : http://matplotlib.sourceforge.net/users/customizing.html


La tracé de ma carte prend du temps !
-------------------------------------

Les tracés de carte qui prennent du temps sont ceux d'une zone réduite
qui font appel à un trait de côte à haute résolution ("f" du GSHHS ou "s" pour Histolitt du SHOM).
Dans le cas des traîts du GSHHS, un système de mise en cache des cartes 
tracées avec :mod:`~mpl_toolkits.basemap` est utilisé 
(voir :mod:`vacumm.misc.grid.basemap`) : si une carte utilisant le même zoom et la même résolution
de traît de côte doit être retracée, celle-ci est chargée à partir du cache.
Pour info, les cartes en cache sont stockées dans le répertoire donné par
la fonction :func:`~vacumm.misc.grid.basemap.get_map_dir` (a priori
:file:`~/.matplotlib/basemap/cached_maps`).

.. note::
    
    Ce système de cache peut générer des erreurs au moment de l'utilisation
    d'une carte sauvegardée. Voir :ref:`user.faq.graphics.err_basemap`.
    
En outre, il se peut qu'une carte veuille s'afficher à l'écran
alors que vous faites un tracé à distance, ce qui peut prendre du temps.
Si vous désirez désactiver cet affichage, lisez :
:ref:`user.faq.graphics.window`.



Comment tuner un graphique à tracer ?
-------------------------------------

Les routines graphiques de VACUMM transfèrent des mots clés vers 
la plupart des fonctions de matplotlib utilisées.
On spécifie un mot-clé à transférer un préfixant celui-ci du nom de
la fonction suffixé avec un "_" (on utilise :func:`~vacumm.misc.misc.kwfilter` pour
filtrer les mots clés).

Par exemple, il est possible de changer la taille du titre d'un graphique,
ainsi que l'épaisseur des contours en procédant ainsi :
    
    >>> map2(sst, title_size=20, contour_linewidths=3)

Dans ce cas précis, pour avoir la liste des autres options possibles
on se réfère à la documentation des fonctions :func:`matplotlib.pyplot.title`
et :func:`matplotlib.pyplot.contour`.


Comment placer soit-même une colorbar ?
---------------------------------------

Désactivez le tracé automatique de la colorbar et l'affichage de la figure, 
puis tracez-la en spécifiant sa position (cf. :func:`matplotlib.pyplot.axes`) :
    
    >>> m = map2(data, colorbar=False, right=.9, show=False)
    >>> m.colorbar(cax=[.91, .2, .3, .6])
    >>> m.show()
    
ou directement

    >>> map2(data, colorbar_cax==[.91, .2, .3, .6], right=.9, show=False)
    
Voir :func:`matplotlib.pyplot.colorbar` pour les options de tracé.


.. _user.faq.graphics.window:

Ma figure ne s'affiche pas / désactiver l'affichage des figures
---------------------------------------------------------------

La créations des figures nécessite de faire appel à une couche logicielle graphique.

Ces couches logicielles sont appellées backends. Le minimum pour la
sauvegarde des figures est d'utiliser le backend ``Agg``, cela permet
aussi de contourner les erreurs lors de l'utilisation de scripts en batch
lorsque l'affichage n'est pas disponible.

Les autres backends nécesitent un affichage, le backend ``TkAgg``
est presque systématiquemnt disponible.
Le backend ``qt4agg`` est également disponible avec UVCDAT.

Si vous voulez faire un affichage de vos figures, utilisez donc ``qt4agg``,
si vous utiliser ``Agg``, aucun affichage ne sera fait.

Il y a deux manières de faire pour modifier le backend:
    
    1) Celle permanente consiste à modifier la ligne correspondante
       dans le fichier de configuration de Matplotlib ::
           
           backend      : Agg
           
    2) La deuxième manière consiste à le faire à la volée,
       au début d'un script (avant tout chargement de :mod:`matplotlib`
       et :mod:`vacumm`) ::
           
           from matplotlib import use
           use('Agg')

Comment désactiver le zoom automatique ?
----------------------------------------

Par exemple quand vous tracez une carte avec :func:`~vacumm.misc.plot.map2`,
si les données non masquées ne couvrent pas tout le domaine, 
un zoom automatique évitera de tracer les zones non couvertes.

Vous pouvez éviter cela en utilisant les options ``xmasked=False``,
``ymasked=False``, voir ``xymasked=False`` pour désactiver le zoom sur les deux axes
(dans le cas d'un plot 2D).

    >>> map(sst, xymasked=False)

Comment ajouter des effets de type ombre ou glow ?
--------------------------------------------------

La plupart du temps, des mots clés sont prévus à cet effet pour les
fonctions de tracer. Par exemple, si vous tracez des contours, vous
pouvez procéder ainsi pour ajouter une ombre à ces derniers et un effets
glow aux labels :
    
    >>> map2(data, contour_shadow=True, clabel_glow=True)

Vous pouvez passer des mots clés pour modifier les effets (voir
:func:`~vacumm.misc.core_plot.add_shadow` and :func:`~vacumm.misc.core_plot.add_glow`
pour connaitre les options) :
    
    >>> map2(data, contour_shadow_width=4, contour_shadow_xoffset=3)
    
Vous pouvez procéder aussi à la main pour un tracé externe :
    
    >>> m = map2(data, show=False)
    >>> m.add_shadow(m.axes.plot(x,y)) # m.axes.plot ou pylab.plot
    
voire :
    
    >>> from vacumm.misc.core_plot import add_shadow
    >>> add_shadow(P.plot(x,y))


.. _user.faq.graphics.err_basemap:

J'ai une erreur liée à la classe :class:`~mpl_toolkits.basemap.Basemap`
-----------------------------------------------------------------------

Il se peut que l'utilisation d'une carte (typiquement lors d'une tracé avec
:func:`~vacumm.misc.plot.map2`) produise une erreur (par exemple un attribut manquant,
tel que :attr:`celestial`` ou :attr:`_mapboundarydrawn`).
C'est probablement lié à la mise en cache automatique des cartes déjà tracées.
Ce procédé permet de gagner du temps lors du tracé d'une carte ayant exactement les mêmes
caractéristiques (domaine, projection, trait de côte).
Or, la carte mise en cache n'est peut-être pas compatible avec la version actuelle de 
:mod:`~mpl_toolkits.basemap`, suite à une mise à jour des paquets python.
Pour régler le problème radicalement, vous pouvez supprimer toutes les cartes en cache.
Elles se trouvent par défaut dans le répertoire donné par la fonction :func:`~vacumm.misc.grid.basemap.get_map_dir`.
Il existe une fonction pour effectuer cette opération :
    
    >>> from vacumm.misc.grid.basemap import reset_cache
    >>> reset_cache()
