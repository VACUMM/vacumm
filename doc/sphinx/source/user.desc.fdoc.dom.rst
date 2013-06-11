Extension fortran à :mod:`sphinx.domains`
=========================================

Le module :mod:`~vacumm.sphinxext.fortran_domain` est
une extension à Sphinx introduisant le "domaine" fortran (voir 
la :sphinx:`documentation de sphinx<domains.html>` sur les domaines syntaxiques),
et donc permettant de documenter du code fortran.

.. highlight:: rst

Configurer Sphinx
-----------------

Il suffit de d'ajouter :mod:`vacumm.sphinxext.fortran_domain`
à la liste décrite par la variable ``extension`` dans votre fichier
de configuration :file:`conf.py` (on suppose que le paquet python :mod:`vacumm` est accessible).

Syntaxe
-------

Cette extension fournit des "directives" pour déclarer (décrire et référencer) 
des entités
fortran (programme, module, type, function, subroutine et variable),
ainsi que de "roles" pour faire référence à des entité déclarées.

Les directives
~~~~~~~~~~~~~~


Toutes les directives sont préfixées par ``f:`` pour spécifier le domaine fortran.

.. note::
    
    Par la suite, les crochets ``[]`` désignent un argument optionel.
    
Toutes les directives acceptent un contenu permettant de décrire
l'entité, et qui sera interpété par sphinx.
Cette description peut en outre tirer partie des :ref:`docfields <docfields>`
pour décrire les arguments des fonctions, subroutines et programmes.

.. rst:directive:: .. f:program:: progname

    Description d'un programme fortran.
    Cette directive accepte des :ref:`docfields <docfields>`.   
    
    Example : ::
        
        .. f:program:: main
            
            Il s'agit du programme principal.
    
    
.. rst:directive:: .. f:module:: modname

    Cette crée une référence à un module
    et ne produit aucune sortie.
    Elle accepte l'option ``:synopsis:`` qui permet de décrire
    brièvement le module dans l'index des modules.
    Elle définit en outre le module courant
    (à l'instar de :rst:dir:`f:currentmodule`).
    
    Example : ::
        
        .. f:module:: monmodule
            :synopsis: Module de statistiques
                
    
.. rst:directive:: .. f:currentmodule:: [modname]

    Cette directive ne produit aucune sortie ; 
    elle fait de ``modname`` le module courant :
    les fonctions, subroutines, types et variables décrits
    par les suite seront considérés comme appartenant à ce module.
    Si ``modname`` est vide ou égal à ``None``, il n'y a plus
    de module courant.
    
    Example : ::
        
        .. f:currentmodule:: mymodule
        
        .. f:variable:: float myvar
        
        .. f:currentmodule::
            
            
    
.. rst:directive:: .. f:type:: [~][modname][/]typename

    Cette directive décrit un type dérivée dans un module.
    Elle accepte le docfield spécial ``:f ...:`` pour
    décrire les champs du module.
    
    Example : ::
        
        .. f:type:: mymod/mytype
        
            :f integer var1: Variable 1
            :f float var2: Variable 2
    
   
.. rst:directive:: .. f:variable:: [type] [~][modname][/]varname[(shape)]

    Cette directive décrit une variable d'un module.
    Elle accèpte les options suivantes :
    
        - ``:type:``: Type de la variable (``float``, ``mytype``, etc). 
          Si présent, un lien est créé vers ce type. 
          Le type peut aussi être indiqué avant le nom de la variable.
        - ``:shape:``: Dimension de la forme ``nx,ny``,
          qui peut aussi être déclarée après le nom (entre parenthèses).
          Une référence à toutes les variables trouvées est créée.
        - ``:attrs:``: Attributs supplémentaires (``in``, ``parameter``, etc).
        
    Exemple : ::
        
        .. f:function:: float myvar
            :shape: nx,ny
            :attrs: in
            
            Description de ma variable.
            
    
.. rst:directive:: .. f:function:: [type] [~][modname][/]funcname(signature)

    Cette directive décrit une function appartenant ou non à un module.
    Elle accèpte l'option `:type:` et utilise les :ref:`docfields <docfields>`
    pour décrire ses arguments, ses appels et ses modules utilisés.
    
    Exemple : ::
        
        .. f:function:: myfunc(a [,b])
        
            Ceci est ma fonction principale.
            
            :p float a: Mon premier argument.
            :o integer b [optional]: Mon deuxième argument.
            :from: :f:subr:`mysub`
   
    
.. rst:directive:: .. f:subroutine:: [~][modname][/]subrname[(signature)]

     Cette directive décrit une subroutine à l'instar de la directive
     :rst:dir:`f:function`.

    Exemple : ::
        
        .. f:subroutine:: mysub(a)
        
            Description.
        
            :param a: Mon paramètre.
            :to: :f:func:`myfunc`
    
Les roles
~~~~~~~~~

Les roles permettent en langage rst de faire référence à des entités (programme, function, etc).
On les utilise avec une syntaxe du type : ::
    
    :role:`cible`
    :role:`nom <cible>` # avec nom alternatif
    
Plusieurs ont été définis relativement aux directives fortran présentées ci-dessus.

.. rst:role:: f:prog

    Référence à un programme déclaré par :rst:dir:`f:program`.

.. rst:role:: f:mod:
    
    Référence à un module déclaré par :rst:dir:`f:module`.
    
.. rst:role:: f:type:
    
    Référence à un type dérivé déclaré par :rst:dir:`f:type`.
    
.. rst:role:: f:var:
    
    Référence à une variable déclaré par :rst:dir:`f:variable`.

.. rst:role:: f:func

    Référence à une fonction ou une subroutine déclarés respectivement par :rst:dir:`f:function` et :rst:dir:`f:subroutine`.

.. rst:role:: f:subr:
    
    Alias pour :rst:role:`f:func`.
    

Il est possible de faire féférence à des types dérivés, variables, 
fonctions et subroutines appartenant à un module, avec la syntaxe typique suivante : ::
    
    :f:func:`monmodule/mafunction`
    
Si une fonction locale et celle d'un module ont le même nom,
il est possible d'utiliser la syntaxe suivante si le module courant est celui de la fonction : ::
    
    :f:func:`/mafunction`
    
Si le `"/"` est omis, la référence portera sur la fonction locale et non celle du module courant.

Si le module est précisé dans la référence, il est possible de ne pas le faire
afficher en préfixant le nom par ``"~"`` : ::
    
    :f:func:`~monmodule/mafunction`
    
.. _docfields: 

Les *docfields*
~~~~~~~~~~~~~~~

Les docfields sont des balises rst de type 
:rstdoc:`field list <restructuredtext.html#field-lists>`,
qui sont interprétées dans le contenu de certaines directives
afin de décrire des paramètres, options et autres champs particuliers.

Le domaine fortran en permet une utilisation assez proche de :sphinx:`celle
implémentée <domains.html#info-field-lists>` pour le :sphinx:`domaine python <domains.html#the-python-domain>`.

Il existe deux familles de *docfields* fortran : celle des fonctions et subroutines, et celle des types dérivés.


.. rubric:: Famille des fonctions et subroutines

Pour cette famille, les *docfields* permettent de 
décrire les arguments obligatoires et optionnels, les modules utilisés,
les programmes, fonction et subroutines qui appellent l'entité courantes,
et les fonctions et subroutines appelées par cette entité.
Certains d'entre eux possèdent des alias.

- ``param`` (ou ``p``, ``parameter``, ``a``, ``arg``, ``argument``) : Argument obligatoire. ::
    
    :param myvar: Description.
    
  Il est possible de spécifier le type, les dimensions et les attributs spéciaux
  lors de la déclaration en suivant l'exemple ci-dessous : ::
      
      :param mytype myvar(nx,ny) [in]: Description.
       
- ``type`` (ou ``paramtype``, ``ptype``) : Type du paramètre (exemple : float). Une référence est faite à ce paramètre s'il est présent. ::
    
    :type: float
    
- ``shape`` (out ``pshape``) : Dimensions du paramètre (séparées par des virgules).

    :shape: nx, ny
    
- ``attrs`` (ou ``pattrs``, ``attr``) : Attributs particuliers (séparés par des virgules). ::
    
    :attr: in/out
    
- ``option`` (ou ``o``, ``optional``, ``keyword``) : Déclaration d'un paramètre optionnel, avec une syntaxe similaire à celle des paramètres obligatoires (``param``).
- ``otype`` (ou ``optparamtype``) : Son type..
- ``oshape`` : Ses dimensios.
- ``oattrs``` (ou ``oattrs``) : Ses attributs.
- ``return`` (ou ``r``, ``returns``) : Variable retournée par la fonction.
- ``rtype`` (ou ``returntype``) : Son type.
- ``rshape`` : Ses dimensions.
- ``rattrs` (ou ``rattrs``) : Ses attributs.
- ``calledfrom`` (ou ``from``) : Fonctions, subroutines ou programmes appelant la routines courante.
- ``callto`` (ou ``to``) : Fonctions ou subroutines que la routines courante appelle.

Les **docfiels** de paramètres sont fusionnés en une liste pour ceux obligatoires et une pour ceux optionnels, et leur **docfiels** associés (type, dimensions et attributs) sont supprimés et leur contenus insérés dans la déclaration du paramètre.

.. note::
    En complément de ces *docfiels* interprétés, vous pouvez les accompagner
    d'autres de votre choix. Par exemple : ::
        
        :actions: Cette fonction procède à
            
            #) Une initialisation.
            #) Un calcul.
            #) Un plot.
            
        :p float myvar [in]: Variable à tracer.
        :use: Fait appel au module :f:mod:`mymod`.


.. rubric:: Famille des types dérivés.

Ces *docfields* décrivent les champs des types dérivés.
Il sont à insérér dans l'entête d'une directive :rst:dir:`f:type`.

    
- ``ftype`` (ou ``f``, typef, typefield) : Champs d'un type dérivé, avec une la syntaxe similaire à celle des paramètres obligatoires (``param``).
- ``ftype`` (ou ``fieldtype``) : Son type..
- ``fshape`` : Ses dimensions.
- ``fattrs``` (ou ``fattrs``) : Ses attributs.

Exemple
-------


.. literalinclude:: user.desc.fdoc.dom.sample.txt
    :language: rst

 
*Ce qui donne :*
    
    .. include:: user.desc.fdoc.dom.sample.txt
    
    
.. note :: 
    Les modules déclarés sont listés dans leur :f:ref:`index <f-modindex>`, et les autres entités fortrans sont aussi accesibles depuis l':ref:`index général <genindex>`.

