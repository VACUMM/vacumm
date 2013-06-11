.. _appendix.conventions:

Coding rules
************

Coding norm : PEP 8
===================

Les conventions de programmation choisies suivent notamment les normes
spécifiées par la :pep:`8`.
Parmi ces normes, citons les points importants suivants :
    
    - Les indentations sont composées de quatre espaces (pas de tabulation).
    - Les noms de classes et exceptions ont des majuscules en début de mot,
      et aucun caractère n'est utilisé pour séparer les mots dans le cas
      des noms composés. Exemple : ``MySuperClass``.
    - Les autres noms (fonctions, méthodes, attributs) sont entièrement en
      minuscules et le caractère underscore ``_`` est utilisé pour séparer les mots.
      Exemple : ``my_super_function``.
    - L'encodage doit être en UTF-8.
    
Les points suivants sont d'importance secondaire :
    
    - Les lignes ne doivent pas dépasser 79 caractères.
    - Des espaces doivent être utilisés pour séparer les opérateurs.
     
     
Docstrings
==========

Les docstings sont les commentaires sous forme de chaîne de caractères situés 
immédiatement après une déclaration.
Elles servent à renseigner sur l'objet déclaré (but, paramètres, sorties...), 
et ces informations vont apparaître donc la documentation de la librairie.

Les docstrings doivent être rédigées en langage de programmation `rst <http://docutils.sourceforge.net/rst.html>`_,
utilisé par la générateur de documentation `Sphinx <http://sphinx.pocoo.org>`_ .
La langue utilisée est l'anglais.

Nous proposons la forme suivante : ::
    
    def myfunc(a, b, c=None, d=5, **kwargs):
        """Short description of myfunc
        
        More information about myfunc.
        You can even include images.
        
        .. note:: A simple note.
        
        :Params:
            
            - **a**: First mandatory parameter.
              Its description continues here.
            - **b**: The second mandatory one.
            - **c**, optional: First optional keyword parameter.
            - **d**, optional: Second one.
            - Other parameters.
            
        :Return:
            
            - *x*: First return variable.
            - *y*: Second one.
            
        :Example:
            
            >>> x, y = myfunc(a, b, d=4, g='name')
        """
        ...
            
