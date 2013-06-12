import numpy as N, MV2
MA = N.ma

# Creation d'un tableau numerique pur
# - plein de zeros (pareil avec ones)
z = N.zeros((3,4))
print z.shape
#  -> (3, 4)
# - avec des valeurs choisies
a = N.array([1.,2,3]) # un des element est reel, donc le tableau aussi
print a
#  -> [ 1., 2., 3.,]

# On cree un tableau maske
# - pas de masquage
a2 = MA.array(a)
# - on maske exactement les 2.
b = MA.masked_object(a,2.)
# - pareil mais avec une marge possible (voir help(masked_values))
print MA.masked_values(a,2.)
# > [1.0 -- 3.0 ]
# - autre selction
print MA.masked_outside(a,0.,2.)
#  -> [1.0 2.0 -- ]

# Masks : 1 = valeur masquee !
# - standard
print b.mask
#  -> [False, True, False] 
# - quand une variable n'a pas de valeur masquee
print a2.mask
#  -> False
print a2.mask is N.ma.nomask, a2.mask is MV2.nomask
#  -> True, True
# - mais si on veut quand meme un tabeau
print MA.getmaskarray(a2)
#  -> [False False False]

# Copy or not copy ? 
#  Si on a copy=0, on fait seulement un "lien" des
#  valeurs numeriques => 
#   - ON UTILISE MOINS DE MEMOIRE
#   - attention aux suprises si le tableau d'origine est modifie
d = MA.masked_object(a,2.,copy=0)
a[0] = 10
print d
#  -> [10.0 -- 3.0]

# Une variable avec des axes
# - creation avec lien
m = MV2.array(b,copy=0,id='yoman')
# - modification des valeurs SANS CREER UN NOUVEL OBJECT !
m[:] = MA.masked_object(m,10.,copy=0)
print repr(m)
# > yoman
# > array(data = 
# >  [ 1.0 -- 3.0],
# >       mask = 
# >  [False True False],
# >      fill_value=2.)
print m
#  -> [-- -- 3.0]
m.info()
#  -> bla bla sur les attributes de la variable et ses axes

