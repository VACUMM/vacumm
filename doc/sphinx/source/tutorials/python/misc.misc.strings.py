# -*- coding: utf8 -*-
# Opérations
print 'aaa'+"bbb"
#   -> aaabbb
print '-'*50
#   -> --------------------------------------------------

# Substitutions
# - ordonné
print '%i' % 6
#   -> 6
print '%05i %f' % (5, 7)
#   -> 00005 7.000000
# - par mots clés
print '[%(toto)s] %(result)i !' % {'toto':'ok', 'result':20}
#   -> [ok] 20 !
# - par variables grace a vars()
adjectif = 'pratique'
print "C'est %(adjectif)s !"%vars()
#   -> C'est pratique !

# Divers
# - centrage
title = (' %s '%'centrage'.title()).center(40,'#')
print '#'*len(title)+'\n'+title+'\n'+'#'*len(title)
#  -> ########################################
#  -> ############### Centrage ###############
#  -> ########################################
# - tests
print 'haha'.startswith('h')
#  -> True
print '01'.isdigit()
#  -> True
# - split-join
print ' '.join('a-b-c d-e'.split('-'))
#  -> a b c d e
