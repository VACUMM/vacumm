# -*- coding: utf8 -*-
# On initialise le document Excel
from xlwt import Workbook, Formula
w = Workbook()

# Définition des styles
from vacumm.misc import xls_style
from copy import deepcopy as cp
style_header = xls_style(b=1, ha='center', bottom=2, top=2, c=2)
style_left = xls_style(ha='left', i=1)
style_num = xls_style(fmt='0.00')
style_bottom = xls_style( top=2)

# Première feuille
# - nom
wc = w.add_sheet('Courants')
# - données
data1 = {
    'Four':[1.15646541, 0.25641],
    'Pointe du Raz':[2.56421, 0.1556122],
}
# - taille
wc.col(0).width = int(wc.col(0).width * 1.5)
wc.col(1).width = int(wc.col(0).width * 0.7)
# - entêtes
for icol, title in enumerate(('','STD [m/s]', 'RMS [%]')):
    wc.write(0, icol, title, style_header)
# - colonne de gauche et données
for irow, title in enumerate(data1.keys()):
    # titre
    wc.write(irow+1, 0, title, style_left)
    # données
    for icol, value in enumerate(data1[title]):
        # RMS en %
        if icol == 1:
            col_style = xls_style(cp(style_num), fmt='0.00%')
        else:
            col_style = style_num
        # écriture
        wc.write(irow+1, icol+1, str(value), col_style)
# - bottom
for icol in xrange(3): wc.write(irow+2, icol, '', style_bottom)

# Deuxième feuille
# - nom
wn = w.add_sheet('Testouille')
# - gros titre
wn.write_merge(0, 1, 0, 5, u"Ça c'est un gros titre", xls_style(o=1, s=500))
# - hyperlien
uri = "http://relay.vacumm.fr/~raynaud/pydoc"
urn = "Python Actimar"
wn.write_merge(3,3,1,10,Formula('HYPERLINK("%s";"%s")'%(uri, urn)),xls_style(c=4))
# - groupes de niveaux
wn.write(4, 1, 'level1')
wn.write(5, 1, 'level2')
wn.write(6, 1, 'level2')
wn.write(7, 1, 'level1')
wn.row(4).level = wn.row(7).level = 1
wn.row(5).level = wn.row(6).level = 2
# - dates
from datetime import datetime
wn.write(8, 0, 'Dates :')
wn.write(8, 1, datetime.now(), xls_style(fmt='h:mm:ss AM/PM'))
wn.write(8, 2, datetime.now(), xls_style(fmt='MMM-YY'))


# Écriture du document
w.save('misc.io.xls.xls')


